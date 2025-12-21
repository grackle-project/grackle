//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines functionality for querying rate data.
///
//===----------------------------------------------------------------------===//

#include <cstring>  // std::strcmp, std::strlen, std::strcpy
#include <limits>
#include "grackle.h"
#include "internal_types.hpp"  // CollisionalRxnRateCollection
#include "LUT.hpp"             // CollisionalRxnLUT
#include "opaque_storage.hpp"  // gr_opaque_storage
#include "ratequery.hpp"
#include "status_reporting.h"

#include <algorithm>

// In comparison to the dynamic API for accessing elements of chemistry_data,
// we have explicitly opted NOT to make use of offsetof to access arbitrary
// values within a struct.
// -> while `offsetof` may lead to somewhat simpler code, the validity of using
//    it with pointer-arithmetic to access members of structs is murky at best.
//    - Things become murky when you consider that the C standard models an
//      abstract machine and doesn't place the strongest (or clearest)
//      constraints on a compiler when you start arbitrarily messing with
//      memory
//    - C's object model (and maybe also pointer provenance) is relevant
// -> An extended discussion can be found in this stackoverflow answer and in
//    the comments about the answer (https://stackoverflow.com/a/69936600).
//    - the most compelling point is that it would be useless for the standard
//      to define the `offsetof` if these semantics didn't work
// -> In any case, I'm much more skeptical that valid C++ code can make use of
//    offsetof in this fashion (the object model is stricter and C++ provides a
//    other machinery to accomplish similar things). Since we want to support
//    compilation with a C++ compiler, this is the most important point
//    offsetof in this fashion
//

namespace grackle::impl::ratequery {
// we have reserved the right to change this value at any time
enum {
  UNDEFINED_RATE_ID_ = std::numeric_limits<grunstable_rateid_type>::max()
};

// Create machinery to lookup Other Miscellaneous Rates
// ----------------------------------------------------
// -> ideally this should get relocated to a more sensible location (but not
//    sure where that is...

enum MiscRxnRateKind_ {
  MiscRxn_k24,
  MiscRxn_k25,
  MiscRxn_k26,
  MiscRxn_k27,
  MiscRxn_k28,
  MiscRxn_k29,
  MiscRxn_k30,
  MiscRxn_k31,

  MiscRxn_NRATES  // <- will hold the number of reactions
};

static Entry get_MiscRxn_Entry(chemistry_data_storage* my_rates, int i) {
  if (my_rates == nullptr) {  // <- shouldn't come up
    return mk_invalid_Entry();
  }
  switch (i) {
    // Radiative rates for 6-species (for external field):
    case MiscRxn_k24:
      return new_Entry(&my_rates->k24, "k24");
    case MiscRxn_k25:
      return new_Entry(&my_rates->k25, "k25");
    case MiscRxn_k26:
      return new_Entry(&my_rates->k26, "k26");
      // Radiative rates for 9-species
    case MiscRxn_k27:
      return new_Entry(&my_rates->k27, "k27");
    case MiscRxn_k28:
      return new_Entry(&my_rates->k28, "k28");
    case MiscRxn_k29:
      return new_Entry(&my_rates->k29, "k29");
    case MiscRxn_k30:
      return new_Entry(&my_rates->k30, "k30");
    case MiscRxn_k31:
      return new_Entry(&my_rates->k31, "k31");
    default: {
      return mk_invalid_Entry();
    }
  }
}
}  // namespace grackle::impl::ratequery

int grackle::impl::ratequery::RegBuilder_misc_recipies(
    RegBuilder* ptr, const chemistry_data* my_chemistry) {
  if (my_chemistry->primordial_chemistry != 0) {
    return RegBuilder_recipe_scalar(ptr, MiscRxn_NRATES, &get_MiscRxn_Entry);
  }
  return GR_SUCCESS;
}

namespace grackle::impl::ratequery {

/// calls delete or delete[] on the specified pointer (depending on the value
/// of is_scalar).
///
/// @note
/// This was originally function-like macro. We can go back to that if people
/// dislike this usage of a template function.
/// function-like macro
template <typename T>
static void careful_delete_(T* ptr, bool is_scalar) {
  (is_scalar) ? delete ptr : delete[] ptr;
}

/// This deallocates data within a list of owned Entries (i.e. all pointers
/// within an Entry are deleted)
///
/// Importantly, this does **NOT** call delete[] on entry_list
static void drop_owned_Entry_list_contents(Entry* entry_list, int n_entries) {
  for (int entry_idx = 0; entry_idx < n_entries; entry_idx++) {
    Entry* cur_entry = &entry_list[entry_idx];

    // invariant: each Entry within an embedded-list should be valid and
    //            fully initialized
    GR_INTERNAL_REQUIRE(cur_entry->name != nullptr, "sanity check!");
    GR_INTERNAL_REQUIRE(!cur_entry->data.is_null(), "sanity check!");

    delete[] cur_entry->name;

    bool is_scalar = cur_entry->props.ndim == 0;
    switch (cur_entry->data.tag()) {
      case PtrKind::const_f64: {
        careful_delete_(cur_entry->data.const_f64(), is_scalar);
        break;
      }
      case PtrKind::mutable_f64: {
        careful_delete_(cur_entry->data.mutable_f64(), is_scalar);
        break;
      }
      case PtrKind::const_str: {
        const char* const* ptr = cur_entry->data.const_str();

        // get the number of strings referenced by cur_entry->data
        int n_strings = 1;
        for (int i = 0; i < cur_entry->props.ndim; i++) {
          n_strings *= cur_entry->props.shape[i];
        }

        // delete each string within ptr
        for (int i = 0; i < n_strings; i++) {
          delete[] ptr[i];
        }

        // now actually delete ptr
        careful_delete_(ptr, is_scalar);
        break;
      }
      default:
        GR_INTERNAL_UNREACHABLE_ERROR();
    }
  }
}

/// resets the instance to the initial (empty) state.
///
/// the skip_dealloc argument will only be true when the data is transferred
/// to a Registry
static void RegBuilder_reset_to_empty(RegBuilder* ptr, bool skip_dealloc) {
  if ((ptr->capacity > 0) && !skip_dealloc) {
    delete[] ptr->sets;
  }
  (*ptr) = new_RegBuilder();
}

static int RegBuilder_recipe_(RegBuilder* ptr, fetch_Entry_recipe_fn* recipe_fn,
                              int n_entries, EntryProps common_props) {
  if (recipe_fn == nullptr) {
    return GrPrintAndReturnErr("recipe_fn is a nullptr");
  } else if (n_entries <= 0) {
    return GrPrintAndReturnErr("n_entries is not positive");
  } else if (!EntryProps_is_valid(common_props)) {
    return GrPrintAndReturnErr("common_props isn't valid");
  }

  if (ptr->capacity == 0) {
    ptr->capacity = 5;
    ptr->sets = new EntrySet[ptr->capacity];
  } else if (ptr->len == ptr->capacity) {
    // consider making this resizable in the future...
    return GrPrintAndReturnErr("out of capacity");
  }

  ptr->sets[ptr->len++] = EntrySet{n_entries, nullptr, recipe_fn, common_props};
  return GR_SUCCESS;
}

}  // namespace grackle::impl::ratequery

void grackle::impl::ratequery::drop_RegBuilder(RegBuilder* ptr) {
  RegBuilder_reset_to_empty(ptr, false);
}

int grackle::impl::ratequery::RegBuilder_recipe_scalar(
    RegBuilder* ptr, int n_entries, fetch_Entry_recipe_fn* recipe_fn) {
  EntryProps common_props = mk_invalid_EntryProps();
  common_props.ndim = 0;
  return RegBuilder_recipe_(ptr, recipe_fn, n_entries, common_props);
}

int grackle::impl::ratequery::RegBuilder_recipe_1d(
    RegBuilder* ptr, int n_entries, fetch_Entry_recipe_fn* recipe_fn,
    int common_len) {
  EntryProps common_props = mk_invalid_EntryProps();
  common_props.ndim = 1;
  common_props.shape[0] = common_len;
  return RegBuilder_recipe_(ptr, recipe_fn, n_entries, common_props);
}

/// build a new Registry.
///
/// In the process, the current Registry is consumed; it's effectively reset to
/// the state immediately after it was initialized. (This lets us avoid
/// reallocating lots of memory)
grackle::impl::ratequery::Registry
grackle::impl::ratequery::RegBuilder_consume_and_build(RegBuilder* ptr) {
  int n_sets = ptr->len;
  if (n_sets == 0) {
    drop_RegBuilder(ptr);
    return Registry{0, n_sets, nullptr, nullptr};
  } else {
    // set up id_offsets and determine the total number of entries
    int* id_offsets = new int[n_sets];
    int tot_entry_count = 0;
    for (int i = 0; i < n_sets; i++) {
      id_offsets[i] = tot_entry_count;
      tot_entry_count += ptr->sets[i].len;
    }
    Registry out{tot_entry_count, n_sets, id_offsets, ptr->sets};
    // reset to ptr to initial state (but don't deallocate since the ptr was
    // transferred to out
    RegBuilder_reset_to_empty(ptr, true);
    return out;
  }
}

grackle::impl::ratequery::Entry grackle::impl::ratequery::EntrySet_access(
    const EntrySet* entry_set, chemistry_data_storage* my_rates, int i) {
  if (i > entry_set->len) {
    return mk_invalid_Entry();
  } else if (entry_set->embedded_list == nullptr) {  // in recipe-mode
    Entry out = (entry_set->recipe_fn)(my_rates, i);
    out.props = entry_set->common_recipe_props;
    return out;
  } else {  // in embedded-list mode
    return entry_set->embedded_list[i];
  }
}

void grackle::impl::ratequery::drop_EntrySet(EntrySet* ptr) {
  if (ptr->embedded_list == nullptr) {
    // nothing to deallocate in recipe-mode
  } else {
    // in embedded-list-mode, we need to deallocate each Entry in the list
    // and the deallocate the actual list-pointer
    drop_owned_Entry_list_contents(ptr->embedded_list, ptr->len);

    // deallocate the memory associated with embedded_list
    delete[] ptr->embedded_list;
    ptr->embedded_list = nullptr;
  }
}

void grackle::impl::ratequery::drop_Registry(
    grackle::impl::ratequery::Registry* ptr) {
  if (ptr->sets != nullptr) {
    delete[] ptr->id_offsets;
    ptr->id_offsets = nullptr;
    for (int i = 0; i < ptr->n_sets; i++) {
      drop_EntrySet(&ptr->sets[i]);
    }
    delete[] ptr->sets;
    ptr->sets = nullptr;
  }
}

namespace grackle::impl::ratequery {

struct ratequery_rslt_ {
  grunstable_rateid_type rate_id;
  Entry entry;
};

static ratequery_rslt_ invalid_rslt_() {
  return {UNDEFINED_RATE_ID_, mk_invalid_Entry()};
}

static const Registry* get_registry(const chemistry_data_storage* my_rates) {
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr)) {
    return nullptr;
  }
  return my_rates->opaque_storage->registry;
}

/// internal function to search for the rate description i
///
/// @note
/// While it would be nice to enforce that rate data in a
/// `const chemistry_data_storage` can't be mutated, that is something that
/// should be done at the public API level.
static ratequery_rslt_ query_Entry(chemistry_data_storage* my_rates,
                                   long long i) {
  const Registry* registry = get_registry(my_rates);

  if (registry == nullptr || i >= registry->n_entries) {
    return invalid_rslt_();
  }

  for (int set_idx = 0; set_idx < registry->n_sets; set_idx++) {
    EntrySet* cur_set = &registry->sets[set_idx];
    int tmp = i - registry->id_offsets[set_idx];
    if ((tmp >= 0) && (tmp < cur_set->len)) {
      return ratequery_rslt_{i, EntrySet_access(cur_set, my_rates, tmp)};
    }
  }
  return invalid_rslt_();
}

/*
/// this only exists for debugging purposes
static void show_Entry(const Entry* entry) {
  std::printf(
    "{.data=%p, .name=\"%s\", .props={.ndim=%d, .shape={",
    reinterpret_cast<void*>(entry->data), entry->name, entry->props.ndim);
  // we should strongly consider factoring out logic for printing arrays
  for (int i = 0; i < entry->props.ndim; i++) {
    if (i > 0) {
      std::printf(", ");
    }
    std::printf("%d", entry->props.shape[i]);
  }
  std::printf("}}\n");
}
*/

/// compute the number of items in an Entry described by @p props
static long long get_n_items(EntryProps props) {
  GR_INTERNAL_REQUIRE(props.ndim >= 0, "sanity check!");

  if (props.ndim == 0) {
    return 1LL;  // a scalar always consists of 1 item
  }
  long long n_items = 1LL;
  for (int i = 0; i < props.ndim; i++) {
    n_items *= static_cast<long long>(props.shape[i]);
  }
  return n_items;
}

}  // namespace grackle::impl::ratequery

// here we implement the public API
// --------------------------------

extern "C" grunstable_rateid_type grunstable_ratequery_id(
    const chemistry_data_storage* my_rates, const char* name) {
  namespace rate_q = grackle::impl::ratequery;

  const rate_q::Registry* registry = rate_q::get_registry(my_rates);
  if ((name == nullptr) || (registry == nullptr)) {
    return rate_q::UNDEFINED_RATE_ID_;
  }

  for (int set_idx = 0; set_idx < registry->n_sets; set_idx++) {
    rate_q::EntrySet* set = &registry->sets[set_idx];
    int set_len = set->len;
    for (int i = 0; i < set_len; i++) {
      rate_q::Entry entry = rate_q::EntrySet_access(
          set, const_cast<chemistry_data_storage*>(my_rates), i);
      if (std::strcmp(name, entry.name) == 0) {
        return registry->id_offsets[set_idx] + i;
      }
    }
  }
  return rate_q::UNDEFINED_RATE_ID_;
}

extern "C" int grunstable_ratequery_get_f64(chemistry_data_storage* my_rates,
                                            grunstable_rateid_type rate_id,
                                            double* buf) {
  namespace rate_q = grackle::impl::ratequery;
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id).entry;

  if (entry.data.is_null()) {  // in this case, the query failed
    return GR_FAIL;
  }

  const double* src;
  switch (entry.data.tag()) {
    case rate_q::PtrKind::const_f64:
      src = entry.data.const_f64();
      break;
    case rate_q::PtrKind::mutable_f64:
      src = const_cast<const double*>(entry.data.mutable_f64());
      break;
    default:
      return GR_FAIL;
  }

  long long n_items = rate_q::get_n_items(entry.props);
  for (long long i = 0; i < n_items; i++) {
    buf[i] = src[i];
  }
  return GR_SUCCESS;
}

extern "C" int grunstable_ratequery_set_f64(chemistry_data_storage* my_rates,
                                            grunstable_rateid_type rate_id,
                                            const double* buf) {
  namespace rate_q = grackle::impl::ratequery;
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id).entry;
  if (entry.data.is_null()) {  // in this case, the query failed
    return GR_FAIL;
  }

  long long n_items = rate_q::get_n_items(entry.props);
  if (entry.data.tag() == rate_q::PtrKind::mutable_f64) {
    double* dst = entry.data.mutable_f64();
    for (long long i = 0; i < n_items; i++) {
      dst[i] = buf[i];
    }
    return GR_SUCCESS;
  } else {  // the retrieved pointer is either immutable or the wrong type
    return GR_FAIL;
  }
}

extern "C" int grunstable_ratequery_get_str(chemistry_data_storage* my_rates,
                                            grunstable_rateid_type rate_id,
                                            char* const* buf) {
  namespace rate_q = grackle::impl::ratequery;
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id).entry;

  if (entry.data.is_null()) {  // in this case, the query failed
    return GR_FAIL;
  }

  if (entry.data.tag() != rate_q::PtrKind::const_str) {
    return GR_FAIL;
  }
  const char* const* src = entry.data.const_str();
  long long n_items = rate_q::get_n_items(entry.props);
  for (long long i = 0; i < n_items; i++) {
    std::strcpy(buf[i], src[i]);
  }
  return GR_SUCCESS;
}

extern "C" int grunstable_ratequery_prop(
    const chemistry_data_storage* my_rates, grunstable_rateid_type rate_id,
    enum grunstable_ratequery_prop_kind prop_kind, long long* ptr) {
  namespace rate_q = grackle::impl::ratequery;

  // short-term hack! (it's bad practice to "cast away the const")
  rate_q::Entry entry =
      rate_q::query_Entry(const_cast<chemistry_data_storage*>(my_rates),
                          rate_id)
          .entry;

  const rate_q::EntryProps& props = entry.props;
  if ((entry.name == nullptr) || !rate_q::EntryProps_is_valid(props)) {
    return GR_FAIL;
  }

  switch (prop_kind) {
    case GRUNSTABLE_QPROP_NDIM: {
      *ptr = static_cast<long long>(props.ndim);
      return GR_SUCCESS;
    }
    case GRUNSTABLE_QPROP_SHAPE: {
      for (int i = 0; i < props.ndim; i++) {
        ptr[i] = static_cast<long long>(props.shape[i]);
      }
      return GR_SUCCESS;
    }
    case GRUNSTABLE_QPROP_MAXITEMSIZE: {
      switch (entry.data.tag()) {
        case rate_q::PtrKind::const_f64:
        case rate_q::PtrKind::mutable_f64: {
          *ptr = static_cast<long long>(sizeof(double));
          return GR_SUCCESS;
        }
        case rate_q::PtrKind::const_str: {
          long long n_items = rate_q::get_n_items(entry.props);
          const char* const* str_list = entry.data.const_str();

          std::size_t max_size = 1;
          for (long long i = 0; i < n_items; i++) {
            // max_size holds the max number of bytes per element
            max_size = std::max(max_size, std::strlen(str_list[i]) + 1);
          }
          *ptr = static_cast<long long>(max_size);
          return GR_SUCCESS;
        }
      }
      // if we reach here then we didn't handle a possible PtrKind
      GR_INTERNAL_UNREACHABLE_ERROR();
    }
    case GRUNSTABLE_QPROP_WRITABLE: {
      *ptr = static_cast<long long>(entry.data.is_const_ptr());
      return GR_SUCCESS;
    }
    case GRUNSTABLE_QPROP_DTYPE: {
      switch (entry.data.tag()) {
        case rate_q::PtrKind::const_f64:
        case rate_q::PtrKind::mutable_f64: {
          *ptr = static_cast<long long>(GRUNSTABLE_TYPE_F64);
          return GR_SUCCESS;
        }
        case rate_q::PtrKind::const_str: {
          *ptr = static_cast<long long>(GRUNSTABLE_TYPE_STR);
          return GR_SUCCESS;
        }
      }
      // if we reach here then we didn't handle a possible PtrKind
      GR_INTERNAL_UNREACHABLE_ERROR();
    }
  }
  return GrPrintAndReturnErr("received an unknown prop_kind");
}

extern "C" const char* grunstable_ith_rate(
    const chemistry_data_storage* my_rates, unsigned long long i,
    grunstable_rateid_type* out_rate_id) {
  namespace rate_q = grackle::impl::ratequery;

  const long long sanitized_i =
      (i < std::numeric_limits<long long>::max()) ? (long long)i : -1;
  // short-term hack! (it's bad practice to "cast away the const")
  rate_q::ratequery_rslt_ tmp = rate_q::query_Entry(
      const_cast<chemistry_data_storage*>(my_rates), sanitized_i);
  if (out_rate_id != nullptr) {
    *out_rate_id = tmp.rate_id;
  }
  return tmp.entry.name;
}

extern "C" unsigned long long grunstable_ratequery_nrates(
    const chemistry_data_storage* my_rates) {
  namespace rate_q = grackle::impl::ratequery;
  const rate_q::Registry* registry = my_rates->opaque_storage->registry;
  return registry->n_entries;
}
