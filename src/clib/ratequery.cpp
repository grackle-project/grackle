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

#include <cstring>  // strcmp
#include <climits>  // LLONG_MAX
#include "grackle.h"
#include "LUT.hpp"             // PhotoRxnLUT
#include "opaque_storage.hpp"  // gr_opaque_storage
#include "ratequery.hpp"
#include "status_reporting.h"

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
enum { UNDEFINED_RATE_ID_ = 0 };

// Create machinery to lookup Other Miscellaneous Rates
// ----------------------------------------------------
// -> ideally this should get relocated to a more sensible location (but not
//    sure where that is...)

static Entry get_PhotoRxn_Entry(chemistry_data_storage* my_rates, int i) {
  if (my_rates == nullptr) {  // <- shouldn't come up
    return mk_invalid_Entry();
  }
  switch (i) {
    case PhotoRxnLUT::k24:
      return new_Entry(&my_rates->k24, "k24");
    case PhotoRxnLUT::k25:
      return new_Entry(&my_rates->k25, "k25");
    case PhotoRxnLUT::k26:
      return new_Entry(&my_rates->k26, "k26");
    case PhotoRxnLUT::k27:
      return new_Entry(&my_rates->k27, "k27");
    case PhotoRxnLUT::k28:
      return new_Entry(&my_rates->k28, "k28");
    case PhotoRxnLUT::k29:
      return new_Entry(&my_rates->k29, "k29");
    case PhotoRxnLUT::k30:
      return new_Entry(&my_rates->k30, "k30");
    case PhotoRxnLUT::k31:
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
    return RegBuilder_recipe_scalar(ptr, PhotoRxnLUT::NUM_ENTRIES,
                                    &get_PhotoRxn_Entry);
  }
  return GR_SUCCESS;
}

namespace grackle::impl::ratequery {

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

  ptr->sets[ptr->len++] = EntrySet{n_entries, recipe_fn, common_props};
  return GR_SUCCESS;
}

}  // namespace grackle::impl::ratequery

void grackle::impl::ratequery::drop_RegBuilder(RegBuilder* ptr) {
  RegBuilder_reset_to_empty(ptr, false);
}

namespace grackle::impl::ratequery {}  // namespace grackle::impl::ratequery

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

void grackle::impl::ratequery::drop_Registry(
    grackle::impl::ratequery::Registry* ptr) {
  if (ptr->sets != nullptr) {
    delete[] ptr->id_offsets;
    ptr->id_offsets = nullptr;
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
static ratequery_rslt_ query_Entry(chemistry_data_storage* my_rates,
                                   long long i) {
  const Registry* registry = get_registry(my_rates);

  if (registry == nullptr || i >= registry->n_entries) {
    return invalid_rslt_();
  }

  for (int set_idx = 0; set_idx < registry->n_sets; set_idx++) {
    const struct EntrySet cur_set = registry->sets[set_idx];
    int cur_id_offset = registry->id_offsets[set_idx];

    const long long tmp = i - static_cast<long long>(cur_id_offset);
    if ((tmp >= 0) && (tmp < cur_set.len)) {
      ratequery_rslt_ out;
      out.rate_id = i;
      out.entry = cur_set.recipe_fn(my_rates, tmp);
      out.entry.props = cur_set.common_props;
      return out;
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
    const rate_q::EntrySet set = registry->sets[set_idx];
    int set_len = set.len;
    for (int i = 0; i < set_len; i++) {
      // short-term hack! (it's bad practice to "cast away the const")
      rate_q::Entry entry =
          set.recipe_fn(const_cast<chemistry_data_storage*>(my_rates), i);
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

  if (entry.data == nullptr) {  // in this case, the query failed
    return GR_FAIL;
  }

  long long n_items = rate_q::get_n_items(entry.props);
  for (long long i = 0; i < n_items; i++) {
    buf[i] = entry.data[i];
  }
  return GR_SUCCESS;
}

extern "C" int grunstable_ratequery_set_f64(chemistry_data_storage* my_rates,
                                            grunstable_rateid_type rate_id,
                                            const double* buf) {
  namespace rate_q = grackle::impl::ratequery;
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id).entry;
  if (entry.data == nullptr) {  // in this case, the query failed
    return GR_FAIL;
  }

  long long n_items = rate_q::get_n_items(entry.props);
  for (long long i = 0; i < n_items; i++) {
    entry.data[i] = buf[i];
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
      *ptr = static_cast<long long>(sizeof(double));
      return GR_SUCCESS;
    }
    default:
      GR_INTERNAL_UNREACHABLE_ERROR();
  }
}

extern "C" const char* grunstable_ith_rate(
    const chemistry_data_storage* my_rates, unsigned long long i,
    grunstable_rateid_type* out_rate_id) {
  namespace rate_q = grackle::impl::ratequery;

  const long long sanitized_i = (i < LLONG_MAX) ? (long long)i : -1;
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
