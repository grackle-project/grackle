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
#include "internal_types.hpp"  // CollisionalRxnRateCollection
#include "LUT.hpp"             // CollisionalRxnLUT
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

// introduce some basic machinery to help us implement dynamic lookup of rates

static Entry new_Entry_standard_kcol_(chemistry_data_storage* my_rates,
                                      const char* name, int index) {
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr) ||
      (my_rates->opaque_storage->kcol_rate_tables == nullptr)) {
    return new_Entry(nullptr, name);
  } else {
    return new_Entry(my_rates->opaque_storage->kcol_rate_tables->data[index],
                     name);
  }
}

#define MKENTRY_(PTR, NAME, FAMILY)                                            \
  new_Entry(((PTR) == nullptr) ? nullptr : (PTR)->NAME, #NAME, FAMILY)
#define MKENTRY_SCALAR_(PTR, NAME)                                             \
  new_Entry(((PTR) == nullptr) ? nullptr : &((PTR)->NAME), #NAME)
#define MKENTRY_STANDARD_KCOL_(PTR, NAME, INDEX)                               \
  new_Entry_standard_kcol_(PTR, #NAME, INDEX)

// Create machinery to lookup Standard-Form Collisional Reaction Rates
// -------------------------------------------------------------------
// see the next section for other macros

// this fn leverages the following properties of the CollisionalRxnLUT enum:
// - each entry of collisional_rxn_rate_members.def has a corresponding
//   enumeration-constant
// - the very first enumeration constant has a value of 0 (since a value wasn't
//   explicitly specified)
// - the value of each other enumeration constants is 1 larger than the value
//   of the previous value (if a value isn't explicitly specified)
// - CollisionalRxnLUT::NUM_ENTRIES specifies the number of other enumeration
//   constants (excluding CollisionalRxnLUT::NUM_ENTRIES) in the enum
static Entry get_CollisionalRxn_Entry(chemistry_data_storage* my_rates, int i) {
  switch (i) {
#define ENTRY(NAME)                                                            \
  case CollisionalRxnLUT::NAME: {                                              \
    return MKENTRY_STANDARD_KCOL_(my_rates, NAME, CollisionalRxnLUT::NAME);    \
  }
#include "collisional_rxn_rate_members.def"
#undef ENTRY
    default: {
      return mk_invalid_Entry();
    }
  }
}

static Entry get_k13dd_Entry(chemistry_data_storage* my_rates, int i) {
  if (i == 0) {
    double* ptr = (my_rates == nullptr) ? nullptr : my_rates->k13dd;
    return new_Entry(ptr, "k13dd");
  } else {
    return mk_invalid_Entry();
  }
}

// Create machinery to lookup Other Miscellaneous Rates
// ----------------------------------------------------

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
  switch (i) {
    // Radiative rates for 6-species (for external field):
    case MiscRxn_k24:
      return MKENTRY_SCALAR_(my_rates, k24);
    case MiscRxn_k25:
      return MKENTRY_SCALAR_(my_rates, k25);
    case MiscRxn_k26:
      return MKENTRY_SCALAR_(my_rates, k26);
    // Radiative rates for 9-species
    case MiscRxn_k27:
      return MKENTRY_SCALAR_(my_rates, k27);
    case MiscRxn_k28:
      return MKENTRY_SCALAR_(my_rates, k28);
    case MiscRxn_k29:
      return MKENTRY_SCALAR_(my_rates, k29);
    case MiscRxn_k30:
      return MKENTRY_SCALAR_(my_rates, k30);
    case MiscRxn_k31:
      return MKENTRY_SCALAR_(my_rates, k31);
    default: {
      return mk_invalid_Entry();
    }
  }
}

}  // namespace grackle::impl::ratequery

grackle::impl::ratequery::Registry grackle::impl::ratequery::new_Registry(
    const chemistry_data& my_chemistry) {
  if (my_chemistry.primordial_chemistry == 0) {
    return Registry{0, nullptr};
  }
  EntryProps props_LogTLinInterp = mk_invalid_EntryProps();
  props_LogTLinInterp.ndim = 1;
  props_LogTLinInterp.shape[0] = my_chemistry.NumberOfTemperatureBins;

  // maybe k13dd should be considered multi-dimensional?
  EntryProps props_k13dd = mk_invalid_EntryProps();
  props_k13dd.ndim = 1;
  props_k13dd.shape[0] = my_chemistry.NumberOfTemperatureBins * 14;

  EntryProps props_scalar = mk_invalid_EntryProps();
  props_scalar.ndim = 0;

  const RecipeEntrySet standard_sets[] = {
      {1000, CollisionalRxnLUT::NUM_ENTRIES, &get_CollisionalRxn_Entry,
       props_LogTLinInterp},
      {1999, 1, &get_k13dd_Entry, props_k13dd},
      {2000, MiscRxn_NRATES, &get_MiscRxn_Entry, props_scalar}};

  int len = static_cast<int>(sizeof(standard_sets) / sizeof(RecipeEntrySet));
  RecipeEntrySet* sets = new RecipeEntrySet[len];
  for (int i = 0; i < len; i++) {
    sets[i] = standard_sets[i];
  }

  return Registry{len, sets};
}

void grackle::impl::ratequery::drop_Registry(
    grackle::impl::ratequery::Registry* ptr) {
  if (ptr->sets != nullptr) {
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
/// We interpret i as rate_id, when use_rate_id is true. Otherwise, we just
/// look for the ith rate description (we introduce an artificial distinction
/// between the 2 cases because we want to reserve the right to be able to
/// change the relationship if it becomes convenient in the future)
static ratequery_rslt_ query_Entry(chemistry_data_storage* my_rates,
                                   long long i, bool use_rate_id) {
  const Registry* registry = get_registry(my_rates);

  if (registry == nullptr) {
    return invalid_rslt_();
  }

  int total_len = 0;  // <- we increment this as we go through the rates

  for (int set_idx = 0; set_idx < registry->len; set_idx++) {
    const struct RecipeEntrySet cur_set = registry->sets[set_idx];

    const long long tmp = (use_rate_id) ? i - cur_set.id_offset : i - total_len;
    if ((tmp >= 0) && (tmp < cur_set.len)) {
      ratequery_rslt_ out;
      out.rate_id = tmp + cur_set.id_offset;
      out.entry = cur_set.fn(my_rates, tmp);
      out.entry.props = cur_set.common_props;
      return out;
    }
    total_len += cur_set.len;
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

  for (int set_idx = 0; set_idx < registry->len; set_idx++) {
    const rate_q::RecipeEntrySet set = registry->sets[set_idx];
    for (int i = 0; i < set.len; i++) {
      rate_q::Entry entry = set.fn(nullptr, i);
      if (std::strcmp(name, entry.name) == 0) {
        return set.id_offset + i;
      }
    }
  }
  return rate_q::UNDEFINED_RATE_ID_;
}

extern "C" int grunstable_ratequery_get_f64(chemistry_data_storage* my_rates,
                                            grunstable_rateid_type rate_id,
                                            double* buf) {
  namespace rate_q = grackle::impl::ratequery;
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id, true).entry;

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
  rate_q::Entry entry = rate_q::query_Entry(my_rates, rate_id, true).entry;
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
                          rate_id, true)
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
      const_cast<chemistry_data_storage*>(my_rates), sanitized_i, false);
  if (out_rate_id != nullptr) {
    *out_rate_id = tmp.rate_id;
  }
  return tmp.entry.name;
}

extern "C" unsigned long long grunstable_ratequery_nrates(
    const chemistry_data_storage* my_rates) {
  namespace rate_q = grackle::impl::ratequery;
  const rate_q::Registry* registry = my_rates->opaque_storage->registry;

  unsigned long long out = 0;
  for (int i = 0; i < registry->len; i++) {
    out += static_cast<unsigned long long>(registry->sets[i].len);
  }
  return out;
}
