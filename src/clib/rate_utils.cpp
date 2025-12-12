//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines functions that perform basic utilities related to rate data.
///
//===----------------------------------------------------------------------===//

#include <cstring>  // strcmp
#include <climits>  // LLONG_MAX
#include "grackle.h"
#include "internal_types.hpp"  // CollisionalRxnRateCollection
#include "LUT.hpp"             // CollisionalRxnLUT
#include "opaque_storage.hpp"  // gr_opaque_storage

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

struct rateprop_ {
  double* data;
  const char* name;
};

static inline rateprop_ mk_rateprop_(double* rate, const char* name) {
  rateprop_ out;
  out.data = rate;
  out.name = name;
  return out;
}

static inline rateprop_ mk_rateprop_standard_kcol_(
    chemistry_data_storage* my_rates, const char* name, int index) {
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr) ||
      (my_rates->opaque_storage->kcol_rate_tables == nullptr)) {
    return mk_rateprop_(nullptr, name);
  } else {
    return mk_rateprop_(my_rates->opaque_storage->kcol_rate_tables->data[index],
                        name);
  }
}

#define MKPROP_(PTR, NAME)                                                     \
  mk_rateprop_(((PTR) == nullptr) ? nullptr : (PTR)->NAME, #NAME)
#define MKPROP_SCALAR_(PTR, NAME)                                              \
  mk_rateprop_(((PTR) == nullptr) ? nullptr : &((PTR)->NAME), #NAME)
#define MKPROP_STANDARD_KCOL_(PTR, NAME, INDEX)                                \
  mk_rateprop_standard_kcol_(PTR, #NAME, INDEX)

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
static rateprop_ get_CollisionalRxn_rateprop_(chemistry_data_storage* my_rates,
                                              int i) {
  switch (i) {
#define ENTRY(NAME)                                                            \
  case CollisionalRxnLUT::NAME: {                                              \
    return MKPROP_STANDARD_KCOL_(my_rates, NAME, CollisionalRxnLUT::NAME);     \
  }
#include "collisional_rxn_rate_members.def"
#undef ENTRY
    default: {
      rateprop_ out = {nullptr, nullptr};
      return out;
    }
  }
}

// Create machinery to lookup Other Miscellaneous Rates
// ----------------------------------------------------

enum MiscRxnRateKind_ {
  MiscRxn_k13dd,
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

static rateprop_ get_MiscRxn_rateprop_(chemistry_data_storage* my_rates,
                                       int i) {
  switch (i) {
    // density dependent version of k13 (which is a CollisionalRxn)
    case MiscRxn_k13dd:
      return MKPROP_(my_rates, k13dd);
    // Radiative rates for 6-species (for external field):
    case MiscRxn_k24:
      return MKPROP_SCALAR_(my_rates, k24);
    case MiscRxn_k25:
      return MKPROP_SCALAR_(my_rates, k25);
    case MiscRxn_k26:
      return MKPROP_SCALAR_(my_rates, k26);
    // Radiative rates for 9-species
    case MiscRxn_k27:
      return MKPROP_SCALAR_(my_rates, k27);
    case MiscRxn_k28:
      return MKPROP_SCALAR_(my_rates, k28);
    case MiscRxn_k29:
      return MKPROP_SCALAR_(my_rates, k29);
    case MiscRxn_k30:
      return MKPROP_SCALAR_(my_rates, k30);
    case MiscRxn_k31:
      return MKPROP_SCALAR_(my_rates, k31);
    default: {
      rateprop_ out = {nullptr, nullptr};
      return out;
    }
  }
}

// define some additional generic machinery
// ----------------------------------------

#define RATE_SET_COUNT 2
typedef rateprop_ fetch_rateprop_fn(chemistry_data_storage*, int);
struct rateprop_set_ {
  int id_offset;
  int len;
  fetch_rateprop_fn* fn;
};
struct rate_registry_type_ {
  int len;
  struct rateprop_set_ sets[RATE_SET_COUNT];
};

static const struct rate_registry_type_ rate_registry_ = {
    /* len: */ RATE_SET_COUNT,
    /* sets: */ {
        {1000, CollisionalRxnLUT::NUM_ENTRIES, &get_CollisionalRxn_rateprop_},
        {2000, MiscRxn_NRATES, &get_MiscRxn_rateprop_}}};

struct ratequery_rslt_ {
  grunstable_rateid_type rate_id;
  rateprop_ prop;
};

/// internal function to search for the rate_property i
///
/// We interpret i as rate_id, when use_rate_id is true. Otherwise, we just
/// look for the ith rateprop (we introduce an artificial distinction between
/// the 2 cases because we want to reserve the right to be able to change the
/// relationship if it becomes convenient in the future)
static struct ratequery_rslt_ query_rateprop_(chemistry_data_storage* my_rates,
                                              long long i, bool use_rate_id) {
  int total_len = 0;  // <- we increment this as we go through the rates

  for (int set_idx = 0; set_idx < rate_registry_.len; set_idx++) {
    const struct rateprop_set_ cur_set = rate_registry_.sets[set_idx];

    const long long tmp = (use_rate_id) ? i - cur_set.id_offset : i - total_len;
    if ((tmp >= 0) && (tmp < cur_set.len)) {
      struct ratequery_rslt_ out;
      out.rate_id = tmp + cur_set.id_offset;
      out.prop = cur_set.fn(my_rates, tmp);
      return out;
    }
    total_len += cur_set.len;
  }
  struct ratequery_rslt_ out = {UNDEFINED_RATE_ID_, {nullptr, nullptr}};
  return out;
}

}  // namespace grackle::impl::ratequery

// here we implement the public API
// --------------------------------

extern "C" grunstable_rateid_type grunstable_ratequery_id(const char* name) {
  namespace rate_q = grackle::impl::ratequery;

  if (name == nullptr) {
    return rate_q::UNDEFINED_RATE_ID_;
  }

  for (int set_idx = 0; set_idx < rate_q::rate_registry_.len; set_idx++) {
    const rate_q::rateprop_set_ cur_set = rate_q::rate_registry_.sets[set_idx];
    for (int i = 0; i < cur_set.len; i++) {
      rate_q::rateprop_ prop = cur_set.fn(nullptr, i);
      if (std::strcmp(name, prop.name) == 0) {
        return cur_set.id_offset + i;
      }
    }
  }
  return rate_q::UNDEFINED_RATE_ID_;
}

extern "C" double* grunstable_ratequery_get_ptr(
    chemistry_data_storage* my_rates, grunstable_rateid_type rate_id) {
  namespace rate_q = grackle::impl::ratequery;

  return rate_q::query_rateprop_(my_rates, rate_id, true).prop.data;
}

extern "C" const char* grunstable_ith_rate(
    unsigned long long i, grunstable_rateid_type* out_rate_id) {
  namespace rate_q = grackle::impl::ratequery;

  const long long sanitized_i = (i < LLONG_MAX) ? (long long)i : -1;
  rate_q::ratequery_rslt_ tmp =
      rate_q::query_rateprop_(nullptr, sanitized_i, false);
  if (out_rate_id != nullptr) {
    *out_rate_id = tmp.rate_id;
  }
  return tmp.prop.name;
}
