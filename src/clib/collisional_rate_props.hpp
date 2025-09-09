//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines the function to initialize extra collisional rates
///
//===----------------------------------------------------------------------===//

#ifndef COLLISIONAL_RATE_PROPS_HPP
#define COLLISIONAL_RATE_PROPS_HPP

#include "grackle.h"

// this function pointer type is defined inside an extern "C" block because all
// described functions have C linkage
extern "C" {

/// the type of a generic temperature dependent rate function
typedef double (*rate_function)(double, double, chemistry_data*);

}  // extern "C"

namespace grackle::impl {


/// holds properties about a collisional rate.
///
/// In the future, it may be useful to track more information (e.g.
/// reactants, products, the reaction name)
struct KColProp {
  /// holds the LUT (lookup table) index associated with the rate.
  ///
  /// (We could potentially refactor to avoid tracking this)
  int kcol_lut_index;

  /// pointer to the function used for computing the rate
  rate_function fn_ptr;

  /// specifies whether this is a 2-body or 3-body rate
  bool is_2body;
};

/// describes the signature of the callback function that is used to digest each
/// instances of `KColProp`.
///
/// As is standard for C-like callback functions, a callback function accepts
/// a function-specific `void*` context that is used argument to pass extra
/// information to the callback and track information between calls to the
/// callback
typedef void visit_kcol_prop_callback(struct KColProp, void*);

/// calls visit_kcol_prop_callback for each rate relevant to the Chemical
/// Network configuration specifier by `my_chemistry`
///
/// Implementation Notes
/// ====================
/// If we are willing to more fully embrace C++, it would be a lot more
/// idiomatic to eliminate the `context` argument and make `fn` a template
/// argument
int visit_rate_props(const chemistry_data* my_chemistry,
                     visit_kcol_prop_callback* cb,
                     void* ctx);

/// initialize misc primordial_chemistry == 4 and metal chemistry rates
///
/// @todo
/// we should refactor this logic so it is implemented that all "standard"
/// collisional rate initialization logic is implemented in a consistent way
int init_extra_collisional_rates(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units);

} // namespace grackle::impl

#endif /* COLLISIONAL_RATE_PROPS_HPP */

