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

namespace grackle::impl {

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

