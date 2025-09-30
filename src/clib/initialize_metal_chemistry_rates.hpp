//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the metal chemistry rates
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_METAL_CHEMISTRY_RATES_HPP
#define INITIALIZE_METAL_CHEMISTRY_RATES_HPP

#include "grackle.h"

namespace grackle::impl {

int initialize_metal_chemistry_rates(chemistry_data* my_chemistry,
                                     chemistry_data_storage* my_rates,
                                     code_units* my_units);

int free_metal_chemistry_rates(chemistry_data* my_chemistry,
                               chemistry_data_storage* my_rates);
}  // namespace grackle::impl

#endif /* INITIALIZE_METAL_CHEMISTRY_RATES_HPP */
