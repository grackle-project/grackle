//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the primordial rates
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_RATES_HPP
#define INITIALIZE_RATES_HPP

#include "grackle.h"

namespace grackle::impl {

int initialize_rates(chemistry_data* my_chemistry,
                     chemistry_data_storage* my_rates, code_units* my_units,
                     double co_length_unit, double co_density_unit);

}  // namespace grackle::impl

#endif /* INITIALIZE_RATES_HPP */
