//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares a function to initialize miscellaneous species cooling rates
///
//===----------------------------------------------------------------------===//

#ifndef INIT_MISC_SPECIES_COOL_RATES_HPP
#define INIT_MISC_SPECIES_COOL_RATES_HPP

#include "grackle.h"

namespace grackle::impl {

int init_misc_species_cool_rates(chemistry_data* my_chemistry,
                                 chemistry_data_storage* my_rates,
                                 code_units* my_units);

int free_misc_species_cool_rates(chemistry_data* my_chemistry,
                                 chemistry_data_storage* my_rates);
}  // namespace grackle::impl

#endif /* INIT_MISC_SPECIES_COOL_RATES_HPP */
