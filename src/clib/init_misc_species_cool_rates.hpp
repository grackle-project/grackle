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

#include "cool/TabMDSpCooler.hpp"
#include "grackle.h"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<TabMDSpCooler, Error> init_misc_species_cool_rates(
    const chemistry_data& my_chemistry, const code_units& my_units);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* INIT_MISC_SPECIES_COOL_RATES_HPP */
