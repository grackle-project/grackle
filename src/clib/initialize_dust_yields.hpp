//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the dust yields
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_DUST_YIELDS_HPP
#define INITIALIZE_DUST_YIELDS_HPP

#include "grackle.h"

namespace grackle::impl {

int initialize_dust_yields(chemistry_data* my_chemistry,
                           chemistry_data_storage* my_rates,
                           code_units* my_units);

}  // namespace grackle::impl

#endif /* INITIALIZE_DUST_YIELDS_HPP */
