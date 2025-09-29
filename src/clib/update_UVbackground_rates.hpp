//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Forward declare the update_UVbackground_rates function
///
//===----------------------------------------------------------------------===//

#ifndef UPDATE_UVBACKGROUND_RATES_HPP
#define UPDATE_UVBACKGROUND_RATES_HPP

#include "grackle.h"

namespace grackle::impl {

/// Update UV background rates to current redshift
int update_UVbackground_rates(chemistry_data* my_chemistry,
                              chemistry_data_storage* my_rates,
                              photo_rate_storage* my_uvb_rates,
                              code_units* my_units);

}  // namespace grackle::impl

#endif  // UPDATE_UVBACKGROUND_RATES_HPP
