//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines self_shielding_err_check
///
//===----------------------------------------------------------------------===//

#ifndef SELF_SHIELDING_ERR_CHECK_HPP
#define SELF_SHIELDING_ERR_CHECK_HPP

#include <cstdio> // std::fprintf, stderr
#include "grackle.h"

/// Perform an error check related to self-shielding. If the check fails, the
/// function returns FAIL and prints an error message to stderr
///
/// @param my_chemistry Holds configuration of chemistry solver
/// @param fields Specify the field names
/// @param func_name Name of the function that is calling the error check
inline int self_shielding_err_check(const chemistry_data *my_chemistry,
                                    const grackle_field_data *fields,
                                    const char* func_name) {
  if (my_chemistry->H2_self_shielding == 1) {
    if (fields->grid_rank != 3) {
      std::fprintf(stderr, "Error in %s: H2 self-shielding option 1 "
                           "will only work for 3D Cartesian grids. Use option 2 "
                           "to provide an array of shielding lengths with "
                           "H2_self_shielding_length or option 3 to use the "
                           "local Jeans length.",
                   func_name);
      return GR_FAIL;
    } else if (my_chemistry->primordial_chemistry >= 2 &&
               fields->grid_dx <= 0) {
      std::fprintf(stderr, "Error in %s: H2 self-shielding option 1 and primordial "
                           "chemistry options of 2 or more require that grid_dx "
                           "has a positive value.",
                   func_name);
      return GR_FAIL;
    }
  }
  return GR_SUCCESS;
}

#endif  // SELF_SHIELDING_ERR_CHECK_HPP
