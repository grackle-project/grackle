//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_temp_cloudy_g function
///
//===----------------------------------------------------------------------===//

/// @file calc_temp_cloudy_g-cpp.h
/// @brief Declares signature of calc_temp_cloudy_g

// This file was initially generated automatically during conversion of the
// calc_temp_cloudy_g function from FORTRAN to C++

#ifndef CALC_TEMP_CLOUDY_G_H
#define CALC_TEMP_CLOUDY_G_H

#include "grackle.h"  // gr_float
#include "internal_units.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Calculate temperature on a 3D grid using mmw from a cloudy table
///
/// @par History
/// written by: Britton Smith May 2015
/// modified1:  February, 2025 by Matthew Abruzzo; ported to C++
///
/// @param[out] temperature_data Array where computed values are written
/// @param[in]  imetal flag if metal field is active (0 = no, 1 = yes)
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  cloudy_primordia specifies the cloudy table
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
///
/// @note
/// We are specifying way more information than necessary!
void calc_temp_cloudy_g(gr_float* temperature_data_, int imetal,
                        chemistry_data* my_chemistry,
                        cloudy_data cloudy_primordial,
                        grackle_field_data* my_fields,
                        InternalGrUnits internalu);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* CALC_TEMP_CLOUDY_G_H */
