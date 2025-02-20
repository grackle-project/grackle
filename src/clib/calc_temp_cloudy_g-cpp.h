// See LICENSE file for license and copyright information

/// @file calc_temp_cloudy_g-cpp.h
/// @brief Declares signature of calc_temp_cloudy_g

// This file was initially generated automatically during conversion of the
// calc_temp_cloudy_g function from FORTRAN to C++

#ifndef MY_FILE_CPP_H
#define MY_FILE_CPP_H

#include "grackle.h"             // gr_float
#include "internal_units.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Calculate temperature on a 3D grid using mmw from a cloudy table
///
/// @par History
/// written by: Britton Smith May 2015
///
/// @param[out] temperature_data Array where computed values are written
/// @param[in]  flag if metal field is active (0 = no, 1 = yes)
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  my_rates specifies the tables
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
///
/// @note
/// We are specifying way more information than necessary!
void calc_temp_cloudy_g(
  gr_float* temperature_data_, int* imetal, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields,
  InternalGrUnits internalu
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* MY_FILE_CPP_H */
