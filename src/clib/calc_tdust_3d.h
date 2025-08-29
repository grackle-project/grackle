// See LICENSE file for license and copyright information

/// @file calc_tdust_3d.h
/// @brief Declares signature of calc_tdust_3d_g

// This file was initially generated automatically during conversion of the
// calc_tdust_3d_g function from FORTRAN to C++

#ifndef MY_FILE_CPP_H
#define MY_FILE_CPP_H

#include "grackle.h"             // gr_float
#include "internal_units.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Calculate dust heat balance to get the dust temperature.
///
/// @par History
/// written by: Britton Smith July 2011
void calc_tdust_3d_g(
  gr_float* gas_temp_data_, gr_float* dust_temp_data_, int imetal,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, InternalGrUnits internalu
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* MY_FILE_CPP_H */
