// See LICENSE file for license and copyright information

/// @file cool_multi_time_g-cpp.h
/// @brief Declares signature of cool_multi_time_g

// This file was initially generated automatically during conversion of the
// cool_multi_time_g function from FORTRAN to C++

#ifndef MY_FILE_CPP_H
#define MY_FILE_CPP_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Solve the energy cooling equations
///
/// @par history
/// * written by: Yu Zhang, Peter Anninos and Tom Abel
/// * modified1: January, 1996 by Greg Bryan; adapted to KRONOS
/// * modified2: October, 1996 by GB; moved to AMR
/// * modified3: February, 2003 by Robert Harkness; iteration mask
/// * modified4: December, 2024 by Matthew Abruzzo; ported to C++
void cool_multi_time_g(
  gr_float* cooltime_data_, int* imetal, double* utem, double* uxyz,
  double* urho, chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  code_units* my_units, grackle_field_data* my_fields,
  photo_rate_storage* my_uvb_rates
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* MY_FILE_CPP_H */
