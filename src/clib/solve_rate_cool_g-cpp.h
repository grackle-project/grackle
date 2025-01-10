// See LICENSE file for license and copyright information

/// @file solve_rate_cool_g-cpp.h
/// @brief Declares signature of solve_rate_cool_g

// This file was initially generated automatically during conversion of the
// solve_rate_cool_g function from FORTRAN to C++

#ifndef MY_FILE_CPP_H
#define MY_FILE_CPP_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Solve the multi-species rate and cooling equations
///
/// @par History
/// written by: Yu Zhang, Peter Anninos and Tom Abel
/// modified1:  January, 1996 by Greg Bryan; converted to KRONOS
/// modified2:  October, 1996 by GB; adapted to AMR
/// modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
/// modified4:  June,    2005 by GB to solve rate & cool at same time
/// modified5:  April,   2009 by JHW to include radiative transfer
/// modified6:  September, 2009 by BDS to include cloudy cooling
/// modified7:  January, 2025 by Matthew Abruzzo; ported to C++
///
/// @return Returns GR_SUCCESS or GR_FAIL to indicate whether there was an error
int solve_rate_cool_g(
  int imetal, double dt, double utem, double uxyz, double urho,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  code_units* my_units, grackle_field_data* my_fields,
  photo_rate_storage* my_uvb_rates
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* MY_FILE_CPP_H */
