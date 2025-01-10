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

void solve_rate_cool_g(
  int* imetal, double* dt, double* utem, double* uxyz, double* urho, int* ierr,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  code_units* my_units, grackle_field_data* my_fields,
  photo_rate_storage* my_uvb_rates
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* MY_FILE_CPP_H */
