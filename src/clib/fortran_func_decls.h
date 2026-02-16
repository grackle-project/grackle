/// @file fortran_func_decls.h
/// @brief Header file defining forward declaring function signatures of
///        Fortran subroutines have not been ported to C++ yet.

#ifndef FORTRAN_FN_DECLARATIONS_HPP
#define FORTRAN_FN_DECLARATIONS_HPP

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "grackle_macros.h" // FORTRAN_NAME
#include "grackle.h"        // gr_float
#include "phys_constants.h" // physical constants
#include <stdint.h>         // int32_t
typedef int32_t gr_mask_type;
typedef long long gr_i64;

#define MASK_TRUE 1
#define MASK_FALSE 0

void FORTRAN_NAME(calc_grain_size_increment_species_1d)(
  const int* igrgr, const gr_mask_type* itmask, const int* SN0_N, int* in, int* jn, int* kn,
  int* is, int* ie, int* j, int* k, double* dom, gr_float* d_data_ptr, int* nSN,
  const gr_float* dsp_data_ptr, gr_float* SN_metal_data_ptr, double* SN_fsp,
  double* SN_r0sp_data_ptr, double* ssp, double* sgsp, double* alsp_data_ptr,
  const int* gr_N, const int* gr_Size, const double* gr_dT, const double* gr_Td,
  double* SN_kp0sp_data_ptr
);

void FORTRAN_NAME(solve_cubic_equation)(
  double* a, double* b, double* c, double* root
);

void FORTRAN_NAME(calc_tdust_1d_g)(
  double* tdust, double* tgas, double* nh, double* gasgr, double* gamma_isrfa,
  double* isrf, gr_mask_type* itmask, double* trad, int* in, int* is, int* ie,
  int* j, int* k, int* gr_N, int* gr_Size, double* gr_dT, double* gr_Td,
  gr_float* alsp_data_ptr, double* kgr, int* idspecies
);

void FORTRAN_NAME(calc_kappa_gr_g)(
  double* tdust, double* kgr, gr_mask_type* itmask, int* in, int* is, int* ie,
  double* t_subl, int* gr_N, int* gr_Size, double* gr_dT, double* gr_Td,
  gr_float* logalsp_data_ptr, int* idspecies
);

void FORTRAN_NAME(calc_gr_balance_g)(
  double* tdust, double* tgas, double* kgr, double* trad4, double* gasgr,
  double* gamma_isrf, double* nh, gr_mask_type* itmask, double* sol, int* in,
  int* is, int* ie
);

void FORTRAN_NAME(calc_temp1d_cloudy_g)(
  gr_float* d_data_ptr, gr_float* metal_data_ptr, gr_float* e_data_ptr,
  double* rhoH, int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
  double* tgas, double* mmw, double* dom, double* zr, double* temstart,
  double* temend, double* gamma, double* utem, int* imetal,
  long long* clGridRank, long long* clGridDim, double* clPar1, double* clPar2,
  double* clPar3, long long* clDataSize, double* clMMW, gr_mask_type* itmask
);

void FORTRAN_NAME(gaussj_g)(
  int* n, double* a_data_ptr, double* b, int* ierr
);

// in the following interpolate functions, all of the arguments are const
// pointers other than value. (but I have just annotated the subset of
// arguments with const that are necessary to get the desired wrapper function
// signature)

void FORTRAN_NAME(interpolate_1d_g)(
  double* input1, const long long* gridDim, const double* gridPar1,
  double* dgridPar1, long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_2d_g)(
  double* input1, double* input2, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  double* dgridPar2, long long* dataSize, const double* dataField,
  double* value
);

void FORTRAN_NAME(interpolate_3d_g)(
  double* input1, double* input2, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  double* dgridPar2, const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_3dz_g)(
  double* input1, double* input2, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  long long* index2, const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, long long* end_int,
  double* value
);

void FORTRAN_NAME(interpolate_2df3d_g)(
  double* input1, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, long long* index2,
  const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_4d_g)(
  double* input1, double* input2, double* input3, double* input4,
  const long long* gridDim, const double* gridPar1, double* dgridPar1,
  const double* gridPar2, double* dgridPar2, const double* gridPar3,
  double* dgridPar3, const double* gridPar4, double* dgridPar4,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_5d_g)(
  double* input1, double* input2, double* input3, double* input4,
  double* input5, const long long* gridDim, const double* gridPar1,
  double* dgridPar1, const double* gridPar2, double* dgridPar2,
  const double* gridPar3, double* dgridPar3, const double* gridPar4,
  double* dgridPar4, const double* gridPar5, double* dgridPar5,
  long long* dataSize, const double* dataField, double* value
);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* FORTRAN_FN_DECLARATIONS_HPP */
