// See LICENSE file for license and copyright information

/// @file calc_grain_size_increment_1d-cpp.h
/// @brief Declares signature of calc_grain_size_increment_species_1d

// This file was initially generated automatically during conversion of the
// calc_grain_size_increment_species_1d function from FORTRAN to C++

#ifndef CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP
#define CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

namespace grackle::impl {

void calc_grain_size_increment_species_1d_cpp(
  int igrgr, const gr_mask_type* itmask, int SN0_N, int in, int jn, int kn,
  int is, int ie, int j, int k, gr_float* d_data_, int nSN,
  const gr_float* dsp_data_, gr_float* SN_metal_data_, double* SN_fsp,
  double* SN_r0sp_data_, double ssp, double* sgsp, double* alsp_data_,
  int* gr_N, int gr_Size, double* SN_kp0sp_data_
);

}

#endif /* CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP */
