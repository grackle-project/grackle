//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_tdust_1d_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_tdust_1d_g function from FORTRAN to C++

#ifndef CALC_TDUST_1D_G_HPP
#define CALC_TDUST_1D_G_HPP

#include "fortran_func_decls.h"  // gr_mask_int
#include "grackle.h"             // gr_float
#include "index_helper.h"

namespace grackle::impl {

///  Calculate equilibrium dust temperature
///
///  TODO: this docstring, and the docstrings of all helper functions should
///  EXPLICITLY document how the meaning of arguments change based on
///  whether we are using the single-field dust model or the
///  multi-species dust model. The different meaning of each variable
///  gets VERY confusing
///
/// @param[out] tdust 1D array to hold dust temperatures
/// @param[in]  tgas 1D array to hold gas temperatures
/// @param[in]  nh 1D array of hydrogen number densities
/// @param[in]  gasgr Array of gas-grain heat transfer rates
/// @param[in]  gamma_isrfa Heating from interstellar radiation field
/// @param[in]  isrf Interstellar radiation field strength in Habing units
/// @param[in]  itmask  Itaration mask
/// @param[in]  trad CMB ratiation temperature
/// @param[in]  in Length of the 1D slice
/// @param[in]  gr_N Number of temperature points in the grain opacity table
/// @param[in]  gr_dT Temperature spacing of the grain opacity table
/// @param[in]  gr_Td Temperature values of the grain opacity table
/// @param[in]  alsp_data_ Grain opacity table data
/// @param[in]  kgr Array to hold computed grain opacities
/// @param[in]  idspecies Array of grain species IDs (only used in certain
/// configurations)
/// @param[in]  idx_range Index range specifying the portion of the grid to
///     operate on
///
/// @par History
/// written by: Britton Smith, 2011
/// modified: January, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
void calc_tdust_1d_g(double* tdust, double* tgas, double* nh, double* gasgr,
                     const double* gamma_isrfa, const double* isrf, const gr_mask_type* itmask,
                     double trad, int in, int gr_N, double* gr_dT,
                     double* gr_Td, gr_float* alsp_data_, double* kgr,
                     int* idspecies, IndexRange idx_range);

}  // namespace grackle::impl
#endif /* CALC_TDUST_1D_G_HPP */
