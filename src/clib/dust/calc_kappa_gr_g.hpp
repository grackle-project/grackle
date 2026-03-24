//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_kappa_gr_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_kappa_gr_g function from FORTRAN to C++

#ifndef CALC_KAPPA_GR_G_HPP
#define CALC_KAPPA_GR_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"

namespace grackle::impl {

///  Calculate grain plank mean opacity
///
///  TODO: this docstring should EXPLICITLY document how the action of this
///  function changes based on whether we are using the single-field
///  dust model or the multi-species dust model.
///  In the classic single single-field model, the returned opacities
///  have units of cm^2/g and they are measured "per unit grain mass"
///  In the multi-species dust model, I THINK that the returned
///  opacities have units of cm^2/g and they are measured "per unit
///  gas mass."
///
/// @param[in]  tdust 1D array to hold gas temperatures
/// @param[out] kgr Array to hold computed grain opacities
///   configurations)
/// @param[in]  itmask Iteration mask
/// @param[in]  in i dimension of 3D fields
/// @param[in]  idx_range Index range specifying the portion of the grid to
///     operate on
/// @param[in]  t_subl Grain sublimation temperature
/// @param[in]  gr_N Number of temperature points in the grain opacity table
/// @param[in]  gr_Size Number of temperature points in the grain opacity table
/// @param[in]  gr_dT Temperature spacing of the grain opacity table
/// @param[in]  gr_Td Temperature values of the grain opacity table
/// @param[in]  logalsp_data_ Grain opacity table data
/// @param[in]  idspecies Array of grain species IDs (only used in certain
///
/// @par History
/// written by: Britton Smith, September, 2011
/// modified: March, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
void calc_kappa_gr_g(double* tdust, double* kgr, const gr_mask_type* itmask, int in,
                     IndexRange idx_range, const double* t_subl, int* gr_N,
                     const int* gr_Size, double* gr_dT, double* gr_Td,
                     gr_float* logalsp_data_, int idspecies);

}  // namespace grackle::impl
#endif /* CALC_KAPPA_GR_G_HPP */
