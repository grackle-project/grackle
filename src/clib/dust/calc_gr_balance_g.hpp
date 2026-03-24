//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_gr_balance_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_gr_balance_g function from FORTRAN to C++

#ifndef CALC_GR_BALANCE_G_HPP
#define CALC_GR_BALANCE_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"

namespace grackle::impl {

///  Calculate grain heating/cooling balance
///
///  TODO: this docstring should EXPLICITLY document how the action of this
///  function changes based on whether we are using the single-field
///  dust model or the multi-species dust model.
///  In the classic single single-field model, the returned values
///  are emission/absorption rate per unit grain mass with units of
///  erg/s/g.
///  In the multi-species dust model, I **THINK** that the returned
///  values are emission/absorption rate per unit gas mass with units
///  of erg/s/g (But I'm not entirely unsure).
///
/// @param[in]  tdust 1D array to hold dust temperatures
/// @param[in]  tgas 1D array to hold gas temperatures
/// @param[in]  kgr Array to hold computed grain opacities
/// @param[in]  trad4 CMB radiation temperature to the 4th power
/// @param[in]  gasgr Array of gas-grain heat transfer rates
/// @param[in]  gamma_isrf Heating from interstellar radiation field
/// @param[in]  nh 1D array of hydrogen number densities
/// @param[in]  itmask  Iteration mask
/// @param[out] sol Array to hold the computed heating/cooling balance
/// @param[in]  idx_range Index range specifying the portion of the grid to
///     operate on
///
/// @par History
/// written by: Britton Smith, 2019
/// modified: January, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
void calc_gr_balance_g(double* tdust, double* tgas, double* kgr, double trad4,
                       double* gasgr, double* gamma_isrf, double* nh,
                       gr_mask_type* itmask, double* sol, IndexRange idx_range);

}  // namespace grackle::impl

#endif /* CALC_GR_BALANCE_G_HPP */
