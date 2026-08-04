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
#include "support/config.hpp"
#include "support/index_helper.hpp"

#include <cmath>  // std::pow

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Defines constants and logic used for computing dust temperature
namespace Tdust_detail {

/// @brief Four times the Stefan-Boltzmann constant
inline constexpr double sigma_sb_times_4 = 4.0 * sigma_sb_grflt;

/// @brief Calculate grain heating/cooling balance
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
/// @param Tdust dust temperature
/// @param Tgas gas temperature
/// @param kgr grain opacity
/// @param Trad4 CMB radiation temperature to the 4th power
/// @param gasgr gas-grain heat transfer rates
/// @param gamma_isrf Heating from interstellar radiation field
/// @param nH hydrogen number density
///
/// @note
/// If we provided the derivative of @p kgr with respect to Tdust, the
/// analytic derivative of this function is straightforward
///
/// @par History
/// original function written by: Britton Smith, 2019
/// modified: January, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
/// modified: August 2026 by Matthew Abruzzo; converted to scalar calculation
inline double calc_grain_balance(double Tdust, double Tgas, double kgr,
                                 double Trad4, double gasgr, double gamma_isrf,
                                 double nH) {
  return gamma_isrf + sigma_sb_times_4 * kgr * (Trad4 - std::pow(Tdust, 4)) +
         (gasgr * nH * (Tgas - Tdust));
  //  Historically, the following comment was present here:
  //      emission/absorption rate per unit grain mass [erg/s/g]
  //      for Z = Zsun (default)
  //  This comment is **ONLY** correct when the function used as
  //  part of the single-field dust model. See the docstring for
  //  more details.
}

}  // namespace Tdust_detail

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
/// @param[in]  itmask Iteration mask
/// @param[in]  trad CMB ratiation temperature
/// @param[in]  buf_len Length of the 1D slice
/// @param[in]  gr_N Number of temperature points in the grain opacity table
/// @param[in]  gr_dT Temperature spacing of the grain opacity table
/// @param[in]  gr_Td Temperature values of the grain opacity table
/// @param[in]  alsp_data_ Grain opacity table data
/// @param[out] kgr Array to hold computed grain opacities
/// @param[in]  idspecies Flag to solve multiple grain species
/// @param[in]  idx_range Index range specifying the portion of the grid to
///     operate on
///
/// @par History
/// written by: Britton Smith, 2011
/// modified: March, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
///
void calc_tdust_1d_g(double* tdust, const double* tgas, const double* nh,
                     const double* gasgr, const double* gamma_isrfa,
                     const double* isrf, const gr_mask_type* itmask,
                     double trad, int buf_len, int gr_N, double gr_dT,
                     const double* gr_Td, const double* alsp_data_, double* kgr,
                     int idspecies, IndexRange idx_range);

}  // namespace GRIMPL_NAMESPACE_DECL
#endif /* CALC_TDUST_1D_G_HPP */
