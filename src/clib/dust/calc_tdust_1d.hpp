//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_tdust_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_tdust_1d_g function from FORTRAN to C++

#ifndef DUST_CALC_TDUST_1D_HPP
#define DUST_CALC_TDUST_1D_HPP

#include "fortran_func_decls.h"  // gr_mask_int
#include "math/solve_eqn.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "utils-cpp.hpp"

#include <cmath>  // std::pow
#include <vector>

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Defines constants and logic used for computing dust temperature
namespace Tdust_detail {

/// @brief grain sublimation temperature in passive dust model
inline constexpr double passive_dust_model_T_sublimation = 1500.0;

/// @brief Four times the Stefan-Boltzmann constant
inline constexpr double sigma_sb_times_4 = 4.0 * constants::sigma_sb_grflt;

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
/// @param[in]  nominal_Trad CMB ratiation temperature
/// @param[in]  buf_len Length of the 1D slice
/// @param[out] kgr Array to hold computed grain opacities
/// @param[in]  idx_range Index range specifying the portion of the grid to
///     operate on
/// @param[in]  calculator Object of template type ``OpacCalcuator`` that
///     supports a method with the signature given by
///     ``double OpacCalculator::calc_opac(double tdust, int i) const``.
///     This method returns the opacity given for the cell in the index range
///     with an index of ``i`` given a dust temperature ``tdust``
///
/// @par History
/// written by: Britton Smith, 2011
/// modified: March, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
///
/// @todo
/// Consider using a hybrid Newton-Bisection strategy. The section in Numerical
/// Recipes on Newton-Raphson describes this kind of strategy. I would be
/// shocked if a hybrid scheme wasn't superior to running pure Newton-Raphson
/// and then falling back to Bisection.
template <typename OpacCalculator>
void calc_tdust_1d(double* tdust, const double* tgas, const double* nh,
                   const double* gasgr, const double* gamma_isrfa,
                   const double* isrf, const gr_mask_type* itmask,
                   double nominal_Trad, int buf_len, double* kgr,
                   IndexRange idx_range, const OpacCalculator& calculator) {
  const double Trad = std::fmax(1., nominal_Trad);
  const double Trad4 = std::pow(Trad, 4);

  // fill a buffer with the heating rate from the interstellar radiation field
  std::vector<double> gamma_isrf(buf_len);
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      gamma_isrf[i] = isrf[i] * gamma_isrfa[i];
    }
  }

  // define the function that we will perform root finding on
  //
  // It computes the grain opacity and grain's heating/cooling rate given a dust
  // temperature. The root finding hunts for the dust temperature where the
  // grain's heating and cooling are balanced (we are also interested in
  // returning the associated opacity).
  // -> we should think about directly computing analytic derivatives. This
  //    shouldn't be too hard and would significantly improve the robustness
  //    (and speed) of the newton solver pass
  auto fn = [&](double Tdust, int i) -> double {
    double kgr = calculator.calc_opac(Tdust, i);
    double sol = Tdust_detail::calc_grain_balance(
        Tdust, tgas[i], kgr, Trad4, gasgr[i], gamma_isrf[i], nh[i]);
    return sol;
  };

  // grain opacity from Omukai (2000, equation 17) normalized by
  // the local dust-to-gas ratio, which in this work is 0.934e-2.
  const double kgr1 = 4.0e-4 / 0.00934;

  const double bi_tol = 1.e-3;
  const int itmax = 50;
  const int bi_itmax = 30;

  const solve_eqn::FDiffNewtonParams newton_params{.maxiter = itmax,
                                                   .perturb_init = 1.e-3,
                                                   .perturb_min = 1.e-10,
                                                   .tol = 1.e-5};

  std::vector<double> sol(buf_len);
  // holds dust temperature guess for the current root-finding iteration
  std::vector<double> tdustnow(buf_len);
  // relative finite difference step size
  std::vector<double> pert(buf_len);
  std::vector<double> bi_t_high(buf_len);

  // specifies the status of each index (a solver skips a calculation if an
  // index corresponds to SolveStatus::CONVERGED)
  std::vector<solve_eqn::SolveStatus> solvemask(
      buf_len, solve_eqn::SolveStatus::CONVERGED);

  // setup solvemask for use with newton-raphson and set initial guess
  bool at_least_one_bisection = false;
  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      if (Trad >= tgas[i]) {
        tdustnow[i] = Trad;
        solvemask[i] = solve_eqn::SolveStatus::CONVERGED;
      } else if (tgas[i] > Tdust_detail::passive_dust_model_T_sublimation) {
        // we'll be skipping the newton method if T_gas exceeds grain
        // sublimation temperature (we'll wait to solve for Tdust until we
        // we use the bisection method)
        solvemask[i] = solve_eqn::SolveStatus::SKIP_SOLVE;
        at_least_one_bisection = true;
      } else {
        // the following is a guess based on the premise that cooling from
        // thermal radiation is balanced by heating from the interstellar
        // radiation field.
        // -> The constants assume that we are using the classic passive dust
        //    model with Tgrain < 200 Kelvin
        // -> strictly speaking, the exponent should be 1/6
        double isrf_balance_guess = std::pow(
            (gamma_isrf[i] / Tdust_detail::sigma_sb_times_4 / kgr1), 0.17);
        tdustnow[i] = std::fmax(Trad, isrf_balance_guess);
        solvemask[i] = solve_eqn::SolveStatus::UNCONVERGED;
      }

    } else {
      solvemask[i] = solve_eqn::SolveStatus::CONVERGED;
    }
  }

  {  // perform newton's method
    double* x = tdustnow.data();
    double max_x = 3000.0;
    // below this x value, we give up!
    double giveup_small_x_threshold = Trad;

    double* f_vals = sol.data();

    solve_eqn::EqnSolveRslt rslt = unchecked_posdomain_fdiff_newton(
        fn, x, solvemask.data(), f_vals, pert.data(), idx_range.i_start,
        idx_range.i_stop, newton_params, giveup_small_x_threshold, max_x);
    at_least_one_bisection = at_least_one_bisection || !rslt.all_solved;
  }

  if (at_least_one_bisection) {
    // set initial guesses solvemask for bisection
    //
    // There's an implicit assumption here that
    // - fn(tdustnow[i], i).f_val > 0
    // - fn(bi_t_high[i], i).f_val < 0
    //
    // We may want to check that assumption is indeed correct. While its
    // unlikely to be false, the iteration will fail silently (with an
    // absolute garbage result) on the off chance that it comes up
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (solvemask[i] == solve_eqn::SolveStatus::CONVERGED) {
        solvemask[i] = solve_eqn::SolveStatus::CONVERGED;
      } else {  // we want SolveStatus::UNCONVERGED and SolveStatus::SKIP_SOLVE
        solvemask[i] = solve_eqn::SolveStatus::UNCONVERGED;
        tdustnow[i] = Trad;
        // bi_t_high(i) = tgas(i)
        bi_t_high[i] = 3e3;
      }
    }

    double max_initial_guess = Tdust_detail::passive_dust_model_T_sublimation;
    int iter =
        unchecked_bisect(fn, tdustnow.data(), bi_t_high.data(),
                         solvemask.data(), idx_range.i_start, idx_range.i_stop,
                         bi_tol, bi_itmax, max_initial_guess)
            .iterations;

    // If iteration count exceeded with bisection, end of the line.
    if (iter > bi_itmax) {
      int n_unconverged = 0;
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        n_unconverged += (solvemask[i] != solve_eqn::SolveStatus::CONVERGED);
      }
      OMP_PRAGMA_CRITICAL {
        eprintf(
            "CALC_TDUST_1D_G failed to converge after %d iterations for %d "
            "cells.\n",
            iter, n_unconverged);
      }
    }
  }

  // Copy values back to thrown slice
  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Check for bad solutions
      if (tdustnow[i] < 0.) {
        OMP_PRAGMA_CRITICAL {
          eprintf(
              "CALC_TDUST_1D_G Newton method -  T_dust < 0: i =  %d j =  %d k "
              "=  %d nh =  %g t_gas =  %g t_rad =  %g t_dust =  %g\n",
              i, idx_range.jp1, idx_range.kp1, nh[i], tgas[i], Trad,
              tdustnow[i]);
        }
        // ERROR_MESSAGE
      }

      tdust[i] = tdustnow[i];
      kgr[i] = calculator.calc_opac(tdust[i], i);
    }
  }

  return;
}

}  // namespace GRIMPL_NAMESPACE_DECL
#endif  // DUST_CALC_TDUST_1D_HPP
