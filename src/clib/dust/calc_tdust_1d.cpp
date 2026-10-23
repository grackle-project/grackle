//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_tdust_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_tdust_1d_g function from FORTRAN to C++

#include <vector>

#include "dust/multi_grain_species/opac_calculator.hpp"
#include "dust/passive/analytic_opac.hpp"
#include "fortran_func_decls.h"
#include "math/solve_eqn.hpp"
#include "support/config.hpp"
#include "utils-cpp.hpp"

#include "calc_tdust_1d.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// grain sublimation temperature in passive dust model
static constexpr double passive_dust_model_T_sublimation = 1500.0;

/// @brief A helper function that helps implement calc_tdust_1d
///
/// The basic premise is that the particulars of the opacity calculation are
/// handled by the `OpacCalculator` template argument.
///
/// In the future, we may want to entirely replace @ref calc_tdust_1d with this
/// function template (if we do that, this will need to be moved to the header)
///
/// @todo
/// Consider using a hybrid Newton-Bisection strategy. The section in Numerical
/// Recipes on Newton-Raphson describes this kind of strategy. I would be
/// shocked if a hybrid scheme wasn't superior to running pure Newton-Raphson
/// and then falling back to Bisection.
template <typename OpacCalculator>
void calc_tdust_1d_(double* tdust, const double* tgas, const double* nh,
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
      } else if (tgas[i] > passive_dust_model_T_sublimation) {
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

  int iter;
  {  // perform newton's method
    double* x = tdustnow.data();
    double max_x = 3000.0;
    // below this x value, we give up!
    double giveup_small_x_threshold = Trad;

    double* f_vals = sol.data();

    solve_eqn::EqnSolveRslt rslt = unchecked_posdomain_fdiff_newton(
        fn, x, solvemask.data(), f_vals, pert.data(), idx_range.i_start,
        idx_range.i_stop, newton_params, giveup_small_x_threshold, max_x);
    iter = rslt.iterations;
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

    double max_initial_guess = passive_dust_model_T_sublimation;
    iter =
        unchecked_bisect(fn, tdustnow.data(), bi_t_high.data(),
                         solvemask.data(), idx_range.i_start, idx_range.i_stop,
                         bi_tol, bi_itmax, max_initial_guess)
            .iterations;

    // If iteration count exceeded with bisection, end of the line.
    // -> Since `(bi_itmax+1) < itmax`, I'm pretty sure the following
    //    if-statement is never invoked
    if (iter > itmax) {
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

    // if (iter .gt. itmax) then
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

void calc_tdust_1d(double* tdust, const double* tgas, const double* nh,
                   const double* gasgr, const double* gamma_isrfa,
                   const double* isrf, const gr_mask_type* itmask, double trad,
                   int buf_len, int gr_N, double gr_dT, const double* gr_Td,
                   const double* alsp_data_, double* kgr,
                   bool analytic_opac_formula, IndexRange idx_range) {
  if (analytic_opac_formula) {
    AnalyticOpacCalc calculator(passive_dust_model_T_sublimation);
    calc_tdust_1d_(tdust, tgas, nh, gasgr, gamma_isrfa, isrf, itmask, trad,
                   buf_len, kgr, idx_range, calculator);
  } else {
    FortranView<const double**> alsp(alsp_data_, gr_N, buf_len);
    std::vector<double> logalsp_data_(gr_N * buf_len);
    FortranView<double**> logalsp(logalsp_data_.data(), gr_N, buf_len);

    int Td_N = gr_N;
    int Td_Size = gr_N;
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        for (int j = 0; j < gr_N; j++) {
          logalsp(j, i) = std::log10(alsp(j, i));
        }
      }
    }

    MultiGrainGrowthOpacCalc calculator(buf_len, Td_N, Td_Size, gr_dT, gr_Td,
                                        logalsp_data_.data());
    calc_tdust_1d_(tdust, tgas, nh, gasgr, gamma_isrfa, isrf, itmask, trad,
                   buf_len, kgr, idx_range, calculator);
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL