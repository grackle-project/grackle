//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_tdust_1d_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_tdust_1d_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "dust/multi_grain_species/opac_calculator.hpp"
#include "dust/passive/analytic_opac.hpp"
#include "grackle.h"
#include "fortran_func_decls.h"
// TODO: to be removed when transcription is done
#include "fortran_func_wrappers.hpp"
#include "phys_constants.h"
#include "support/config.hpp"
#include "utils-cpp.hpp"

#include "calc_tdust_1d_g.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// grain sublimation temperature in passive dust model
static constexpr double passive_dust_model_T_sublimation = 1500.0;

/// @brief Aggregates results of a function that is being iteratively solved
///
/// @note
/// To more fully support Newton-Raphson, we may want to consider adding a slot
/// for the derivative
struct FnEval {
  /// the calculated value
  double f_val;

  /// @brief an associated intermediate value
  ///
  /// The basic premise is that this value is calculated in the process of
  /// computing @ref f_val and we are interested in its value at the location
  /// where the function has been solved.
  double associated_val;
};

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
  auto fn = [&](double Tdust, int i) -> FnEval {
    double kgr = calculator.calc_opac(Tdust, i);
    double sol = Tdust_detail::calc_grain_balance(
        Tdust, tgas[i], kgr, Trad4, gasgr[i], gamma_isrf[i], nh[i]);
    return {sol, kgr};
  };

  // grain opacity from Omukai (2000, equation 17) normalized by
  // the local dust-to-gas ratio, which in this work is 0.934e-2.
  const double kgr1 = 4.0e-4 / 0.00934;

  const double tol = 1.e-5;
  const double bi_tol = 1.e-3;
  const double minpert = 1.e-10;
  const int itmax = 50;
  const int bi_itmax = 30;

  // Locals

  int iter, c_done, c_total, nm_done;

  // Slice Locals

  std::vector<double> kgrplus(buf_len);
  std::vector<double> sol(buf_len);
  std::vector<double> solplus(buf_len);
  // holds dust temperature guess from the last root-finding iteration
  std::vector<double> tdustold(buf_len);
  // holds dust temperature guess for the current root-finding iteration
  std::vector<double> tdustnow(buf_len);
  // holds a value offset from the current dust temperature (tdustnow)
  // it is used to get a finite difference estimate for the derivative
  std::vector<double> tdplus(buf_len);
  // relative finite difference step size
  std::vector<double> pert(buf_len);
  std::vector<double> bi_t_high(buf_len);
  // iteration mask specifies where we use newton's method with finite
  // differences
  std::vector<gr_mask_type> nm_itmask(buf_len);
  // iteration mask specifies where we use bisection
  std::vector<gr_mask_type> bi_itmask(buf_len);

  double pert_i = 1.e-3;

  // Set total cells for calculation

  c_done = 0;
  nm_done = 0;
  c_total = idx_range.i_stop - idx_range.i_start;

  // Set local iteration mask and initial guess

  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    nm_itmask[i] = itmask[i];
    bi_itmask[i] = itmask[i];
    if (nm_itmask[i] != MASK_FALSE) {
      if (Trad >= tgas[i]) {
        tdustnow[i] = Trad;
        nm_itmask[i] = MASK_FALSE;
        bi_itmask[i] = MASK_FALSE;
        c_done = c_done + 1;
        nm_done = nm_done + 1;
      } else if (tgas[i] > passive_dust_model_T_sublimation) {
        // Use bisection if T_gas > grain sublimation temperature.
        nm_itmask[i] = MASK_FALSE;
        nm_done = nm_done + 1;
      } else {
        // we the following is a guess based on the premise that cooling from
        // thermal radiation is balanced by heating from the interstellar
        // radiation field.
        // -> The constants assume that we are using the classic passive dust
        //    model with Tgrain < 200 Kelvin
        // -> striclty speaking, the exponent should be 1/6
        double isrf_balance_guess = std::pow(
            (gamma_isrf[i] / Tdust_detail::sigma_sb_times_4 / kgr1), 0.17);
        tdustnow[i] = std::fmax(Trad, isrf_balance_guess);
        pert[i] = pert_i;
      }

    } else {
      c_done = c_done + 1;
      nm_done = nm_done + 1;
    }
  }

  // Iterate to convergence with Newton's method

  for (iter = 1; iter <= (itmax); iter++) {
    // Loop over slice

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        tdplus[i] = std::fmax(1.e-3, ((1. + pert[i]) * tdustnow[i]));
      }
    }

    // Calculate grain opacities AND heating/cooling balance
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        FnEval eval_rslt = fn(tdustnow[i], i);
        kgr[i] = eval_rslt.associated_val;
        sol[i] = eval_rslt.f_val;
      }
    }
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        // can we get rid of tdplus?
        FnEval eval_rslt = fn(tdplus[i], i);
        // todo: get rid of kgrplus
        kgrplus[i] = eval_rslt.associated_val;
        solplus[i] = eval_rslt.f_val;
      }
    }

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        // Check if the solution has converged (if not prepare the next guess)

        double slope = (solplus[i] - sol[i]) / (pert[i] * tdustnow[i]);

        tdustold[i] = tdustnow[i];
        // tdustnow(i) = tdustnow(i) - (sol(i) / slope)
        tdustnow[i] = std::fmin(tdustnow[i] - (sol[i] / slope), 3e3);

        pert[i] = std::fmax(
            std::fmin(pert[i], (0.5 * std::fabs(tdustnow[i] - tdustold[i]) /
                                tdustnow[i])),
            minpert);

        // If negative solution calculated, give up and wait for bisection step.
        if (tdustnow[i] < Trad) {
          nm_itmask[i] = MASK_FALSE;
          nm_done = nm_done + 1;
          // Check for convergence of solution
        } else if (std::fabs(sol[i] / solplus[i]) < tol) {
          nm_itmask[i] = MASK_FALSE;
          c_done = c_done + 1;
          bi_itmask[i] = MASK_FALSE;
          nm_done = nm_done + 1;
        }

        // if ( nm_itmask(i) )
      }

      // End loop over slice
    }

    // Check for all cells converged
    if (c_done >= c_total) {
      break;
    }

    // Check for all cells done with Newton method
    // This includes attempts where a negative solution was found
    if (nm_done >= c_total) {
      break;
    }

    // End iteration loop for Newton's method
  }

  // If iteration count exceeded, try once more with bisection
  if (c_done < c_total) {
    // set initial guesses
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (bi_itmask[i] != MASK_FALSE) {
        tdustnow[i] = Trad;
        // bi_t_high(i) = tgas(i)
        bi_t_high[i] = 3e3;
      }
    }

    double max_initial_guess = passive_dust_model_T_sublimation;

    // implicitly assumption that
    //   - sol[i] > 0 for tdustnow[i]
    //   - sol[i] < 0 for bi_t_high
    // TODO: we should probably check this assumption (especially for
    //       multi-species dust grains)

    std::vector<double> x_mid(buf_len);
    for (iter = 1; iter <= (bi_itmax); iter++) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (bi_itmask[i] != MASK_FALSE) {
          x_mid[i] = 0.5 * (tdustnow[i] + bi_t_high[i]);
          if (iter == 1) {
            x_mid[i] = std::fmin(x_mid[i], max_initial_guess);
          }
        }
      }

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (bi_itmask[i] != MASK_FALSE) {
          FnEval eval_rslt = fn(x_mid[i], i);
          kgr[i] = eval_rslt.associated_val;
          sol[i] = eval_rslt.f_val;
        }
      }

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (bi_itmask[i] != MASK_FALSE) {
          if (sol[i] > 0.) {
            tdustnow[i] = x_mid[i];
          } else {
            bi_t_high[i] = x_mid[i];
          }

          if ((std::fabs(bi_t_high[i] - tdustnow[i]) / tdustnow[i]) <= bi_tol) {
            bi_itmask[i] = MASK_FALSE;
            c_done = c_done + 1;
          }

          // Check for all cells converged
          if (c_done >= c_total) {
            break;
          }

          // if ( bi_itmask(i) )
        }

        // End loop over slice
      }

      if (c_done >= c_total) {
        break;
      }
      // End iteration loop for bisection
    }

    // If iteration count exceeded with bisection, end of the line.
    if (iter > itmax) {
      OMP_PRAGMA_CRITICAL {
        eprintf(
            "CALC_TDUST_1D_G failed to converge after %d iterations for %d "
            "cells.\n",
            iter, (c_total - c_done));
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
    }
  }

  return;
}

void calc_tdust_1d_g(double* tdust, const double* tgas, const double* nh,
                     const double* gasgr, const double* gamma_isrfa,
                     const double* isrf, const gr_mask_type* itmask,
                     double trad, int buf_len, int gr_N, double gr_dT,
                     const double* gr_Td, const double* alsp_data_, double* kgr,
                     int idspecies, IndexRange idx_range) {
  if (idspecies == 0) {
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