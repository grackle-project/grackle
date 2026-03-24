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

#include "dust/calc_gr_balance_g.hpp"
#include "grackle.h"
#include "fortran_func_decls.h"
// TODO: to be removed when transcription is done
#include "fortran_func_wrappers.hpp"
#include "utils-cpp.hpp"

#include "calc_tdust_1d_g.hpp"

void grackle::impl::calc_tdust_1d_g(
    double* tdust, double* tgas, double* nh, double* gasgr,
    const double* gamma_isrfa, const double* isrf, const gr_mask_type* itmask,
    double trad, int in, int gr_N, double* gr_dT, double* gr_Td,
    gr_float* alsp_data_, double* kgr, int* idspecies, IndexRange idx_range) {
  // opacity table of a grain species
  //
  // In some configurations gr_N can be 0 while the backing buffer may still be
  // non-null. The View invariant disallows non-null data with a zero leading
  // extent, so pass nullptr for the zero-length case.
  gr_float* alsp_ptr = (gr_N > 0) ? alsp_data_ : nullptr;
  grackle::impl::View<gr_float**> alsp(alsp_ptr, gr_N, in);
  std::vector<gr_float> logalsp_data_(gr_N * in);
  gr_float* logalsp_ptr = (gr_N > 0) ? logalsp_data_.data() : nullptr;
  grackle::impl::View<gr_float**> logalsp(logalsp_ptr, gr_N, in);
  int Td_Size;
  int Td_N;

  // Parameters

  // grain sublimation temperature
  double t_subl = 1.5e3;  // TODO: should be const
  const double radf = 4.;

  // grain opacity from Omukai (2000, equation 17) normalized by
  // the local dust-to-gas ratio, which in this work is 0.934e-2.
  const double kgr1 = 4.0e-4;

  std::vector<double> gamma_isrf(in);
  const double tol = 1.e-5;
  const double bi_tol = 1.e-3;
  const double minpert = 1.e-10;
  const int itmax = 50;
  const int bi_itmax = 30;

  // Locals

  int i, iter, c_done, c_total, nm_done;

  double pert_i, trad4;

  // Slice Locals

  std::vector<double> kgrplus(in);
  std::vector<double> sol(in);
  std::vector<double> solplus(in);
  std::vector<double> slope(in);
  std::vector<double> tdplus(in);
  std::vector<double> tdustnow(in);
  std::vector<double> tdustold(in);
  std::vector<double> pert(in);
  std::vector<double> bi_t_mid(in);
  std::vector<double> bi_t_high(in);
  std::vector<gr_mask_type> nm_itmask(in);
  std::vector<gr_mask_type> bi_itmask(in);

  pert_i = 1.e-3;

  trad = std::fmax(
      1.,
      trad);  // TODO: do we really want to write in a passed by copy input par?
  trad4 = std::pow(trad, 4);

  // \sum rho_SN kappa_SN / \sum rho_SN ndust_SN

  Td_N = gr_N;
  Td_Size = gr_N;
  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      for (int j = 0; j < gr_N; j++) {
        logalsp(j, i) = std::log10(alsp(j, i));
      }
    }
  }

  // Set total cells for calculation

  c_done = 0;
  nm_done = 0;
  c_total = idx_range.i_end - idx_range.i_start + 1;

  // Set local iteration mask and initial guess

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      gamma_isrf[i] = isrf[i] * gamma_isrfa[i];
    }
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    nm_itmask[i] = itmask[i];
    bi_itmask[i] = itmask[i];
    if (nm_itmask[i] != MASK_FALSE) {
      if (trad >= tgas[i]) {
        tdustnow[i] = trad;
        nm_itmask[i] = MASK_FALSE;
        bi_itmask[i] = MASK_FALSE;
        c_done = c_done + 1;
        nm_done = nm_done + 1;
      } else if (tgas[i] > t_subl) {
        // Use bisection if T_gas > grain sublimation temperature.
        nm_itmask[i] = MASK_FALSE;
        nm_done = nm_done + 1;
      } else {
        tdustnow[i] =
            std::fmax(trad, std::pow((gamma_isrf[i] / radf / kgr1), 0.17));
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

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        tdplus[i] = std::fmax(1.e-3, ((1. + pert[i]) * tdustnow[i]));
      }
    }

    // Calculate grain opacities

    FORTRAN_NAME(calc_kappa_gr_g)(tdustnow.data(), kgr, nm_itmask.data(), &in,
                                  &idx_range.i_start, &idx_range.i_end, &t_subl,
                                  &Td_N, &Td_Size, gr_dT, gr_Td, logalsp.data(),
                                  idspecies);

    FORTRAN_NAME(calc_kappa_gr_g)(tdplus.data(), kgrplus.data(),
                                  nm_itmask.data(), &in, &idx_range.i_start,
                                  &idx_range.i_end, &t_subl, &Td_N, &Td_Size,
                                  gr_dT, gr_Td, logalsp.data(), idspecies);

    // Calculate heating/cooling balance

    // TODO: to be removed when transcription is done
    // FORTRAN_NAME(calc_gr_balance_g)(tdustnow.data(), tgas, kgr, &trad4,
    // gasgr,
    //                                 gamma_isrf.data(), nh, nm_itmask.data(),
    //                                 sol.data(), &in, &idx_range.i_start,
    //                                 &idx_range.i_end);
    // grackle::impl::fortran_wrapper::calc_gr_balance_g(tdustnow.data(), tgas,
    // kgr, &trad4,
    //                                 gasgr, gamma_isrf.data(), nh,
    //                                 nm_itmask.data(), sol.data(), &in,
    //                                 &idx_range.i_start, &idx_range.i_end);
    grackle::impl::calc_gr_balance_g(tdustnow.data(), tgas, kgr, trad4, gasgr,
                                     gamma_isrf.data(), nh, nm_itmask.data(),
                                     sol.data(), idx_range);

    // TODO: to be removed when transcription is done
    // FORTRAN_NAME(calc_gr_balance_g)(tdplus.data(), tgas, kgrplus.data(),
    // &trad4,
    //                                 gasgr, gamma_isrf.data(), nh,
    //                                 nm_itmask.data(), solplus.data(), &in,
    //                                 &idx_range.i_start, &idx_range.i_end);
    // grackle::impl::fortran_wrapper::calc_gr_balance_g(tdplus.data(), tgas,
    // kgrplus.data(), &trad4,
    //                                 gasgr, gamma_isrf.data(), nh,
    //                                 nm_itmask.data(), solplus.data(), &in,
    //                                 &idx_range.i_start, &idx_range.i_end);
    grackle::impl::calc_gr_balance_g(
        tdplus.data(), tgas, kgrplus.data(), trad4, gasgr, gamma_isrf.data(),
        nh, nm_itmask.data(), solplus.data(), idx_range);

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (nm_itmask[i] != MASK_FALSE) {
        // Use Newton's method to solve for Tdust

        slope[i] = (solplus[i] - sol[i]) / (pert[i] * tdustnow[i]);

        tdustold[i] = tdustnow[i];
        // tdustnow(i) = tdustnow(i) - (sol(i) / slope(i))
        tdustnow[i] = std::fmin(tdustnow[i] - (sol[i] / slope[i]), 3e3);

        pert[i] = std::fmax(
            std::fmin(pert[i], (0.5 * std::fabs(tdustnow[i] - tdustold[i]) /
                                tdustnow[i])),
            minpert);

        // If negative solution calculated, give up and wait for bisection step.
        if (tdustnow[i] < trad) {
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
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (bi_itmask[i] != MASK_FALSE) {
        tdustnow[i] = trad;
        // bi_t_high(i) = tgas(i)
        bi_t_high[i] = 3e3;
      }
    }

    for (iter = 1; iter <= (bi_itmax); iter++) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (bi_itmask[i] != MASK_FALSE) {
          bi_t_mid[i] = 0.5 * (tdustnow[i] + bi_t_high[i]);
          if (iter == 1) {
            bi_t_mid[i] = std::fmin(bi_t_mid[i], t_subl);
          }
        }
      }

      FORTRAN_NAME(calc_kappa_gr_g)(bi_t_mid.data(), kgr, bi_itmask.data(), &in,
                                    &idx_range.i_start, &idx_range.i_end,
                                    &t_subl, &Td_N, &Td_Size, gr_dT, gr_Td,
                                    logalsp.data(), idspecies);

      // TODO: to be removed when transcription is done
      // FORTRAN_NAME(calc_gr_balance_g)(bi_t_mid.data(), tgas, kgr, &trad4,
      // gasgr,
      //                                 gamma_isrf.data(), nh,
      //                                 bi_itmask.data(), sol.data(), &in,
      //                                 &idx_range.i_start, &idx_range.i_end);
      // grackle::impl::fortran_wrapper::calc_gr_balance_g(bi_t_mid.data(),
      // tgas, kgr, &trad4, gasgr,
      //                                 gamma_isrf.data(), nh,
      //                                 bi_itmask.data(), sol.data(), &in,
      //                                 &idx_range.i_start, &idx_range.i_end);
      grackle::impl::calc_gr_balance_g(bi_t_mid.data(), tgas, kgr, trad4, gasgr,
                                       gamma_isrf.data(), nh, bi_itmask.data(),
                                       sol.data(), idx_range);

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (bi_itmask[i] != MASK_FALSE) {
          if (sol[i] > 0.) {
            tdustnow[i] = bi_t_mid[i];
          } else {
            bi_t_high[i] = bi_t_mid[i];
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
  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Check for bad solutions
      if (tdustnow[i] < 0.) {
        OMP_PRAGMA_CRITICAL {
          eprintf(
              "CALC_TDUST_1D_G Newton method -  T_dust < 0: i =  %d j =  %d k "
              "=  %d nh =  %g t_gas =  %g t_rad =  %g t_dust =  %g\n",
              i, idx_range.jp1, idx_range.kp1, nh[i], tgas[i], trad,
              tdustnow[i]);
        }
        // ERROR_MESSAGE
      }

      tdust[i] = tdustnow[i];
    }
  }

  return;
}
