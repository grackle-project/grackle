//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines/declares generic functionality for integrating ordinary
/// differential equations. The basic premise is that this machinery should
/// generally written in a generic enough form so that we can unit test it.
///
//===----------------------------------------------------------------------===//

#ifndef MATH_INTEGRATE_HPP
#define MATH_INTEGRATE_HPP

#include <cmath>
#include <vector>

#include "grackle.h"
#include "../fortran_func_wrappers.hpp"  // GRIMPL_NS::fortran_wrapper::gaussj_g
#include "../support/config.hpp"
#include "../support/View.hpp"
#include "../time_deriv_0d.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace integrate {

/// attempts to integrate over a timestep dt
template <typename Fn>
GRIMPL_FORCE_INLINE void stiff_newton_raphson(
    double dt, const int*& imp_eng, FortranView<double***>& d,
    const Fn& calc_deriv, int& ierror, int& nsp, int& isp, int& jsp,
    double& dspj, double err_max, std::vector<double>& dsp,
    std::vector<double>& dsp1, std::vector<double>& dspdot,
    std::vector<double>& dspdot1, std::vector<double>& ddsp,
    FortranView<double**>& jacobian, std::vector<int>& idsp,
    FortranView<double**>& mtrx, std::vector<double>& vec, const double& eps,
    int& i, const int j, const int k, bool enforce_positive_non_NaN) {
  // shorten `GRIMPL_NS::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::GRIMPL_NS::fortran_wrapper;

  int itr = 0;

  while (err_max > 1.e-8) {
    if (itr >= 20) {
      ierror = 1;
      return;
    }

    // calc the time derivatives
    calc_deriv(dsp.data(), dspdot.data());

    // fill in the jacobian matrix for the time derivative
    // -> to accomplish this, we use finite differences to estimate
    //    partial derivative for each evolved variable (i.e. the species
    //    densities and possibly the total energy)
    for (jsp = 1; jsp <= (nsp); jsp++) {
      dspj = eps * dsp[idsp[jsp - 1]];
      for (isp = 1; isp <= (nsp); isp++) {
        if (isp == jsp) {
          dsp1[idsp[isp - 1]] = dsp[idsp[isp - 1]] + dspj;
        } else {
          dsp1[idsp[isp - 1]] = dsp[idsp[isp - 1]];
        }
      }

      calc_deriv(dsp1.data(), dspdot1.data());

      for (isp = 1; isp <= (nsp); isp++) {
        if ((dsp[idsp[isp - 1]] == 0.e0) &&
            (dspdot1[idsp[isp - 1]] == dspdot[idsp[isp - 1]])) {
          jacobian(idsp[isp - 1], idsp[jsp - 1]) = 0.e0;
        } else {
          jacobian(idsp[isp - 1], idsp[jsp - 1]) =
              (dspdot1[idsp[isp - 1]] - dspdot[idsp[isp - 1]]) / dspj;
        }
      }
    }

    for (isp = 1; isp <= (nsp); isp++) {
      for (jsp = 1; jsp <= (nsp); jsp++) {
        if (isp == jsp) {
          mtrx(isp - 1, jsp - 1) =
              1.e0 - dt * jacobian(idsp[isp - 1], idsp[jsp - 1]);
        } else {
          mtrx(isp - 1, jsp - 1) = -dt * jacobian(idsp[isp - 1], idsp[jsp - 1]);
        }
      }
    }

    for (isp = 1; isp <= (nsp); isp++) {
      vec[isp - 1] = dspdot[idsp[isp - 1]] * dt - ddsp[idsp[isp - 1]];
    }

    // to get more accuracy
    for (isp = 1; isp <= (nsp); isp++) {
      vec[isp - 1] = vec[isp - 1] / d(i, j, k);
    }

    ierror = f_wrap::gaussj_g(nsp, mtrx.data(), vec.data());
    if (ierror == 1) {
      return;
    }

    // multiply with density again
    for (isp = 1; isp <= (nsp); isp++) {
      vec[isp - 1] = vec[isp - 1] * d(i, j, k);
    }

    for (isp = 1; isp <= (nsp); isp++) {
      ddsp[idsp[isp - 1]] = ddsp[idsp[isp - 1]] + vec[isp - 1];
      dsp[idsp[isp - 1]] = dsp[idsp[isp - 1]] + vec[isp - 1];
    }

    if (enforce_positive_non_NaN) {
      for (isp = 1; isp <= (nsp); isp++) {
        if (std::isnan(dsp[idsp[isp - 1]]) || (dsp[idsp[isp - 1]] <= 0.)) {
          ierror = 1;
          return;
        }
      }
    }

    err_max = 0.e0;
    for (isp = 1; isp <= (nsp); isp++) {
      double err;
      if (dsp[idsp[isp - 1]] > tiny8) {
        err = grackle::impl::dabs(vec[isp - 1] / dsp[idsp[isp - 1]]);
      } else {
        err = 0.e0;
      }
      if (err > err_max) {
        err_max = err;
      }
    }

    itr = itr + 1;
  }
}

}  // namespace integrate
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_INTEGRATE_HPP