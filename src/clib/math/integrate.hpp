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
    const Fn& calc_deriv, int& ierror, int nsp, std::vector<double>& dsp,
    std::vector<double>& dsp1, std::vector<double>& dspdot,
    std::vector<double>& dspdot1, std::vector<double>& ddsp,
    FortranView<double**>& jacobian, std::vector<int>& idsp,
    FortranView<double**>& mtrx, std::vector<double>& vec, const double& eps,
    int& i, const int j, const int k, bool enforce_positive_non_NaN) {
  // shorten `GRIMPL_NS::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::GRIMPL_NS::fortran_wrapper;

  int itr = 0;

  double err_max = huge8;  // <- arbitrary initial value that exceeds threshold
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
    for (int jsp = 0; jsp < nsp; jsp++) {
      double dspj = eps * dsp[idsp[jsp]];
      for (int isp = 0; isp < nsp; isp++) {
        if (isp == jsp) {
          dsp1[idsp[isp]] = dsp[idsp[isp]] + dspj;
        } else {
          dsp1[idsp[isp]] = dsp[idsp[isp]];
        }
      }

      calc_deriv(dsp1.data(), dspdot1.data());

      for (int isp = 0; isp < nsp; isp++) {
        if ((dsp[idsp[isp]] == 0.0) &&
            (dspdot1[idsp[isp]] == dspdot[idsp[isp]])) {
          jacobian(idsp[isp], idsp[jsp]) = 0.0;
        } else {
          jacobian(idsp[isp], idsp[jsp]) =
              (dspdot1[idsp[isp]] - dspdot[idsp[isp]]) / dspj;
        }
      }
    }

    for (int isp = 0; isp < nsp; isp++) {
      for (int jsp = 0; jsp < nsp; jsp++) {
        if (isp == jsp) {
          mtrx(isp, jsp) = 1.0 - dt * jacobian(idsp[isp], idsp[jsp]);
        } else {
          mtrx(isp, jsp) = -dt * jacobian(idsp[isp], idsp[jsp]);
        }
      }
    }

    for (int isp = 0; isp < nsp; isp++) {
      vec[isp] = dspdot[idsp[isp]] * dt - ddsp[idsp[isp]];
    }

    // to get more accuracy
    for (int isp = 0; isp < nsp; isp++) {
      vec[isp] = vec[isp] / d(i, j, k);
    }

    ierror = f_wrap::gaussj_g(nsp, mtrx.data(), vec.data());
    if (ierror == 1) {
      return;
    }

    // multiply with density again
    for (int isp = 0; isp < nsp; isp++) {
      vec[isp] = vec[isp] * d(i, j, k);
    }

    for (int isp = 0; isp < nsp; isp++) {
      ddsp[idsp[isp]] = ddsp[idsp[isp]] + vec[isp];
      dsp[idsp[isp]] = dsp[idsp[isp]] + vec[isp];
    }

    if (enforce_positive_non_NaN) {
      for (int isp = 0; isp < nsp; isp++) {
        if (std::isnan(dsp[idsp[isp]]) || (dsp[idsp[isp]] <= 0.)) {
          ierror = 1;
          return;
        }
      }
    }

    err_max = 0.0;
    for (int isp = 0; isp < nsp; isp++) {
      double cur = dsp[idsp[isp]];
      // todo: double check that our behavior, when cur~0, makes sense
      double err = (cur > tiny8) ? std::fabs(vec[isp] / cur) : 0.0;
      err_max = std::fmax(err, err_max);
    }

    itr = itr + 1;
  }
}

}  // namespace integrate
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_INTEGRATE_HPP