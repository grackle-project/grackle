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

class StiffNewtonRaphson {
  /// the maximum number of variables that can be evolved in the current
  /// Grackle configuration
  int max_evolved_variables_;

public:  // public API
  StiffNewtonRaphson() = delete;

  /// @brief primary constructor
  explicit StiffNewtonRaphson(int max_evolved_variables)
      : max_evolved_variables_{max_evolved_variables} {}

  // upper bound on the scratch space needed in step
  int num_scratch_buf_elements() const noexcept {
    int n_vectors = 3;   // dspdot, vec and ddsp
    int n_matrices = 2;  // jacobian and mtrx
    int entries_per_matrix = max_evolved_variables_ * max_evolved_variables_;
    return n_vectors * max_evolved_variables_ + n_matrices * entries_per_matrix;
  }

  /// attempts to integrate over a timestep dt
  ///
  /// @return GR_SUCCESS indicates that the solution converged. Other values
  ///     denote a problem
  template <typename Fn>
  GRIMPL_FORCE_INLINE int step(double dt, double local_density,
                               const Fn& calc_deriv_and_jacobian, int nsp,
                               std::vector<double>& dsp, double* scratch_ptr,
                               bool enforce_positive_non_NaN) const {
    // shorten `GRIMPL_NS::fortran_wrapper` to `f_wrap` within this function
    namespace f_wrap = ::GRIMPL_NS::fortran_wrapper;

    // I think this is a highly relevant pattern
    int scratch_offset = 0;
    double* ddsp = scratch_ptr + scratch_offset;
    scratch_offset += nsp;
    double* vec = scratch_ptr + scratch_offset;
    scratch_offset += nsp;
    double* dspdot = scratch_ptr + scratch_offset;
    scratch_offset += nsp;
    FortranView<double**> jacobian(scratch_ptr + scratch_offset, nsp, nsp);
    scratch_offset += nsp * nsp;
    FortranView<double**> mtrx(scratch_ptr + scratch_offset, nsp, nsp);

    for (int i = 0; i < nsp; i++) {
      ddsp[i] = 0.0;
    }

    const double max_error_exit_thresh = 1.e-8;
    const int maxiter = 20;
    for (int itr = 0; itr < maxiter; itr++) {
      // calc the time derivatives & the jacobian matrix for the time derivative
      calc_deriv_and_jacobian(dsp.data(), dspdot, jacobian);

      for (int isp = 0; isp < nsp; isp++) {
        for (int jsp = 0; jsp < nsp; jsp++) {
          if (isp == jsp) {
            mtrx(isp, jsp) = 1.0 - dt * jacobian(isp, jsp);
          } else {
            mtrx(isp, jsp) = -dt * jacobian(isp, jsp);
          }
        }
      }

      for (int isp = 0; isp < nsp; isp++) {
        vec[isp] = dspdot[isp] * dt - ddsp[isp];
      }

      // to get more accuracy
      for (int isp = 0; isp < nsp; isp++) {
        vec[isp] = vec[isp] / local_density;
      }

      // todo: consider adjusting gaussj_g's return value so that its more
      //       consistent with GR_SUCCESS
      int ierror = f_wrap::gaussj_g(nsp, mtrx.data(), vec);
      if (ierror == 1) {
        return GR_FAIL;
      }

      // multiply with density again
      for (int isp = 0; isp < nsp; isp++) {
        vec[isp] = vec[isp] * local_density;
      }

      for (int isp = 0; isp < nsp; isp++) {
        ddsp[isp] = ddsp[isp] + vec[isp];
        dsp[isp] = dsp[isp] + vec[isp];
      }

      if (enforce_positive_non_NaN) {
        for (int isp = 0; isp < nsp; isp++) {
          if (std::isnan(dsp[isp]) || (dsp[isp] <= 0.)) {
            return GR_FAIL;
          }
        }
      }

      double err_max = 0.0;
      for (int isp = 0; isp < nsp; isp++) {
        double cur = dsp[isp];
        // todo: double check that our behavior, when cur~0, makes sense
        double err = (cur > tiny8) ? std::fabs(vec[isp] / cur) : 0.0;
        err_max = std::fmax(err, err_max);
      }

      if (err_max <= max_error_exit_thresh) {
        return GR_SUCCESS;
      }
    }

    return GR_FAIL;  // only reached if iterations exceeded maxiter
  }
};

}  // namespace integrate
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_INTEGRATE_HPP