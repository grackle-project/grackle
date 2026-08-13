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
    int n_vectors = 3;   // f, vec and dy
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
                               const Fn& calc_deriv_and_jacobian, int n,
                               std::vector<double>& y, double* scratch_ptr,
                               bool enforce_positive_non_NaN) const {
    // shorten `GRIMPL_NS::fortran_wrapper` to `f_wrap` within this function
    namespace f_wrap = ::GRIMPL_NS::fortran_wrapper;

    // I think this is a highly relevant pattern
    int scratch_offset = 0;
    double* dy = scratch_ptr + scratch_offset;
    scratch_offset += n;
    double* vec = scratch_ptr + scratch_offset;
    scratch_offset += n;
    double* f = scratch_ptr + scratch_offset;
    scratch_offset += n;
    FortranView<double**> jacobian(scratch_ptr + scratch_offset, n, n);
    scratch_offset += n * n;
    FortranView<double**> mtrx(scratch_ptr + scratch_offset, n, n);

    for (int i = 0; i < n; i++) {
      dy[i] = 0.0;
    }

    const double max_error_exit_thresh = 1.e-8;
    const int maxiter = 20;
    for (int itr = 0; itr < maxiter; itr++) {
      // calc the time derivatives & the jacobian matrix for the time derivative
      calc_deriv_and_jacobian(y.data(), f, jacobian);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i == j) {
            mtrx(i, j) = 1.0 - dt * jacobian(i, j);
          } else {
            mtrx(i, j) = -dt * jacobian(i, j);
          }
        }
      }

      for (int i = 0; i < n; i++) {
        vec[i] = f[i] * dt - dy[i];
      }

      // to get more accuracy
      for (int i = 0; i < n; i++) {
        vec[i] = vec[i] / local_density;
      }

      // todo: consider adjusting gaussj_g's return value so that its more
      //       consistent with GR_SUCCESS
      int ierror = f_wrap::gaussj_g(n, mtrx.data(), vec);
      if (ierror == 1) {
        return GR_FAIL;
      }

      // multiply with density again
      for (int i = 0; i < n; i++) {
        vec[i] = vec[i] * local_density;
      }

      for (int i = 0; i < n; i++) {
        dy[i] = dy[i] + vec[i];
        y[i] = y[i] + vec[i];
      }

      if (enforce_positive_non_NaN) {
        for (int i = 0; i < n; i++) {
          if (std::isnan(y[i]) || (y[i] <= 0.)) {
            return GR_FAIL;
          }
        }
      }

      double err_max = 0.0;
      for (int i = 0; i < n; i++) {
        double cur = y[i];
        // todo: double check that our behavior, when cur~0, makes sense
        double err = (cur > tiny8) ? std::fabs(vec[i] / cur) : 0.0;
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