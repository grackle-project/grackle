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

/// @brief A class for integrating stiff systems of differential equations
///
/// This is intended to be a general-purpose tool that can be tested against
/// various systems of differential equations.
class StiffNewtonRaphson {
  /// the maximum number of variables that can be evolved in the current
  /// Grackle configuration
  int max_evolved_variables_;

  /// integration is considered converged when the max relative difference
  /// magnitude of any vector component does not exceed this value
  double max_rtol_;

  /// maximum number of iterations
  int maxiter_;

public:  // public API
  StiffNewtonRaphson() = delete;

  /// @brief primary constructor
  StiffNewtonRaphson(int max_evolved_variables, double max_rtol, int maxiter)
      : max_evolved_variables_{max_evolved_variables},
        max_rtol_{max_rtol},
        maxiter_{maxiter} {}

  // upper bound on the scratch space needed in step
  int num_scratch_buf_elements() const noexcept {
    int n_vectors = 3;   // f, vec and ycur_minus_ystart
    int n_matrices = 2;  // jacobian_f and mtrx
    int entries_per_matrix = max_evolved_variables_ * max_evolved_variables_;
    return n_vectors * max_evolved_variables_ + n_matrices * entries_per_matrix;
  }

  /// @brief performs a single numerical integration step for a system of ODEs
  ///
  /// For the sake of discussion:
  /// - let `y(t)` specify a vector representing the state of the system of
  ///   ODEs at a time `t`
  /// - let `dy/dt = f(y)`, where `f(y)` is a vector valued function that
  ///   effectively specifies the differential equations
  ///
  /// This method uses Newton-Raphson to solve for the value of `y` at the end
  /// of the timestep, given the value at the start of the timestep. We are
  /// specifically solving the equation for updating a system of ODEs according
  /// to the backward Euler method (aka the 1-stage backwards differentiation
  /// formula).
  ///
  /// @param[in,out] y buffer of @p n elements. It specifies the initial state
  ///     for the system of ODEs and will be updated (in-place) to hold the
  ///     state at the end of the timestep.
  /// @param[in]     h the timestep
  /// @param[in]     local_density is used for rescaling. Frankly, this
  ///     parameter is a historical artifact that should be removed (when we're
  ///     prepared to update the gold-standard). All rescaling should occur
  ///     outside of this function.
  /// @param[in]     calc_deriv_and_jacobian function object with the signature
  ///     `void f(const double* y_buf, double* f_buf, double* jacobian_f_buf)`.
  ///     It should treat `y_buf` argument as an input buffer specifyng a state
  ///     for the system of ODEs. When called, it should fill `f_buf` with the
  ///     the value of `f(y)` for the state in `y_buf` and `jacobian_f_buf`
  ///     with jacobian of `f(y)` (evaluated for the state in `y_buf`). The
  ///     function should expect `y_buf` and `f_buf` to have @p n elements
  ///     `jacobian_f_buf` to have space for @p n * @p n elements.
  /// @param[in]     n number of equations in the system of ODE. This must
  ///     **NOT** exceed the `max_evolved_variables` argument passed to the
  ///     constructor of this type.
  /// @param[in]    scratch_ptr buffer of scratch space for in the calculation.
  ///     The expected length is provided by the
  ///     @ref StiffNewtonRaphson::num_scratch_buf_elements method.
  /// @param[in]    enforce_positive_non_NaN when `true`, the calculation will
  ///     perform an explicit check that each guessed state of the differential
  ///     equations only contains positive, non-NaN components.
  ///
  /// @return GR_SUCCESS indicates that the solution converged. Other values
  ///     denote a problem
  ///
  /// High Level Algorithm Overview
  /// -----------------------------
  /// At a high level, we're looking to compute the value of yₑ given by:
  /// >     yₑ = yₛ + h * f(yₑ)
  /// where
  /// - y denotes the (vector-valued) state of a system of a ODEs
  /// - f(y) is a vector-valued function that gives the time derivative of y
  /// - h is the size of a timestep
  /// - yₛ denotes the state at the system at the start of the timestep
  /// - yₑ denotes the state at the end of the timestep
  ///
  /// To solve for yₑ, let's define:
  /// >    g(yₑ) = yₑ - yₛ - h * f(yₑ)
  /// We'll use the Newton-Raphson Method to find the root of g(yₑ), which
  /// coincides with the desired value of yₑ.
  ///
  /// If J_g(y) specifies the Jacobian matrix (evaluated for some value y) for
  /// the function g(y), then the Newton-Raphson methods tells us that
  /// >    yₑ,ₖ₊₁ = yₑ,ₖ - (J_g(yₑ,ₖ))⁻¹ ⋅ g(yₑ,ₖ)
  ///
  /// More Detail
  /// -----------
  /// Let's further define J_f(y) as the Jacobian (evaluated for some value y)
  /// for the function f(y). We can then write
  /// >    J_g(y) = I - h * J_f(y)
  /// where I represents the identity matrix.
  ///
  /// We can use that formula (and the definition for g(y)) to rewrite the
  /// formula for yₑ,ₖ₊₁ as:
  /// >   yₑ,ₖ₊₁ = yₑ,ₖ - (I - h * J_f(yₑ,ₖ))⁻¹ ⋅ (yₑ,ₖ - yₛ - h * f(yₑ,ₖ))
  ///
  /// We can rewrite this as:
  /// >   δₖ = (I - h * J_f(yₑ,ₖ))⁻¹ ⋅ (h * f(yₑ,ₖ) - (yₑ,ₖ - yₛ))
  /// >   yₑ,ₖ₊₁ = yₑ,ₖ + δₖ
  /// This is the equation implemented in this function
  template <typename Fn>
  GRIMPL_FORCE_INLINE int step(double* y, double h, double local_density,
                               const Fn& calc_deriv_and_jacobian, int n,
                               double* scratch_ptr,
                               bool enforce_positive_non_NaN) const {
    // shorten `GRIMPL_NS::fortran_wrapper` to `f_wrap` within this function
    namespace f_wrap = ::GRIMPL_NS::fortran_wrapper;

    // I suspect that this practice of creating local buffers from a scratch
    // buffer is a useful idiom that we'll want to repeat
    int scratch_offset = 0;
    double* ycur_minus_ystart = scratch_ptr + scratch_offset;
    scratch_offset += n;
    double* vec = scratch_ptr + scratch_offset;
    scratch_offset += n;
    double* f = scratch_ptr + scratch_offset;
    scratch_offset += n;
    FortranView<double**> jacobian_f(scratch_ptr + scratch_offset, n, n);
    scratch_offset += n * n;
    FortranView<double**> mtrx(scratch_ptr + scratch_offset, n, n);

    // initialize ycur_minus_ystart (i.e. difference between current estimate
    // for y at the end of the timestep and y at the start of the timestep)
    for (int i = 0; i < n; i++) {
      ycur_minus_ystart[i] = 0.0;
    }

    // as we enter the loop, note that `y` holds `yₑ,₀ = yₛ`
    for (int k = 0; k < maxiter_; k++) {
      // evaluate f(yₑ,ₖ) and the jacobian matrix for f at yₑ,ₖ.
      calc_deriv_and_jacobian(y, f, jacobian_f);

      // store the value of (I - h * J_f(yₑ,ₖ)) inside mtrx
      // -> I is the identity matrix and J_f is the jacobian matrix for f
      //    evaluated at the current guess for yₑ.
      // -> h is the size of the timestep
      // -> as the derivation of the method in the docstring shows, `mtrx`, is
      //    technically the jacobian of the function for which we're performing
      //    root finding
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          if (i == j) {
            mtrx(i, j) = 1.0 - h * jacobian_f(i, j);
          } else {
            mtrx(i, j) = -h * jacobian_f(i, j);
          }
        }
      }

      // fill vec with the value of (h * f(yₑ,ₖ) - (yₑ,ₖ - yₛ))
      for (int i = 0; i < n; i++) {
        vec[i] = f[i] * h - ycur_minus_ystart[i];
      }

      // to get more accuracy  (TODO: shouldn't be part of this fn)
      for (int i = 0; i < n; i++) {
        vec[i] = vec[i] / local_density;
      }

      // compute  δₖ = (I - h * J_f(yₑ,ₖ))⁻¹ ⋅ (h * f(yₑ,ₖ) - (yₑ,ₖ - yₛ))
      // and overwrite elements in the `vec` buffer with the result
      // - to that end, we solve following eqn for `x` (and write it to `vec`):
      //      A⋅x = b
      //      | |   |
      //      | |   ---- = `h * f(yₑ,ₖ) - (yₑ,ₖ - yₛ)`  (currently in vec)
      //      | -------- = `δₖ`                        (will be written to vec)
      //      ---------- = `(I - h * J_f(yₑ,ₖ))`       (currently in matrx)
      // - the current algorithm inverts `mtrx` in place, but that's not
      //   important to us
      //
      // todo: consider adjusting gaussj_g's return value so that its more
      //       consistent with GR_SUCCESS
      int ierror = f_wrap::gaussj_g(n, mtrx.data(), vec);
      if (ierror == 1) {
        return GR_FAIL;
      }

      // multiply with density again  (TODO: shouldn't be part of this fn)
      for (int i = 0; i < n; i++) {
        vec[i] = vec[i] * local_density;
      }

      // update values in ycur_minus_ystart and y for next iteration using:
      //  (yₑ,ₖ₊₁ - yₛ) = (yₑ,ₖ - yₛ) + δₖ
      //   yₑ,ₖ₊₁ = yₑ,ₖ + δₖ
      for (int i = 0; i < n; i++) {
        ycur_minus_ystart[i] = ycur_minus_ystart[i] + vec[i];
        y[i] = y[i] + vec[i];
      }

      if (enforce_positive_non_NaN) {
        for (int i = 0; i < n; i++) {
          if (std::isnan(y[i]) || (y[i] <= 0.)) {
            return GR_FAIL;
          }
        }
      }

      double max_rel_diff_mag = 0.0;
      for (int i = 0; i < n; i++) {
        double cur = y[i];
        // todo: double check that our behavior, when cur ~ 0, makes sense
        // reminder: at this point, vec holds δₖ = yₑ,ₖ₊₁ - yₑ,ₖ
        double err = (cur > tiny8) ? std::fabs(vec[i] / cur) : 0.0;
        max_rel_diff_mag = std::fmax(err, max_rel_diff_mag);
      }

      if (max_rel_diff_mag <= max_rtol_) {
        return GR_SUCCESS;
      }
    }

    return GR_FAIL;  // only reached if iterations exceeded maxiter
  }
};

}  // namespace integrate
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_INTEGRATE_HPP