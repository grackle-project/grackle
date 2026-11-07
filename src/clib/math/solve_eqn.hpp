//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares logic fo solving equations
///
//===----------------------------------------------------------------------===//

#ifndef MATH_SOLVE_EQN_HPP
#define MATH_SOLVE_EQN_HPP

#include <cmath>

#include "./misc.hpp"
#include "support/config.hpp"  // GRIMPL_NAMESPACE_DECL

namespace GRIMPL_NAMESPACE_DECL {
namespace solve_eqn {

///  Find the real root of a cubic equation of the form
///   x^3 + a x^2 + b x + c = 0
///
/// @param[in]  a coefficient of the quadratic term
/// @param[in]  b coefficient of the linear term
/// @param[in]  c constant term
/// @param[out] root the real root of the cubic equation
/// @return     0 on success, 1 if there are three real roots
///
/// @par History
/// modified: March, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
int cubic_eqn_root(double a, double b, double c, double& root);

// numerical-solve logic
// ---------------------
// starting below, we define a bunch logic pertaining to numerically solving
// a scalar-valued function. This logic was factored out of calc_tdust_1d (in
// certain respects, it's still somewhat tailored to that use-case) and the
// longer-term goal is to make this functionality usable in other parts of the
// codebase
//
// A few items to reconsider for the future:
// -> the machinery is designed to let you run a pseudo-newton's method first
//    and then run a bisection method afterwards. It would be MUCH better to
//    adopt a hybrid method
//    - if we don't know derivatives, we could adopt a hybrid secant-bisection
//      method (for context, Brent's method also involves inverse quadratic
//      interpolation)
//    - if we do know derivatives, it's very common to modify newton's method
//      to work with bisection. This helps when derivatives are close to 0
// -> the machinery currently uses finite-differences to approximate derivatives
//    in newton's method. In calc_tdust_1d and other places where this logic
//    could be used (namely, calculating gas temperature in tabulated-mode),
//    the analytic function for computing the derivative is straight-forward.
//    The trick is that due to piecewise functions/linear interpolation, the
//    analytic derivative is not continuous. The trick to using the analytic
//    derivatives is to an initial phase of the calculation using
//    bisection/binary to find the interval containing bounding the solution
//    where the analytic derivative is continuous

/// Encodes a mask value specifying the status at a given index of a numerical
/// root-finding algorithm
enum struct SolveStatus {
  UNCONVERGED,  ///< the solver has not reached a converged solution
  CONVERGED,    ///< the solver has successfully converged
  SKIP_SOLVE    ///< the solver will skip the calculation
};

/// @brief aggregates the result of a solve
///
/// The name is a little more generic that just RootResult since we may want
/// to use it with fixed point iteration
struct EqnSolveRslt {
  /// Indicates if all indices initially assocated with a
  /// @ref SolveStatus::UNCONVERGED status have been solved
  bool all_solved;
  /// The number of completed iterations
  int iterations;
};

// this could definitely be better generalized
/// @brief Find root of a function within an interval at a number of locations
///
/// Uses the bisection method to find an array of roots for @p fn (a
/// scalar-valued function) for an array of bracktting intervals.
///
/// @param[in]    fn Function object representing the models the mathematical
///    function for which roots are computed. This should have a signature
///    of the form `double fn(double x, int i)` where `i` refers to an
///    associated array index
/// @param[inout] x_a, x_b Arrays of interval bounds bracketting the root.
///    Elements are updated in place as the bounds around the root are
///    narrowed. Historically, we have called the final value in x_a the
///    root. See the note below for additional requirements.
/// @param[inout] solvemask Array specifying a @ref SolveStatus value at each
///    index. The solver will only look for roots at indices that initially have
///    values of @ref SolveStatus::UNCONVERGED. Elements will be updated to hold
///    @ref SolveStatus::CONVERGED, as convergence is achieved.
/// @param[in]    i_start, i_stop The range of indices for which root finding
///    is performed.
/// @param[in]    rtol Convergence is achieved when at an index `i`, when
///    `fabs((x_b[i] - x_a[i]) / x_a[i]) <= rtol`.
/// @param[in]    max_iter The maximum number of iterations
/// @param[in]    max_initial_guess The maximum value of the very first
///     midpoint.
///
/// @returns An object that describes the solution. You can identify
///     unconverged locations by checking the values of @p solvemask
///
/// @todo
/// We should get rid of the @p max_initial_guess argument. It primarily exists
/// for historical consistency. We could get rid of it by providing better
/// initial x_a and x_b arguments.
///
/// @note
/// For an index `i`, the current implementation requires (but does not
/// explicitly enforce) that `fn(x_a[i], i).f_val > 0` and
/// `fn(x_b[i], i).f_val < 0`. If this condition is not satisfied, then the
/// results are undefined. The fact that these conditons aren't checked is
/// reflected by the name of this function template.
///
/// @warning
/// The current implementation for comparing to `rtol` introduces an implicit
/// assumption that `x_a` and `x_b` both initially hold positive values.
template <typename Fn>
static EqnSolveRslt unchecked_bisect(
    const Fn& fn, double* x_a, double* x_b, SolveStatus* solvemask, int i_start,
    int i_stop, double rtol, int max_iter,
    double max_initial_guess = std::numeric_limits<double>::max()) {
  int n_unconverged = 0;
  for (int i = i_start; i < i_stop; i++) {
    n_unconverged += (solvemask[i] == SolveStatus::UNCONVERGED);
  }

  // implicitly assumption that
  //   - fn(x_a[i], i).f_val > 0
  //   - fn(x_b[i], i).f_val < 0
  // TODO: we should probably check this assumption (especially for
  //       multi-species dust grains)

  int iter;
  for (iter = 0; n_unconverged > 0 && iter < max_iter; iter++) {
    for (int i = i_start; i < i_stop; i++) {
      if (solvemask[i] == SolveStatus::UNCONVERGED) {
        double x_mid = 0.5 * (x_a[i] + x_b[i]);
        if (iter == 0) {
          x_mid = std::fmin(x_mid, max_initial_guess);
        }

        double f_val = fn(x_mid, i);

        // TODO: consider implementing common bisection strategy that tracks
        //       x_a and current bracket width. This is advantageous because
        //       we can eliminate the conditional update of x_b (we know that
        //       bracket width is always cut in half each iteration)
        //
        // if x_a, x_b, or x_mid isn't finite, branchless choice will act
        // weird (but we will be pretty doomed in that scenario, anyway)
        bool update_a = f_val > 0.0;
        x_a[i] = update_a ? x_mid : x_a[i];
        x_b[i] = update_a ? x_b[i] : x_mid;

        if (std::fabs(x_b[i] - x_a[i]) <= std::fabs(x_a[i]) * rtol) {
          solvemask[i] = SolveStatus::CONVERGED;
          n_unconverged--;
        }
      }
    }
  }
  return EqnSolveRslt{n_unconverged == 0, iter + 1};
}

/// @brief Parameters for the positive-domain, finite difference newton method
struct FDiffNewtonParams {
  /// maximum number of iterations
  int maxiter;

  /// @brief the initial perturbation for computing finite differences
  ///
  /// A finite difference involves evaluating a function at a value ``x`` and
  /// a value ``(1+perturb)*x``. This is the initial perturbation. The size of
  /// the perturbation gets smaller with subsequent iterations.
  double perturb_init;

  /// @brief the minimum perturbation for computing finite differences
  ///
  /// @note
  /// If the relative difference between 2 successive guesses exceeds this
  /// value, we should probably consider the result converged
  double perturb_min;

  /// @brief convergence tolerance
  ///
  /// This is a somewhat atypical definition for tolerance. A guess ``x`` is
  /// considered converged if ``|f(x)| < |tol*f(x*(1+perturb))|``
  double tol;
};

/// @brief Find root of a function within an interval at a number of locations
///
/// Use a variant of newton's method to find an array of roots for @p fn (a
/// scalar-valued function that does NOT directly compute the derivative).
/// Derivatives are estimated via finite differences.
///
/// This function assumes that the root of @p fn occurs for a positive x value.
/// It is also designed to help the caller fall back to an alternative method
/// (e.g. bisection) upon failure.
///
/// @param[in]    fn Function object representing the models the mathematical
///    function for which roots are computed. This should have a signature
///    of the form `double fn(double x, int i)` where `i` refers to an
///    associated array index
/// @param[inout] x Array of initial guesses for the root. Elements are updated
///    in place as the root is searched for.
///    narrowed. Historically, we have called the final value in x_a the
///    root. See the note below for additional requirements.
/// @param[inout] solvemask Array specifying a @ref SolveStatus value at each
///    index. The solver will only look for roots at indices that initially have
///    values of @ref SolveStatus::UNCONVERGED. Elements will be updated to hold
///    @ref SolveStatus::CONVERGED, as convergence is achieved or to hold
///    @ref SolveStatus::SKIP_SOLVE at locations where the guess drops below
///    @p giveup_small_x_threshold
/// @param[inout] f_vals, pert Scratch arrays employed by the function.
/// @param[in]    i_start, i_stop The range of indices for which root finding
///    is performed.
/// @param[in]    param Carries a number of assorted configuration options
/// @param[in]    giveup_small_x_threshold A positive value. Below this
///    threshold, the solver gives up
/// @param[in]    max_x The maximum x value in the domain that can be guessed
///    as a root. Be mindful that for a guessed root of ``max_x``, a derivative
///    will be estimated by evaluating the function with an x value that may be
///    as large as ``max_x * (1.0 + perturb_init)``, where ``perturb_init`` is
///    specified in the @p param argument
///
/// @note
/// For an index `i`, the current implementation requires (but does not
/// explicitly enforce) that ``giveup_small_x_threshold <= x[i]`` and
/// ``x[i] <= max_x``. If this condition is not satisfied, then the
/// results are undefined. The fact that these conditons aren't checked is
/// reflected by the name of this function template.
template <typename Fn>
static EqnSolveRslt unchecked_posdomain_fdiff_newton(
    const Fn& fn, double* x, SolveStatus* solvemask, double* f_vals,
    double* pert, int i_start, int i_stop, const FDiffNewtonParams& param,
    double giveup_small_x_threshold, double max_x) {
  int n_to_solve = 0;

  for (int i = i_start; i < i_stop; i++) {
    n_to_solve += solvemask[i] == SolveStatus::UNCONVERGED;
    pert[i] = param.perturb_init;
  }
  // Iterate to convergence with Newton's method
  bool any_giveups = false;
  int iter;
  for (iter = 0; (iter < param.maxiter) && (n_to_solve > 0); iter++) {
    for (int i = i_start; i < i_stop; i++) {
      if (solvemask[i] == SolveStatus::UNCONVERGED) {
        f_vals[i] = fn(x[i], i);
      }
    }
    for (int i = i_start; i < i_stop; i++) {
      if (solvemask[i] == SolveStatus::UNCONVERGED) {
        // (outside of the very first iteration), the combination of using
        // giveup_small_x_threshold and applying fmax to xplus doesn't make a
        // lot of sense. Given that pert is guaranteed to be positive:
        // - when giveup_small_x_threshold > 1.0e-3, this operation is totally
        //   unnecessary
        // - 1.0e-3 > when giveup_small_x_threshold > 0, then the value of
        //   (x_plus-x[i])/x[i] may get huge!
        // thus, instead of applying fmax to x_plus, it would be more robust
        // to apply fmax to giveup_small_x_threshold and add an initial check
        // to comparing x and giveup_small_x_threshold
        double x_plus = std::fmax(1.e-3, (1. + pert[i]) * x[i]);
        double fplus_val = fn(x_plus, i);
        // it would be slightly faster to compute 1.0/slope
        double slope = (fplus_val - f_vals[i]) / (pert[i] * x[i]);

        double x_old = x[i];
        // note: next line is erroneous if fplus_val - f_vals is exactly 0
        x[i] = std::fmin(x[i] - (f_vals[i] / slope), max_x);

        // try to ensure next x_plus is closer to x than x_old
        pert[i] = fclamp(0.5 * std::fabs(x[i] - x_old) / x[i],
                         param.perturb_min, pert[i]);

        if (x[i] < giveup_small_x_threshold) {
          solvemask[i] = SolveStatus::SKIP_SOLVE;
          n_to_solve--;
          any_giveups = true;
        } else if (std::fabs(f_vals[i]) < std::fabs(fplus_val * param.tol)) {
          solvemask[i] = SolveStatus::CONVERGED;
          n_to_solve--;
        }
      }
    }
  }
  return {(n_to_solve == 0) && !any_giveups, iter + 1};
}
}  // namespace solve_eqn
}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* MATH_SOLVE_EQN_HPP */
