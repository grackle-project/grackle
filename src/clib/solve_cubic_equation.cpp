//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the solve_cubic_equation function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// solve_cubic_equation function from FORTRAN to C++

#include <cmath>
#include <cstdio>

#include "phys_constants.h"  // pi_fortran_val
#include "solve_cubic_equation.hpp"

int grackle::impl::solve_cubic_equation_cpp(double a, double b, double c,
                                            double& root) {
  double root1, root2, root3;
  double q, r, m;
  double th;
  double s, t;
  const double pi_local_var = pi_fortran_val;

  q = (a * a - 3.e0 * b) / 9.e0;
  r = (2.e0 * a * a * a - 9.e0 * a * b + 27.e0 * c) / 54.e0;
  m = r * r - q * q * q;

  if (m < 0.e0)  // three real roots
  {
    th = acos(r / std::sqrt(q * q * q));
    root1 = -(2.e0 * std::sqrt(q) * std::cos(th / 3.e0)) - a / 3.e0;
    root2 =
        -(2.e0 * std::sqrt(q) * std::cos((th + 2.e0 * pi) / 3.e0)) - a / 3.e0;
    root3 =
        -(2.e0 * std::sqrt(q) * std::cos((th - 2.e0 * pi) / 3.e0)) - a / 3.e0;
    fprintf(stderr, "three real roots %g %g %g\n", root1, root2, root3);
    return 1;
  } else {  // one real root
    if (r > 0.e0) {
      s = -std::pow((r + std::sqrt(m)), (1.e0 / 3.e0));
    } else {
      s = std::pow((-r + std::sqrt(m)), (1.e0 / 3.e0));
    }
    t = q / s;
    root = s + t - a / 3.e0;
  }

  return 0;
}
