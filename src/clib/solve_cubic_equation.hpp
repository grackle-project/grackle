//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the solve_cubic_equation function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// solve_cubic_equation function from FORTRAN to C++

#ifndef SOLVE_CUBIC_EQUATION_HPP
#define SOLVE_CUBIC_EQUATION_HPP

namespace grackle::impl {

///  Find the real root of a cubic equation of the form
///
/// @param[in]  a coefficient of the quadratic term
/// @param[in]  b coefficient of the linear term
/// @param[in]  c constant term
/// @param[out] root the real root of the cubic equation
/// @return     0 on success, 1 if there are three real roots
///
/// @par History
/// modified: January, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
int solve_cubic_equation_cpp(
  double a, double b, double c, double& root
);

} // namespace grackle::impl

#endif /* SOLVE_CUBIC_EQUATION_HPP */
