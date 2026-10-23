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
}  // namespace solve_eqn
}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* MATH_SOLVE_EQN_HPP */
