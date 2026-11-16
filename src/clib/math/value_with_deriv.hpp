//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares @ref ValueWithDeriv
///
//===----------------------------------------------------------------------===//

#ifndef MATH_VALUE_WITH_DERIV_HPP
#define MATH_VALUE_WITH_DERIV_HPP

#include "support/config.hpp"  // GRIMPL_NAMESPACE_DECL

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Holds the value of a scalar function and its derivative, both
///     evaluated at the same point
///
/// Callers are encouraged to give the members meaningful names with a
/// structured binding, e.g. ``auto [kgr, dkgr_dT] = calc_opac_with_deriv(...)``
struct ValueWithDeriv {
  /// the function value
  double value;
  /// the derivative of the function with respect to its argument
  double deriv;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_VALUE_WITH_DERIV_HPP
