//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares @ref OpacWithDeriv
///
//===----------------------------------------------------------------------===//

#ifndef DUST_OPAC_WITH_DERIV_HPP
#define DUST_OPAC_WITH_DERIV_HPP

#include "../support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Holds an opacity and its derivative with respect to dust temperature
///
/// This is returned by the ``calc_opac_with_deriv`` method of the opacity
/// calculators
struct OpacWithDeriv {
  // Opacity value
  double kgr;
  // Derivative of opacity with respect to dust temperature
  double dkgr_dT;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_OPAC_WITH_DERIV_HPP
