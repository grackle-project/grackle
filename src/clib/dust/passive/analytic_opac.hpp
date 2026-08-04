//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/defines @ref AnalyticOpacCalc
///
//===----------------------------------------------------------------------===//

#ifndef DUST_PASSIVE_OPAC_CALC_HPP
#define DUST_PASSIVE_OPAC_CALC_HPP

#include "../../support/config.hpp"
#include "../../grackle_macros.h"  // tiny_fortran_val

#include <cmath>

namespace GRIMPL_NAMESPACE_DECL {

/// @brief class for computing opacity
///
/// This class is temporarily constructed to facillitate repeated calculations
/// of the opacity
class AnalyticOpacCalc {
  double t_subl;

  /// grain opacity from Omukai (2000, equation 17) normalized by
  /// the local dust-to-gas ratio, which in this work is 0.934e-2.
  ///
  /// We apply a Td^2 scaling up to 200 K following Dopcke et al. (2011).
  /// However, note that:
  /// - Dopcke et al. (2011) adopt a different normalization (i.e., for kgr1)
  /// - Omukai (2000) says the Td^2 proportionality is only valid
  ///   for Td < 50 K. We use this scaling for Td > 50 K, anyways, because
  ///   Omukai (2000) does not suggest what should be done for Td > 50 K.
  static constexpr double kgr1 = 4.0e-4 / 0.00934;

  /// the coefficient used above 200 K
  static constexpr double kgr200 = 16.0 / 0.00934;

public:
  /// @brief Primary constructor
  ///
  /// @param t_subl Grain sublimation temperature
  explicit AnalyticOpacCalc(double t_subl) : t_subl{t_subl} {}

  /// @brief Calculate mean plank opacity for the dust grain
  ///
  /// the returned opacities have units of cm^2/g and they are measured "per
  /// unit grain mass"
  ///
  /// @par History
  /// This logic is descended from logic originally written by Britton Smith in
  /// September 2011. That logic was subsequently ported to C++ in March 2026
  /// by Christopher Bignamini & Matthew Abruzzo
  double calc_opac(double tdust) const noexcept {
    // Temperature dependence from Dopcke et al. (2011).
    // Normalized to Omukai (2000).
    // See comment above for note about Td dependence for kgr.
    if (tdust < 200.) {
      return kgr1 * std::pow(tdust, 2);
    } else if (tdust < t_subl) {
      return kgr200;
    } else {
      return std::fmax(tiny_fortran_val,
                       (kgr200 * std::pow(tdust / 1.5e3, -12)));
    }
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_PASSIVE_OPAC_CALC_HPP