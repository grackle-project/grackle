//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/defines @ref MultiGrainGrowthOpacCalc
///
//===----------------------------------------------------------------------===//

#ifndef DUST_MULTI_GRAIN_SPECIES_OPAC_CALC_HPP
#define DUST_MULTI_GRAIN_SPECIES_OPAC_CALC_HPP

#include "../../interpolate.hpp"
#include "../../support/config.hpp"
#include "../../support/View.hpp"

#include <cstdint>

namespace GRIMPL_NAMESPACE_DECL {

/// @brief class for computing opacity
///
/// This class is temporarily constructed to facillitate repeated calculations
/// of the opacity
class MultiGrainGrowthOpacCalc {
  int64_t gr_N_;
  double dTdust_;
  const double* Tdust_vals_;
  View<const double**> logalsp;

public:
  /// @brief Construct a new instance
  ///
  /// @param[in]  in i dimension of 3D fields
  /// @param[in]  gr_N Number of temperature points in the grain opacity table
  /// @param[in]  gr_dT Temperature spacing of the grain opacity table
  /// @param[in]  gr_Td Temperature values of the grain opacity table
  /// @param[in]  logalsp_data Grain opacity table data
  MultiGrainGrowthOpacCalc(int in, int gr_N, double gr_dT, const double* gr_Td,
                           const double* logalsp_data)
      : gr_N_{static_cast<int64_t>(gr_N)}, dTdust_{gr_dT}, Tdust_vals_{gr_Td} {
    logalsp = View<const double**>(logalsp_data, in, gr_N);
  }

  /// @brief Calculate mean plank opacity for the dust grain
  ///
  /// I **THINK** the returned opacity has units of cm^2/g and it's measured
  /// "per unit gas mass" (this contrasts with @ref AnalyticOpac::calc_opac
  /// which returns a value measured "per unit grain mass")
  double calc_opac(double tdust, int i) const noexcept {
    double logkgr = interpolate_1d(std::log10(tdust), &gr_N_, Tdust_vals_,
                                   dTdust_, gr_N_, logalsp.contig1d_ptr(i));
    return std::pow(10., logkgr);
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_MULTI_GRAIN_SPECIES_OPAC_CALC_HPP