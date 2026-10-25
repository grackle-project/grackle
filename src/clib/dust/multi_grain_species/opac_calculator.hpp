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
  double dlog10Tdust_;
  const double* log10Tdust_vals_;
  View<const double**> log10kappa;

public:
  /// @brief Construct a new instance
  ///
  /// @param[in]  in i dimension of 3D fields
  /// @param[in]  gr_N Number of temperature points in the grain opacity table
  /// @param[in]  dlog10Tdust Temperature spacing of the grain opacity table
  /// @param[in]  log10Tdust_vals Temperature values of the grain opacity table
  /// @param[in]  log10kappa_data Grain opacity table data
  MultiGrainGrowthOpacCalc(int in, int gr_N, double dlog10Tdust,
                           const double* log10Tdust_vals,
                           const double* log10kappa_data)
      : gr_N_{static_cast<int64_t>(gr_N)},
        dlog10Tdust_{dlog10Tdust},
        log10Tdust_vals_{log10Tdust_vals} {
    log10kappa = View<const double**>(log10kappa_data, in, gr_N);
  }

  /// @brief Calculate mean plank opacity for the dust grain
  ///
  /// I **THINK** the returned opacity has units of cm^2/g and it's measured
  /// "per unit gas mass" (this contrasts with @ref AnalyticOpac::calc_opac
  /// which returns a value measured "per unit grain mass")
  double calc_opac(double tdust, int i) const noexcept {
    double logkgr =
        interpolate_1d(std::log10(tdust), &gr_N_, log10Tdust_vals_,
                       dlog10Tdust_, gr_N_, log10kappa.contig1d_ptr(i));
    return std::pow(10., logkgr);
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_MULTI_GRAIN_SPECIES_OPAC_CALC_HPP