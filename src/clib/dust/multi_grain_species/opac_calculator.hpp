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

#include <memory>  // std::unique_ptr, std::make_unique

namespace GRIMPL_NAMESPACE_DECL {

/// @brief class for computing opacity
///
/// This class is temporarily constructed to facillitate repeated calculations
/// of the opacity
///
/// @note
/// The future goal is to get rid of intermediate buffers so that we can
/// directly compute opacity from the injection pathway data
class MultiGrainGrowthOpacCalc {
  long long gr_N_;
  int gr_Size_;
  double dTdust_;
  const double* Tdust_vals_;
  FortranView<const double**> logalsp;
  std::unique_ptr<double[]> logalsp_tmp;

public:
  /// @brief Construct a new instance
  ///
  /// @param[in]  in i dimension of 3D fields
  /// @param[in]  gr_N Number of temperature points in the grain opacity table
  /// @param[in]  gr_Size Number of temperature points in the grain opacity
  ///     table
  /// @param[in]  gr_dT Temperature spacing of the grain opacity table
  /// @param[in]  gr_Td Temperature values of the grain opacity table
  /// @param[in]  logalsp_data Grain opacity table data
  MultiGrainGrowthOpacCalc(int in, int gr_N, int gr_Size, double gr_dT,
                           const double* gr_Td, const double* logalsp_data)
      : gr_N_{static_cast<long long>(gr_N)},
        gr_Size_{gr_Size},
        dTdust_{gr_dT},
        Tdust_vals_{gr_Td} {
    logalsp = FortranView<const double**>(logalsp_data, gr_N, in);
    logalsp_tmp = std::make_unique<double[]>(gr_Size);
  }

  /// @brief Calculate mean plank opacity for the dust grain
  ///
  /// I **THINK** the returned opacity has units of cm^2/g and it's measured
  /// "per unit gas mass" (this contrasts with @ref AnalyticOpac::calc_opac
  /// which returns a value measured "per unit grain mass")
  double calc_opac(double tdust, int i) const noexcept {
    for (int j = 0; j < gr_Size_; j++) {
      logalsp_tmp[j] = logalsp(j, i);
    }
    double log10tdust = std::log10(tdust);

    double logkgr = interpolate_1d(log10tdust, &gr_N_, Tdust_vals_, dTdust_,
                                   gr_N_, logalsp_tmp.get());

    return std::pow(10., logkgr);
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_MULTI_GRAIN_SPECIES_OPAC_CALC_HPP