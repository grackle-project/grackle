//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_kappa_grain function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_kappa_grain function from FORTRAN to C++

#include <vector>

#include "interpolate.hpp"
#include "support/View.hpp"
#include "passive/analytic_opac.hpp"

#include "calc_kappa_grain.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void calc_kappa_grain(const double* tdust, double* kgr,
                      const gr_mask_type* itmask, int in, IndexRange idx_range,
                      double t_subl, int gr_N, int gr_Size, double gr_dT,
                      const double* gr_Td, const double* logalsp_data_,
                      int idspecies) {
  if (idspecies == 0) {
    AnalyticOpacCalc calculator(t_subl);
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        kgr[i] = calculator.calc_opac(tdust[i]);
      }
    }
  } else {
    // Opacity table
    FortranView<const double**> logalsp(logalsp_data_, gr_N, in);
    std::vector<double> logalsp1(gr_Size);
    long long gr_N_i64 = static_cast<long long>(gr_N);

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        for (int j = 0; j < gr_Size; j++) {
          logalsp1[j] = logalsp(j, i);
        }
        double log10tdust = std::log10(tdust[i]);

        double logkgr = interpolate_1d(log10tdust, &gr_N_i64, gr_Td, gr_dT,
                                       gr_N_i64, logalsp1.data());

        kgr[i] = std::pow(10., logkgr);
      }
    }
  }
  return;
}

}  // namespace GRIMPL_NAMESPACE_DECL