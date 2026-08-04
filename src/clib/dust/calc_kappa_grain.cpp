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
#include "multi_grain_species/opac_calculator.hpp"
#include "passive/analytic_opac.hpp"
#include "support/View.hpp"

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
        kgr[i] = calculator.calc_opac(tdust[i], i);
      }
    }
  } else {
    MultiGrainGrowthOpacCalc calculator(in, gr_N, gr_Size, gr_dT, gr_Td,
                                        logalsp_data_);
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        kgr[i] = calculator.calc_opac(tdust[i], i);
      }
    }
  }
  return;
}

}  // namespace GRIMPL_NAMESPACE_DECL