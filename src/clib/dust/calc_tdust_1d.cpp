//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_tdust_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_tdust_1d_g function from FORTRAN to C++

#include <vector>

#include "dust/multi_grain_species/opac_calculator.hpp"
#include "dust/passive/analytic_opac.hpp"
#include "fortran_func_decls.h"
#include "support/config.hpp"

#include "calc_tdust_1d.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void calc_tdust_1d(double* tdust, const double* tgas, const double* nh,
                   const double* gasgr, const double* gamma_isrfa,
                   const double* isrf, const gr_mask_type* itmask, double trad,
                   int buf_len, int gr_N, double gr_dT, const double* gr_Td,
                   const double* alsp_data_, double* kgr,
                   bool analytic_opac_formula, IndexRange idx_range) {
  if (analytic_opac_formula) {
    AnalyticOpacCalc calculator(Tdust_detail::passive_dust_model_T_sublimation);
    calc_tdust_1d_(tdust, tgas, nh, gasgr, gamma_isrfa, isrf, itmask, trad,
                   buf_len, kgr, idx_range, calculator);
  } else {
    FortranView<const double**> alsp(alsp_data_, gr_N, buf_len);
    std::vector<double> logalsp_data_(gr_N * buf_len);
    FortranView<double**> logalsp(logalsp_data_.data(), gr_N, buf_len);

    int Td_N = gr_N;
    int Td_Size = gr_N;
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        for (int j = 0; j < gr_N; j++) {
          logalsp(j, i) = std::log10(alsp(j, i));
        }
      }
    }

    MultiGrainGrowthOpacCalc calculator(buf_len, Td_N, Td_Size, gr_dT, gr_Td,
                                        logalsp_data_.data());
    calc_tdust_1d_(tdust, tgas, nh, gasgr, gamma_isrfa, isrf, itmask, trad,
                   buf_len, kgr, idx_range, calculator);
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL