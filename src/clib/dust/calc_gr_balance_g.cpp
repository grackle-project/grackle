//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_gr_balance_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_gr_balance_g function from FORTRAN to C++

#include "fortran_func_decls.h"

#include "calc_tdust_1d_g.hpp"
#include "calc_gr_balance_g.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void calc_gr_balance_g(const double* tdust, const double* tgas,
                       const double* kgr, double trad4, const double* gasgr,
                       const double* gamma_isrf, const double* nh,
                       const gr_mask_type* itmask, double* sol,
                       IndexRange idx_range) {
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      sol[i] = Tdust_detail::calc_grain_balance(
          tdust[i], tgas[i], kgr[i], trad4, gasgr[i], gamma_isrf[i], nh[i]);
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL