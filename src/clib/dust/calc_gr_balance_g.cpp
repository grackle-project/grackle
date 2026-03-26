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

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "phys_constants.h"
#include "utils-cpp.hpp"

#include "calc_gr_balance_g.hpp"

void grackle::impl::calc_gr_balance_g(double* tdust, const double* tgas, const double* kgr,
                                      double trad4, const double* gasgr,
                                      const double* gamma_isrf, const double* nh,
                                      const gr_mask_type* itmask, double* sol,
                                      IndexRange idx_range) {
  // Parameters

  const double radf = 4. * sigma_sb_grflt;

  // Locals

  int i;

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      sol[i] = gamma_isrf[i] + radf * kgr[i] * (trad4 - std::pow(tdust[i], 4)) +
               (gasgr[i] * nh[i] * (tgas[i] - tdust[i]));
      // Historically, the following comment was present here:
      //     emission/absorption rate per unit grain mass [erg/s/g]
      //     for Z = Zsun (default)
      // This comment is **ONLY** correct when the function used as
      // part of the single-field dust model. See the docstring for
      // more details.
    }
  }

  return;
}
