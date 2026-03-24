//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_kappa_gr_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_kappa_gr_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_wrappers.hpp"
#include "utils-cpp.hpp"

#include "calc_kappa_gr_g.hpp"

void grackle::impl::calc_kappa_gr_g(double* tdust, double* kgr,
                                    const gr_mask_type* itmask, int in,
                                    IndexRange idx_range, const double* t_subl,
                                    int* gr_N, const int* gr_Size, double* gr_dT,
                                    double* gr_Td, gr_float* logalsp_data_,
                                    int idspecies) {
  // Parameters

  // grain opacity from Omukai (2000, equation 17) normalized by
  // the local dust-to-gas ratio, which in this work is 0.934e-2.
  const double kgr1 = 4.0e-4 / 0.00934;

  // We then apply a Td^2 scaling up to 200 K following Dopcke et al. (2011).
  // However, note that Dopcke et al. (2011) adopt a different
  // normalization (i.e., for kgr1).
  // Also note that Omukai (2000) says the Td^2 proportionality is only valid
  // for Td < 50 K. We go with this because Omukai (2000) does not suggest
  // what should be done for Td > 50 K.
  const double kgr200 = 16.0 / 0.00934;

  // Opacity table
  grackle::impl::View<gr_float**> logalsp(logalsp_data_, gr_N[1 - 1], in);

  // Locals

  int i;
  gr_float logkgr;
  std::vector<gr_float> logalsp1((*gr_Size));
  long long gr_N_i64;
  double log10tdust;

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================
  if (idspecies != 0) {
    gr_N_i64 = (long long)(gr_N[0]);
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      if (idspecies == 0) {
        // Temperature dependence from Dopcke et al. (2011).
        // Normalized to Omukai (2000).
        // See comment above for note about Td dependence for kgr.
        if (tdust[i] < 200.) {
          kgr[i] = kgr1 * std::pow(tdust[i], 2);
        } else if (tdust[i] < (*t_subl)) {
          kgr[i] = kgr200;
        } else {
          kgr[i] = std::fmax(tiny_fortran_val,
                             (kgr200 * std::pow((tdust[i] / 1.5e3), (-12))));
        }
      } else {
        for (int j = 0; j < *gr_Size; j++) {
          logalsp1[j] = logalsp(j, i);
        }
        log10tdust = std::log10(tdust[i]);

        logkgr = grackle::impl::fortran_wrapper::interpolate_1d_g(log10tdust, &gr_N_i64, gr_Td,
          *gr_dT, gr_N_i64, logalsp1.data());

        kgr[i] = std::pow(10., logkgr);
      }
    }
  }
  return;
}
