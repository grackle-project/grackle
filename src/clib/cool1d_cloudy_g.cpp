//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the cool1d_cloudy_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool1d_cloudy_g function from FORTRAN to C++

#include <cmath>
#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "cool1d_cloudy_g.hpp"

void grackle::impl::cool1d_cloudy_g(
    double* rhoH, double* metallicity, double* logtem, double* edot,
    double comp2, double* dom, double* zr, int* icmbTfloor, int* iClHeat,
    int* iZscale, long long* clGridRank, long long* clGridDim, double* clPar1,
    double* clPar2, double* clPar3, long long* clDataSize, double* clCooling,
    double* clHeating, gr_mask_type* itmask, grackle_field_data* my_fields,
    IndexRange idx_range) {
  // Locals

  int i, get_heat;
  long long zindex, zmidpt, zhighpt;
  double inv_log10, log10_tCMB;
  std::vector<double> dclPar((*clGridRank));
  long long end_int;

  // Slice locals

  std::vector<double> log_n_h(my_fields->grid_dimension[0]);
  std::vector<double> log_cool(my_fields->grid_dimension[0]);
  std::vector<double> log_cool_cmb(my_fields->grid_dimension[0]);
  std::vector<double> log_heat(my_fields->grid_dimension[0]);
  std::vector<double> edot_met(my_fields->grid_dimension[0]);
  std::vector<double> log10tem(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  end_int = 0;
  get_heat = (*iClHeat);

  inv_log10 = 1. / std::log(10.);
  log10_tCMB = std::log10(comp2);

  // Calculate parameter value slopes

  dclPar[1 - 1] = (clPar1[clGridDim[1 - 1] - 1] - clPar1[1 - 1]) /
                  (double)(clGridDim[1 - 1] - 1);
  if ((*clGridRank) > 1) {
    dclPar[2 - 1] = (clPar2[clGridDim[2 - 1] - 1] - clPar2[1 - 1]) /
                    (double)(clGridDim[2 - 1] - 1);
  }
  if ((*clGridRank) > 2) {
    dclPar[3 - 1] = (clPar3[clGridDim[3 - 1] - 1] - clPar3[1 - 1]) /
                    (double)(clGridDim[3 - 1] - 1);
  }

  for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
    if (itmask[i - 1] != MASK_FALSE) {
      log10tem[i - 1] = logtem[i - 1] * inv_log10;

      // Calculate proper log(n_H)

      log_n_h[i - 1] = std::log10(rhoH[i - 1] * (*dom));

      // Calculate index for redshift dimension

      if ((*clGridRank) > 2) {
        // Get index for redshift dimension via bisection

        if ((*zr) <= clPar2[1 - 1]) {
          zindex = 1;
        } else if ((*zr) >= clPar2[clGridDim[2 - 1] - 1 - 1]) {
          zindex = clGridDim[2 - 1];
          end_int = 1;
          get_heat = 0;
        } else if ((*zr) >= clPar2[clGridDim[2 - 1] - 2 - 1]) {
          zindex = clGridDim[2 - 1] - 2;
        } else {
          zindex = 1;
          zhighpt = clGridDim[2 - 1] - 2;
          while ((zhighpt - zindex) > 1) {
            zmidpt = int((zhighpt + zindex) / 2);
            if ((*zr) >= clPar2[zmidpt - 1]) {
              zindex = zmidpt;
            } else {
              zhighpt = zmidpt;
            }
          }
        }
      }

      // Call interpolation functions to get heating/cooling

      // Interpolate over temperature.
      if ((*clGridRank) == 1) {
        FORTRAN_NAME(interpolate_1d_g)(&log10tem[i - 1], clGridDim, clPar1,
                                       &dclPar[1 - 1], clDataSize, clCooling,
                                       &log_cool[i - 1]);
        edot_met[i - 1] = -std::pow(10., log_cool[i - 1]);

        // Ignore CMB term if T >> T_CMB
        if (((*icmbTfloor) == 1) && ((log10tem[i - 1] - log10_tCMB) < 2.)) {
          FORTRAN_NAME(interpolate_1d_g)(&log10_tCMB, clGridDim, clPar1,
                                         &dclPar[1 - 1], clDataSize, clCooling,
                                         &log_cool_cmb[i - 1]);
          edot_met[i - 1] =
              edot_met[i - 1] + std::pow(10., log_cool_cmb[i - 1]);
        }

        if (get_heat == 1) {
          FORTRAN_NAME(interpolate_1d_g)(&log10tem[i - 1], clGridDim, clPar1,
                                         &dclPar[1 - 1], clDataSize, clHeating,
                                         &log_heat[i - 1]);
          edot_met[i - 1] = edot_met[i - 1] + std::pow(10., log_heat[i - 1]);
        }

        // Interpolate over density and temperature.
      } else if ((*clGridRank) == 2) {
        FORTRAN_NAME(interpolate_2d_g)(&log_n_h[i - 1], &log10tem[i - 1],
                                       clGridDim, clPar1, &dclPar[1 - 1],
                                       clPar2, &dclPar[2 - 1], clDataSize,
                                       clCooling, &log_cool[i - 1]);
        edot_met[i - 1] = -std::pow(10., log_cool[i - 1]);

        // Ignore CMB term if T >> T_CMB
        if (((*icmbTfloor) == 1) && ((log10tem[i - 1] - log10_tCMB) < 2.)) {
          FORTRAN_NAME(interpolate_2d_g)(&log_n_h[i - 1], &log10_tCMB,
                                         clGridDim, clPar1, &dclPar[1 - 1],
                                         clPar2, &dclPar[2 - 1], clDataSize,
                                         clCooling, &log_cool_cmb[i - 1]);
          edot_met[i - 1] =
              edot_met[i - 1] + std::pow(10., log_cool_cmb[i - 1]);
        }

        if (get_heat == 1) {
          FORTRAN_NAME(interpolate_2d_g)(&log_n_h[i - 1], &log10tem[i - 1],
                                         clGridDim, clPar1, &dclPar[1 - 1],
                                         clPar2, &dclPar[2 - 1], clDataSize,
                                         clHeating, &log_heat[i - 1]);
          edot_met[i - 1] = edot_met[i - 1] + std::pow(10., log_heat[i - 1]);
        }

        // Interpolate over density, redshift, and temperature.
      } else if ((*clGridRank) == 3) {
        FORTRAN_NAME(interpolate_3dz_g)(
            &log_n_h[i - 1], zr, &log10tem[i - 1], clGridDim, clPar1,
            &dclPar[1 - 1], clPar2, &zindex, clPar3, &dclPar[3 - 1], clDataSize,
            clCooling, &end_int, &log_cool[i - 1]);
        edot_met[i - 1] = -std::pow(10., log_cool[i - 1]);

        // Ignore CMB term if T >> T_CMB
        if (((*icmbTfloor) == 1) && ((log10tem[i - 1] - log10_tCMB) < 2.)) {
          FORTRAN_NAME(interpolate_3dz_g)(
              &log_n_h[i - 1], zr, &log10_tCMB, clGridDim, clPar1,
              &dclPar[1 - 1], clPar2, &zindex, clPar3, &dclPar[3 - 1],
              clDataSize, clCooling, &end_int, &log_cool_cmb[i - 1]);
          edot_met[i - 1] =
              edot_met[i - 1] + std::pow(10., log_cool_cmb[i - 1]);
        }

        if (get_heat == 1) {
          FORTRAN_NAME(interpolate_3dz_g)(
              &log_n_h[i - 1], zr, &log10tem[i - 1], clGridDim, clPar1,
              &dclPar[1 - 1], clPar2, &zindex, clPar3, &dclPar[3 - 1],
              clDataSize, clHeating, &end_int, &log_heat[i - 1]);
          edot_met[i - 1] = edot_met[i - 1] + std::pow(10., log_heat[i - 1]);
        }

      } else {
        OMP_PRAGMA_CRITICAL {
          printf("Maximum cooling data grid rank is 3!\n");
        }
        return;
      }

      // Scale cooling by metallicity.

      if ((*iZscale) == 1) {
        edot_met[i - 1] = edot_met[i - 1] * metallicity[i - 1];
      }

      edot[i - 1] = edot[i - 1] + (edot_met[i - 1] * rhoH[i - 1] * rhoH[i - 1]);
    }
  }

  return;
}
