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
#include "fortran_func_wrappers.hpp"
#include "utils-cpp.hpp"

#include "cool1d_cloudy_g.hpp"

void grackle::impl::cool1d_cloudy_g(
    const double* rhoH, const double* metallicity, const double* logtem,
    double* edot, double comp2, double dom, double zr, int icmbTfloor,
    int iClHeat, int iZscale, const gr_mask_type* itmask,
    cloudy_data cloudy_table, IndexRange idx_range) {
  // Locals

  int i, get_heat;
  long long zindex, zmidpt, zhighpt;
  double inv_log10, log10_tCMB;
  double dclPar[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION] = {};
  long long end_int;

  // Slice locals

  std::vector<double> log_n_h(idx_range.i_stop);
  std::vector<double> log_cool(idx_range.i_stop);
  std::vector<double> log_cool_cmb(idx_range.i_stop);
  std::vector<double> log_heat(idx_range.i_stop);
  std::vector<double> edot_met(idx_range.i_stop);
  std::vector<double> log10tem(idx_range.i_stop);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  end_int = 0;
  get_heat = iClHeat;

  inv_log10 = 1. / std::log(10.);
  log10_tCMB = std::log10(comp2);

  // Calculate parameter value slopes

  dclPar[0] =
      (cloudy_table.grid_parameters[0][cloudy_table.grid_dimension[0] - 1] -
       cloudy_table.grid_parameters[0][0]) /
      (double)(cloudy_table.grid_dimension[0] - 1);
  if (cloudy_table.grid_rank > 1) {
    dclPar[1] =
        (cloudy_table.grid_parameters[1][cloudy_table.grid_dimension[1] - 1] -
         cloudy_table.grid_parameters[1][0]) /
        (double)(cloudy_table.grid_dimension[1] - 1);
  }
  if (cloudy_table.grid_rank > 2) {
    dclPar[2] =
        (cloudy_table.grid_parameters[2][cloudy_table.grid_dimension[2] - 1] -
         cloudy_table.grid_parameters[2][0]) /
        (double)(cloudy_table.grid_dimension[2] - 1);
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      log10tem[i] = logtem[i] * inv_log10;

      // Calculate proper log(n_H)

      log_n_h[i] = std::log10(rhoH[i] * dom);

      // Calculate index for redshift dimension

      if (cloudy_table.grid_rank > 2) {
        // Get index for redshift dimension via bisection

        if (zr <= cloudy_table.grid_parameters[1][0]) {
          zindex = 1;
        } else if (zr >=
                   cloudy_table
                       .grid_parameters[1][cloudy_table.grid_dimension[1] - 1 -
                                           1]) {
          zindex = cloudy_table.grid_dimension[1];
          end_int = 1;
          get_heat = 0;
        } else if (zr >=
                   cloudy_table
                       .grid_parameters[1][cloudy_table.grid_dimension[1] - 2 -
                                           1]) {
          zindex = cloudy_table.grid_dimension[1] - 2;
        } else {
          zindex = 1;
          zhighpt = cloudy_table.grid_dimension[1] - 2;
          while ((zhighpt - zindex) > 1) {
            zmidpt = int((zhighpt + zindex) / 2);
            if (zr >= cloudy_table.grid_parameters[1][zmidpt - 1]) {
              zindex = zmidpt;
            } else {
              zhighpt = zmidpt;
            }
          }
        }
      }

      // Call interpolation functions to get heating/cooling

      // Interpolate over temperature.
      if (cloudy_table.grid_rank == 1) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
            log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (get_heat == 1) {
          log_cool[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density and temperature.
      } else if (cloudy_table.grid_rank == 2) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
            log_n_h[i], log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0],
            cloudy_table.grid_parameters[1], dclPar[1], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log_n_h[i], log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (get_heat == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log_n_h[i], log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density, redshift, and temperature.
      } else if (cloudy_table.grid_rank == 3) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_3dz_g(
            log_n_h[i], zr, log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0],
            cloudy_table.grid_parameters[1], zindex,
            cloudy_table.grid_parameters[2], dclPar[2], cloudy_table.data_size,
            cloudy_table.cooling_data, end_int);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((icmbTfloor == 1) && ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_3dz_g(
              log_n_h[i], zr, log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], zindex,
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.cooling_data, end_int);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (get_heat == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_3dz_g(
              log_n_h[i], zr, log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], zindex,
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.heating_data, end_int);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

      } else {
        OMP_PRAGMA_CRITICAL {
          printf("Maximum cooling data grid rank is 3!\n");
        }
        return;
      }

      // Scale cooling by metallicity.

      if (iZscale == 1) {
        edot_met[i] = edot_met[i] * metallicity[i];
      }

      edot[i] = edot[i] + (edot_met[i] * rhoH[i] * rhoH[i]);
    }
  }

  return;
}
