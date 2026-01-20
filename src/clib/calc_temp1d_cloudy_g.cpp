//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the the calc_temp1d_cloudy_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_temp1d_cloudy_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "grackle.h"
#include "utils-cpp.hpp"

#include "calc_temp1d_cloudy_g.hpp"

void grackle::impl::calc_temp1d_cloudy_g(
    const double* rhoH, double* tgas, double* mmw, double dom, double zr,
    int imetal, const gr_mask_type* itmask, chemistry_data* my_chemistry,
    cloudy_data cloudy_table, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range) {
  // General Arguments

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(
      my_fields->internal_energy, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  std::vector<double> logtem(my_fields->grid_dimension[0]);

  // Parameters

  // approx. mean molecular weight of metals
  const double mu_metal = 16.;
  const int ti_max = 20;

  // Locals

  int i, ti;
  long long zindex, zmidpt, zhighpt;
  double inv_log10, muold, munew;
  double dclPar[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
  long long end_int;

  // Slice locals

  std::vector<double> log_n_h(my_fields->grid_dimension[0]);
  std::vector<double> log10tem(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  end_int = 0;

  inv_log10 = 1. / std::log(10.);

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

  // Calculate index for redshift dimension
  zindex = 1;
  if (cloudy_table.grid_rank > 2) {
    // Get index for redshift dimension via bisection

    if (zr <= cloudy_table.grid_parameters[1][0]) {
      zindex = 1;
    } else if (zr >=
               cloudy_table
                   .grid_parameters[1][cloudy_table.grid_dimension[1] - 2]) {
      zindex = cloudy_table.grid_dimension[1];
      end_int = 1;
    } else if (zr >=
               cloudy_table
                   .grid_parameters[1][cloudy_table.grid_dimension[1] - 3]) {
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

  for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Calculate proper log(n_H)
      log_n_h[i] = std::log10(rhoH[i] * dom);
    }
  }

  for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
    if (itmask[i] != MASK_FALSE) {
      munew = 1.;
      bool skip_mmw_update = false;
      for (ti = 1; ti <= (ti_max); ti++) {
        muold = munew;

        tgas[i] = std::fmax((my_chemistry->Gamma - 1.) *
                                e(i, idx_range.j, idx_range.k) * munew *
                                internalu.utem,
                            my_chemistry->TemperatureStart);
        logtem[i] = std::log(tgas[i]);

        log10tem[i] = logtem[i] * inv_log10;

        // Call interpolation functions to get heating/cooling

        // Interpolate over temperature.
        if (cloudy_table.grid_rank == 1) {
          munew = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.mmw_data);

          // Interpolate over density and temperature.
        } else if (cloudy_table.grid_rank == 2) {
          munew = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log_n_h[i], log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.mmw_data);

          // Interpolate over density, redshift, and temperature.
        } else if (cloudy_table.grid_rank == 3) {
          munew = grackle::impl::fortran_wrapper::interpolate_3dz_g(
              log_n_h[i], zr, log10tem[i],
              cloudy_table.grid_dimension,  // 3 elements
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], zindex,
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.mmw_data, end_int);
        } else {
          printf("Maximum mmw data grid rank is 3!\n");
          return;
        }

        munew = 0.5f * (munew + muold);
        tgas[i] = tgas[i] * munew / muold;

        if (std::fabs((munew / muold) - 1.) <= 1.e-2) {
          muold = munew;

          // Add metal species to mean molecular weight

          if (imetal == 1) {
            munew = d(i, idx_range.j, idx_range.k) /
                    ((d(i, idx_range.j, idx_range.k) -
                      metal(i, idx_range.j, idx_range.k)) /
                         munew +
                     metal(i, idx_range.j, idx_range.k) / mu_metal);
            tgas[i] = tgas[i] * munew / muold;
          }

          mmw[i] = munew;
          skip_mmw_update = true;
          break;
        }
      }

      if (!skip_mmw_update) {
        mmw[i] = munew;
        printf("Mean molecular weight not converged! %e %e %e\n", munew,
                muold, std::fabs((munew / muold) - 1.));
      }
    }
  }

  return;
}
