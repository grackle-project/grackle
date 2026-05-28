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

#include "../fortran_func_decls.h"
#include "../fortran_func_wrappers.hpp"
#include "grackle.h"
#include "../utils-cpp.hpp"
#include "./common.hpp"

#include "calc_temp1d_cloudy.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void calc_temp1d_cloudy(const double* rhoH, double* tgas, double* mmw,
                        double dom, double zr, int imetal,
                        const gr_mask_type* itmask,
                        const chemistry_data* my_chemistry,
                        cloudy_data cloudy_table,
                        const grackle_field_data* my_fields,
                        InternalGrUnits internalu, IndexRange idx_range) {
  // General Arguments

  View<const gr_float***> d(my_fields->density, my_fields->grid_dimension[0],
                            my_fields->grid_dimension[1],
                            my_fields->grid_dimension[2]);
  View<const gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  View<const gr_float***> e(
      my_fields->internal_energy, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  std::vector<double> logtem(my_fields->grid_dimension[0]);

  // Parameters

  // approx. mean molecular weight of metals
  const double mu_metal = 16.;
  const int ti_max = 20;
  // todo: maybe factor this constant out and have it be precomputed
  //       (and rename it to INV_LN10). Alternatively, if we directly compute
  //       log10(T) we can drop this entirely
  const double inv_log10 = 1. / std::log(10.);

  // Slice locals
  std::vector<double> log10tem(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Calculate parameter value slopes
  const std::array<double, tabulated_detail::MAX_RANK> dclPar =
      tabulated_detail::param_deltas(cloudy_table);

  // Calculate index for redshift dimension
  const tabulated_detail::FindZIndexRslt zindex_pair =
      tabulated_detail::find_zindex(zr, cloudy_table);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Calculate proper log10(n_H)
      // -> we **may** want to precompute log(rhoH) and log(dom) since it seems
      //    like log(rhoH) is used in a fair number of places
      double log10_nH = std::log10(rhoH[i] * dom);

      // begin the iterative solve for tgas and munew
      double munew = 1.;
      double muold;
      bool skip_mmw_update = true;
      for (int ti = 1; ti <= ti_max; ti++) {
        muold = munew;

        tgas[i] = std::fmax((my_chemistry->Gamma - 1.) *
                                e(i, idx_range.j, idx_range.k) * munew *
                                internalu.utem,
                            my_chemistry->TemperatureStart);
        logtem[i] = std::log(tgas[i]);

        log10tem[i] = logtem[i] * inv_log10;

        // Call interpolation functions to get mmw

        // Interpolate over temperature.
        if (cloudy_table.grid_rank == 1) {
          munew = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.mmw_data);

          // Interpolate over density and temperature.
        } else if (cloudy_table.grid_rank == 2) {
          munew = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log10_nH, log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.mmw_data);

          // Interpolate over density, redshift, and temperature.
        } else if (cloudy_table.grid_rank == 3) {
          munew = grackle::impl::fortran_wrapper::interpolate_3dz_g(
              log10_nH, zr, log10tem[i],
              cloudy_table.grid_dimension,  // 3 elements
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], zindex_pair.zindex,
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.mmw_data,
              zindex_pair.end_int);
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
          skip_mmw_update = false;
          break;
        }
      }

      if (skip_mmw_update) {
        mmw[i] = munew;
        printf("Mean molecular weight not converged! %e %e %e\n", munew, muold,
               std::fabs((munew / muold) - 1.));
      }
    }
  }

  return;
}

}  // namespace GRIMPL_NAMESPACE_DECL