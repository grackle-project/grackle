//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_temp_cloudy_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_temp_cloudy_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "index_helper.h"
#include "scale_fields.hpp"
#include "utils-cpp.hpp"

#include "calc_temp_cloudy_g.h"
#include "calc_temp1d_cloudy_g.hpp"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void calc_temp_cloudy_g(gr_float* temperature_data_, int imetal,
                        chemistry_data* my_chemistry,
                        cloudy_data cloudy_primordial,
                        grackle_field_data* my_fields,
                        InternalGrUnits internalu) {

  // Calc quantities using values specified by internalu
  const double dom = internalu_calc_dom_(internalu);
  const double zr = 1. / (internalu.a_value * internalu.a_units) - 1.;

  // Convert densities from comoving to proper

  if (internalu.extfields_in_comoving == 1) {
    double factor = std::pow(internalu.a_value, -3);
    grackle::impl::scale_fields_table(my_fields, factor);
  }

  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  OMP_PRAGMA("omp parallel") {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

    grackle::impl::View<gr_float***> d(
        my_fields->density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    grackle::impl::View<gr_float***> metal;

    if (imetal == 1) {
      metal = grackle::impl::View<gr_float***>(
          my_fields->metal_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    }

    grackle::impl::View<gr_float***> temperature(
        temperature_data_, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    // these are used to temporarily hold values from each idx_range
    std::vector<double> tgas(my_fields->grid_dimension[0]);
    std::vector<double> rhoH(my_fields->grid_dimension[0]);
    std::vector<double> mmw(my_fields->grid_dimension[0]);
    std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);

    // The following for-loop is a flattened loop over every k,j combination.
    // OpenMP divides this loop between all threads. Within the loop, we
    // complete calculations for the constructed index-range construct
    // (an index range corresponds to an "i-slice")
    OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);

      // Initialize iteration mask to true for all cells and compute the mass
      // density of Hydrogen
      const double f_H = my_chemistry->HydrogenFractionByMass;
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;

        if (imetal == 1) {
          gr_float metal_free_density = (d(i, idx_range.j, idx_range.k) -
                                         metal(i, idx_range.j, idx_range.k));
          rhoH[i] = f_H * metal_free_density;
        } else {
          rhoH[i] = f_H * d(i, idx_range.j, idx_range.k);
        }
      }

      // Calculate temperature and mean molecular weight
      grackle::impl::calc_temp1d_cloudy_g(
          rhoH.data(), tgas.data(), mmw.data(), dom, zr, imetal, itmask.data(),
          my_chemistry, cloudy_primordial, my_fields, internalu, idx_range);

      // Record the computed temperature values in the output array
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        temperature(i, idx_range.j, idx_range.k) = tgas[i];
      }
    }
  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from proper

  if (internalu.extfields_in_comoving == 1) {
    double factor = std::pow(internalu.a_value, 3);
    grackle::impl::scale_fields_table(my_fields, factor);
  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
