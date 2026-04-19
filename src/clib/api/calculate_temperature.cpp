/***********************************************************************
/
/ Calculate temperature field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cstdio>
#include <vector>

#include "gas_props.hpp"
#include "grackle.h"
#include "index_helper.h"
#include "internal_units.hpp"
#include "scale_fields.hpp"
#include "support/config.hpp"
#include "utils-cpp.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// Set the mean molecular mass for metals
// -> TODO: this should really be defined by a (internal) header
// -> currently, it's also defined by cool1d_multi_g and calc_temp1d_cloudy_g
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0

namespace GRIMPL_NAMESPACE_DECL {

/// Calculate temperature on a 3D grid using mmw from a cloudy table
///
/// @par History
/// written by: Britton Smith May 2015
/// modified1:  February, 2025 by Matthew Abruzzo; ported to C++
///
/// @param[out] temperature_data Array where computed values are written
/// @param[in]  imetal flag if metal field is active (0 = no, 1 = yes)
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  cloudy_primordia specifies the cloudy table
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
static void calc_temp_cloudy(gr_float* temperature_data_, int imetal,
                             chemistry_data* my_chemistry,
                             cloudy_data cloudy_primordial,
                             grackle_field_data* my_fields,
                             InternalGrUnits internalu)
{
    // this assertion is a hint to clang-analyzer about the relationship between
  // `imetal` and whether `metal_density` is a nullptr
  // -> for context, `scale_fields_table` is inlined into this function. Rather
  //    than look at `imetal`, it checks if `metal_density` is a nullptr
  // -> without this assertion clang-tidy will infer that since there's an
  //    explicit check for whether `metal_density` is null, regardless of
  //    `imetal`'s value, then it must be possible for `metal_density` to be
  //    null when `imetal` is 1. Then, it would report an error
  GR_INTERNAL_REQUIRE((imetal != 1) || (my_fields->metal_density != nullptr),
                      "imetal has an incorrect value");

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

      // Initialize iteration mask to true for all cells
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;
      }

      // calculate the basic gas properties (tgas, mmw, rhoH)
      GRIMPL_NS::basic_gas_props(
          tgas.data(), mmw.data(), rhoH.data(), imetal, itmask.data(),
          my_chemistry, &cloudy_primordial, my_fields, internalu, idx_range);

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

}  // namespace GRIMPL_NAMESPACE_DECL
 
extern "C" int local_calculate_temperature(chemistry_data *my_chemistry,
                                           chemistry_data_storage *my_rates,
                                           code_units *my_units,
                                           grackle_field_data *my_fields,
                                           gr_float *temperature)
{
  if (!my_chemistry->use_grackle) { return GR_SUCCESS; }

  const int imetal = (my_fields->metal_density != NULL) ? 1 : 0;

  // we have special handling for tabulated-chemistry-mode
  if (my_chemistry->primordial_chemistry == 0) {
    GRIMPL_NS::calc_temp_cloudy(temperature, imetal, my_chemistry,
                                my_rates->cloudy_primordial, my_fields,
                                GRIMPL_NS::new_internalu_(my_units));
    return GR_SUCCESS;
  };


  /* Compute the pressure first. */
  if (local_calculate_pressure(my_chemistry, my_rates, my_units,
                               my_fields, temperature) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in calculate_pressure.\n");
    return GR_FAIL;
  }

  // Calculate temperature units and fetch some constants

  const double temperature_units = get_temperature_units(my_units);
  const double tiny_number = 1.e-20;
  const double inv_metal_mol = 1.0 / MU_METAL;

  /* Compute properties used to index the field. */
  const grackle_index_helper ind_helper = build_index_helper_(my_fields);

  /* Compute temperature with mu calculated directly. */

  /* parallelize the k and j loops with OpenMP
   * (these loops are flattened them for better parallelism) */
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime )
# endif
  for (int outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const field_flat_index_range range = inner_flat_range_(outer_ind,
                                                           &ind_helper);

    for (int index = range.start; index <= range.end; index++) {

      // we will only be in this loop if primordial_chemistry > 0
      double number_density =
        0.25 * (my_fields->HeI_density[index] +
		    my_fields->HeII_density[index] +
        my_fields->HeIII_density[index]) +
        my_fields->HI_density[index] + my_fields->HII_density[index] +
        my_fields->e_density[index];

      /* Add in H2. */
 
      if (my_chemistry->primordial_chemistry > 1) {
	number_density += my_fields->HM_density[index] +
	  0.5 * (my_fields->H2I_density[index] +
		 my_fields->H2II_density[index]);
      }

      if (imetal) {
	number_density += my_fields->metal_density[index] * inv_metal_mol;
      }
 
      /* Ignore deuterium. */
 
      temperature[index] *= temperature_units / fmax(number_density,
						    tiny_number);
      temperature[index] = fmax(temperature[index], MINIMUM_TEMPERATURE);
    } // end: loop over i
  } // end: loop over outer_ind

  return GR_SUCCESS;
}


extern "C" int calculate_temperature(code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *temperature)
{
  if (local_calculate_temperature(grackle_data, &grackle_rates, my_units,
                                  my_fields, temperature) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_calculate_temperature.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}
