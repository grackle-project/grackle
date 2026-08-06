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
#include "internal_units.hpp"
#include "scale_fields.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "utils-cpp.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Fill temperature on a 3D grid.
///
/// @par History
/// written by: Britton Smith May 2015
/// modified1:  February, 2025 by Matthew Abruzzo; ported to C++
///
/// @param[out] temperature_data Array where computed values are written
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  cloudy_primordia specifies the cloudy table
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
template<typename Fn>
static int calc_T_related_(Fn callback,
                   const chemistry_data* my_chemistry,
                   cloudy_data cloudy_primordial, grackle_field_data* my_fields,
                   InternalGrUnits internalu)
{
  const int imetal = (my_fields->metal_density != NULL) ? 1 : 0;

  const DensityUnitKind du_kind = (internalu.extfields_in_comoving == 1) ?
    DensityUnitKind::COMOVING : DensityUnitKind::PROPER;

  const IndexHelper idx_helper = build_index_helper_(my_fields);

  OMP_PRAGMA("omp parallel") {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

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
      basic_gas_props(
          tgas.data(), mmw.data(), rhoH.data(), imetal, itmask.data(),
          my_chemistry, &cloudy_primordial, my_fields, internalu, idx_range,
          du_kind);

      // Record the computed temperature values in the output array
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        callback(tgas[i], mmw[i], i, idx_range.j, idx_range.k);
      }
    }
  }  // OMP_PRAGMA("omp parallel")

  return GR_SUCCESS;
}

}  // namespace GRIMPL_NAMESPACE_DECL
 
extern "C" int local_calculate_temperature(chemistry_data *my_chemistry,
                                           chemistry_data_storage *my_rates,
                                           code_units *my_units,
                                           grackle_field_data *my_fields,
                                           gr_float *temperature)
{
  if (!my_chemistry->use_grackle) { return GR_SUCCESS; }

  // define a callback function that records T values to temperature
  GRIMPL_NS::FortranView<gr_float***> T_view(
    temperature, my_fields->grid_dimension[0],
    my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  auto callback = [T_view](double T_gas, double mmw, int i, int j, int k) {
    T_view(i, j, k) = static_cast<gr_float>(T_gas);
  };

  return GRIMPL_NS::calc_T_related_(callback, my_chemistry,
                                    my_rates->cloudy_primordial, my_fields,
                                    GRIMPL_NS::new_internalu_(my_units));
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

extern "C" int local_calculate_gamma(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *my_gamma)
{

  if (!my_chemistry->use_grackle)
    return GR_SUCCESS;
 
  const GRIMPL_NS::IndexHelper ind_helper
      = GRIMPL_NS::build_index_helper_(my_fields);
  int outer_ind, index;
  
  /* If molecular hydrogen is not being used, just use monotonic.
     (this should not really be called, but provide it just in case). */

  for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const GRIMPL_NS::FieldFlatIndexRange range = GRIMPL_NS::inner_flat_range_
        (outer_ind, &ind_helper);

    for (index = range.start; index <= range.end; index++) {
      my_gamma[index] = my_chemistry->Gamma;
    }
  }
 
  if (my_chemistry->primordial_chemistry > 1) {

    /* Compute the temperature first. */
 
    if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                    my_fields, my_gamma) != GR_SUCCESS) {
      std::fprintf(stderr, "Error in local_calculate_temperature.\n");
      return GR_FAIL;
    }

    /* Compute Gamma with molecular Hydrogen formula from Omukau \& Nishi
       astro-ph/9811308. */
 
    double x, nH2, number_density, GammaH2Inverse, 
      GammaInverse = 1 / (my_chemistry->Gamma - 1.0);

    /* parallelize the k and j loops with OpenMP
     * (these loops are flattened them for better parallelism) */
#   ifdef _OPENMP
#   pragma omp parallel for schedule( runtime ) \
    private( outer_ind, index, x, nH2, number_density, GammaH2Inverse )
#   endif
    for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

      const GRIMPL_NS::FieldFlatIndexRange range = GRIMPL_NS::inner_flat_range_
          (outer_ind, &ind_helper);

      for (index = range.start; index <= range.end; index++) {
 
	/* Compute relative number abundence of molecular hydrogen. */
 
	number_density =
	  0.25 * (my_fields->HeI_density[index] +
		  my_fields->HeII_density[index] +
		  my_fields->HeIII_density[index]) +
	  my_fields->HI_density[index] + my_fields->HII_density[index] +
	  my_fields->HM_density[index] + my_fields->e_density[index];
 
	nH2 = 0.5 * (my_fields->H2I_density[index] +
		     my_fields->H2II_density[index]);
 
	/* Only do full computation if there is a reasonable amount of H2.
	   The second term in GammaH2Inverse accounts for the vibrational
	   degrees of freedom. */
 
	GammaH2Inverse = 0.5*5.0;
	if (nH2 / number_density > 1e-3) {
	  x = 6100.0 / my_gamma[index];
	  if (x < 10.0)
	    GammaH2Inverse = 0.5*(5 + 2.0 * x*x * std::exp(x)/std::pow(std::exp(x)-1.0,2.));
	}

	/* Add in H2. */

	my_gamma[index] = 1.0 + (nH2 + number_density) /
	  (nH2 * GammaH2Inverse + number_density * GammaInverse);
 
      } // end: loop over index
    } // end: loop over outer_ind
 
  } // end: if (my_chemistry->primordial_chemistry > 1)

  return GR_SUCCESS;
}

extern "C" int calculate_gamma(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *my_gamma)
{
  if (local_calculate_gamma(grackle_data, &grackle_rates, my_units,
                            my_fields, my_gamma) != GR_SUCCESS) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}
 