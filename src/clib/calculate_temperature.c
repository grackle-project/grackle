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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "calc_temp_cloudy_g-cpp.h"
#include "grackle.h"
#include "grackle_macros.h"
#include "index_helper.h"
#include "internal_units.h"
#include "phys_constants.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Set the mean molecular mass for metals
// -> TODO: this should really be defined by a (internal) header
// -> currently, it's also defined by cool1d_multi_g and calc_temp1d_cloudy_g
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature);
 
int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Compute the pressure first. */
 
  if (my_chemistry->primordial_chemistry > 0) {
    if (local_calculate_pressure(my_chemistry, my_rates, my_units,
                                 my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in calculate_pressure.\n");
      return FAIL;
    }
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  double number_density, tiny_number = 1.-20;
  double inv_metal_mol = 1.0 / MU_METAL;
  
  if (my_chemistry->primordial_chemistry == 0) {
    if (local_calculate_temperature_table(my_chemistry, my_rates, my_units,
                                          my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in local_calculcate_temperature_table.\n");
      return FAIL;
    }
    return SUCCESS;
  }

  /* Compute properties used to index the field. */
  const grackle_index_helper ind_helper = build_index_helper_(my_fields);
  int outer_ind, index;

  /* Compute temperature with mu calculated directly. */

  /* parallelize the k and j loops with OpenMP
   * (these loops are flattened them for better parallelism) */
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) \
  private( outer_ind, index, number_density )
# endif
  for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const field_flat_index_range range = inner_flat_range_(outer_ind,
                                                           &ind_helper);

    for (index = range.start; index <= range.end; index++) {
 
      if (my_chemistry->primordial_chemistry > 0) {
	number_density =
	  0.25 * (my_fields->HeI_density[index] +
		  my_fields->HeII_density[index] +
		  my_fields->HeIII_density[index]) +
	  my_fields->HI_density[index] + my_fields->HII_density[index] +
	  my_fields->e_density[index];
      }

      /* Add in H2. */
 
      if (my_chemistry->primordial_chemistry > 1) {
	number_density += my_fields->HM_density[index] +
	  0.5 * (my_fields->H2I_density[index] +
		 my_fields->H2II_density[index]);
      }

      if (my_fields->metal_density != NULL) {
	number_density += my_fields->metal_density[index] * inv_metal_mol;
      }
 
      /* Ignore deuterium. */
 
      temperature[index] *= temperature_units / max(number_density,
						    tiny_number);
      temperature[index] = max(temperature[index], MINIMUM_TEMPERATURE);
    } // end: loop over i
  } // end: loop over outer_ind

  return SUCCESS;
}

int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }

  InternalGrUnits internalu = new_internalu_(my_units);

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  calc_temp_cloudy_g(
    temperature, metal_field_present, my_chemistry, my_rates->cloudy_primordial,
    my_fields, internalu
  );
  return SUCCESS;
}

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature)
{
  if (local_calculate_temperature(grackle_data, &grackle_rates, my_units,
                                  my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
