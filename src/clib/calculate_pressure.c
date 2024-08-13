/***********************************************************************
/
/ Calculate pressure field
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
#include "grackle.h"
#include "grackle_macros.h"
#include "phys_constants.h"
#include "index_helper.h"
#include "unit_handling.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* do unit-handling */
  code_units units = determine_code_units(my_units, my_rates,
                                          my_fields->current_a_value,
                                          my_chemistry->unit_handling,
                                          "calculate_pressure");
  if (units.a_units < 0) {
    return FAIL;
  } else {
    my_units = &units;
  }


  double tiny_number = 1.e-20;
  const grackle_index_helper ind_helper = _build_index_helper(my_fields);
  int outer_ind, index;

  /* parallelize the k and j loops with OpenMP
   * (these loops are flattened them for better parallelism) */
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) private( outer_ind, index )
# endif
  for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const grackle_index_range range = _inner_range(outer_ind, &ind_helper);

    for (index = range.start; index <= range.end; index++) {

      pressure[index] = ((my_chemistry->Gamma - 1.0) *
			 my_fields->density[index] *
			 my_fields->internal_energy[index]);
 
      if (pressure[index] < tiny_number)
        pressure[index] = tiny_number;
    } // end: loop over i
  } // end: loop over outer_ind

  /* Correct for Gamma from H2. */

  if (my_chemistry->primordial_chemistry > 1) {
 
    /* Calculate temperature units. */

    double temperature_units = get_temperature_units(my_units);

    double number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(my_chemistry->Gamma-1.0), x, Gamma1, temp;
  
#   ifdef _OPENMP
#   pragma omp parallel for schedule( runtime ) \
    private( outer_ind, index, \
             number_density, nH2, GammaH2Inverse, x, Gamma1, temp )
#   endif
    for (int outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

      const grackle_index_range range = _inner_range(outer_ind, &ind_helper);

      for (index = range.start; index <= range.end; index++) {

        number_density =
          0.25 * (my_fields->HeI_density[index] +
		  my_fields->HeII_density[index] +
                  my_fields->HeIII_density[index]) +
          my_fields->HI_density[index] + my_fields->HII_density[index] +
          my_fields->HM_density[index] + my_fields->e_density[index];

        nH2 = 0.5 * (my_fields->H2I_density[index] +
		     my_fields->H2II_density[index]);

        /* First, approximate temperature. */

        if (number_density == 0)
          number_density = tiny_number;
        temp = max(temperature_units * pressure[index] / (number_density + nH2),
		   1);

        /* Only do full computation if there is a reasonable amount of H2.
	   The second term in GammaH2Inverse accounts for the vibrational
	   degrees of freedom. */

        GammaH2Inverse = 0.5*5.0;
        if (nH2 / number_density > 1e-3) {
          x = 6100.0 / temp;
	  if (x < 10.0)
	    GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
        }

	Gamma1 = 1.0 + (nH2 + number_density) /
	               (nH2 * GammaH2Inverse + number_density * GammaInverse);
	
	/* Correct pressure with improved Gamma. */
 
	pressure[index] *= (Gamma1 - 1.0) / (my_chemistry->Gamma - 1.0);
 
      } // end: loop over i
    } // end: loop over outer_ind
 
  } // end: if (my_chemistry->primordial_chemistry > 1)
 
  return SUCCESS;
}

int calculate_pressure(code_units *my_units,
                       grackle_field_data *my_fields,
                       gr_float *pressure)
{
  if (local_calculate_pressure(grackle_data, &grackle_rates, my_units,
                               my_fields, pressure) == FAIL) {
    fprintf(stderr, "Error in local_calculate_pressure.\n");
    return FAIL;
  }
  return SUCCESS;
}
