/***********************************************************************
/
/ Calculate gamma (ratio of specific heats) field
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
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#include "index_helper.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature);

grackle_index_helper _build_index_helper(const grackle_field_data *my_fields);

int local_calculate_gamma(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *my_gamma)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;
 
  const grackle_index_helper ind_helper = _build_index_helper(my_fields);
  int outer_ind, index;
  
  /* If molecular hydrogen is not being used, just use monotonic.
     (this should not really be called, but provide it just in case). */

  for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const grackle_index_range range = _inner_range(outer_ind, &ind_helper);

    for (index = range.start; index <= range.end; index++) {
      my_gamma[index] = my_chemistry->Gamma;
    }
  }
 
  if (my_chemistry->primordial_chemistry > 1) {

    /* Compute the temperature first. */
 
    if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                    my_fields, my_gamma) == FAIL) {
      fprintf(stderr, "Error in local_calculate_temperature.\n");
      return FAIL;
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

      const grackle_index_range range = _inner_range(outer_ind, &ind_helper);

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
	    GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
	}

	/* Add in H2. */

	my_gamma[index] = 1.0 + (nH2 + number_density) /
	  (nH2 * GammaH2Inverse + number_density * GammaInverse);
 
      } // end: loop over index
    } // end: loop over outer_ind
 
  } // end: if (my_chemistry->primordial_chemistry > 1)

  return SUCCESS;
}

int _calculate_gamma(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *my_gamma)
{

  grackle_field_data my_fields;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;

  if (local_calculate_gamma(my_chemistry, my_rates, my_units,
                       &my_fields, my_gamma) == FAIL) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_gamma(code_units *my_units,
                    grackle_field_data *my_fields,
                    gr_float *my_gamma)
{
  if (local_calculate_gamma(grackle_data, &grackle_rates, my_units,
                            my_fields, my_gamma) == FAIL) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}
