/***********************************************************************
/
/ Calculate dust temperature field
/
/
/ Copyright (c) Grackle Development Team. All rights reserved.
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
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

int _calculate_temperature(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           code_units *my_units,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density,
                           gr_float *HM_density, gr_float *HeI_density,
                           gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density,
                           gr_float *HDI_density, gr_float *e_density,
                           gr_float *metal_density, gr_float *temperature);

extern void FORTRAN_NAME(calc_tdust_3d_g)(
	gr_float *d, gr_float *de, gr_float *HI, gr_float *HII,
	gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
	int *in, int *jn, int *kn,
	int *nratec, int *iexpand,
	int *ispecies, int *idim,
	int *is, int *js, int *ks,
	int *ie, int *je, int *ke,
	double *aye, double *temstart, double *temend,
	double *fgr, double *gasgra,
	double *utem, double *uxyz, double *uaye,
	double *urho, double *utim,
	gr_float *gas_temp, gr_float *dust_temp);

int _calculate_dust_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                int grid_rank, int *grid_dimension,
                                int *grid_start, int *grid_end,
                                gr_float *density, gr_float *internal_energy,
                                gr_float *HI_density, gr_float *HII_density,
                                gr_float *HM_density, gr_float *HeI_density,
                                gr_float *HeII_density, gr_float *HeIII_density,
                                gr_float *H2I_density, gr_float *H2II_density,
                                gr_float *DI_density, gr_float *DII_density,
                                gr_float *HDI_density, gr_float *e_density,
                                gr_float *metal_density, gr_float *dust_temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry < 1)
    return SUCCESS;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }
  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  gr_float *temperature;
  temperature = malloc(size * sizeof(gr_float));
  if (_calculate_temperature(my_chemistry, my_rates, my_units,
                             grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, internal_energy,
                             HI_density, HII_density,
                             HM_density, HeI_density,
                             HeII_density, HeIII_density,
                             H2I_density, H2II_density,
                             DI_density, DII_density,
                             HDI_density, e_density,
                             metal_density, temperature) == FAIL) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return FAIL;
  }

  /* Call the appropriate FORTRAN routine to do the work. */

  FORTRAN_NAME(calc_tdust_3d_g)(
       density, e_density, HI_density, HII_density,
       HeI_density, HeII_density, HeIII_density,
       HM_density, H2I_density, H2II_density,
       grid_dimension, grid_dimension+1, grid_dimension+2,
       &my_chemistry->NumberOfTemperatureBins, &my_units->comoving_coordinates,
       &my_chemistry->primordial_chemistry, &grid_rank,
       grid_start, grid_start+1, grid_start+2,
       grid_end, grid_end+1, grid_end+2,
       &my_units->a_value,
       &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
       &my_chemistry->local_dust_fraction_by_mass,
       my_rates->gas_grain,
       &temperature_units, &co_length_units, &my_units->a_units, 
       &co_density_units, &my_units->time_units,
       temperature, dust_temperature);

  free(temperature);

  return SUCCESS;
}

int local_calculate_dust_temperature(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *dust_temperature)
{
  if (_calculate_dust_temperature(
          my_chemistry, my_rates, my_units,
          my_fields->grid_rank, my_fields->grid_dimension,
          my_fields->grid_start, my_fields->grid_end,
          my_fields->density, my_fields->internal_energy,
          my_fields->HI_density, my_fields->HII_density,
          my_fields->HM_density, my_fields->HeI_density,
          my_fields->HeII_density, my_fields->HeIII_density,
          my_fields->H2I_density, my_fields->H2II_density,
          my_fields->DI_density, my_fields->DII_density,
          my_fields->HDI_density, my_fields->e_density,
          my_fields->metal_density, dust_temperature) == FAIL) {
    fprintf(stderr, "Error in _calculate_dust_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature)
{
  if (local_calculate_dust_temperature(
          grackle_data, &grackle_rates, my_units,
          my_fields, dust_temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_dust_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
