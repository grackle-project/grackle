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
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */ 

extern void FORTRAN_NAME(calc_temp_cloudy_g)(
        gr_float *d, gr_float *e, gr_float *metal, gr_float *temperature,
	int *in, int *jn, int *kn, int *iexpand, int *imetal,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3, 
 	long long *priDataSize, double *priMMW);

int _calculate_pressure(chemistry_data *my_chemistry,
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
                        gr_float *pressure);

int _calculate_temperature_table(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units,
                                 int grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *metal_density,
                                 gr_float *temperature);
 
int _calculate_temperature(chemistry_data *my_chemistry,
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
                           gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Compute the pressure first. */
 
  if (my_chemistry->primordial_chemistry > 0) {
    if (_calculate_pressure(my_chemistry, my_rates, my_units,
                            grid_rank, grid_dimension,
                            grid_start, grid_end,
                            density, internal_energy,
                            HI_density, HII_density, HM_density,
                            HeI_density, HeII_density, HeIII_density,
                            H2I_density, H2II_density,
                            DI_density, DII_density, HDI_density,
                            e_density, metal_density,
                            temperature) == FAIL) {
      fprintf(stderr, "Error in calculate_pressure.\n");
      return FAIL;
    }
  }
 
  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  /* Calculate temperature units. */

  double temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

  double number_density, tiny_number = 1.-20;
  double inv_metal_mol = 1.0 / MU_METAL;
  
  if (my_chemistry->primordial_chemistry == 0) {
    if (_calculate_temperature_table(my_chemistry, my_rates, my_units,
                                     grid_rank, grid_dimension,
                                     grid_start, grid_end,
                                     density, internal_energy,
                                     metal_density,
                                     temperature) == FAIL) {
      fprintf(stderr, "Error in calculcate_temperature_table.\n");
      return FAIL;
    }
    return SUCCESS;
  }

 /* Compute temperature with mu calculated directly. */
 
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) private( i, number_density )
# endif
  for (i = 0; i < size; i++) {
 
    if (my_chemistry->primordial_chemistry > 0) {
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] +  HeIII_density[i]) +
        HI_density[i] + HII_density[i] + e_density[i];
    }

    /* Add in H2. */
 
    if (my_chemistry->primordial_chemistry > 1) {
      number_density += HM_density[i] + 
        0.5 * (H2I_density[i] + H2II_density[i]);
    }

    if (metal_density != NULL) {
      number_density += metal_density[i] * inv_metal_mol;
    }
 
    /* Ignore deuterium. */
 
    temperature[i] *= temperature_units / max(number_density, tiny_number);
    temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
  }
 
  return SUCCESS;
}

int _calculate_temperature_table(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units,
                                 int grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *metal_density,
                                 gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }
 
  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (metal_density == NULL)
    metal_field_present = FALSE;

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

  /* Calculate temperature units. */

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  FORTRAN_NAME(calc_temp_cloudy_g)(
        density, internal_energy, metal_density, temperature,
        grid_dimension, grid_dimension+1, grid_dimension+2,
        &my_units->comoving_coordinates, &metal_field_present,
        grid_start, grid_start+1, grid_start+2,
        grid_end, grid_end+1, grid_end+2,
        &my_units->a_value,
        &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
        &temperature_units, &co_length_units, &my_units->a_units, 
        &co_density_units, &my_units->time_units,
        &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
        &my_rates->cloudy_primordial.grid_rank,
        my_rates->cloudy_primordial.grid_dimension,
        my_rates->cloudy_primordial.grid_parameters[0],
        my_rates->cloudy_primordial.grid_parameters[1],
        my_rates->cloudy_primordial.grid_parameters[2],
        &my_rates->cloudy_primordial.data_size,
        my_rates->cloudy_primordial.mmw_data);

  return SUCCESS;
}

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature)
{
  if (_calculate_temperature(grackle_data, &grackle_rates, my_units,
                             my_fields->grid_rank, my_fields->grid_dimension,
                             my_fields->grid_start, my_fields->grid_end,
                             my_fields->density, my_fields->internal_energy,
                             my_fields->HI_density, my_fields->HII_density,
                             my_fields->HM_density,
                             my_fields->HeI_density, my_fields->HeII_density,
                             my_fields->HeIII_density,
                             my_fields->H2I_density, my_fields->H2II_density,
                             my_fields->DI_density, my_fields->DII_density,
                             my_fields->HDI_density,
                             my_fields->e_density, my_fields->metal_density,
                             temperature) == FAIL) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
