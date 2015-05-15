/***********************************************************************
/
/ Calculate temperature field (tabulated cooling function)
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
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"
#include "fortran.def"

extern chemistry_data grackle_data;

/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
  
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
 
int _calculate_temperature_table(chemistry_data *my_chemistry,
                                 code_units *my_units, double a_value,
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
      a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(a_value * my_units->a_units, 3);
  }

  /* Calculate temperature units. */

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  FORTRAN_NAME(calc_temp_cloudy_g)(
        density, internal_energy, metal_density, temperature,
        grid_dimension, grid_dimension+1, grid_dimension+2,
        &my_units->comoving_coordinates, &metal_field_present,
        grid_start, grid_start+1, grid_start+2,
        grid_end, grid_end+1, grid_end+2,
        &a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
        &temperature_units, &co_length_units, &my_units->a_units, 
        &co_density_units, &my_units->time_units,
        &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
        &my_chemistry->cloudy_primordial.grid_rank,
        my_chemistry->cloudy_primordial.grid_dimension,
        my_chemistry->cloudy_primordial.grid_parameters[0],
        my_chemistry->cloudy_primordial.grid_parameters[1],
        my_chemistry->cloudy_primordial.grid_parameters[2],
        &my_chemistry->cloudy_primordial.data_size,
        my_chemistry->cloudy_primordial.mmw_data);

  return SUCCESS;
}

int calculate_temperature_table(code_units *my_units, double a_value,
                                int grid_rank, int *grid_dimension,
                                int *grid_start, int *grid_end,
                                gr_float *density, gr_float *internal_energy,
                                gr_float *metal_density,
                                gr_float *temperature)
{
  if (_calculate_temperature_table(&grackle_data,
                                   my_units, a_value,
                                   grid_rank, grid_dimension,
                                   grid_start, grid_end,
                                   density, internal_energy,
                                   metal_density,
                                   temperature) == FAIL) {
    fprintf(stderr, "Error in _calculate_temperature_table.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_temperature_table_(int *comoving_coordinates,
                                 double *density_units, double *length_units,
                                 double *time_units, double *velocity_units,
                                 double *a_units, double *a_value,
                                 int *grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *metal_density,
                                 gr_float *temperature)
{

  code_units my_units;
  my_units.comoving_coordinates = *comoving_coordinates;
  my_units.density_units = *density_units;
  my_units.length_units = *length_units;
  my_units.time_units = *time_units;
  my_units.velocity_units = *velocity_units;
  my_units.a_units = *a_units;

  int rval;
  rval = calculate_temperature_table(&my_units, *a_value,
                                     *grid_rank, grid_dimension,
                                     grid_start, grid_end,
                                     density, internal_energy,
                                     metal_density,
                                     temperature);

  return rval;

}
