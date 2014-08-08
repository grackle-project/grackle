/***********************************************************************
/
/ Calculate pressure field (tabulated cooling function)
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

extern chemistry_data grackle_data;

int _calculate_pressure_table(chemistry_data *my_chemistry,
                              code_units *my_units,
                              gr_int grid_rank, gr_int *grid_dimension,
                              gr_float *density, gr_float *internal_energy,
                              gr_float *pressure)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }

  gr_float tiny_number = 1.e-20;
  gr_int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  for (i = 0; i < size; i++) {
 
    pressure[i] = (my_chemistry->Gamma - 1.0) * density[i] * internal_energy[i];
 
    if (pressure[i] < tiny_number)
      pressure[i] = tiny_number;
  }  
 
  return SUCCESS;
}

int calculate_pressure_table(code_units *my_units,
                             gr_int grid_rank, gr_int *grid_dimension,
                             gr_float *density, gr_float *internal_energy,
                             gr_float *pressure)
{
  if (_calculate_pressure_table(&grackle_data,
                                my_units,
                                grid_rank, grid_dimension,
                                density, internal_energy,
                                pressure) == FAIL) {
    fprintf(stderr, "Error in _calculate_pressure_table.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_pressure_table_(gr_int *comoving_coordinates,
                              gr_float *density_units, gr_float *length_units,
                              gr_float *time_units, gr_float *velocity_units,
                              gr_float *a_units,
                              gr_int *grid_rank, gr_int *grid_dimension,
                              gr_float *density, gr_float *internal_energy,
                              gr_float *pressure)
{

  code_units my_units;
  my_units.comoving_coordinates = *comoving_coordinates;
  my_units.density_units = *density_units;
  my_units.length_units = *length_units;
  my_units.time_units = *time_units;
  my_units.velocity_units = *velocity_units;
  my_units.a_units = *a_units;

  int rval;
  rval = calculate_pressure_table(&my_units,
                                  *grid_rank, grid_dimension,
                                  density, internal_energy,
                                  pressure);

  return rval;

}
