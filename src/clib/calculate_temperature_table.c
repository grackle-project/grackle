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

extern chemistry_data grackle_data;

/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
  
/* function prototypes */ 
 
int _calculate_temperature_table(chemistry_data *my_chemistry,
                                 code_units *my_units,
                                 int grid_rank, int *grid_dimension,
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

  /* Calculate temperature units. */

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;
  double munew, muold;
  int ti, ti_max, index;
  ti_max = 20;

  double logtem0 = log(my_chemistry->TemperatureStart);
  double logtem9 = log(my_chemistry->TemperatureEnd);
  double dlogtem = (log(my_chemistry->TemperatureEnd) - 
                    log(my_chemistry->TemperatureStart)) / 
    (my_chemistry->NumberOfTemperatureBins - 1);
  double logtem, t1, t2, tdef;

  /* Compute temperature with mu calculated directly. */
 
  for (i = 0; i < size; i++) {

    munew = 1.0;

    for (ti = 0; ti < ti_max; ti++) {

      muold = munew;
      temperature[i] = max((my_chemistry->Gamma - 1.) * 
                           internal_energy[i] *
                           munew * temperature_units,
                           my_chemistry->TemperatureStart);
      logtem = log(temperature[i]);
      logtem = max(logtem, logtem0);
      logtem = min(logtem, logtem9);

      index = min(my_chemistry->NumberOfTemperatureBins - 2,
                  max(0, (int) ((logtem-logtem0)/dlogtem)));
      t1 = (logtem0 + (index)     * dlogtem);
      t2 = (logtem0 + (index + 1) * dlogtem);
      tdef = (logtem - t1) / (t2 - t1);
      munew = my_chemistry->mu[index] + tdef
        * (my_chemistry->mu[index+1] - my_chemistry->mu[index]);

      temperature[i] = temperature[i] * munew / muold;

      if (fabs((munew/muold) - 1.) <= 1.e-2) {
        muold = munew;

        // Add metal species to mean molecular weight
          
        munew = density[i] / (density[i] / munew +
                              metal_density[i] / MU_METAL);
        temperature[i] = temperature[i] * munew / muold;
        break;

      }

    }

    if (ti >= ti_max) {
      fprintf(stderr, "Warning: mean molecular weight failed to converge!\n");
    }

  }

  return SUCCESS;
}

int calculate_temperature_table(code_units *my_units,
                                int grid_rank, int *grid_dimension,
                                gr_float *density, gr_float *internal_energy,
                                gr_float *metal_density,
                                gr_float *temperature)
{
  if (_calculate_temperature_table(&grackle_data,
                                   my_units,
                                   grid_rank, grid_dimension,
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
                                 double *a_units,
                                 int *grid_rank, int *grid_dimension,
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
  rval = calculate_temperature_table(&my_units,
                                     *grid_rank, grid_dimension,
                                     density, internal_energy,
                                     metal_density,
                                     temperature);

  return rval;

}
