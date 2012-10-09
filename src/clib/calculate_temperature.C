/***********************************************************************
/
/  GRID CLASS (COMPUTE THE TEMPERATURE FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"
 
/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */ 

int calculate_pressure(chemistry_data &my_chemistry,
                       code_units &my_units,
                       int grid_rank, int *grid_dimension,
                       float *density, float *internal_energy,
                       float *HI_density, float *HII_density, float *HM_density,
                       float *HeI_density, float *HeII_density, float *HeIII_density,
                       float *H2I_density, float *H2II_density,
                       float *DI_density, float *DII_density, float *HDI_density,
                       float *e_density, float *metal_density,
                       float *pressure);
 
int calculate_temperature(chemistry_data &my_chemistry,
                          code_units &my_units,
                          int grid_rank, int *grid_dimension,
                          float *density, float *internal_energy,
                          float *HI_density, float *HII_density, float *HM_density,
                          float *HeI_density, float *HeII_density, float *HeIII_density,
                          float *H2I_density, float *H2II_density,
                          float *DI_density, float *DII_density, float *HDI_density,
                          float *e_density, float *metal_density,
                          float *temperature)
{

  if (!my_chemistry.use_chemistry) {
    return SUCCESS;
  }

  /* Compute the pressure first. */
 
  if (calculate_pressure(my_chemistry, my_units,
                         grid_rank, grid_dimension,
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
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  float number_density, tiny_number = 1.-20;
  float TemperatureUnits =  mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;
  float inv_metal_mol = 1.0 / MU_METAL;
  
  /* Compute temperature with mu calculated directly. */
 
  for (i = 0; i < size; i++) {
 
    number_density =
      0.25 * (HeI_density[i] + HeII_density[i] +  HeIII_density[i]) +
      HI_density[i] + HII_density[i] + e_density[i];

    /* Add in H2. */
 
    if (my_chemistry.primordial_chemistry > 1)
      number_density += HM_density[i] + 
        0.5 * (H2I_density[i] + H2II_density[i]);

    if (my_chemistry.metal_cooling)
      number_density += metal_density[i] * inv_metal_mol;
 
    /* Ignore deuterium. */
 
    temperature[i] *= TemperatureUnits / max(number_density, tiny_number);
    temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
  }
 
  return SUCCESS;
}
