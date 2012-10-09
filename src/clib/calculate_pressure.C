/***********************************************************************
/
/  GRID CLASS (COMPUTE THE PRESSURE FIELD AT THE GIVEN TIME)
/
/  written by: Greg Bryan
/  date:       November, 1994
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

int calculate_pressure(chemistry_data &my_chemistry,
                       code_units &my_units,
                       int grid_rank, int *grid_dimension,
                       float *density, float *internal_energy,
                       float *HI_density, float *HII_density, float *HM_density,
                       float *HeI_density, float *HeII_density, float *HeIII_density,
                       float *H2I_density, float *H2II_density,
                       float *DI_density, float *DII_density, float *HDI_density,
                       float *e_density, float *metal_density,
                       float *pressure)
{

  float tiny_number = 1.e-20;
  int i, size = 1;
  for (int dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

  for (i = 0; i < size; i++) {
 
    pressure[i] = (my_chemistry.Gamma - 1.0) * density[i] * internal_energy[i];
 
    if (pressure[i] < tiny_number)
      pressure[i] = tiny_number;
  }
  
  /* Correct for Gamma from H2. */
 
  if (my_chemistry.primordial_chemistry > 1) {
 
    float TemperatureUnits =  mh*POW(my_units.length_units/
                                     my_units.time_units,2)/kboltz;

    float number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(my_chemistry.Gamma-1.0), x, Gamma1, temp;
  
    for (i = 0; i < size; i++) {
 
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
        HI_density[i] + HII_density[i] + HM_density[i] +
        e_density[i];
 
      nH2 = 0.5 * (H2I_density[i] + H2II_density[i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = tiny_number;
      temp = max(TemperatureUnits * pressure[i] / (number_density + nH2), 1);
 
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
 
      pressure[i] *= (Gamma1 - 1.0) / (my_chemistry.Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (my_chemistry.primordial_chemistry > 1)
 
  return SUCCESS;
}
