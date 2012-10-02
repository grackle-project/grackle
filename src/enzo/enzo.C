/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       August 12th 2006
/              May 13th 2008
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/    then calls the appropriate routines depending on whether we are
/    doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
 
#define DEFINE_STORAGE
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"

//  ENZO Main Program

#ifdef SHARED_LIBRARY
#define MAIN_NAME enzo_main
#else
#define MAIN_NAME main
#endif

int InitializeRateData(chemistry_data &my_chemistry,
                       code_units &my_units, float a_value);
int set_default_chemistry_data(chemistry_data &my_chemistry);
int RadiationFieldCalculateRates(chemistry_data &my_chemistry,
                                 code_units &my_units, float a_value);

Eint32 MAIN_NAME(Eint32 argc, char *argv[])
{

  chemistry_data my_chemistry;
  if (set_default_chemistry_data(my_chemistry) == FAIL) {
    fprintf(stderr, "Error in set_default_chemistry_data.\n");
    return FAIL;
  }

  code_units my_units;
  my_units.comoving_coordinates = 0;
  my_units.density_units = 1.0;
  my_units.length_units = 1.0;
  my_units.velocity_units = 1.0;
  my_units.time_units = 1.0;
  my_units.a_units = 1.0;

  float a_value = 1.0;

  if (InitializeRateData(my_chemistry, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in InitializeRateData.\n");
    return FAIL;
  }

  return SUCCESS;
 
 }

 
