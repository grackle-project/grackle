/***********************************************************************
/
/ Initialize chemistry and cooling rate data
/
/
/ Copyright (c) Enzo/Grackle Development Team. All rights reserved.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"

double get_velocity_units(const code_units *my_units)
{
  double velocity_units = my_units->length_units / my_units->time_units;
  if (my_units->comoving_coordinates == 1) {
    velocity_units /= my_units->a_value;
  }
  return velocity_units;
}

void set_velocity_units(code_units *my_units)
{
  my_units->velocity_units = get_velocity_units(my_units);
}

double get_temperature_units(const code_units *my_units)
{
  double velocity_units = get_velocity_units(my_units);
  return mh * POW(velocity_units, 2) / kboltz;
}
