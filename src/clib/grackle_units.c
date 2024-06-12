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
#include "grackle.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"

double get_velocity_units(code_units *my_units)
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

static double get_temperature_units_(double velocity_units)
{
  return mh * POW(velocity_units, 2) / kboltz;
}


double get_temperature_units(code_units *my_units)
{
  return get_temperature_units_(get_velocity_units(my_units));
}

/// helper-function used to build help gr_required_units implement a
/// quick-and-dirty calculations of velocity_units
static double required_velocity_units_(const chemistry_data_storage * my_rates,
                                       double current_a_value)
{
  // assume: (my_rates != NULL) AND (current_a_value > 0.0)

  // quick-and-dirty solution, build fully-configured temporary code_units
  // instance for the desired a_value, then pass to get_velocity_units
  // (this can definitely be optimized!)
  code_units tmp;
  tmp.a_value = current_a_value;
  tmp.a_units = my_rates->initial_units.a_units;
  tmp.comoving_coordinates = my_rates->initial_units.comoving_coordinates;
  tmp.length_units = gr_required_units(my_rates, "length_units",
                                       current_a_value);
  tmp.density_units = gr_required_units(my_rates, "density_units",
                                        current_a_value);
  tmp.time_units = gr_required_units(my_rates, "time_units", current_a_value); 

  return get_velocity_units(&tmp);
}

#define _ERR_UNIT_RETURN -1.0

double gr_required_units(const chemistry_data_storage * my_rates,
                         const char* units_name, double current_a_value)
{
  if (my_rates == NULL) {
    fprintf(stderr, "my_rates argument is NULL\n");
    return _ERR_UNIT_RETURN; // maybe we should abort?
  }
  const code_units* initial_units = &(my_rates->initial_units);

  // tracks whether the redshift has changed from the initial value
  int is_unchanged = (current_a_value == initial_units->a_value);

  // sanitize the current_a_value arg
  if (current_a_value == -1.0) {
    current_a_value = initial_units->a_value;
    is_unchanged = 1;
  } else if ((initial_units->comoving_coordinates == 0) && !is_unchanged) {
    fprintf(stderr, ("for non-comoving coordinates, current_a_value must be "
                     "-1.0 or it must EXACTLY match the initial a_value\n"));
    return _ERR_UNIT_RETURN;

  } else if (current_a_value <= 0.0) {
    fprintf(stderr, "current_a_value must be -1 or a positive value\n");
    return _ERR_UNIT_RETURN;

  }

  // now handle check the units name
  if (units_name == NULL) {
    fprintf(stderr, "units_name argument is NULL\n");
    return _ERR_UNIT_RETURN; // maybe we should abort?

  } else if (strcmp(units_name, "density_units") == 0) {
    double init_d_u = initial_units->density_units;
    return (is_unchanged) ? init_d_u
      : init_d_u * pow( (initial_units->a_value / current_a_value), 3);

  } else if (strcmp(units_name, "length_units") == 0) {
    double init_l_u = initial_units->length_units;
    return (is_unchanged) ? init_l_u
      : init_l_u * (current_a_value / initial_units->a_value);

  } else if (strcmp(units_name, "time_units") == 0) {
    return initial_units->time_units;

  } else if (strcmp(units_name, "velocity_units") == 0) {
    return required_velocity_units_(my_rates, current_a_value);

  } else if (strcmp(units_name, "temperature_units") == 0) {
    return get_temperature_units_(required_velocity_units_(my_rates,
                                                           current_a_value));

  } else if (strcmp(units_name, "a_value") == 0) {
    return current_a_value;

  } else if (strcmp(units_name, "a_units") == 0) {
    return initial_units->a_units;

  } else {
    fprintf(stderr, "unknown units name: \"%s\"\n", units_name);
    return _ERR_UNIT_RETURN;
  }
}
