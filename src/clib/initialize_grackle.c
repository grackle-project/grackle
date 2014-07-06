/***********************************************************************
/
/ Grackle initialization wrapper.
/
/
/ Copyright (c) 2014, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

extern chemistry_data my_chemistry;

int set_default_chemistry_parameters();
 
int initialize_chemistry_data(code_units *my_units, gr_float a_value);

int initialize_grackle(gr_int comoving_coordinates,
                       gr_float density_units, gr_float length_units,
                       gr_float time_units, gr_float velocity_units,
                       gr_float a_units, gr_float a_value,
                       gr_int use_grackle, gr_int with_radiative_cooling,
                       char *grackle_data_file,
                       gr_int primordial_chemistry, gr_int metal_cooling,
                       gr_int h2_on_dust, gr_int cmb_temperature_floor)
{

  code_units my_units;
  my_units.comoving_coordinates = comoving_coordinates;
  my_units.density_units = density_units;
  my_units.length_units = length_units;
  my_units.time_units = time_units;
  my_units.velocity_units = velocity_units;
  my_units.a_units = a_units;

  if (set_default_chemistry_parameters() == FAIL) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return FAIL;
  }

  my_chemistry.use_grackle = use_grackle;
  my_chemistry.with_radiative_cooling = with_radiative_cooling;
  my_chemistry.grackle_data_file = grackle_data_file;
  my_chemistry.primordial_chemistry = primordial_chemistry;
  my_chemistry.metal_cooling = metal_cooling;
  my_chemistry.h2_on_dust = h2_on_dust;
  my_chemistry.cmb_temperature_floor = cmb_temperature_floor;

  if (initialize_chemistry_data(&my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return FAIL;
  }

  return SUCCESS;
}

int initialize_grackle_(gr_int *comoving_coordinates,
                        gr_float *density_units, gr_float *length_units,
                        gr_float *time_units, gr_float *velocity_units,
                        gr_float *a_units, gr_float *a_value,
                        gr_int *use_grackle, gr_int *with_radiative_cooling,
                        char *grackle_file,
                        gr_int *primordial_chemistry, gr_int *metal_cooling,
                        gr_int *h2_on_dust, gr_int *cmb_temperature_floor,
                        int n1)
{

  int i;
  char *grackle_data_file = malloc(n1 * sizeof(char));
  for (i = 0; i < n1; i++) {
    grackle_data_file[i] = grackle_file[i];
  }

  return initialize_grackle(*comoving_coordinates,
                            *density_units, *length_units,
                            *time_units, *velocity_units,
                            *a_units, *a_value,
                            *use_grackle, *with_radiative_cooling,
                            grackle_data_file,
                            *primordial_chemistry, *metal_cooling,
                            *h2_on_dust, *cmb_temperature_floor);

}
