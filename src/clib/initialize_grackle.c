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

#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

extern chemistry_data grackle_data;

int set_default_chemistry_parameters();
 
int initialize_chemistry_data(code_units *my_units, double a_value);

int initialize_grackle(int comoving_coordinates,
                       double density_units, double length_units,
                       double time_units, double velocity_units,
                       double a_units, double a_value,
                       int use_grackle, int with_radiative_cooling,
                       char *grackle_data_file,
                       int primordial_chemistry, int metal_cooling,
                       int UVbackground, int h2_on_dust,
                       int cmb_temperature_floor, double gamma)
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

  grackle_data.use_grackle = use_grackle;
  grackle_data.with_radiative_cooling = with_radiative_cooling;
  grackle_data.grackle_data_file = grackle_data_file;
  grackle_data.primordial_chemistry = primordial_chemistry;
  grackle_data.metal_cooling = metal_cooling;
  grackle_data.UVbackground = UVbackground;
  grackle_data.h2_on_dust = h2_on_dust;
  grackle_data.cmb_temperature_floor = cmb_temperature_floor;
  grackle_data.Gamma = gamma;

  if (initialize_chemistry_data(&my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return FAIL;
  }

  return SUCCESS;
}

int initialize_grackle_(int *comoving_coordinates,
                        double *density_units, double *length_units,
                        double *time_units, double *velocity_units,
                        double *a_units, double *a_value,
                        int *use_grackle, int *with_radiative_cooling,
                        char *grackle_file,
                        int *primordial_chemistry, int *metal_cooling,
                        int *UVbackground, int *h2_on_dust,
                        int *cmb_temperature_floor, double *gamma,
                        int n1)
{

  int i;
  char *grackle_data_file = malloc((n1+1) * sizeof(char));
  for (i = 0; i < n1; i++) {
    grackle_data_file[i] = grackle_file[i];
  }
  grackle_data_file[n1] = NULL; // make NULL-terminated

  int rval;
  rval = initialize_grackle(*comoving_coordinates,
                            *density_units, *length_units,
                            *time_units, *velocity_units,
                            *a_units, *a_value,
                            *use_grackle, *with_radiative_cooling,
                            grackle_data_file,
                            *primordial_chemistry, *metal_cooling,
                            *UVbackground, *h2_on_dust,
                            *cmb_temperature_floor, *gamma);
  return rval;

}
