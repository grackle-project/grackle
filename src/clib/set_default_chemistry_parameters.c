/***********************************************************************
/
/ Set default parameter values
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
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
#include "grackle_chemistry_data.h"

int grackle_verbose = FALSE;

chemistry_data *grackle_data = NULL;
chemistry_data_storage grackle_rates;

int local_initialize_chemistry_parameters(chemistry_data *my_chemistry)
{
  if (my_chemistry == NULL){
    return FAIL;
  }
  // assign the default value to each field of my_chemistry
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) my_chemistry->FIELD = DEFAULT_VAL;
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
  return SUCCESS;
}

int set_default_chemistry_parameters(chemistry_data *my_grackle)
{
  grackle_data = my_grackle;
  return local_initialize_chemistry_parameters(my_grackle);
}

int gr_initialize_field_data(grackle_field_data *my_fields)
{
  if (my_fields == NULL) {
    fprintf(stderr, "gr_initial_field_data was passed a NULL pointer\n");
    return FAIL;
  }

  my_fields->grid_rank = -1;
  my_fields->grid_dimension = NULL;
  my_fields->grid_start = NULL;
  my_fields->grid_end = NULL;
  my_fields->grid_dx = -1.0;

  // below, we use X-Macros to modify all members of my_fields, which are used
  // to hold datafields, to have values of NULL

  // part 1: modify species-field members that grackle can evolve
  #define ENTRY(SPECIES_NAME) my_fields->SPECIES_NAME ## _density = NULL;
  #include "field_data_evolved_species.def"
  #undef ENTRY

  // part 2: modify all other field members
  #define ENTRY(MEMBER_NAME) my_fields->MEMBER_NAME = NULL;
  #include "field_data_misc_fdatamembers.def"
  #undef ENTRY

  // Part 3: modify inject pathway density field slots
  for (int i = 0; i < GRIMPL_MAX_INJ_PATHWAYS; i++) {
    my_fields->inject_pathway_metal_density[i] = NULL;
  }

  return SUCCESS;
}
