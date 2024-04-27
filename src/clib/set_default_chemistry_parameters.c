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

#include "grackle_unstable.h" // prototype for grunstable_initialize_field_data

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

int grunstable_initialize_field_data(grackle_field_data *my_fields,
                                     int only_init_datafields)
{
  if (my_fields == NULL) {
    fprintf(stderr, "gr_initial_field_data was passed a NULL pointer\n");
    return FAIL;
  }

  // branching logic based on information encoded by only_init_datafields
  if (only_init_datafields == 0) { // initialize all members
    my_fields->grid_rank = -1;
    my_fields->grid_dimension = NULL;
    my_fields->grid_start = NULL;
    my_fields->grid_end = NULL;
    my_fields->grid_dx = -1.0;
  } else if (only_init_datafields != 1) {
    fprintf(stderr, "gr_initial_field_data received an invalid arg\n");
    return FAIL;
  }

  // now, modify all members holding datafields to have values of NULL
  // (we use X-Macros to do this)
  #define ENTRY(MEMBER_NAME, _1) my_fields->MEMBER_NAME = NULL;
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY

  return SUCCESS;
}
