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
