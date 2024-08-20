/***********************************************************************
/
/ Implement utility functions used internally by Grackle (across routines)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include "utils.h"
#include <stdio.h> // fprintf, stderr
#include "grackle_macros.h"

int self_shielding_err_check(const chemistry_data *my_chemistry,
                             const grackle_field_data *fields,
                             const char* func_name) {
  if (my_chemistry->H2_self_shielding == 1) {
    if (fields->grid_rank != 3) {
      fprintf(stderr, "Error in %s: H2 self-shielding option 1 "
                      "will only work for 3D Cartesian grids. Use option 2 "
                      "to provide an array of shielding lengths with "
                      "H2_self_shielding_length or option 3 to use the "
                      "local Jeans length.",
              func_name);
      return FAIL;
    } else if (my_chemistry->primordial_chemistry >= 2 &&
               fields->grid_dx <= 0) {
      fprintf(stderr, "Error in %s: H2 self-shielding option 1 and primordial "
                      "chemistry options of 2 or more require that grid_dx "
                      "has a positive value.",
              func_name);
      return FAIL;
    }
  }
  return SUCCESS;
}
