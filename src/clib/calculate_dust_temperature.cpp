/***********************************************************************
/
/ Calculate dust temperature field
/
/
/ Copyright (c) Grackle Development Team. All rights reserved.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cstdio>
#include "calc_tdust_3d.h"
#include "grackle.h"
#include "internal_units.h"

extern "C" int local_calculate_dust_temperature(
    chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
    code_units *my_units, grackle_field_data *my_fields,
    gr_float *dust_temperature)
{

  if (!my_chemistry->use_grackle)
    return GR_SUCCESS;

  if (my_chemistry->dust_chemistry < 1 && my_chemistry->h2_on_dust < 1)
    return GR_SUCCESS;

  InternalGrUnits internalu = new_internalu_(my_units);

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  /* Compute the size of the fields. */
 
  int size = 1;
  for (int dim = 0; dim < my_fields->grid_rank; dim++) {
    size *= my_fields->grid_dimension[dim];
  }

  gr_float *temperature = new gr_float[size];
  if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                  my_fields, temperature) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_calculate_temperature.\n");
    return GR_FAIL;
  }

  calc_tdust_3d_g(
    temperature, dust_temperature, metal_field_present, my_chemistry, my_rates,
    my_fields, internalu
  );
  delete[] temperature;

  return GR_SUCCESS;
}

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature)
{
  if (local_calculate_dust_temperature(
          grackle_data, &grackle_rates, my_units,
          my_fields, dust_temperature) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_calculate_dust_temperature.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}
