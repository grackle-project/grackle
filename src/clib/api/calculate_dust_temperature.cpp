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
#include "dust/calc_tdust_3d.hpp"
#include "grackle.h"
#include "internal_units.hpp"

extern "C" int local_calculate_dust_temperature(
    chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
    code_units *my_units, grackle_field_data *my_fields,
    gr_float *dust_temperature)
{

  if (!my_chemistry->use_grackle)
    return GR_SUCCESS;

  if (my_chemistry->dust_chemistry < 1 && my_chemistry->h2_on_dust < 1)
    return GR_SUCCESS;

  GRIMPL_NS::InternalGrUnits internalu = GRIMPL_NS::new_internalu_(my_units);

  int has_metal_field = (my_fields->metal_density == nullptr) ? FALSE : TRUE;
  GRIMPL_NS::calc_tdust_3d(
    dust_temperature, has_metal_field, my_chemistry, my_rates,
    my_fields, internalu
  );

  return GR_SUCCESS;
}

extern "C" int calculate_dust_temperature(code_units *my_units,
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
