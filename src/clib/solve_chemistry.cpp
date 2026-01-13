/***********************************************************************
/
/ Solve the chemistry and cooling
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cstdio>
#include "grackle.h"
#include "internal_units.h"
#include "solve_rate_cool_g-cpp.h"
#include "update_UVbackground_rates.hpp"
#include "utils.h"

extern "C" int local_solve_chemistry(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     double dt_value)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return GR_SUCCESS;

  /* Update UV background rates. */
  photo_rate_storage my_uvb_rates;

  my_uvb_rates.k24 = my_uvb_rates.k25 = my_uvb_rates.k26 =
    my_uvb_rates.k27 = my_uvb_rates.k28 = my_uvb_rates.k29 =
    my_uvb_rates.k30 = my_uvb_rates.k31 = my_uvb_rates.piHI =
    my_uvb_rates.piHeI = my_uvb_rates.piHeII = my_uvb_rates.crsHI =
    my_uvb_rates.crsHeI = my_uvb_rates.crsHeII =
    my_uvb_rates.comp_xray = my_uvb_rates.temp_xray = 0.;

  if (my_chemistry->UVbackground == 1) {
    if (grackle::impl::update_UVbackground_rates(
          my_chemistry, my_rates, &my_uvb_rates, my_units) != GR_SUCCESS) {
      std::fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return GR_FAIL;
    }
  }
  else {
    my_uvb_rates.k24       = my_rates->k24;
    my_uvb_rates.k25       = my_rates->k25;
    my_uvb_rates.k26       = my_rates->k26;
    my_uvb_rates.k27       = my_rates->k27;
    my_uvb_rates.k28       = my_rates->k28;
    my_uvb_rates.k29       = my_rates->k29;
    my_uvb_rates.k30       = my_rates->k30;
    my_uvb_rates.k31       = my_rates->k31;
    my_uvb_rates.piHI      = my_rates->piHI;
    my_uvb_rates.piHeI     = my_rates->piHeI;
    my_uvb_rates.piHeII    = my_rates->piHeII;
    my_uvb_rates.crsHI     = my_rates->crsHI;
    my_uvb_rates.crsHeI    = my_rates->crsHeI;
    my_uvb_rates.crsHeII   = my_rates->crsHeII;
    my_uvb_rates.comp_xray = my_rates->comp_xray;
    my_uvb_rates.temp_xray = my_rates->temp_xray;
  }

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  InternalGrUnits internalu = new_internalu_(my_units);

  /* Error checking for H2 shielding approximation */
  if (self_shielding_err_check(my_chemistry, my_fields,
                               "local_solve_chemistry") != GR_SUCCESS) {
    return GR_SUCCESS;
  }

  /* Call the routine to solve cooling equations. */

  int ierr = solve_rate_cool_g(
    metal_field_present, dt_value, internalu,
    my_chemistry, my_rates, my_fields, &my_uvb_rates
  );

  if (ierr != GR_SUCCESS) {
    std::fprintf(stderr, "Error in solve_rate_cool_g.\n");
  }

  return ierr;

}

extern "C" int solve_chemistry(code_units *my_units,
                               grackle_field_data *my_fields,
                               double dt_value)
{
  if (local_solve_chemistry(grackle_data, &grackle_rates,
                            my_units, my_fields, dt_value) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_solve_chemistry.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}
