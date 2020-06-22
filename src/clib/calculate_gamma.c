/***********************************************************************
/
/ Calculate gamma (ratio of specific heats) field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature);

int local_calculate_gamma(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *my_gamma)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;
 
  int i, dim, size = 1;
  for (dim = 0; dim < my_fields->grid_rank; dim++)
    size *= my_fields->grid_dimension[dim];
  
  /* If molecular hydrogen is not being used, just use monotonic.
     (this should not really be called, but provide it just in case). */
 
  for (i = 0; i < size; i++) {
    my_gamma[i] = my_chemistry->Gamma;
  }
 
  if (my_chemistry->primordial_chemistry > 1) {

    /* Compute the temperature first. */
 
    if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                    my_fields, my_gamma) == FAIL) {
      fprintf(stderr, "Error in local_calculate_temperature.\n");
      return FAIL;
    }

    /* Compute Gamma with molecular Hydrogen formula from Omukau \& Nishi
       astro-ph/9811308. */
 
    double x, nH2, number_density, GammaH2Inverse, 
      GammaInverse = 1 / (my_chemistry->Gamma - 1.0);

#   ifdef _OPENMP
#   pragma omp parallel for schedule( runtime ) \
    private( x, nH2, number_density, GammaH2Inverse )
#   endif
    for (i = 0; i < size; i++) {
 
      /* Compute relative number abundence of molecular hydrogen. */
 
      number_density =
        0.25 * (my_fields->HeI_density[i] + my_fields->HeII_density[i] +
                my_fields->HeIII_density[i]) +
        my_fields->HI_density[i] + my_fields->HII_density[i] +
        my_fields->HM_density[i] + my_fields->e_density[i];
 
      nH2 = 0.5 * (my_fields->H2I_density[i] + my_fields->H2II_density[i]);
 
      /* Only do full computation if there is a reasonable amount of H2.
         The second term in GammaH2Inverse accounts for the vibrational
         degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2 / number_density > 1e-3) {
	x = 6100.0 / my_gamma[i];
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      /* Add in H2. */
 
      my_gamma[i] = 1.0 + (nH2 + number_density) /
        (nH2 * GammaH2Inverse + number_density * GammaInverse);
 
    } // end: loop over i
 
  } // end: if (my_chemistry->primordial_chemistry > 1)

  return SUCCESS;
}

int _calculate_gamma(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
//#ifdef GRACKLE_MD
                     gr_float *DM_density, gr_float *HDII_density, gr_float *HeHII_density,
                     gr_float *CI_density,  gr_float *CII_density,  gr_float *CO_density,  gr_float *CO2_density,
                     gr_float *OI_density,  gr_float *OH_density,  gr_float *H2O_density,  gr_float *O2_density,
                     gr_float *SiI_density,  gr_float *SiOI_density,  gr_float *SiO2I_density,
                     gr_float *CH_density,  gr_float *CH2_density,  gr_float *COII_density,  gr_float *OII_density,
                     gr_float *OHII_density,  gr_float *H2OII_density,  gr_float *H3OII_density,  gr_float *O2II_density,
                     gr_float *Mg_density,  gr_float *Al_density,  gr_float *S_density,  gr_float *Fe_density,
                     gr_float *SiM_density,  gr_float *FeM_density,  gr_float *Mg2SiO4_density,  gr_float *MgSiO3_density,  gr_float *Fe3O4_density,
                     gr_float *AC_density,  gr_float *SiO2D_density,  gr_float *MgO_density,  gr_float *FeS_density,  gr_float *Al2O3_density,
//#endif
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *my_gamma)
{

  grackle_field_data my_fields;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
#ifdef GRACKLE_MD
  my_fields.DM_density               = DM_density;
  my_fields.HDII_density             = HDII_density;
  my_fields.HeHII_density            = HeHII_density;
  my_fields.CI_density               = CI_density;
  my_fields.CII_density              = CII_density;
  my_fields.CO_density               = CO_density;
  my_fields.CO2_density              = CO2_density;
  my_fields.OI_density               = OI_density;
  my_fields.OH_density               = OH_density;
  my_fields.H2O_density              = H2O_density;
  my_fields.O2_density               = O2_density;
  my_fields.SiI_density              = SiI_density;
  my_fields.SiOI_density             = SiOI_density;
  my_fields.SiO2I_density            = SiO2I_density;
  my_fields.CH_density               = CH_density;
  my_fields.CH2_density              = CH2_density;
  my_fields.COII_density             = COII_density;
  my_fields.OII_density              = OII_density;
  my_fields.OHII_density             = OHII_density;
  my_fields.H2OII_density            = H2OII_density;
  my_fields.H3OII_density            = H3OII_density;
  my_fields.O2II_density             = O2II_density;
  my_fields.Mg_density               = Mg_density;
  my_fields.Al_density               = Al_density;
  my_fields.S_density                = S_density;
  my_fields.Fe_density               = Fe_density;
  my_fields.SiM_density              = SiM_density;
  my_fields.FeM_density              = FeM_density;
  my_fields.Mg2SiO4_density          = Mg2SiO4_density;
  my_fields.MgSiO3_density           = MgSiO3_density;
  my_fields.Fe3O4_density            = Fe3O4_density;
  my_fields.AC_density               = AC_density;
  my_fields.SiO2D_density            = SiO2D_density;
  my_fields.MgO_density              = MgO_density;
  my_fields.FeS_density              = FeS_density;
  my_fields.Al2O3_density            = Al2O3_density;
#endif
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;

  if (local_calculate_gamma(my_chemistry, my_rates, my_units,
                       &my_fields, my_gamma) == FAIL) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_gamma(code_units *my_units,
                    grackle_field_data *my_fields,
                    gr_float *my_gamma)
{
  if (local_calculate_gamma(grackle_data, &grackle_rates, my_units,
                            my_fields, my_gamma) == FAIL) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}
