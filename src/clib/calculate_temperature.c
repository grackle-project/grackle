/***********************************************************************
/
/ Calculate temperature field
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
#include "index_helper.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */ 

extern void FORTRAN_NAME(calc_temp_cloudy_g)(
        gr_float *d, gr_float *e, gr_float *metal, gr_float *temperature,
	int *in, int *jn, int *kn, int *iexpand, int *imetal,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
	double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3, 
 	long long *priDataSize, double *priMMW);

double get_temperature_units(code_units *my_units);

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure);

int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature);
 
int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Compute the pressure first. */
 
  if (my_chemistry->primordial_chemistry > 0) {
    if (local_calculate_pressure(my_chemistry, my_rates, my_units,
                                 my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in calculate_pressure.\n");
      return FAIL;
    }
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  double number_density, tiny_number = 1.-20;
  double inv_metal_mol = 1.0 / MU_METAL;
  
  if (my_chemistry->primordial_chemistry == 0) {
    if (local_calculate_temperature_table(my_chemistry, my_rates, my_units,
                                          my_fields, temperature) == FAIL) {
      fprintf(stderr, "Error in local_calculcate_temperature_table.\n");
      return FAIL;
    }
    return SUCCESS;
  }

  /* Compute properties used to index the field. */
  const grackle_index_helper ind_helper = _build_index_helper(my_fields);
  int outer_ind, index;

  /* Compute temperature with mu calculated directly. */

  /* parallelize the k and j loops with OpenMP
   * (these loops are flattened them for better parallelism) */
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) \
  private( outer_ind, index, number_density )
# endif
  for (outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){

    const grackle_index_range range = _inner_range(outer_ind, &ind_helper);

    for (index = range.start; index <= range.end; index++) {
 
      if (my_chemistry->primordial_chemistry > 0) {
	number_density =
	  0.25 * (my_fields->HeI_density[index] +
		  my_fields->HeII_density[index] +
		  my_fields->HeIII_density[index]) +
	  my_fields->HI_density[index] + my_fields->HII_density[index] +
	  my_fields->e_density[index];
      }

      /* Add in H2. */
 
      if (my_chemistry->primordial_chemistry > 1) {
	number_density += my_fields->HM_density[index] +
	  0.5 * (my_fields->H2I_density[index] +
		 my_fields->H2II_density[index]);
      }

      if (my_fields->metal_density != NULL) {
	number_density += my_fields->metal_density[index] * inv_metal_mol;
      }
 
      /* Ignore deuterium. */
 
      temperature[index] *= temperature_units / max(number_density,
						    tiny_number);
      temperature[index] = max(temperature[index], MINIMUM_TEMPERATURE);
    } // end: loop over i
  } // end: loop over outer_ind

  return SUCCESS;
}

int local_calculate_temperature_table(chemistry_data *my_chemistry,
                                      chemistry_data_storage *my_rates,
                                      code_units *my_units,
                                      grackle_field_data *my_fields,
                                      gr_float *temperature)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  FORTRAN_NAME(calc_temp_cloudy_g)(
        my_fields->density,
        my_fields->internal_energy,
        my_fields->metal_density,
        temperature,
        my_fields->grid_dimension,
        my_fields->grid_dimension+1,
        my_fields->grid_dimension+2,
        &my_units->comoving_coordinates,
        &metal_field_present,
        my_fields->grid_start,
        my_fields->grid_start+1,
        my_fields->grid_start+2,
        my_fields->grid_end,
        my_fields->grid_end+1,
        my_fields->grid_end+2,
        &my_units->a_value,
        &my_chemistry->TemperatureStart,
        &my_chemistry->TemperatureEnd,
        &temperature_units,
        &co_length_units,
        &my_units->a_units,
        &co_density_units,
        &my_units->time_units,
        &my_chemistry->Gamma,
        &my_chemistry->HydrogenFractionByMass,
        &my_rates->cloudy_primordial.grid_rank,
        my_rates->cloudy_primordial.grid_dimension,
        my_rates->cloudy_primordial.grid_parameters[0],
        my_rates->cloudy_primordial.grid_parameters[1],
        my_rates->cloudy_primordial.grid_parameters[2],
        &my_rates->cloudy_primordial.data_size,
        my_rates->cloudy_primordial.mmw_data);

  return SUCCESS;
}

int _calculate_temperature(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           code_units *my_units,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                           gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                           gr_float *DM_density, gr_float *HDII_density, gr_float *HeHII_density,
                           gr_float *CI_density,  gr_float *CII_density,  gr_float *CO_density,  gr_float *CO2_density,
                           gr_float *OI_density,  gr_float *OH_density,  gr_float *H2O_density,  gr_float *O2_density,
                           gr_float *SiI_density,  gr_float *SiOI_density,  gr_float *SiO2I_density,
                           gr_float *CH_density,  gr_float *CH2_density,  gr_float *COII_density,  gr_float *OII_density,
                           gr_float *OHII_density,  gr_float *H2OII_density,  gr_float *H3OII_density,  gr_float *O2II_density,
                           gr_float *Mg_density,  gr_float *Al_density,  gr_float *S_density,  gr_float *Fe_density,
                           gr_float *SiM_density,  gr_float *FeM_density,  gr_float *Mg2SiO4_density,  gr_float *MgSiO3_density,  gr_float *Fe3O4_density,
                           gr_float *AC_density,  gr_float *SiO2D_density,  gr_float *MgO_density,  gr_float *FeS_density,  gr_float *Al2O3_density,
                           gr_float *reforg_density, gr_float *volorg_density, gr_float *H2Oice_density,
                           gr_float *e_density, gr_float *metal_density,
                           gr_float *temperature)
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
  my_fields.reforg_density           = reforg_density;
  my_fields.volorg_density           = volorg_density;
  my_fields.H2Oice_density           = H2Oice_density;
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;

  if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                  &my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature)
{
  if (local_calculate_temperature(grackle_data, &grackle_rates, my_units,
                                  my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
