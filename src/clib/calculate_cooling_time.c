/***********************************************************************
/
/ Calculate cooling time field
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

extern chemistry_data grackle_data;

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              code_units *my_units);
 
extern void FORTRAN_NAME(cool_multi_time_g)(
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	gr_float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, 
        int *ih2co, int *ipiht, int *igammah,
	double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *z_solar,
	double *ceHIa, double *ceHeIa, double *ceHeIIa, 
        double *ciHIa, double *ciHeIa,
	double *ciHeISa, double *ciHeIIa, double *reHIIa, double *reHeII1a,
	double *reHeII2a, double *reHeIIIa, double *brema, double *compa, 
        double *gammaha, double *comp_xraya, double *comp_temp, 
        double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II, 
        gr_float *DI, gr_float *DII, gr_float *HDI, gr_float *metal,
	double *hyd01ka, double *h2k01a, double *vibha, 
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra,
        int *iradtrans, double *photogamma, // AJE-RT
	int *ih2optical, int *iciecool, double *ciecoa,
 	int *icmbTfloor, int *iClHeat, 
        long long *priGridRank, long long *priGridDim,
 	double *priPar1, double *priPar2, double *priPar3, 
 	long long *priDataSize, double *priCooling, double *priHeating,
        double *priMMW,
        long long *metGridRank, long long *metGridDim,
 	double *metPar1, double *metPar2, double *metPar3, 
 	long long *metDataSize, double *metCooling, double *metHeating,
        int *iVheat, int *iMheat, gr_float *Vheat, gr_float *Mheat);

int _calculate_cooling_time(chemistry_data *my_chemistry,
                            code_units *my_units,
                            int grid_rank, int *grid_dimension,
                            int *grid_start, int *grid_end,
                            gr_float *density, gr_float *internal_energy,
                            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                            gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                            gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                            gr_float *H2I_density, gr_float *H2II_density,
                            gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                            gr_float *e_density, gr_float *metal_density,
                            gr_float *cooling_time, gr_float *gammaNum,
                            gr_float *volumetric_heating_rate, gr_float *specific_heating_rate)
{
 
  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_units) == FAIL) {
      fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return FAIL;
    }
  }

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (metal_density == NULL)
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

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Call the fortran routine to solve cooling equations. */

    FORTRAN_NAME(cool_multi_time_g)(
       density, internal_energy, x_velocity, y_velocity, z_velocity,
       e_density, HI_density, HII_density,
       HeI_density, HeII_density, HeIII_density,
       cooling_time,
       grid_dimension, grid_dimension+1, grid_dimension+2,
       &my_chemistry->NumberOfTemperatureBins, &my_units->comoving_coordinates,
       &my_chemistry->primordial_chemistry, &metal_field_present, &my_chemistry->metal_cooling, 
       &my_chemistry->h2_on_dust,
       &grid_rank, grid_start, grid_start+1, grid_start+2,
       grid_end, grid_end+1, grid_end+2,
       &my_chemistry->ih2co, &my_chemistry->ipiht, &my_chemistry->photoelectric_heating,
       &my_units->a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
       &temperature_units, &co_length_units, &my_units->a_units, 
       &co_density_units, &my_units->time_units,
       &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
       &my_chemistry->SolarMetalFractionByMass,
       my_chemistry->ceHI, my_chemistry->ceHeI, my_chemistry->ceHeII, my_chemistry->ciHI,
       my_chemistry->ciHeI, my_chemistry->ciHeIS, my_chemistry->ciHeII, my_chemistry->reHII,
       my_chemistry->reHeII1, my_chemistry->reHeII2, my_chemistry->reHeIII, 
       my_chemistry->brem, &my_chemistry->comp, &my_chemistry->gammah,
       &my_chemistry->comp_xray, &my_chemistry->temp_xray,
       &my_chemistry->piHI, &my_chemistry->piHeI, &my_chemistry->piHeII,
       HM_density, H2I_density, H2II_density,
       DI_density, DII_density, HDI_density, metal_density,
       my_chemistry->hyd01k, my_chemistry->h2k01, my_chemistry->vibh,
       my_chemistry->roth, my_chemistry->rotl,
       my_chemistry->GP99LowDensityLimit, my_chemistry->GP99HighDensityLimit,
       my_chemistry->HDlte, my_chemistry->HDlow,
       my_chemistry->GAHI, my_chemistry->GAH2, my_chemistry->GAHe, my_chemistry->GAHp,
       my_chemistry->GAel, my_chemistry->H2LTE, my_chemistry->gas_grain,
       &my_chemistry->use_radiative_transfer, gammaNum, // AJE-RT
       &my_chemistry->h2_optical_depth_approximation, 
       &my_chemistry->cie_cooling, my_chemistry->cieco,
       &my_chemistry->cmb_temperature_floor,
       &my_chemistry->UVbackground,
       &my_chemistry->cloudy_primordial.grid_rank,
       my_chemistry->cloudy_primordial.grid_dimension,
       my_chemistry->cloudy_primordial.grid_parameters[0],
       my_chemistry->cloudy_primordial.grid_parameters[1],
       my_chemistry->cloudy_primordial.grid_parameters[2],
       &my_chemistry->cloudy_primordial.data_size,
       my_chemistry->cloudy_primordial.cooling_data, 
       my_chemistry->cloudy_primordial.heating_data,
       my_chemistry->cloudy_primordial.mmw_data,
       &my_chemistry->cloudy_metal.grid_rank,
       my_chemistry->cloudy_metal.grid_dimension,
       my_chemistry->cloudy_metal.grid_parameters[0],
       my_chemistry->cloudy_metal.grid_parameters[1],
       my_chemistry->cloudy_metal.grid_parameters[2],
       &my_chemistry->cloudy_metal.data_size,
       my_chemistry->cloudy_metal.cooling_data,
       my_chemistry->cloudy_metal.heating_data,
       &my_chemistry->use_volumetric_heating_rate,
       &my_chemistry->use_specific_heating_rate,
       volumetric_heating_rate, specific_heating_rate);
 
  return SUCCESS;
}

int calculate_cooling_time(code_units *my_units,
                           grackle_field_data *my_fields,
                           gr_float *cooling_time)
{

  if (_calculate_cooling_time(&grackle_data,
                              my_units,
                              my_fields->grid_rank, my_fields->grid_dimension,
                              my_fields->grid_start, my_fields->grid_end,
                              my_fields->density, my_fields->internal_energy,
                              my_fields->x_velocity, my_fields->y_velocity,
                              my_fields->z_velocity,
                              my_fields->HI_density, my_fields->HII_density,
                              my_fields->HM_density,
                              my_fields->HeI_density, my_fields->HeII_density,
                              my_fields->HeIII_density,
                              my_fields->H2I_density, my_fields->H2II_density,
                              my_fields->DI_density, my_fields->DII_density,
                              my_fields->HDI_density,
                              my_fields->e_density, my_fields->metal_density,
                              cooling_time, my_fields->gammaNum,
                              my_fields->volumetric_heating_rate,
                              my_fields->specific_heating_rate) == FAIL) {
    fprintf(stderr, "Error in _calculate_cooling_time.\n");
    return FAIL;
  }
  return SUCCESS;
}
