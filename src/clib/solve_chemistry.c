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

#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              code_units *my_units);

extern void FORTRAN_NAME(solve_rate_cool_g)(
        int *icool,
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
        int *ih2co, int *ipiht, int *igammah,
	double *dx, double *dt, double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *dtoh, double *z_solar,
	double *k1a, double *k2a, double *k3a, double *k4a, double *k5a,
	double *k6a, double *k7a, double *k8a, double *k9a, double *k10a,
	double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a,
	double *k15a, double *k16a, double *k17a, double *k18a, double *k19a,
        double *k22a,	double *k24, double *k25, double *k26, double *k27,
        double *k28, double *k29, double *k30, double *k31,
	double *k50a, double *k51a, double *k52a, double *k53a, double *k54a,
	double *k55a, double *k56a, double *k57a, double *k58a,
	int *ndratec, double *dtemstart, double *dtemend, double *h2dusta,
	double *ncrna, double *ncrd1a, double *ncrd2a,
	double *ceHIa, double *ceHeIa, double *ceHeIIa, double *ciHIa,
	double *ciHeIa, double *ciHeISa, double *ciHeIIa,
        double *reHIIa, double *reHeII1a, double *reHeII2a, double *reHeIIIa,
        double *brema, double *compa, double *gammaha,
	double *comp_xraya, double *comp_temp,
	double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
        gr_float *DI, gr_float *DII, gr_float *HDI, gr_float *metal,
	double *hyd01ka, double *h2k01a, double *vibha,
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra, int *iH2shield,
        int *iradshield, double *avgsighi, double *avgsighei, double *avgsigheii,
        int *iradtrans, int *iradcoupled, int *iradstep, int *irt_honly,
        double *kphHI, double *kphHeI, double *kphHeII, double *kdissH2I,
        double *photogamma,
	int *ierr,
	int *ih2optical, int *iciecool, int *ithreebody, double *ciecoa,
 	int *icmbTfloor, int *iClHeat, double *clEleFra,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3,
        double *priPar4, double *priPar5,
 	long long *priDataSize, double *priCooling,
        double *priHeating, double *priMMW,
        long long *metGridRank, long long *metGridDim,
 	double *metPar1, double *metPar2, double *metPar3,
        double *metPar4, double *metPar5,
 	long long *metDataSize, double *metCooling,
        double *metHeating, int *clnew,
        int *iVheat, int *iMheat, gr_float *Vheat, gr_float *Mheat);

int _solve_chemistry(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units, double dt_value, double dx_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *volumetric_heating_rate, gr_float *specific_heating_rate,
                     gr_float *RT_heating_rate, gr_float *RT_HI_ionization_rate, gr_float *RT_HeI_ionization_rate,
                     gr_float *RT_HeII_ionization_rate, gr_float *RT_H2_dissociation_rate)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_rates, my_units) == FAIL) {
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

  /* Error checking for H2 shielding approximation */
  if (my_chemistry->H2_self_shielding && grid_rank != 3){
    fprintf(stderr, "Error in solve_chemistry: H2 self shielding approximation "
                    "is turned on, yet is only valid for 3D Cartesian grids. "
                    "grid_rank is currently %i \n.", grid_rank);
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  FORTRAN_NAME(solve_rate_cool_g)(
    &my_chemistry->with_radiative_cooling,
    density, internal_energy, x_velocity, y_velocity, z_velocity,
    e_density, HI_density, HII_density,
    HeI_density, HeII_density, HeIII_density,
    grid_dimension, grid_dimension+1, grid_dimension+2,
    &my_chemistry->NumberOfTemperatureBins, &my_units->comoving_coordinates,
    &my_chemistry->primordial_chemistry, &metal_field_present, &my_chemistry->metal_cooling,
    &my_chemistry->h2_on_dust, &grid_rank, grid_start, grid_start+1, grid_start+2,
    grid_end, grid_end+1, grid_end+2,
    &my_chemistry->ih2co, &my_chemistry->ipiht, &my_chemistry->photoelectric_heating,
    &dx_value, &dt_value, &my_units->a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
    &temperature_units, &co_length_units, &my_units->a_units,
    &co_density_units, &my_units->time_units, &my_chemistry->Gamma,
    &my_chemistry->HydrogenFractionByMass, &my_chemistry->DeuteriumToHydrogenRatio,
    &my_chemistry->SolarMetalFractionByMass,
    my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5,
    my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
    my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd,
    my_rates->k14, my_rates->k15, my_rates->k16,
    my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
    &my_rates->k24, &my_rates->k25, &my_rates->k26, &my_rates->k27,
    &my_rates->k28, &my_rates->k29, &my_rates->k30, &my_rates->k31,
    my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53,
    my_rates->k54, my_rates->k55, my_rates->k56,
    my_rates->k57, my_rates->k58,
    &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart,
    &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
    my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
    my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI,
    my_rates->ciHeI, my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII,
    my_rates->reHeII1, my_rates->reHeII2, my_rates->reHeIII, my_rates->brem,
    &my_rates->comp, &my_rates->gammah,
    &my_rates->comp_xray, &my_rates->temp_xray,
    &my_rates->piHI, &my_rates->piHeI, &my_rates->piHeII,
    HM_density, H2I_density, H2II_density,
    DI_density, DII_density, HDI_density, metal_density,
    my_rates->hyd01k, my_rates->h2k01, my_rates->vibh,
    my_rates->roth, my_rates->rotl,
    my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit,
    my_rates->HDlte, my_rates->HDlow,
    my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp,
    my_rates->GAel, my_rates->H2LTE, my_rates->gas_grain,
    &my_chemistry->H2_self_shielding,
    &my_chemistry->self_shielding_method, &my_rates->hi_avg_crs,
    &my_rates->hei_avg_crs, &my_rates->heii_avg_crs,
    &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_coupled_rate_solver,
    &my_chemistry->radiative_transfer_intermediate_step, &my_chemistry->radiative_transfer_hydrogen_only,
    RT_HI_ionization_rate, RT_HeI_ionization_rate, RT_HeII_ionization_rate,
    RT_H2_dissociation_rate, RT_heating_rate,
    &ierr,
    &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling,
    &my_chemistry->three_body_rate, my_rates->cieco,
    &my_chemistry->cmb_temperature_floor,
    &my_chemistry->UVbackground,
    &my_chemistry->cloudy_electron_fraction_factor,
    &my_rates->cloudy_primordial.grid_rank,
    my_rates->cloudy_primordial.grid_dimension,
    my_rates->cloudy_primordial.grid_parameters[0],
    my_rates->cloudy_primordial.grid_parameters[1],
    my_rates->cloudy_primordial.grid_parameters[2],
    my_rates->cloudy_primordial.grid_parameters[3],
    my_rates->cloudy_primordial.grid_parameters[4],
    &my_rates->cloudy_primordial.data_size,
    my_rates->cloudy_primordial.cooling_data,
    my_rates->cloudy_primordial.heating_data,
    my_rates->cloudy_primordial.mmw_data,
    &my_rates->cloudy_metal.grid_rank,
    my_rates->cloudy_metal.grid_dimension,
    my_rates->cloudy_metal.grid_parameters[0],
    my_rates->cloudy_metal.grid_parameters[1],
    my_rates->cloudy_metal.grid_parameters[2],
    my_rates->cloudy_metal.grid_parameters[3],
    my_rates->cloudy_metal.grid_parameters[4],
    &my_rates->cloudy_metal.data_size,
    my_rates->cloudy_metal.cooling_data,
    my_rates->cloudy_metal.heating_data,
    &my_rates->cloudy_data_new,
    &my_chemistry->use_volumetric_heating_rate,
    &my_chemistry->use_specific_heating_rate,
    volumetric_heating_rate, specific_heating_rate);

  return SUCCESS;

}

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value)
{
  if (_solve_chemistry(grackle_data, &grackle_rates,
                       my_units, dt_value, my_fields->grid_dx,
                       my_fields->grid_rank,   my_fields->grid_dimension,
                       my_fields->grid_start,  my_fields->grid_end,
                       my_fields->density,     my_fields->internal_energy,
                       my_fields->x_velocity,  my_fields->y_velocity,
                       my_fields->z_velocity,
                       my_fields->HI_density,  my_fields->HII_density,
                       my_fields->HM_density,
                       my_fields->HeI_density, my_fields->HeII_density,
                       my_fields->HeIII_density,
                       my_fields->H2I_density, my_fields->H2II_density,
                       my_fields->DI_density,  my_fields->DII_density,
                       my_fields->HDI_density,
                       my_fields->e_density,   my_fields->metal_density,
                       my_fields->volumetric_heating_rate,
                       my_fields->specific_heating_rate,
                       my_fields->RT_heating_rate, my_fields->RT_HI_ionization_rate,
                       my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate,
                       my_fields->RT_H2_dissociation_rate) == FAIL) {
    fprintf(stderr, "Error in _solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}

