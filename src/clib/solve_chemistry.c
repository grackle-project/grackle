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
#include "grackle.h"
#include "grackle_macros.h"
#include "phys_constants.h"
#include "unit_handling.h"
#include "utils.h"

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units);

extern void FORTRAN_NAME(solve_rate_cool_g)(
        int *icool,
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idustall,
        int *idustfield, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
        int *ih2co, int *ipiht, int *idustrec, int *igammah,
	double *dx, double *dt, double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *dtoh, double *z_solar, double *fgr,
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
        double *brema, double *compa, double *gammaha, double *isrf,
        double *regra, double *gamma_isrfa, double *comp_xraya, double *comp_temp,
	double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
        gr_float *DI, gr_float *DII, gr_float *HDI,
        gr_float *metal, gr_float *dust,
	double *hyd01ka, double *h2k01a, double *vibha,
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra, int *iH2shield,
        int *iradshield, double *avgsighi, double *avgsighei, double *avgsigheii,
        int *iradtrans, int *iradcoupled, int *iradstep, int *irt_honly,
        gr_float *kphHI, gr_float *kphHeI, gr_float *kphHeII, gr_float *kdissH2I,
        gr_float *photogamma, gr_float *xH2shield,
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
        int *iVheat, int *iMheat, gr_float *Vheat, gr_float *Mheat,
        int *iTfloor, gr_float *Tfloor_scalar, gr_float *Tfloor,
        int *iisrffield, gr_float* isrf_habing, 
        int *iH2shieldcustom, gr_float* f_shield_custom,
        int *itmax, int *exititmax);

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* do unit-handling */
  code_units units = determine_code_units(my_units, my_rates,
                                          my_fields->current_a_value,
                                          my_chemistry->unit_handling,
                                          "local_solve_chemistry");
  if (units.a_units < 0) {
    return FAIL;
  } else {
    my_units = &units;
  }

  /* Update UV background rates. */
  photo_rate_storage my_uvb_rates;

  my_uvb_rates.k24 = my_uvb_rates.k25 = my_uvb_rates.k26 =
    my_uvb_rates.k27 = my_uvb_rates.k28 = my_uvb_rates.k29 =
    my_uvb_rates.k30 = my_uvb_rates.k31 = my_uvb_rates.piHI =
    my_uvb_rates.piHeI = my_uvb_rates.piHeII = my_uvb_rates.crsHI =
    my_uvb_rates.crsHeI = my_uvb_rates.crsHeII =
    my_uvb_rates.comp_xray = my_uvb_rates.temp_xray = 0.;

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_rates,
                                  &my_uvb_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return FAIL;
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
  if (self_shielding_err_check(my_chemistry, my_fields,
                               "local_solve_chemistry") == FAIL) {
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  FORTRAN_NAME(solve_rate_cool_g)(
    &my_chemistry->with_radiative_cooling,
    my_fields->density,
    my_fields->internal_energy,
    my_fields->x_velocity,
    my_fields->y_velocity,
    my_fields->z_velocity,
    my_fields->e_density,
    my_fields->HI_density,
    my_fields->HII_density,
    my_fields->HeI_density,
    my_fields->HeII_density,
    my_fields->HeIII_density,
    my_fields->grid_dimension,
    my_fields->grid_dimension+1,
    my_fields->grid_dimension+2,
    &my_chemistry->NumberOfTemperatureBins,
    &my_units->comoving_coordinates,
    &my_chemistry->primordial_chemistry,
    &metal_field_present,
    &my_chemistry->metal_cooling,
    &my_chemistry->h2_on_dust,
    &my_chemistry->dust_chemistry,
    &my_chemistry->use_dust_density_field,
    &(my_fields->grid_rank),
    my_fields->grid_start,
    my_fields->grid_start+1,
    my_fields->grid_start+2,
    my_fields->grid_end,
    my_fields->grid_end+1,
    my_fields->grid_end+2,
    &my_chemistry->ih2co,
    &my_chemistry->ipiht,
    &my_chemistry->dust_recombination_cooling,
    &my_chemistry->photoelectric_heating,
    &(my_fields->grid_dx),
    &dt_value,
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
    &my_chemistry->DeuteriumToHydrogenRatio,
    &my_chemistry->SolarMetalFractionByMass,
    &my_chemistry->local_dust_to_gas_ratio,
    my_rates->k1,
    my_rates->k2,
    my_rates->k3,
    my_rates->k4,
    my_rates->k5,
    my_rates->k6,
    my_rates->k7,
    my_rates->k8,
    my_rates->k9,
    my_rates->k10,
    my_rates->k11,
    my_rates->k12,
    my_rates->k13,
    my_rates->k13dd,
    my_rates->k14,
    my_rates->k15,
    my_rates->k16,
    my_rates->k17,
    my_rates->k18,
    my_rates->k19,
    my_rates->k22,
    &my_uvb_rates.k24,
    &my_uvb_rates.k25,
    &my_uvb_rates.k26,
    &my_uvb_rates.k27,
    &my_uvb_rates.k28,
    &my_uvb_rates.k29,
    &my_uvb_rates.k30,
    &my_uvb_rates.k31,
    my_rates->k50,
    my_rates->k51,
    my_rates->k52,
    my_rates->k53,
    my_rates->k54,
    my_rates->k55,
    my_rates->k56,
    my_rates->k57,
    my_rates->k58,
    &my_chemistry->NumberOfDustTemperatureBins,
    &my_chemistry->DustTemperatureStart,
    &my_chemistry->DustTemperatureEnd,
    my_rates->h2dust,
    my_rates->n_cr_n,
    my_rates->n_cr_d1,
    my_rates->n_cr_d2,
    my_rates->ceHI,
    my_rates->ceHeI,
    my_rates->ceHeII,
    my_rates->ciHI,
    my_rates->ciHeI,
    my_rates->ciHeIS,
    my_rates->ciHeII,
    my_rates->reHII,
    my_rates->reHeII1,
    my_rates->reHeII2,
    my_rates->reHeIII,
    my_rates->brem,
    &my_rates->comp,
    &my_rates->gammah,
    &my_chemistry->interstellar_radiation_field,
    my_rates->regr,
    &my_rates->gamma_isrf,
    &my_uvb_rates.comp_xray,
    &my_uvb_rates.temp_xray,
    &my_uvb_rates.piHI,
    &my_uvb_rates.piHeI,
    &my_uvb_rates.piHeII,
    my_fields->HM_density,
    my_fields->H2I_density,
    my_fields->H2II_density,
    my_fields->DI_density,
    my_fields->DII_density,
    my_fields->HDI_density,
    my_fields->metal_density,
    my_fields->dust_density,
    my_rates->hyd01k,
    my_rates->h2k01,
    my_rates->vibh,
    my_rates->roth,
    my_rates->rotl,
    my_rates->GP99LowDensityLimit,
    my_rates->GP99HighDensityLimit,
    my_rates->HDlte,
    my_rates->HDlow,
    my_rates->GAHI,
    my_rates->GAH2,
    my_rates->GAHe,
    my_rates->GAHp,
    my_rates->GAel,
    my_rates->H2LTE,
    my_rates->gas_grain,
    &my_chemistry->H2_self_shielding,
    &my_chemistry->self_shielding_method,
    &my_uvb_rates.crsHI,
    &my_uvb_rates.crsHeI,
    &my_uvb_rates.crsHeII,
    &my_chemistry->use_radiative_transfer,
    &my_chemistry->radiative_transfer_coupled_rate_solver,
    &my_chemistry->radiative_transfer_intermediate_step,
    &my_chemistry->radiative_transfer_hydrogen_only,
    my_fields->RT_HI_ionization_rate,
    my_fields->RT_HeI_ionization_rate,
    my_fields->RT_HeII_ionization_rate,
    my_fields->RT_H2_dissociation_rate,
    my_fields->RT_heating_rate,
    my_fields-> H2_self_shielding_length,
    &ierr,
    &my_chemistry->h2_optical_depth_approximation,
    &my_chemistry->cie_cooling,
    &my_chemistry->three_body_rate,
    my_rates->cieco,
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
    my_fields->volumetric_heating_rate,
    my_fields->specific_heating_rate,
    &my_chemistry->use_temperature_floor,
    &my_chemistry->temperature_floor_scalar,
    my_fields->temperature_floor,
    &my_chemistry->use_isrf_field,
    my_fields->isrf_habing,
    &my_chemistry->H2_custom_shielding,
    my_fields->H2_custom_shielding_factor,
    &my_chemistry->max_iterations,
    &my_chemistry->exit_after_iterations_exceeded);

  if (ierr == FAIL) {
    fprintf(stderr, "Error in solve_rate_cool_g.\n");
  }

  return ierr;

}

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value)
{
  if (local_solve_chemistry(grackle_data, &grackle_rates,
                            my_units, my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}
