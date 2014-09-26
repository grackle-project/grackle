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
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"
#include "fortran.def"

extern chemistry_data grackle_data;

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              code_units *my_units, double a_value);

extern void FORTRAN_NAME(solve_rate_cool_g)(
        int *icool,
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand, 
        int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
        int *ih2co, int *ipiht, int *igammah,
	double *dt, double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *dtoh, double *z_solar,
	double *k1a, double *k2a, double *k3a, double *k4a, double *k5a, 
	double *k6a, double *k7a, double *k8a, double *k9a, double *k10a,
	double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a, 
	double *k15a, double *k16a, double *k17a, double *k18a, double *k19a, 
        double *k22a,	double *k24, double *k25, double *k26, double *k27, 
        double *k28, double *k29, double *k30, double *k31,
	double *k50a, double *k51a, double *k52a, double *k53a, double *k54a,
	double *k55a, double *k56a,
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
	double *gasgra,
	int *ierr,
	int *ih2optical, int *iciecool, int *ithreebody, double *ciecoa,
 	int *icmbTfloor, int *iClHeat,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3, 
 	long long *priDataSize, double *priCooling, double *priHeating,
        long long *metGridRank, long long *metGridDim,
 	double *metPar1, double *metPar2, double *metPar3, 
 	long long *metDataSize, double *metCooling, double *metHeating,
        double *mutaba);

int _solve_chemistry(chemistry_data *my_chemistry,
                     code_units *my_units,
                     double a_value, double dt_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_units, a_value) == FAIL) {
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
      a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(a_value * my_units->a_units, 3);
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
    &dt_value, &a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
    &temperature_units, &co_length_units, &my_units->a_units, 
    &co_density_units, &my_units->time_units, &my_chemistry->Gamma,
    &my_chemistry->HydrogenFractionByMass, &my_chemistry->DeuteriumToHydrogenRatio,
    &my_chemistry->SolarMetalFractionByMass,
    my_chemistry->k1, my_chemistry->k2, my_chemistry->k3, my_chemistry->k4, my_chemistry->k5, 
    my_chemistry->k6, my_chemistry->k7, my_chemistry->k8, my_chemistry->k9, my_chemistry->k10,
    my_chemistry->k11, my_chemistry->k12, my_chemistry->k13, my_chemistry->k13dd, 
    my_chemistry->k14, my_chemistry->k15, my_chemistry->k16,
    my_chemistry->k17, my_chemistry->k18, my_chemistry->k19, my_chemistry->k22,
    &my_chemistry->k24, &my_chemistry->k25, &my_chemistry->k26, &my_chemistry->k27,
    &my_chemistry->k28, &my_chemistry->k29, &my_chemistry->k30, &my_chemistry->k31,
    my_chemistry->k50, my_chemistry->k51, my_chemistry->k52, my_chemistry->k53,
    my_chemistry->k54, my_chemistry->k55, my_chemistry->k56,
    &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, 
    &my_chemistry->DustTemperatureEnd, my_chemistry->h2dust, 
    my_chemistry->n_cr_n, my_chemistry->n_cr_d1, my_chemistry->n_cr_d2,
    my_chemistry->ceHI, my_chemistry->ceHeI, my_chemistry->ceHeII, my_chemistry->ciHI,
    my_chemistry->ciHeI, my_chemistry->ciHeIS, my_chemistry->ciHeII, my_chemistry->reHII, 
    my_chemistry->reHeII1, my_chemistry->reHeII2, my_chemistry->reHeIII, my_chemistry->brem, 
    &my_chemistry->comp, &my_chemistry->gammah,
    &my_chemistry->comp_xray, &my_chemistry->temp_xray,
    &my_chemistry->piHI, &my_chemistry->piHeI, &my_chemistry->piHeII,
    HM_density, H2I_density, H2II_density,
    DI_density, DII_density, HDI_density, metal_density,
    my_chemistry->hyd01k, my_chemistry->h2k01, my_chemistry->vibh, 
    my_chemistry->roth, my_chemistry->rotl,
    my_chemistry->GP99LowDensityLimit, my_chemistry->GP99HighDensityLimit, 
    my_chemistry->HDlte, my_chemistry->HDlow,
    my_chemistry->GAHI, my_chemistry->GAH2, my_chemistry->GAHe, my_chemistry->GAHp,
    my_chemistry->GAel, my_chemistry->gas_grain,
    &ierr,
    &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, 
    &my_chemistry->three_body_rate, my_chemistry->cieco,
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
    &my_chemistry->cloudy_metal.grid_rank,
    my_chemistry->cloudy_metal.grid_dimension,
    my_chemistry->cloudy_metal.grid_parameters[0],
    my_chemistry->cloudy_metal.grid_parameters[1],
    my_chemistry->cloudy_metal.grid_parameters[2],
    &my_chemistry->cloudy_metal.data_size,
    my_chemistry->cloudy_metal.cooling_data, 
    my_chemistry->cloudy_metal.heating_data,
    my_chemistry->mu);

  return SUCCESS;

}

int solve_chemistry(code_units *my_units,
                    double a_value, double dt_value,
                    int grid_rank, int *grid_dimension,
                    int *grid_start, int *grid_end,
                    gr_float *density, gr_float *internal_energy,
                    gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                    gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                    gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                    gr_float *H2I_density, gr_float *H2II_density,
                    gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                    gr_float *e_density, gr_float *metal_density)
{
  if (_solve_chemistry(&grackle_data,
                       my_units,
                       a_value, dt_value,
                       grid_rank, grid_dimension,
                       grid_start, grid_end,
                       density, internal_energy,
                       x_velocity, y_velocity, z_velocity,
                       HI_density, HII_density, HM_density,
                       HeI_density, HeII_density, HeIII_density,
                       H2I_density, H2II_density,
                       DI_density, DII_density, HDI_density,
                       e_density, metal_density) == FAIL) {
    fprintf(stderr, "Error in _solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}

int solve_chemistry_(int *comoving_coordinates,
                     double *density_units, double *length_units,
                     double *time_units, double *velocity_units,
                     double *a_units, double *a_value, double *dt_value,
                     int *grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density)
{

  code_units my_units;
  my_units.comoving_coordinates = *comoving_coordinates;
  my_units.density_units = *density_units;
  my_units.length_units = *length_units;
  my_units.time_units = *time_units;
  my_units.velocity_units = *velocity_units;
  my_units.a_units = *a_units;

  int rval;
  rval = solve_chemistry(&my_units,
                         *a_value, *dt_value,
                         *grid_rank, grid_dimension,
                         grid_start, grid_end,
                         density, internal_energy,
                         x_velocity, y_velocity, z_velocity,
                         HI_density, HII_density, HM_density,
                         HeI_density, HeII_density, HeIII_density,
                         H2I_density, H2II_density,
                         DI_density, DII_density, HDI_density,
                         e_density, metal_density);
  return rval;

}
