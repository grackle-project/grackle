/***********************************************************************
/
/ Initialize chemistry and cooling rate data
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

extern chemistry_data grackle_data;

int initialize_cloudy_data(chemistry_data *my_chemistry,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units, gr_float a_value,
                           gr_int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry);

extern void FORTRAN_NAME(calc_rates_g)(
     gr_int *ispecies,
     gr_int *nratec, gr_float *aye, gr_float *temstart, gr_float *temend, 
     gr_int *casebrates, gr_int *threebody,
     gr_float *utem, gr_float *uxyz, gr_float *uaye, gr_float *urho, gr_float *utim,
     gr_float *ceHIa, gr_float *ceHeIa, gr_float *ceHeIIa, gr_float *ciHIa, gr_float *ciHeIa,
     gr_float *ciHeISa, gr_float *ciHeIIa, gr_float *reHIIa, gr_float *reHeII1a,
     gr_float *reHeII2a, gr_float *reHeIIIa, gr_float *brema, gr_float *compa, 
     gr_float *gammahacgs, gr_float *gammaha,
     gr_float *hyd01ka, gr_float *h2k01a, gr_float *vibha, gr_float *rotha, gr_float *rotla,
     gr_float *gpldl, gr_float *gphdl, gr_float *hdlte, gr_float *hdlow, gr_float *hdcool, gr_float *cieco,
     gr_float *gaHIa, gr_float *gaH2a, gr_float *gaHea, gr_float *gaHpa, gr_float *gaela, gr_float *gasgr, 
     gr_float *k1a, gr_float *k2a, gr_float *k3a, gr_float *k4a, gr_float *k5a, gr_float *k6a,
        gr_float *k7a, gr_float *k8a, gr_float *k9a, gr_float *k10a,
     gr_float *k11a, gr_float *k12a, gr_float *k13a, gr_float *k13dda, gr_float *k14a,
        gr_float *k15a, gr_float *k16a, gr_float *k17a, gr_float *k18a,
     gr_float *k19a, gr_float *k20a, gr_float *k21a, gr_float *k22, gr_float *k23,
     gr_float *k50, gr_float *k51, gr_float *k52, gr_float *k53, gr_float *k54, gr_float *k55,
        gr_float *k56, gr_int *ndratec, gr_float *dtemstart, gr_float *dtemend, gr_float *h2dusta, 
     gr_float *ncrca, gr_float *ncrd1a, gr_float *ncrd2a, 
     gr_float *mutab, gr_int *ioutput);

int _initialize_chemistry_data(chemistry_data *my_chemistry, 
                               code_units *my_units, gr_float a_value)
{

  fprintf(stderr, "Initializing chemistry data.\n");

  /* Only allow a units to be one with proper coordinates. */
  if (my_units->comoving_coordinates == FALSE && 
      my_units->a_units != 1.0) {
    fprintf(stderr, "ERROR: a_units must be 1.0 if comoving_coordinates is 0.\n");
    return FAIL;
  }

  /* Allocate CoolData space for rates. */

  if (my_chemistry->primordial_chemistry == 0) {

    my_chemistry->mu    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));

  }
  else {
 
    my_chemistry->ceHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ceHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ceHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ciHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ciHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ciHeIS  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->ciHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->reHII   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->reHeII1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->reHeII2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->reHeIII = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->brem    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->hyd01k  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->h2k01   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->vibh    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->roth    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->rotl    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GP99LowDensityLimit  = malloc(my_chemistry->NumberOfTemperatureBins *
                                               sizeof(gr_float));
    my_chemistry->GP99HighDensityLimit = malloc(my_chemistry->NumberOfTemperatureBins * 
                                               sizeof(gr_float));
    my_chemistry->HDlte   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->HDlow   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->HDcool  = malloc(5 * my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->cieco   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GAHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GAH2    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GAHe    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GAHp    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->GAel    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->gas_grain = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));

  /* Allocate space in my_chemistry for rates. */
 
    my_chemistry->k1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k3 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k4 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k5 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k6 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k7 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k8 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k9 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k10 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k11 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k12 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k13 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k13dd = malloc(7 * my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k14 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k15 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k16 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k17 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k18 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k19 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k20 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k21 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k22 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k23 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k50 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k51 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k52 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k53 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k54 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k55 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->k56 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->h2dust = malloc(my_chemistry->NumberOfTemperatureBins *
                                 my_chemistry->NumberOfDustTemperatureBins * sizeof(gr_float));
    my_chemistry->n_cr_n = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->n_cr_d1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));
    my_chemistry->n_cr_d2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(gr_float));

    my_chemistry->k24 = 0;
    my_chemistry->k25 = 0;
    my_chemistry->k26 = 0;
    my_chemistry->k27 = 0;
    my_chemistry->k28 = 0;
    my_chemistry->k29 = 0;
    my_chemistry->k30 = 0;
    my_chemistry->k31 = 0;
    my_chemistry->piHI = 0;
    my_chemistry->piHeII = 0;
    my_chemistry->piHeI = 0;

  }

  gr_int ioutput = 1;

  gr_float co_length_units, co_density_units;
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

  gr_float temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates_g)(
     &my_chemistry->primordial_chemistry,
     &my_chemistry->NumberOfTemperatureBins, &a_value, &my_chemistry->TemperatureStart,
        &my_chemistry->TemperatureEnd,
        &my_chemistry->CaseBRecombination, &my_chemistry->three_body_rate,
     &temperature_units, &co_length_units, &my_units->a_units, 
     &co_density_units, &my_units->time_units,
     my_chemistry->ceHI, my_chemistry->ceHeI, my_chemistry->ceHeII, my_chemistry->ciHI,
        my_chemistry->ciHeI,
     my_chemistry->ciHeIS, my_chemistry->ciHeII, my_chemistry->reHII,
        my_chemistry->reHeII1,
     my_chemistry->reHeII2, my_chemistry->reHeIII, my_chemistry->brem, &my_chemistry->comp, 
     &my_chemistry->photoelectric_heating_rate, &my_chemistry->gammah,
     my_chemistry->hyd01k, my_chemistry->h2k01, my_chemistry->vibh, my_chemistry->roth,
        my_chemistry->rotl,
     my_chemistry->GP99LowDensityLimit, my_chemistry->GP99HighDensityLimit,
        my_chemistry->HDlte, my_chemistry->HDlow, my_chemistry->HDcool, my_chemistry->cieco,
     my_chemistry->GAHI, my_chemistry->GAH2, my_chemistry->GAHe, my_chemistry->GAHp,
        my_chemistry->GAel, my_chemistry->gas_grain, 
     my_chemistry->k1, my_chemistry->k2, my_chemistry->k3, my_chemistry->k4, my_chemistry->k5,
        my_chemistry->k6, my_chemistry->k7, my_chemistry->k8, my_chemistry->k9, my_chemistry->k10,
     my_chemistry->k11, my_chemistry->k12, my_chemistry->k13, my_chemistry->k13dd, my_chemistry->k14,
        my_chemistry->k15, my_chemistry->k16, my_chemistry->k17, my_chemistry->k18,
     my_chemistry->k19, my_chemistry->k20, my_chemistry->k21, my_chemistry->k22, my_chemistry->k23,
     my_chemistry->k50, my_chemistry->k51, my_chemistry->k52, my_chemistry->k53, my_chemistry->k54,
        my_chemistry->k55, my_chemistry->k56, 
     &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, 
     &my_chemistry->DustTemperatureEnd, my_chemistry->h2dust, 
     my_chemistry->n_cr_n, my_chemistry->n_cr_d1, my_chemistry->n_cr_d2, 
     my_chemistry->mu, &ioutput);

  /* Initialize Cloudy cooling. */
  gr_int read_data;

  /* Primordial tables. */
  read_data = my_chemistry->primordial_chemistry == 0;
  if (initialize_cloudy_data(my_chemistry,
                             &my_chemistry->cloudy_primordial,
                             "Primordial",
                             my_units, a_value, read_data) == FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return FAIL;
  }

  /* Metal tables. */
  read_data = my_chemistry->metal_cooling == TRUE;
  if (initialize_cloudy_data(my_chemistry,
                             &my_chemistry->cloudy_metal,
                             "Metals",
                             my_units, a_value, read_data) == FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return FAIL;
  }

  /* Initialize UV Background data. */
  if (initialize_UVbackground_data(my_chemistry) == FAIL) {
    fprintf(stderr, "Error in initialize_UVbackground_data.\n");
    return FAIL;
  }

  return SUCCESS;
}

int initialize_chemistry_data(code_units *my_units, gr_float a_value)
{
  if (_initialize_chemistry_data(&grackle_data, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in _initialize_chemistry_data.\n");
    return FAIL;
  }
  return SUCCESS;
}
