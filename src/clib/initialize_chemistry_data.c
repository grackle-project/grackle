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

#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

extern int grackle_verbose;

extern chemistry_data grackle_data;

void auto_show_config(FILE *fp);
void auto_show_flags(FILE *fp);
void auto_show_version(FILE *fp);

int initialize_cloudy_data(chemistry_data *my_chemistry,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units, double a_value,
                           int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry);

extern void FORTRAN_NAME(calc_rates_g)(
     int *ispecies,
     int *nratec, double *aye, double *temstart, double *temend, 
     int *casebrates, int *threebody,
     double *uxyz, double *uaye, double *urho, double *utim,
     double *ceHIa, double *ceHeIa, double *ceHeIIa, double *ciHIa, double *ciHeIa,
     double *ciHeISa, double *ciHeIIa, double *reHIIa, double *reHeII1a,
     double *reHeII2a, double *reHeIIIa, double *brema, double *compa, 
     double *gammahacgs, double *gammaha,
     double *hyd01ka, double *h2k01a, double *vibha, double *rotha, double *rotla,
     double *gpldl, double *gphdl, double *hdlte, double *hdlow, double *cieco,
     double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela, double *gasgr, 
     double *k1a, double *k2a, double *k3a, double *k4a, double *k5a, double *k6a,
     double *k7a, double *k8a, double *k9a, double *k10a,
     double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a,
     double *k15a, double *k16a, double *k17a, double *k18a,
     double *k19a, double *k20a, double *k21a, double *k22, double *k23,
     double *k50, double *k51, double *k52, double *k53, double *k54, double *k55,
     double *k56, int *ndratec, double *dtemstart, double *dtemend, double *h2dusta, 
     double *ncrca, double *ncrd1a, double *ncrd2a, 
     int *ioutput);

int _initialize_chemistry_data(chemistry_data *my_chemistry, 
                               code_units *my_units, double a_value)
{

  if (grackle_verbose) {
    FILE *fptr = fopen("GRACKLE_INFO", "w");
    auto_show_version(fptr);
    fprintf(fptr, "Grackle build options:\n");
    auto_show_config(fptr);
    fprintf(fptr, "Grackle build flags:\n");
    auto_show_flags(fptr);
    fclose(fptr);

    auto_show_version(stderr);
    fprintf(stderr, "Initializing grackle data.\n");
    fprintf(stderr, "with_radiative_cooling: %d.\n", my_chemistry->with_radiative_cooling);
    fprintf(stderr, "primordial_chemistry: %d.\n", my_chemistry->primordial_chemistry);
    fprintf(stderr, "metal_cooling: %d.\n", my_chemistry->metal_cooling);
    fprintf(stderr, "UVbackground: %d.\n", my_chemistry->UVbackground);
  }

  /* Only allow a units to be one with proper coordinates. */
  if (my_units->comoving_coordinates == FALSE && 
      my_units->a_units != 1.0) {
    fprintf(stderr, "ERROR: a_units must be 1.0 if comoving_coordinates is 0.\n");
    return FAIL;
  }

  /* Allocate CoolData space for rates. */

  if (my_chemistry->primordial_chemistry > 0) {
 
    my_chemistry->ceHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ceHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ceHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ciHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ciHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ciHeIS  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->ciHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->reHII   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->reHeII1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->reHeII2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->reHeIII = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->brem    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->hyd01k  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->h2k01   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->vibh    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->roth    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->rotl    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GP99LowDensityLimit  = malloc(my_chemistry->NumberOfTemperatureBins *
                                               sizeof(double));
    my_chemistry->GP99HighDensityLimit = malloc(my_chemistry->NumberOfTemperatureBins * 
                                               sizeof(double));
    my_chemistry->HDlte   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->HDlow   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->cieco   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GAHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GAH2    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GAHe    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GAHp    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->GAel    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->gas_grain = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

  /* Allocate space in my_chemistry for rates. */
 
    my_chemistry->k1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k3 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k4 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k5 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k6 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k7 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k8 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k9 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k10 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k11 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k12 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k13 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k13dd = malloc(7 * my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k14 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k15 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k16 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k17 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k18 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k19 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k20 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k21 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k22 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k23 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k50 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k51 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k52 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k53 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k54 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k55 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->k56 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->h2dust = malloc(my_chemistry->NumberOfTemperatureBins *
                                 my_chemistry->NumberOfDustTemperatureBins * sizeof(double));
    my_chemistry->n_cr_n = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->n_cr_d1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_chemistry->n_cr_d2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

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

  int ioutput = 0;

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

  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates_g)(
     &my_chemistry->primordial_chemistry,
     &my_chemistry->NumberOfTemperatureBins, &a_value, &my_chemistry->TemperatureStart,
        &my_chemistry->TemperatureEnd,
        &my_chemistry->CaseBRecombination, &my_chemistry->three_body_rate,
     &co_length_units, &my_units->a_units, 
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
        my_chemistry->HDlte, my_chemistry->HDlow, my_chemistry->cieco,
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
     &ioutput);

  /* Initialize Cloudy cooling. */
  int read_data;

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

int initialize_chemistry_data(code_units *my_units, double a_value)
{
  if (_initialize_chemistry_data(&grackle_data, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in _initialize_chemistry_data.\n");
    return FAIL;
  }
  return SUCCESS;
}
