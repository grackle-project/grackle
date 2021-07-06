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
#include <time.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern int grackle_verbose;

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

void auto_show_config(FILE *fp);
void auto_show_flags(FILE *fp);
void auto_show_version(FILE *fp);
void show_parameters(FILE *fp, chemistry_data *my_chemistry);
int calc_rates_metal(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units);
int calc_rates_dust(chemistry_data *my_chemistry,
                    chemistry_data_storage *my_rates,
                    code_units *my_units);
int _free_cloudy_data(cloudy_data *my_cloudy, chemistry_data *my_chemistry, int primordial);
int initialize_cloudy_data(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units,
                           int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates);

extern void FORTRAN_NAME(calc_rates_g)(
     int *ispecies, int *igammah, int *idust, int *idustall,
     int *nratec, double *aye, double *temstart, double *temend, 
     int *casebrates, int *threebody,
     double *uxyz, double *uaye, double *urho, double *utim,
     double *ceHIa, double *ceHeIa, double *ceHeIIa, double *ciHIa, double *ciHeIa,
     double *ciHeISa, double *ciHeIIa, double *reHIIa, double *reHeII1a,
     double *reHeII2a, double *reHeIIIa, double *brema, double *compa, 
     double *gammahacgs, double *gammaha, double *regra, double *gamma_isrfa,
     double *hyd01ka, double *h2k01a, double *vibha, double *rotha, double *rotla,
     double *gpldl, double *gphdl, double *hdlte, double *hdlow, double *cieco,
     double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela, 
     double *h2ltea, double *gasgra,
     double *k1a, double *k2a, double *k3a, double *k4a, double *k5a, double *k6a,
     double *k7a, double *k8a, double *k9a, double *k10a,
     double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a,
     double *k15a, double *k16a, double *k17a, double *k18a,
     double *k19a, double *k20a, double *k21a, double *k22, double *k23,
     double *k50, double *k51, double *k52, double *k53, double *k54, double *k55,
     double *k56, double *k57, double *k58, int *ndratec, double *dtemstart, 
     double *dtemend, double *h2dusta, double *ncrca, double *ncrd1a, double *ncrd2a, 
     int *ioutput
   , double *h2dustSa, double *h2dustCa, double *gasgr2a, double *gamma_isrf2a, double *grogra);

int _initialize_chemistry_data(chemistry_data *my_chemistry,
                               chemistry_data_storage *my_rates,
                               code_units *my_units)
{

  if (grackle_verbose) {
    auto_show_version(stdout);
    fprintf(stdout, "Initializing grackle data.\n");
  }

  // Activate dust chemistry machinery.
  if (my_chemistry->dust_chemistry > 0) {

    if (my_chemistry->metal_cooling < 1) {
      fprintf(stderr, "ERROR: dust_chemistry > 0 requires metal_cooling > 0.\n");
      return FAIL;
    }

    if (my_chemistry->photoelectric_heating < 0) {
      my_chemistry->photoelectric_heating = 2;
      if (grackle_verbose) {
        fprintf(stdout, "Dust chemistry enabled, setting photoelectric_heating to 2.\n");
      }
    }

    if (my_chemistry->primordial_chemistry > 1 &&
        my_chemistry->h2_on_dust == 0) {
      my_chemistry->h2_on_dust = 1;
      if (grackle_verbose) {
        fprintf(stdout, "Dust chemistry enabled, setting h2_on_dust to 1.\n");
      }
    }

  }

  // Default photo-electric heating to off if unset.
  if (my_chemistry->photoelectric_heating < 0) {
    my_chemistry->photoelectric_heating = 0;
  }

//initialize OpenMP
# ifdef _OPENMP
//number of threads
  omp_set_num_threads( my_chemistry->omp_nthreads );

//schedule
//const int chunk_size = -1;  // determined by default
  const int chunk_size = 1;

//omp_set_schedule( omp_sched_static,  chunk_size );
//omp_set_schedule( omp_sched_dynamic, chunk_size );
  omp_set_schedule( omp_sched_guided,  chunk_size );
//omp_set_schedule( omp_sched_auto,    chunk_size );
# endif

  /* Only allow a units to be one with proper coordinates. */
  if (my_units->comoving_coordinates == FALSE && 
      my_units->a_units != 1.0) {
    fprintf(stderr, "ERROR: a_units must be 1.0 if comoving_coordinates is 0.\n");
    return FAIL;
  }

  if (my_chemistry->primordial_chemistry == 0) {
    /* In fully tabulated mode, set H mass fraction according to
       the abundances in Cloudy, which assumes n_He / n_H = 0.1.
       This gives a value of about 0.716. Using the default value
       of 0.76 will result in negative electron densities at low
       temperature. Below, we set X = 1 / (1 + m_He * n_He / n_H). */
    my_chemistry->HydrogenFractionByMass = 1. / (1. + 0.1 * 3.971);
  }

  /* Allocate CoolData space for rates. */

  if (my_chemistry->primordial_chemistry > 0) {
 
    my_rates->ceHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ceHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ceHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ciHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ciHeI   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ciHeIS  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->ciHeII  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->reHII   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->reHeII1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->reHeII2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->reHeIII = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->brem    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->hyd01k  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->h2k01   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->vibh    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->roth    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->rotl    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GP99LowDensityLimit  = malloc(my_chemistry->NumberOfTemperatureBins *
                                            sizeof(double));
    my_rates->GP99HighDensityLimit = malloc(my_chemistry->NumberOfTemperatureBins * 
                                            sizeof(double));
    my_rates->HDlte   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->HDlow   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->cieco   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GAHI    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GAH2    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GAHe    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GAHp    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->GAel    = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->H2LTE   =  malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

  /* Allocate space in my_rates for rates. */
 
    my_rates->k1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k3 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k4 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k5 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k6 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k7 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k8 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k9 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k10 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k11 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k12 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k13 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k13dd = malloc(14 * my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k14 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k15 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k16 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k17 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k18 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k19 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k20 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k21 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k22 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k23 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k50 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k51 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k52 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k53 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k54 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k55 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k56 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k57 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k58 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->k125 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k129 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k130 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k131 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k132 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k133 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k134 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k135 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k136 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k137 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k148 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k149 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k150 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k151 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k152 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k153 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->kz15 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz16 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz17 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz18 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz19 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz20 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz21 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz22 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz23 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz24 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz25 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz26 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz27 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz28 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz29 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz30 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz31 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz32 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz33 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz34 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz35 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz36 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz37 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz38 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz39 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz40 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz41 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz42 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz43 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz44 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz45 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz46 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz47 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz48 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz49 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz50 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz51 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz52 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz53 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz54 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->cieY06  = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->h2dust = malloc(my_chemistry->NumberOfTemperatureBins *
                              my_chemistry->NumberOfDustTemperatureBins * sizeof(double));
    my_rates->n_cr_n = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->n_cr_d1 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->n_cr_d2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->h2dustS = malloc(my_chemistry->NumberOfTemperatureBins *
                      my_chemistry->NumberOfDustTemperatureBins * sizeof(double));
    my_rates->h2dustC = malloc(my_chemistry->NumberOfTemperatureBins *
                              my_chemistry->NumberOfDustTemperatureBins * sizeof(double));
    my_rates->grogr   = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->k24 = 0;
    my_rates->k25 = 0;
    my_rates->k26 = 0;
    my_rates->k27 = 0;
    my_rates->k28 = 0;
    my_rates->k29 = 0;
    my_rates->k30 = 0;
    my_rates->k31 = 0;
    my_rates->piHI = 0;
    my_rates->piHeII = 0;
    my_rates->piHeI = 0;
    my_rates->crsHI = 0;
    my_rates->crsHeI = 0;
    my_rates->crsHeII = 0;
    my_rates->comp_xray = 0;
    my_rates->temp_xray = 0;

  }

  if (my_chemistry->h2_on_dust > 0 || my_chemistry->dust_chemistry > 0) {
    my_rates->gas_grain = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->regr      = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->gas_grain2 = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
  }

  int ioutput = 0;

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

  /* Call FORTRAN routine to do the hard work. */

  FORTRAN_NAME(calc_rates_g)(
     &my_chemistry->primordial_chemistry, &my_chemistry->photoelectric_heating,
     &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry,
     &my_chemistry->NumberOfTemperatureBins, &my_units->a_value,
     &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
     &my_chemistry->CaseBRecombination, &my_chemistry->three_body_rate,
     &co_length_units, &my_units->a_units, 
     &co_density_units, &my_units->time_units,
     my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI,
        my_rates->ciHeI,
     my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII,
        my_rates->reHeII1,
     my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, 
     &my_chemistry->photoelectric_heating_rate, &my_rates->gammah, my_rates->regr,
     &my_rates->gamma_isrf,
     my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth,
        my_rates->rotl,
     my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit,
        my_rates->HDlte, my_rates->HDlow, my_rates->cieco,
     my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp,
     my_rates->GAel, my_rates->H2LTE, my_rates->gas_grain,
     my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5,
        my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
     my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14,
        my_rates->k15, my_rates->k16, my_rates->k17, my_rates->k18,
     my_rates->k19, my_rates->k20, my_rates->k21, my_rates->k22, my_rates->k23,
     my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54,
        my_rates->k55, my_rates->k56, my_rates->k57, my_rates->k58,
     &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, 
     &my_chemistry->DustTemperatureEnd, my_rates->h2dust, 
     my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2, 
     &ioutput
   , my_rates->h2dustS, my_rates->h2dustC, my_rates->gas_grain2, &my_rates->gamma_isrf2, my_rates->grogr);

  /* Metal chemistry rates */
  if (calc_rates_metal(my_chemistry, my_rates, my_units) == FAIL) {
    fprintf(stderr, "Error in calc_rates_metal.\n");
    return FAIL;
  }
  /* Dust rates */
  if (calc_rates_dust(my_chemistry, my_rates, my_units) == FAIL) {
    fprintf(stderr, "Error in calc_rates_dust.\n");
    return FAIL;
  }

  /* Initialize Cloudy cooling. */
  my_rates->cloudy_data_new = 1;
  int read_data;

  /* Primordial tables. */
  read_data = my_chemistry->primordial_chemistry == 0;
  if (initialize_cloudy_data(my_chemistry, my_rates,
                             &my_rates->cloudy_primordial,
                             "Primordial", my_units, read_data) == FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return FAIL;
  }

  /* Metal tables. */
  read_data = my_chemistry->metal_cooling == TRUE;
  if (initialize_cloudy_data(my_chemistry, my_rates,
                             &my_rates->cloudy_metal,
                             "Metals", my_units, read_data) == FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return FAIL;
  }

  /* Initialize UV Background data. */
  my_rates->UVbackground_table.Nz     = 0;
  my_rates->UVbackground_table.z      = NULL;
  my_rates->UVbackground_table.k24    = NULL;
  my_rates->UVbackground_table.k25    = NULL;
  my_rates->UVbackground_table.k26    = NULL;
  my_rates->UVbackground_table.k27    = NULL;
  my_rates->UVbackground_table.k28    = NULL;
  my_rates->UVbackground_table.k29    = NULL;
  my_rates->UVbackground_table.k30    = NULL;
  my_rates->UVbackground_table.k31    = NULL;
  my_rates->UVbackground_table.piHI   = NULL;
  my_rates->UVbackground_table.piHeII = NULL;
  my_rates->UVbackground_table.piHeI  = NULL;
  my_rates->UVbackground_table.crsHI  = NULL;
  my_rates->UVbackground_table.crsHeII = NULL;
  my_rates->UVbackground_table.crsHeI = NULL;

  if (initialize_UVbackground_data(my_chemistry, my_rates) == FAIL) {
    fprintf(stderr, "Error in initialize_UVbackground_data.\n");
    return FAIL;
  }

  if (grackle_verbose) {
    time_t timer;
    char tstr[80];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(tstr, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    FILE *fptr = fopen("GRACKLE_INFO", "w");
    fprintf(fptr, "%s\n", tstr);
    auto_show_version(fptr);
    fprintf(fptr, "Grackle build options:\n");
    auto_show_config(fptr);
    fprintf(fptr, "Grackle build flags:\n");
    auto_show_flags(fptr);
    fprintf(fptr, "Grackle run-time parameters:\n");
    show_parameters(fptr, my_chemistry);
    fclose(fptr);

    fprintf(stdout, "Grackle run-time parameters:\n");
    show_parameters(stdout, my_chemistry);

#   ifdef _OPENMP
    int omp_nthread, omp_chunk_size;
    omp_sched_t omp_schedule;

    omp_get_schedule( &omp_schedule, &omp_chunk_size );
#   pragma omp parallel
#   pragma omp master
    { omp_nthread = omp_get_num_threads(); }

    fprintf( stdout, "OpenMP: on\n" );
    fprintf( stdout, "  num_threads: %d\n", omp_nthread );
    fprintf( stdout, "  schedule: %s\n", ( omp_schedule == omp_sched_static  ) ? "static"  :
                                         ( omp_schedule == omp_sched_dynamic ) ? "dynamic" :
                                         ( omp_schedule == omp_sched_guided  ) ? "guided"  :
                                         ( omp_schedule == omp_sched_auto    ) ? "auto"    :
                                                                                 "unknown" );
    fprintf( stdout, "  chunk size: %d\n", omp_chunk_size );
#   else
    fprintf( stdout, "OpenMP: off\n" );
#   endif
  }

  return SUCCESS;
}

int initialize_chemistry_data(code_units *my_units)
{
  if (_initialize_chemistry_data(grackle_data, &grackle_rates,
                                 my_units) == FAIL) {
    fprintf(stderr, "Error in _initialize_chemistry_data.\n");
    return FAIL;
  }
  return SUCCESS;
}

void show_parameters(FILE *fp, chemistry_data *my_chemistry)
{
  fprintf(fp, "use_grackle                       = %d\n",
          my_chemistry->use_grackle);
  fprintf(fp, "with_radiative_cooling            = %d\n",
          my_chemistry->with_radiative_cooling);
  fprintf(fp, "primordial_chemistry              = %d\n",
          my_chemistry->primordial_chemistry);
  fprintf(fp, "dust_chemistry                    = %d\n",
          my_chemistry->dust_chemistry);
  fprintf(fp, "metal_cooling                     = %d\n",
          my_chemistry->metal_cooling);
  fprintf(fp, "UVbackground                      = %d\n",
          my_chemistry->UVbackground);
  fprintf(fp, "grackle_data_file                 = %s\n",
          my_chemistry->grackle_data_file);
  fprintf(fp, "cmb_temperature_floor             = %d\n",
          my_chemistry->cmb_temperature_floor);
  fprintf(fp, "Gamma                             = %g\n",
          my_chemistry->Gamma);
  fprintf(fp, "h2_on_dust                        = %d\n",
          my_chemistry->h2_on_dust);
  fprintf(fp, "use_dust_density_field            = %d\n",
          my_chemistry->use_dust_density_field);
  fprintf(fp, "metal_chemistry                   = %d\n",
          my_chemistry->metal_chemistry);
  fprintf(fp, "multi_metals                      = %d\n",
          my_chemistry->multi_metals);
  fprintf(fp, "metal_abundances                  = %d\n",
          my_chemistry->metal_abundances);
  fprintf(fp, "dust_species                      = %d\n",
          my_chemistry->dust_species);
  fprintf(fp, "dust_temperature_multi            = %d\n",
          my_chemistry->dust_temperature_multi);
  fprintf(fp, "dust_sublimation                  = %d\n",
          my_chemistry->dust_sublimation);
  fprintf(fp, "grain_growth                      = %d\n",
          my_chemistry->grain_growth);
  fprintf(fp, "photoelectric_heating             = %d\n",
          my_chemistry->photoelectric_heating);
  fprintf(fp, "photoelectric_heating_rate        = %g\n",
          my_chemistry->photoelectric_heating_rate);
  fprintf(fp, "use_isrf_field                    = %d\n",
          my_chemistry->use_isrf_field);
  fprintf(fp, "interstellar_radiation_field      = %g\n",
          my_chemistry->interstellar_radiation_field);
  fprintf(fp, "use_volumetric_heating_rate       = %d\n",
          my_chemistry->use_volumetric_heating_rate);
  fprintf(fp, "use_specific_heating_rate         = %d\n",
          my_chemistry->use_specific_heating_rate);
  fprintf(fp, "three_body_rate                   = %d\n",
          my_chemistry->three_body_rate);
  fprintf(fp, "cie_cooling                       = %d\n",
          my_chemistry->cie_cooling);
  fprintf(fp, "h2_optical_depth_approximation    = %d\n",
          my_chemistry->h2_optical_depth_approximation);
  fprintf(fp, "ih2co                             = %d\n",
          my_chemistry->ih2co);
  fprintf(fp, "ipiht                             = %d\n",
          my_chemistry->ipiht);
  fprintf(fp, "HydrogenFractionByMass            = %g\n",
          my_chemistry->HydrogenFractionByMass);
  fprintf(fp, "DeuteriumToHydrogenRatio          = %g\n",
          my_chemistry->DeuteriumToHydrogenRatio);
  fprintf(fp, "SolarMetalFractionByMass          = %g\n",
          my_chemistry->SolarMetalFractionByMass);
  fprintf(fp, "local_dust_to_gas_ratio           = %g\n",
          my_chemistry->local_dust_to_gas_ratio);
  fprintf(fp, "NumberOfTemperatureBins           = %d\n",
          my_chemistry->NumberOfTemperatureBins);
  fprintf(fp, "CaseBRecombination                = %d\n",
          my_chemistry->CaseBRecombination);
  fprintf(fp, "TemperatureStart                  = %g\n",
          my_chemistry->TemperatureStart);
  fprintf(fp, "TemperatureEnd                    = %g\n",
          my_chemistry->TemperatureEnd);
  fprintf(fp, "NumberOfDustTemperatureBins       = %d\n",
          my_chemistry->NumberOfDustTemperatureBins);
  fprintf(fp, "DustTemperatureStart              = %g\n",
          my_chemistry->DustTemperatureStart);
  fprintf(fp, "DustTemperatureEnd                = %g\n",
          my_chemistry->DustTemperatureEnd);
  fprintf(fp, "Compton_xray_heating              = %d\n",
          my_chemistry->Compton_xray_heating);
  fprintf(fp, "LWbackground_sawtooth_suppression = %d\n",
          my_chemistry->LWbackground_sawtooth_suppression);
  fprintf(fp, "LWbackground_intensity            = %g\n",
          my_chemistry->LWbackground_intensity);
  fprintf(fp, "UVbackground_redshift_on          = %g\n",
          my_chemistry->UVbackground_redshift_on);
  fprintf(fp, "UVbackground_redshift_off         = %g\n",
          my_chemistry->UVbackground_redshift_off);
  fprintf(fp, "UVbackground_redshift_fullon      = %g\n",
          my_chemistry->UVbackground_redshift_fullon);
  fprintf(fp, "UVbackground_redshift_drop        = %g\n",
          my_chemistry->UVbackground_redshift_drop);
  fprintf(fp, "cloudy_electron_fraction_factor   = %g\n",
          my_chemistry->cloudy_electron_fraction_factor);
  fprintf(fp, "use_radiative_transfer            = %d\n",
          my_chemistry->use_radiative_transfer);
  fprintf(fp, "radiative_transfer_coupled_rate_solver = %d\n",
          my_chemistry->radiative_transfer_coupled_rate_solver);
  fprintf(fp, "radiative_transfer_intermediate_step = %d\n",
          my_chemistry->radiative_transfer_intermediate_step);
  fprintf(fp, "radiative_transfer_hydrogen_only  = %d\n",
          my_chemistry->radiative_transfer_hydrogen_only);
  fprintf(fp, "self_shielding_method             = %d\n",
          my_chemistry->self_shielding_method);
  fprintf(fp, "H2_self_shielding                 = %d\n",
          my_chemistry->H2_self_shielding);
  fprintf(fp, "radiative_transfer_H2II_diss      = %d\n",
          my_chemistry->radiative_transfer_H2II_diss);
  fprintf(fp, "radiative_transfer_HDI_diss       = %d\n",
          my_chemistry->radiative_transfer_HDI_diss);
  fprintf(fp, "radiative_transfer_metal_ion      = %d\n",
          my_chemistry->radiative_transfer_metal_ion);
  fprintf(fp, "radiative_transfer_metal_diss     = %d\n",
          my_chemistry->radiative_transfer_metal_diss);

# ifdef _OPENMP
  fprintf(fp, "omp_nthreads                      = %d\n",
          my_chemistry->omp_nthreads);
# endif
}


int _free_chemistry_data(chemistry_data *my_chemistry,
			 chemistry_data_storage *my_rates) {
  if (my_chemistry->primordial_chemistry > 0) {
    GRACKLE_FREE(my_rates->ceHI);
    GRACKLE_FREE(my_rates->ceHeI);
    GRACKLE_FREE(my_rates->ceHeII);
    GRACKLE_FREE(my_rates->ciHI);
    GRACKLE_FREE(my_rates->ciHeI);
    GRACKLE_FREE(my_rates->ciHeIS);
    GRACKLE_FREE(my_rates->ciHeII);
    GRACKLE_FREE(my_rates->reHII);
    GRACKLE_FREE(my_rates->reHeII1);
    GRACKLE_FREE(my_rates->reHeII2);
    GRACKLE_FREE(my_rates->reHeIII);
    GRACKLE_FREE(my_rates->brem);
    GRACKLE_FREE(my_rates->hyd01k);
    GRACKLE_FREE(my_rates->h2k01);
    GRACKLE_FREE(my_rates->vibh);
    GRACKLE_FREE(my_rates->roth);
    GRACKLE_FREE(my_rates->rotl);
    GRACKLE_FREE(my_rates->GP99LowDensityLimit);
    GRACKLE_FREE(my_rates->GP99HighDensityLimit);

    GRACKLE_FREE(my_rates->HDlte);
    GRACKLE_FREE(my_rates->HDlow);
    GRACKLE_FREE(my_rates->cieco);
    GRACKLE_FREE(my_rates->GAHI);
    GRACKLE_FREE(my_rates->GAH2);
    GRACKLE_FREE(my_rates->GAHe);
    GRACKLE_FREE(my_rates->GAHp);
    GRACKLE_FREE(my_rates->GAel);
    GRACKLE_FREE(my_rates->H2LTE);
    GRACKLE_FREE(my_rates->gas_grain);
    GRACKLE_FREE(my_rates->gas_grain2);

    GRACKLE_FREE(my_rates->k1);
    GRACKLE_FREE(my_rates->k2);
    GRACKLE_FREE(my_rates->k3);
    GRACKLE_FREE(my_rates->k4);
    GRACKLE_FREE(my_rates->k5);
    GRACKLE_FREE(my_rates->k6);
    GRACKLE_FREE(my_rates->k7);
    GRACKLE_FREE(my_rates->k8);
    GRACKLE_FREE(my_rates->k9);
    GRACKLE_FREE(my_rates->k10);
    GRACKLE_FREE(my_rates->k11);
    GRACKLE_FREE(my_rates->k12);
    GRACKLE_FREE(my_rates->k13);
    GRACKLE_FREE(my_rates->k13dd);
    GRACKLE_FREE(my_rates->k14);
    GRACKLE_FREE(my_rates->k15);
    GRACKLE_FREE(my_rates->k16);
    GRACKLE_FREE(my_rates->k17);
    GRACKLE_FREE(my_rates->k18);
    GRACKLE_FREE(my_rates->k19);
    GRACKLE_FREE(my_rates->k20);
    GRACKLE_FREE(my_rates->k21);
    GRACKLE_FREE(my_rates->k22);
    GRACKLE_FREE(my_rates->k23);
    GRACKLE_FREE(my_rates->k50);
    GRACKLE_FREE(my_rates->k51);
    GRACKLE_FREE(my_rates->k52);
    GRACKLE_FREE(my_rates->k53);
    GRACKLE_FREE(my_rates->k54);
    GRACKLE_FREE(my_rates->k55);
    GRACKLE_FREE(my_rates->k56);
    GRACKLE_FREE(my_rates->k57);
    GRACKLE_FREE(my_rates->k58);
    GRACKLE_FREE(my_rates->h2dust);
    GRACKLE_FREE(my_rates->n_cr_n);
    GRACKLE_FREE(my_rates->n_cr_d1);
    GRACKLE_FREE(my_rates->n_cr_d2);
    GRACKLE_FREE(my_rates->h2dustS);
    GRACKLE_FREE(my_rates->h2dustC);
    GRACKLE_FREE(my_rates->grogr);
  }



  _free_cloudy_data(&my_rates->cloudy_primordial, my_chemistry, /* primordial */ 1);
  _free_cloudy_data(&my_rates->cloudy_metal, my_chemistry, /* primordial */ 0);
  
  GRACKLE_FREE(my_rates->UVbackground_table.z);
  GRACKLE_FREE(my_rates->UVbackground_table.k24);
  GRACKLE_FREE(my_rates->UVbackground_table.k25);
  GRACKLE_FREE(my_rates->UVbackground_table.k26);

  if (my_chemistry->primordial_chemistry > 1) {
    GRACKLE_FREE(my_rates->UVbackground_table.k27);
    GRACKLE_FREE(my_rates->UVbackground_table.k28);
    GRACKLE_FREE(my_rates->UVbackground_table.k29);
    GRACKLE_FREE(my_rates->UVbackground_table.k30);
    GRACKLE_FREE(my_rates->UVbackground_table.k31);
  }

  GRACKLE_FREE(my_rates->UVbackground_table.piHI);
  GRACKLE_FREE(my_rates->UVbackground_table.piHeII);
  GRACKLE_FREE(my_rates->UVbackground_table.piHeI);

  if (my_chemistry->self_shielding_method > 0){    
    GRACKLE_FREE(my_rates->UVbackground_table.crsHI);
    GRACKLE_FREE(my_rates->UVbackground_table.crsHeII);
    GRACKLE_FREE(my_rates->UVbackground_table.crsHeI);
  }

}
