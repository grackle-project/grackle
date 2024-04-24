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
#include "state_dump_utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern int grackle_verbose;

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

void auto_show_config(FILE *fp);
void auto_show_flags(FILE *fp);
grackle_version get_grackle_version();

int _free_cloudy_data(cloudy_data *my_cloudy, chemistry_data *my_chemistry, int primordial);
int initialize_cloudy_data(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units,
                           int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates);

int local_free_chemistry_data(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

int initialize_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units,
                double co_length_units, double co_density_units);

void initialize_empty_UVBtable_struct(UVBtable *table);

/**
 * Initializes an empty #chemistry_data_storage struct with zeros and NULLs.
 */
void initialize_empty_chemistry_data_storage_struct(chemistry_data_storage *my_rates)
{
  my_rates->k1 = NULL;
  my_rates->k2 = NULL;
  my_rates->k3 = NULL;
  my_rates->k4 = NULL;
  my_rates->k5 = NULL;
  my_rates->k6 = NULL;
  my_rates->k7 = NULL;
  my_rates->k8 = NULL;
  my_rates->k9 = NULL;
  my_rates->k10 = NULL;
  my_rates->k11 = NULL;
  my_rates->k12 = NULL;
  my_rates->k13 = NULL;
  my_rates->k14 = NULL;
  my_rates->k15 = NULL;
  my_rates->k16 = NULL;
  my_rates->k17 = NULL;
  my_rates->k18 = NULL;
  my_rates->k19 = NULL;
  my_rates->k20 = NULL;
  my_rates->k21 = NULL;
  my_rates->k22 = NULL;
  my_rates->k23 = NULL;
  my_rates->k13dd = NULL;
  my_rates->k24 = 0.;
  my_rates->k25 = 0.;
  my_rates->k26 = 0.;
  my_rates->k27 = 0.;
  my_rates->k28 = 0.;
  my_rates->k29 = 0.;
  my_rates->k30 = 0.;
  my_rates->k31 = 0.;
  my_rates->k50 = NULL;
  my_rates->k51 = NULL;
  my_rates->k52 = NULL;
  my_rates->k53 = NULL;
  my_rates->k54 = NULL;
  my_rates->k55 = NULL;
  my_rates->k56 = NULL;
  my_rates->k57 = NULL;
  my_rates->k58 = NULL;
  my_rates->h2dust = NULL;
  my_rates->n_cr_n = NULL;
  my_rates->n_cr_d1 = NULL;
  my_rates->n_cr_d2 = NULL;
  my_rates->ceHI = NULL;
  my_rates->ceHeI = NULL;
  my_rates->ceHeII = NULL;
  my_rates->ciHI = NULL;
  my_rates->ciHeI = NULL;
  my_rates->ciHeIS = NULL;
  my_rates->ciHeII = NULL;
  my_rates->reHII = NULL;
  my_rates->reHeII1 = NULL;
  my_rates->reHeII2 = NULL;
  my_rates->reHeIII = NULL;
  my_rates->brem = NULL;
  my_rates->comp = 0.;
  my_rates->comp_xray = 0.;
  my_rates->temp_xray = 0.;
  my_rates->piHI = 0.;
  my_rates->piHeI = 0.;
  my_rates->piHeII = 0.;
  my_rates->crsHI = 0.;
  my_rates->crsHeI = 0.;
  my_rates->crsHeII = 0.;
  my_rates->hyd01k = NULL;
  my_rates->h2k01 = NULL;
  my_rates->vibh = NULL;
  my_rates->roth = NULL;
  my_rates->rotl = NULL;
  my_rates->GP99LowDensityLimit = NULL;
  my_rates->GP99HighDensityLimit = NULL;
  my_rates->GAHI = NULL;
  my_rates->GAH2 = NULL;
  my_rates->GAHe = NULL;
  my_rates->GAHp = NULL;
  my_rates->GAel = NULL;
  my_rates->H2LTE = NULL;
  my_rates->HDlte = NULL;
  my_rates->HDlow = NULL;
  my_rates->cieco = NULL;
  my_rates->gammah = 0.;
  my_rates->regr = NULL;
  my_rates->gamma_isrf = 0.;
  my_rates->gas_grain = NULL;
  my_rates->cloudy_data_new = -1;
}

int local_initialize_chemistry_data(chemistry_data *my_chemistry,
                                    chemistry_data_storage *my_rates,
                                    code_units *my_units)
{

  /* Better safe than sorry: Initialize everything to NULL/0 */
  initialize_empty_chemistry_data_storage_struct(my_rates);

  if (grackle_verbose) {
    grackle_version gversion = get_grackle_version();
    show_version_(stdout, &gversion);
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

    if (my_chemistry->dust_recombination_cooling < 0) {
      my_chemistry->dust_recombination_cooling = 1;
      if (grackle_verbose) {
        fprintf(stdout, "Dust chemistry enabled, setting dust_recombination_cooling to 1.\n");
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
# ifndef _OPENMP
  if (my_chemistry->omp_nthreads > 1) {
    fprintf(stdout,
            "omp_nthreads can't be set when Grackle isn't compiled with "
            "OPENMP\n");
    return FAIL;
  }
# else _OPENMP
  if (my_chemistry->omp_nthreads < 1) {
    // this is the default behavior (unless the user intervenes)
    my_chemistry->omp_nthreads = omp_get_max_threads();
  }
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

  //* Call initialise_rates to compute rate tables.
  initialize_rates(my_chemistry, my_rates, my_units, co_length_units, co_density_units);

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
  initialize_empty_UVBtable_struct(&(my_rates->UVbackground_table));
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
    grackle_version gversion = get_grackle_version();
    show_version_(fptr, &gversion);
    fprintf(fptr, "Grackle build options:\n");
    auto_show_config(fptr);
    fprintf(fptr, "Grackle build flags:\n");
    auto_show_flags(fptr);
    fprintf(fptr, "Grackle run-time parameters:\n");
    show_parameters_(fptr, my_chemistry);
    fclose(fptr);

    fprintf(stdout, "Grackle run-time parameters:\n");
    show_parameters_(stdout, my_chemistry);

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
  if (local_initialize_chemistry_data(grackle_data, &grackle_rates,
                                      my_units) == FAIL) {
    fprintf(stderr, "Error in local_initialize_chemistry_data.\n");
    return FAIL;
  }
  return SUCCESS;
}

int free_chemistry_data(void){
  if (local_free_chemistry_data(grackle_data, &grackle_rates) == FAIL) {
    fprintf(stderr, "Error in local_free_chemistry_data.\n");
    return FAIL;
  }
  return SUCCESS;
}

int local_free_chemistry_data(chemistry_data *my_chemistry,
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

  return SUCCESS;
}
