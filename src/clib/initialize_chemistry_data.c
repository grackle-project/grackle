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

int _free_cloudy_data(cloudy_data *my_cloudy, chemistry_data *my_chemistry, int primordial);
int initialize_cloudy_data(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units,
                           int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates);


int initialise_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units, 
                double co_length_units, double co_density_units);

int writeRates(char language[50], chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

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

  if (my_chemistry->h2_on_dust > 0 || my_chemistry->dust_chemistry > 0) {
    my_rates->gas_grain = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->regr      = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
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
  initialise_rates(my_chemistry, my_rates, my_units, co_length_units, co_density_units);

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

}
