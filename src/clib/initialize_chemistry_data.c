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
#include "grackle.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#if (GR_SUCCESS != SUCCESS) || (GR_FAIL != FAIL)
#error "Sanity check failure: GR_SUCCESS must be consistent with SUCCESS and GR_FAIL must be consistent with FAIL"
#endif

extern int grackle_verbose;

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

void auto_show_config(FILE *fp);
void auto_show_flags(FILE *fp);
grackle_version get_grackle_version();
void show_parameters(FILE *fp, chemistry_data *my_chemistry);
int _free_cloudy_data(cloudy_data *my_cloudy, chemistry_data *my_chemistry, int primordial);
int initialize_cloudy_data(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           cloudy_data *my_cloudy, char *group_name,
                           code_units *my_units,
                           int read_data);

int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates);

int local_free_chemistry_data(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

int local_free_metal_chemistry_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int local_free_dust_yields(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

int initialize_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, code_units *my_units,
                double co_length_units, double co_density_units);

void initialize_empty_UVBtable_struct(UVBtable *table);

static void show_version(FILE *fp)
{
  grackle_version gversion = get_grackle_version();
  fprintf (fp, "\n");
  fprintf (fp, "The Grackle Version %s\n", gversion.version);
  fprintf (fp, "Git Branch   %s\n", gversion.branch);
  fprintf (fp, "Git Revision %s\n", gversion.revision);
  fprintf (fp, "\n");
}

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

  my_rates->k24 = 0;
  my_rates->k25 = 0;
  my_rates->k26 = 0;

  my_rates->k27 = 0;
  my_rates->k28 = 0;
  my_rates->k29 = 0;
  my_rates->k30 = 0;
  my_rates->k31 = 0;

  my_rates->k50 = NULL;
  my_rates->k51 = NULL;
  my_rates->k52 = NULL;
  my_rates->k53 = NULL;
  my_rates->k54 = NULL;
  my_rates->k55 = NULL;
  my_rates->k56 = NULL;

  my_rates->k57 = NULL;
  my_rates->k58 = NULL;

  my_rates->k125 = NULL;
  my_rates->k129 = NULL;
  my_rates->k130 = NULL;
  my_rates->k131 = NULL;
  my_rates->k132 = NULL;
  my_rates->k133 = NULL;
  my_rates->k134 = NULL;
  my_rates->k135 = NULL;
  my_rates->k136 = NULL;
  my_rates->k137 = NULL;
  my_rates->k148 = NULL;
  my_rates->k149 = NULL;
  my_rates->k150 = NULL;
  my_rates->k151 = NULL;
  my_rates->k152 = NULL;
  my_rates->k153 = NULL;

  my_rates->kz15 = NULL;
  my_rates->kz16 = NULL;
  my_rates->kz17 = NULL;
  my_rates->kz18 = NULL;
  my_rates->kz19 = NULL;
  my_rates->kz20 = NULL;
  my_rates->kz21 = NULL;
  my_rates->kz22 = NULL;
  my_rates->kz23 = NULL;
  my_rates->kz24 = NULL;
  my_rates->kz25 = NULL;
  my_rates->kz26 = NULL;
  my_rates->kz27 = NULL;
  my_rates->kz28 = NULL;
  my_rates->kz29 = NULL;
  my_rates->kz30 = NULL;
  my_rates->kz31 = NULL;
  my_rates->kz32 = NULL;
  my_rates->kz33 = NULL;
  my_rates->kz34 = NULL;
  my_rates->kz35 = NULL;
  my_rates->kz36 = NULL;
  my_rates->kz37 = NULL;
  my_rates->kz38 = NULL;
  my_rates->kz39 = NULL;
  my_rates->kz40 = NULL;
  my_rates->kz41 = NULL;
  my_rates->kz42 = NULL;
  my_rates->kz43 = NULL;
  my_rates->kz44 = NULL;
  my_rates->kz45 = NULL;
  my_rates->kz46 = NULL;
  my_rates->kz47 = NULL;
  my_rates->kz48 = NULL;
  my_rates->kz49 = NULL;
  my_rates->kz50 = NULL;
  my_rates->kz51 = NULL;
  my_rates->kz52 = NULL;
  my_rates->kz53 = NULL;
  my_rates->kz54 = NULL;

  my_rates->h2dust = NULL;
  my_rates->h2dustS = NULL;
  my_rates->h2dustC = NULL;

  my_rates->grain_growth_rate = NULL;

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
  my_rates->comp = 0;
  my_rates->comp_xray = 0;
  my_rates->temp_xray = 0;

  my_rates->piHI = 0;
  my_rates->piHeI = 0;
  my_rates->piHeII = 0;

  my_rates->crsHI = 0;
  my_rates->crsHeI = 0;
  my_rates->crsHeII = 0;

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

  my_rates->gammah = 0;

  my_rates->regr = NULL;

  my_rates->gamma_isrf = 0;
  my_rates->gamma_isrf2 = 0;

  my_rates->gas_grain = NULL;
  my_rates->gas_grain2 = NULL;

  my_rates->cieY06 = NULL;

  my_rates->LH2_N = NULL;
  my_rates->LH2_Size = 0;
  my_rates->LH2_D = NULL;
  my_rates->LH2_T = NULL;
  my_rates->LH2_H = NULL;
  my_rates->LH2_dD = 0.;
  my_rates->LH2_dT = 0.;
  my_rates->LH2_dH = 0.;
  my_rates->LH2_L = NULL;
  my_rates->LHD_N = NULL;
  my_rates->LHD_Size = 0;
  my_rates->LHD_D = NULL;
  my_rates->LHD_T = NULL;
  my_rates->LHD_H = NULL;
  my_rates->LHD_dD = 0.;
  my_rates->LHD_dT = 0.;
  my_rates->LHD_dH = 0.;
  my_rates->LHD_L = NULL;

  my_rates->LCI_N = NULL;
  my_rates->LCI_Size = 0;
  my_rates->LCI_D = NULL;
  my_rates->LCI_T = NULL;
  my_rates->LCI_H = NULL;
  my_rates->LCI_dD = 0.;
  my_rates->LCI_dT = 0.;
  my_rates->LCI_dH = 0.;
  my_rates->LCI_L = NULL;
  my_rates->LCII_N = NULL;
  my_rates->LCII_Size = 0;
  my_rates->LCII_D = NULL;
  my_rates->LCII_T = NULL;
  my_rates->LCII_H = NULL;
  my_rates->LCII_dD = 0.;
  my_rates->LCII_dT = 0.;
  my_rates->LCII_dH = 0.;
  my_rates->LCII_L = NULL;
  my_rates->LOI_N = NULL;
  my_rates->LOI_Size = 0;
  my_rates->LOI_D = NULL;
  my_rates->LOI_T = NULL;
  my_rates->LOI_H = NULL;
  my_rates->LOI_dD = 0.;
  my_rates->LOI_dT = 0.;
  my_rates->LOI_dH = 0.;
  my_rates->LOI_L = NULL;

  my_rates->LCO_N = NULL;
  my_rates->LCO_Size = 0;
  my_rates->LCO_D = NULL;
  my_rates->LCO_T = NULL;
  my_rates->LCO_H = NULL;
  my_rates->LCO_dD = 0.;
  my_rates->LCO_dT = 0.;
  my_rates->LCO_dH = 0.;
  my_rates->LCO_L = NULL;
  my_rates->LOH_N = NULL;
  my_rates->LOH_Size = 0;
  my_rates->LOH_D = NULL;
  my_rates->LOH_T = NULL;
  my_rates->LOH_H = NULL;
  my_rates->LOH_dD = 0.;
  my_rates->LOH_dT = 0.;
  my_rates->LOH_dH = 0.;
  my_rates->LOH_L = NULL;
  my_rates->LH2O_N = NULL;
  my_rates->LH2O_Size = 0;
  my_rates->LH2O_D = NULL;
  my_rates->LH2O_T = NULL;
  my_rates->LH2O_H = NULL;
  my_rates->LH2O_dD = 0.;
  my_rates->LH2O_dT = 0.;
  my_rates->LH2O_dH = 0.;
  my_rates->LH2O_L = NULL;

  my_rates->alphap_N = NULL;
  my_rates->alphap_Size = 0;
  my_rates->alphap_D = NULL;
  my_rates->alphap_T = NULL;
  my_rates->alphap_dD = 0.;
  my_rates->alphap_dT;
  my_rates->alphap_Data;

  my_rates->gr_N = NULL;
  my_rates->gr_Size = 0;
  my_rates->gr_dT = 0.;
  my_rates->gr_Td = NULL;
  my_rates->SN0_N = 0;
  my_rates->SN0_XC = NULL;
  my_rates->SN0_XO = NULL;
  my_rates->SN0_XMg = NULL;
  my_rates->SN0_XAl = NULL;
  my_rates->SN0_XSi = NULL;
  my_rates->SN0_XS = NULL;
  my_rates->SN0_XFe = NULL;
  my_rates->SN0_fC = NULL;
  my_rates->SN0_fO = NULL;
  my_rates->SN0_fMg = NULL;
  my_rates->SN0_fAl = NULL;
  my_rates->SN0_fSi = NULL;
  my_rates->SN0_fS = NULL;
  my_rates->SN0_fFe = NULL;
  my_rates->SN0_fSiM = NULL;
  my_rates->SN0_fFeM = NULL;
  my_rates->SN0_fMg2SiO4 = NULL;
  my_rates->SN0_fMgSiO3 = NULL;
  my_rates->SN0_fFe3O4 = NULL;
  my_rates->SN0_fAC = NULL;
  my_rates->SN0_fSiO2D = NULL;
  my_rates->SN0_fMgO = NULL;
  my_rates->SN0_fFeS = NULL;
  my_rates->SN0_fAl2O3 = NULL;
  my_rates->SN0_freforg = NULL;
  my_rates->SN0_fvolorg = NULL;
  my_rates->SN0_fH2Oice = NULL;
  my_rates->SN0_r0SiM = NULL;
  my_rates->SN0_r0FeM = NULL;
  my_rates->SN0_r0Mg2SiO4 = NULL;
  my_rates->SN0_r0MgSiO3 = NULL;
  my_rates->SN0_r0Fe3O4 = NULL;
  my_rates->SN0_r0AC = NULL;
  my_rates->SN0_r0SiO2D = NULL;
  my_rates->SN0_r0MgO = NULL;
  my_rates->SN0_r0FeS = NULL;
  my_rates->SN0_r0Al2O3 = NULL;
  my_rates->SN0_r0reforg = NULL;
  my_rates->SN0_r0volorg = NULL;
  my_rates->SN0_r0H2Oice = NULL;
  my_rates->SN0_kpSiM = NULL;
  my_rates->SN0_kpFeM = NULL;
  my_rates->SN0_kpMg2SiO4 = NULL;
  my_rates->SN0_kpMgSiO3 = NULL;
  my_rates->SN0_kpFe3O4 = NULL;
  my_rates->SN0_kpAC = NULL;
  my_rates->SN0_kpSiO2D = NULL;
  my_rates->SN0_kpMgO = NULL;
  my_rates->SN0_kpFeS = NULL;
  my_rates->SN0_kpAl2O3 = NULL;
  my_rates->SN0_kpreforg = NULL;
  my_rates->SN0_kpvolorg = NULL;
  my_rates->SN0_kpH2Oice = NULL;

  my_rates->cloudy_data_new = -1;
}

int local_initialize_chemistry_data(chemistry_data *my_chemistry,
                                    chemistry_data_storage *my_rates,
                                    code_units *my_units)
{

  /* Better safe than sorry: Initialize everything to NULL/0 */
  initialize_empty_chemistry_data_storage_struct(my_rates);

  if (grackle_verbose) {
    show_version(stdout);
    fprintf(stdout, "Initializing grackle data.\n");
  }

  // Activate dust chemistry machinery.
  if (my_chemistry->dust_chemistry > 0) {

    if (my_chemistry->metal_cooling < 1) {
      fprintf(stderr, "ERROR: dust_chemistry > 0 requires metal_cooling > 0.\n");
      return GR_FAIL;
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
    return GR_FAIL;
  }
# else /* _OPENMP */
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
    return GR_FAIL;
  }

  // deal with my_chemistry->HydrogenFractionByMass
  // ==============================================
  // - set_default_chemistry_parameters & local_initialize_chemistry_parameters
  //   historically configured my_chemistry such that this member had a value
  //   of 0.76 (the default when primordial_chemistry >= 1). If
  //   primordial_chemistry == 0, this would always silently overwrite the
  //   value in fully tabulated mode. This led to cases where downstream users
  //   might mistakenly think they modified this parameter
  //
  // - New behavior: the member is now initialized to an undefined value
  //   - if the user hasn't changed it, we now set it to the appropriate
  //     default value
  //   - if the user has changed it, the precise handling depends on the choice
  //     of primordial_chemistry

  if (my_chemistry->HydrogenFractionByMass > 1) {
    fprintf(stderr, "ERROR: HydrogenFractionByMass cannot exceed 1.0\n");
    return GR_FAIL;
  } else if (my_chemistry->primordial_chemistry == 0) {
    /* In fully tabulated mode, set H mass fraction according to
       the abundances in Cloudy, which assumes n_He / n_H = 0.1.
       This gives a value of about 0.716. Using the other default value
       of 0.76 will result in negative electron densities at low
       temperature. Below, use X = 1 / (1 + m_He * n_He / n_H). */

    // we precomputed this value (in enhanced precision) rather than compute
    // it on the fly to allow users to directly specify this value (without
    // worrying about round-off issues)
    const double default_Hfrac = 0.715768377353088514;
    // TODO: at some point in the future, we should store the above value as
    //       an attribute of the HDF5 table and read it in directly.

    if (my_chemistry->HydrogenFractionByMass < 0) {
      my_chemistry->HydrogenFractionByMass = default_Hfrac;
    } else if (my_chemistry->HydrogenFractionByMass != default_Hfrac) {
      fprintf(stderr,
              "ERROR: Invalid HydrogenFractionByMass value is specified.\n"
              " -> when primordial_chemistry == 0, the allowed values are\n"
              "    are strictly enforced. You have 2 options: \n"
              "      1. leave the value unset or set it to a negative number\n"
              "         to have it calculated internally as %.18f\n"
              "      2. set the value exactly to the above number\n"
              " -> NOTE: for primordial_chemistry == 0, prior versions of\n"
              "    Grackle would silently overwrite the value of\n"
              "    HydrogenFractionByMass instead of reporting this error\n",
              default_Hfrac);
      return GR_FAIL;
    }
  } else {
    const double default_Hfrac = 0.76;
    if (my_chemistry->HydrogenFractionByMass < 0) {
      my_chemistry->HydrogenFractionByMass = default_Hfrac;
    }
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
                             "Primordial", my_units, read_data) == GR_FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return GR_FAIL;
  }

  /* Metal tables. */
  read_data = my_chemistry->metal_cooling == TRUE;
  if (initialize_cloudy_data(my_chemistry, my_rates,
                             &my_rates->cloudy_metal,
                             "Metals", my_units, read_data) == GR_FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.\n");
    return GR_FAIL;
  }

  /* Initialize UV Background data. */
  initialize_empty_UVBtable_struct(&(my_rates->UVbackground_table));
  if (initialize_UVbackground_data(my_chemistry, my_rates) == GR_FAIL) {
    fprintf(stderr, "Error in initialize_UVbackground_data.\n");
    return GR_FAIL;
  }

  /* store a copy of the initial units */
  my_rates->initial_units = *my_units;

  if (grackle_verbose) {
    time_t timer;
    char tstr[80];
    struct tm* tm_info;
    time(&timer);
    tm_info = localtime(&timer);
    strftime(tstr, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    FILE *fptr = fopen("GRACKLE_INFO", "w");
    fprintf(fptr, "%s\n", tstr);
    show_version(fptr);
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

  return GR_SUCCESS;
}

int initialize_chemistry_data(code_units *my_units)
{
  if (local_initialize_chemistry_data(grackle_data, &grackle_rates,
                                      my_units) == GR_FAIL) {
    fprintf(stderr, "Error in local_initialize_chemistry_data.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}

// Define helpers for the show_parameters function
// NOTE: it's okay that these functions all begin with an underscore since they
//       each have internal linkage (i.e. they are each declared static)
static void _show_field_INT(FILE *fp, const char* field, int val)
{ fprintf(fp, "%-33s = %d\n", field, val); }
static void _show_field_DOUBLE(FILE *fp, const char* field, double val)
{ fprintf(fp, "%-33s = %g\n", field, val); }
static void _show_field_STRING(FILE *fp, const char* field, const char* val)
{ fprintf(fp, "%-33s = %s\n", field, val); }

// this function writes each field of my_chemistry to fp
void show_parameters(FILE *fp, chemistry_data *my_chemistry){
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) \
    _show_field_ ## TYPE (fp, #FIELD, my_chemistry->FIELD);
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
}

int free_chemistry_data(void){
  if (local_free_chemistry_data(grackle_data, &grackle_rates) == GR_FAIL) {
    fprintf(stderr, "Error in local_free_chemistry_data.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
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
    GRACKLE_FREE(my_rates->gas_grain2);

    GRACKLE_FREE(my_rates->LH2_N);
    GRACKLE_FREE(my_rates->LH2_D);
    GRACKLE_FREE(my_rates->LH2_T);
    GRACKLE_FREE(my_rates->LH2_H);
    GRACKLE_FREE(my_rates->LH2_L);
    GRACKLE_FREE(my_rates->LHD_N);
    GRACKLE_FREE(my_rates->LHD_D);
    GRACKLE_FREE(my_rates->LHD_T);
    GRACKLE_FREE(my_rates->LHD_H);
    GRACKLE_FREE(my_rates->LHD_L);
    GRACKLE_FREE(my_rates->alphap_N);
    GRACKLE_FREE(my_rates->alphap_D);
    GRACKLE_FREE(my_rates->alphap_T);
    GRACKLE_FREE(my_rates->alphap_Data);

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
    GRACKLE_FREE(my_rates->grain_growth_rate);
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

  if (local_free_metal_chemistry_rates(my_chemistry, my_rates) == FAIL) {
    fprintf(stderr, "Error in local_free_metal_chemistry_rates.\n");
    return FAIL;
  }

  if (local_free_dust_yields(my_chemistry, my_rates) == FAIL) {
    fprintf(stderr, "Error in local_free_dust_yields.\n");
    return FAIL;
  }

  return GR_SUCCESS;
}
