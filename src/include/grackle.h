/***********************************************************************
/
/ Grackle function prototypes
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_H__
#define __GRACKLE_H__

// this logic must occur before including Grackle's headers
#ifndef GRIMPL_PUBLIC_INCLUDE
  #define GRIMPL_PUBLIC_INCLUDE 0
#elif GRIMPL_PUBLIC_INCLUDE == 1
  // we temporarily allow this case for internal use within the Grackle lib to
  // avoid merge conflicts (external Grackle users shouldn't depend on this!)
#else
  #error "it is an error for GRIMPL_PUBLIC_INCLUDE to be defined"
#endif /* GRIMPL_PUBLIC_INCLUDE */


// include the private headers
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// standard error codes returned by the grackle functions
#define GR_SUCCESS 1
#define GR_FAIL 0

#define GR_SPECIFY_INITIAL_A_VALUE -1

extern int grackle_verbose;

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

double get_velocity_units(const code_units *my_units);

void set_velocity_units(code_units *my_units);

double get_temperature_units(const code_units *my_units);

double gr_query_units(const chemistry_data_storage * my_rates,
                      const char* units_name, double current_a_value);

int set_default_chemistry_parameters(chemistry_data *my_grackle);

int local_initialize_chemistry_parameters(chemistry_data *my_chemistry);

int initialize_chemistry_data(code_units *my_units);

int local_initialize_chemistry_data(chemistry_data *my_chemistry, 
                                    chemistry_data_storage *my_rates,
                                    code_units *my_units);

int* local_chemistry_data_access_int(chemistry_data* my_chemistry,
                                     const char* param_name);
double* local_chemistry_data_access_double(chemistry_data* my_chemistry,
                                           const char* param_name);
char** local_chemistry_data_access_string(chemistry_data* my_chemistry,
                                          const char* param_name);

const char* param_name_int(unsigned int i);
const char* param_name_double(unsigned int i);
const char* param_name_string(unsigned int i);

unsigned int grackle_num_params(const char* type_name);

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value);

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value);

int calculate_cooling_time(code_units *my_units,
                           grackle_field_data *my_fields,
                           gr_float *cooling_time);

int local_calculate_cooling_time(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units,
                                 grackle_field_data *my_fields,
                                 gr_float *cooling_time);

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature);

int local_calculate_dust_temperature(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *dust_temperature);

int calculate_gamma(code_units *my_units,
                    grackle_field_data *my_fields,
                    gr_float *my_gamma);

int local_calculate_gamma(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *my_gamma);

int calculate_pressure(code_units *my_units,
                       grackle_field_data *my_fields,
                       gr_float *pressure);

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure);

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature);

int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature);

int free_chemistry_data(void);

int local_free_chemistry_data(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

grackle_version get_grackle_version(void);

int gr_initialize_field_data(grackle_field_data *my_fields);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */


// this is important for letting us diagnose improper use of private headrs
#if GRIMPL_PUBLIC_INCLUDE == 0
  #undef GRIMPL_PUBLIC_INCLUDE
#endif

#endif /* __GRACKLE_H__ */
