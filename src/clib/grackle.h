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

#include "grackle_types.h"
#include "code_units.h"
#include "chemistry_data.h"

extern int grackle_verbose;

extern chemistry_data grackle_data;

int set_default_chemistry_parameters();

chemistry_data _set_default_chemistry_parameters();

int initialize_chemistry_data(code_units *my_units, double a_value);

int _initialize_chemistry_data(chemistry_data *my_chemistry, 
                               code_units *my_units, double a_value);

int initialize_grackle(int comoving_coordinates,
                       double density_units, double length_units,
                       double time_units, double velocity_units,
                       double a_units, double a_value,
                       int use_grackle, int with_radiative_cooling,
                       char *grackle_data_file,
                       int primordial_chemistry, int metal_cooling,
                       int UVbackground, int h2_on_dust,
                       int cmb_temperature_floor, double gamma);

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
                    gr_float *e_density, gr_float *metal_density);

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
                     gr_float *e_density, gr_float *metal_density);

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
                     gr_float *e_density, gr_float *metal_density);

int calculate_cooling_time(code_units *my_units, double a_value,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                           gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                           gr_float *e_density, gr_float *metal_density,
                           gr_float *cooling_time);

int _calculate_cooling_time(chemistry_data *my_chemistry,
                            code_units *my_units, double a_value,
                            int grid_rank, int *grid_dimension,
                            int *grid_start, int *grid_end,
                            gr_float *density, gr_float *internal_energy,
                            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                            gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                            gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                            gr_float *H2I_density, gr_float *H2II_density,
                            gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                            gr_float *e_density, gr_float *metal_density,
                            gr_float *cooling_time);

int calculate_cooling_time_(int *comoving_coordinates,
                            double *density_units, double *length_units,
                            double *time_units, double *velocity_units,
                            double *a_units, double *a_value,
                            int *grid_rank, int *grid_dimension,
                            int *grid_start, int *grid_end,
                            gr_float *density, gr_float *internal_energy,
                            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                            gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                            gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                            gr_float *H2I_density, gr_float *H2II_density,
                            gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                            gr_float *e_density, gr_float *metal_density,
                            gr_float *cooling_time);

int calculate_gamma(code_units *my_units, double a_value,
                    int grid_rank, int *grid_dimension,
                    int *grid_start, int *grid_end,
                    gr_float *density, gr_float *internal_energy,
                    gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                    gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                    gr_float *H2I_density, gr_float *H2II_density,
                    gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                    gr_float *e_density, gr_float *metal_density,
                    gr_float *my_gamma);

int _calculate_gamma(chemistry_data *my_chemistry,
                     code_units *my_units, double a_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *my_gamma);

int calculate_gamma_(int *comoving_coordinates,
                     double *density_units, double *length_units,
                     double *time_units, double *velocity_units,
                     double *a_units, double *a_value,
                     int *grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *e_density, gr_float *metal_density,
                     gr_float *my_gamma);

int calculate_pressure(code_units *my_units,
                       int grid_rank, int *grid_dimension,
                       gr_float *density, gr_float *internal_energy,
                       gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                       gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                       gr_float *H2I_density, gr_float *H2II_density,
                       gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                       gr_float *e_density, gr_float *metal_density,
                       gr_float *pressure);

int _calculate_pressure(chemistry_data *my_chemistry,
                        code_units *my_units,
                        int grid_rank, int *grid_dimension,
                        gr_float *density, gr_float *internal_energy,
                        gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                        gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                        gr_float *H2I_density, gr_float *H2II_density,
                        gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                        gr_float *e_density, gr_float *metal_density,
                        gr_float *pressure);

int calculate_pressure_(int *comoving_coordinates,
                        double *density_units, double *length_units,
                        double *time_units, double *velocity_units,
                        double *a_units,
                        int *grid_rank, int *grid_dimension,
                        gr_float *density, gr_float *internal_energy,
                        gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                        gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                        gr_float *H2I_density, gr_float *H2II_density,
                        gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                        gr_float *e_density, gr_float *metal_density,
                        gr_float *pressure);

int calculate_temperature(code_units *my_units, double a_value,
                          int grid_rank, int *grid_dimension,
                          int *grid_start, int *grid_end,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                          gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                          gr_float *H2I_density, gr_float *H2II_density,
                          gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                          gr_float *e_density, gr_float *metal_density,
                          gr_float *temperature);

int _calculate_temperature(chemistry_data *my_chemistry,
                           code_units *my_units, double a_value,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                           gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                           gr_float *e_density, gr_float *metal_density,
                           gr_float *temperature);

int calculate_temperature_(int *comoving_coordinates,
                           double *density_units, double *length_units,
                           double *time_units, double *velocity_units,
                           double *a_units, double *a_value,
                           int *grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                           gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                           gr_float *e_density, gr_float *metal_density,
                           gr_float *temperature);

/* Tabular-only functions. */

int solve_chemistry_table(code_units *my_units,
                          double a_value, double dt_value,
                          int grid_rank, int *grid_dimension,
                          int *grid_start, int *grid_end,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                          gr_float *metal_density);

int _solve_chemistry_table(chemistry_data *my_chemistry,
                           code_units *my_units,
                           double a_value, double dt_value,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *metal_density);

int solve_chemistry_table_(int *comoving_coordinates,
                           double *density_units, double *length_units,
                           double *time_units, double *velocity_units,
                           double *a_units, double *a_value, double *dt_value,
                           int *grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *metal_density);

int calculate_cooling_time_table(code_units *my_units, double a_value,
                                 int grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                                 gr_float *metal_density,
                                 gr_float *cooling_time);

int _calculate_cooling_time_table(chemistry_data *my_chemistry,
                                  code_units *my_units, double a_value,
                                  int grid_rank, int *grid_dimension,
                                  int *grid_start, int *grid_end,
                                  gr_float *density, gr_float *internal_energy,
                                  gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                                  gr_float *metal_density,
                                  gr_float *cooling_time);

int calculate_cooling_time_table_(int *comoving_coordinates,
                                  double *density_units, double *length_units,
                                  double *time_units, double *velocity_units,
                                  double *a_units, double *a_value,
                                  int *grid_rank, int *grid_dimension,
                                  int *grid_start, int *grid_end,
                                  gr_float *density, gr_float *internal_energy,
                                  gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                                  gr_float *metal_density,
                                  gr_float *cooling_time);

int calculate_pressure_table(code_units *my_units,
                             int grid_rank, int *grid_dimension,
                             gr_float *density, gr_float *internal_energy,
                             gr_float *pressure);

int _calculate_pressure_table(chemistry_data *my_chemistry,
                              code_units *my_units,
                              int grid_rank, int *grid_dimension,
                              gr_float *density, gr_float *internal_energy,
                              gr_float *pressure);

int calculate_pressure_table_(int *comoving_coordinates,
                              double *density_units, double *length_units,
                              double *time_units, double *velocity_units,
                              double *a_units,
                              int *grid_rank, int *grid_dimension,
                              gr_float *density, gr_float *internal_energy,
                              gr_float *pressure);

int calculate_temperature_table(code_units *my_units, double a_value,
                                int grid_rank, int *grid_dimension,
                                int *grid_start, int *grid_end,
                                gr_float *density, gr_float *internal_energy,
                                gr_float *metal_density,
                                gr_float *temperature);

int _calculate_temperature_table(chemistry_data *my_chemistry,
                                 code_units *my_units, double a_value,
                                 int grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *metal_density,
                                 gr_float *temperature);

int calculate_temperature_table_(int *comoving_coordinates,
                                 double *density_units, double *length_units,
                                 double *time_units, double *velocity_units,
                                 double *a_units, double *a_value,
                                 int *grid_rank, int *grid_dimension,
                                 int *grid_start, int *grid_end,
                                 gr_float *density, gr_float *internal_energy,
                                 gr_float *metal_density,
                                 gr_float *temperature);

#endif
