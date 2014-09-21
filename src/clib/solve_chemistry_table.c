/***********************************************************************
/
/ Solve cooling (tabulated cooling function)
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

int _solve_chemistry_table(chemistry_data *my_chemistry,
                           code_units *my_units,
                           double a_value, double dt_value,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *metal_density)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->primordial_chemistry != 0) {
    fprintf(stderr, "ERROR: this function requires primordial_chemistry set to 0.\n");
    return FAIL;
  }

  gr_float *HI_density, *HII_density, *HM_density,
    *HeI_density, *HeII_density, *HeIII_density,
    *H2I_density, *H2II_density,
    *DI_density, *DII_density, *HDI_density,
    *e_density;

  HI_density = HII_density = HM_density =
    HeI_density = HeII_density = HeIII_density =
    H2I_density = H2II_density =
    DI_density = DII_density = HDI_density =
    e_density;

  if (solve_chemistry(my_units,
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
      fprintf(stderr, "Error in solve_chemistry.\n");
      return FAIL;
    }

  return SUCCESS;

}

int solve_chemistry_table(code_units *my_units,
                          double a_value, double dt_value,
                          int grid_rank, int *grid_dimension,
                          int *grid_start, int *grid_end,
                          gr_float *density, gr_float *internal_energy,
                          gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                          gr_float *metal_density)
{
  if (_solve_chemistry_table(&grackle_data,
                             my_units,
                             a_value, dt_value,
                             grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, internal_energy,
                             x_velocity, y_velocity, z_velocity,
                             metal_density) == FAIL) {
    fprintf(stderr, "Error in _solve_chemistry_table.\n");
    return FAIL;
  }
  return SUCCESS;
}

int solve_chemistry_table_(int *comoving_coordinates,
                           double *density_units, double *length_units,
                           double *time_units, double *velocity_units,
                           double *a_units, double *a_value, double *dt_value,
                           int *grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                           gr_float *metal_density)
{
  
  code_units my_units;
  my_units.comoving_coordinates = *comoving_coordinates;
  my_units.density_units = *density_units;
  my_units.length_units = *length_units;
  my_units.time_units = *time_units;
  my_units.velocity_units = *velocity_units;
  my_units.a_units = *a_units;

  int rval;

  rval = solve_chemistry_table(&my_units,
                               *a_value, *dt_value,
                               *grid_rank, grid_dimension,
                               grid_start, grid_end,
                               density, internal_energy,
                               x_velocity, y_velocity, z_velocity,
                               metal_density);

  return rval;

}
