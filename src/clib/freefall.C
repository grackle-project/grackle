/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       August 12th 2006
/              May 13th 2008
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/    then calls the appropriate routines depending on whether we are
/    doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
 
#define DEFINE_STORAGE
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"

//  ENZO Main Program

#ifdef SHARED_LIBRARY
#define MAIN_NAME freefall_main
#else
#define MAIN_NAME main
#endif

int set_default_chemistry_parameters(chemistry_data &my_chemistry);
int initialize_chemistry_data(chemistry_data &my_chemistry,
                              code_units &my_units, float a_value);

int solve_chemistry(chemistry_data &my_chemistry,
                    code_units &my_units,
                    float a_value, float dt_value,
                    int grid_rank, int *grid_dimension,
                    int *grid_start, int *grid_end,
                    float *density, float *internal_energy,
                    float *x_velocity, float *y_velocity, float  *z_velocity,
                    float *HI_density, float *HII_density, float *HM_density,
                    float *HeI_density, float *HeII_density, float *HeIII_density,
                    float *H2I_density, float *H2II_density,
                    float *DI_density, float *DII_density, float *HD_density,
                    float *e_density, float *metal_density);

Eint32 MAIN_NAME(Eint32 argc, char *argv[])
{

  chemistry_data my_chemistry;
  if (set_default_chemistry_parameters(my_chemistry) == FAIL) {
    fprintf(stderr, "Error in set_default_chemistry_data.\n");
    return FAIL;
  }
  my_chemistry.use_chemistry = 1;
  my_chemistry.primordial_chemistry = 2;
  my_chemistry.metal_cooling = 1;
  my_chemistry.cloudy_table_file = (char*) "solar_2008_3D_metals.h5";

  code_units my_units;
  my_units.comoving_coordinates = 0;
  my_units.density_units = mh;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.a_units = 1.0;

  float energy_units = POW((my_units.length_units / my_units.time_units), 2.0);

  float gravitational_constant = 4.0 * 3.1415926 * 6.6726e-8 * 
    my_units.density_units * POW(my_units.time_units, 2);

  float a_value = 1.0;

  if (initialize_chemistry_data(my_chemistry, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return FAIL;
  }

  float *density, *energy, *x_velocity, *y_velocity, *z_velocity;
  float *HI_density, *HII_density, *HM_density,
    *HeI_density, *HeII_density, *HeIII_density,
    *H2I_density, *H2II_density,
    *DI_density, *DII_density, *HDI_density,
    *e_density, *metal_density;
  float tiny_number = 1.e-20;

  int my_size = 1;
  density = new float[my_size];
  energy = new float[my_size];
  x_velocity = new float[my_size];
  y_velocity = new float[my_size];
  z_velocity = new float[my_size];
  HI_density = new float[my_size];
  HII_density = new float[my_size];
  HM_density = new float[my_size];
  HeI_density = new float[my_size];
  HeII_density = new float[my_size];
  HeIII_density = new float[my_size];
  H2I_density = new float[my_size];
  H2II_density = new float[my_size];
  DI_density = new float[my_size];
  DII_density = new float[my_size];
  HDI_density = new float[my_size];
  e_density = new float[my_size];
  metal_density = new float[my_size];

  density[0] = 1.0;
  HI_density[0] = my_chemistry.HydrogenFractionByMass * density[0];
  HII_density[0] = tiny_number * density[0];
  HM_density[0] = tiny_number * density[0];
  HeI_density[0] = (1.0 - my_chemistry.HydrogenFractionByMass) * density[0];
  HeII_density[0] = tiny_number * density[0];
  HeIII_density[0] = tiny_number * density[0];
  H2I_density[0] = tiny_number * density[0];
  H2II_density[0] = tiny_number * density[0];
  DI_density[0] = 2.0 * 3.4e-5 * density[0];
  DII_density[0] = tiny_number * density[0];
  HDI_density[0] = tiny_number * density[0];
  e_density[0] = tiny_number * density[0];
  metal_density[0] = 1.e-5 * density[0];

  float temperature_units = mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;

  energy[0] = 1000. / temperature_units;
  x_velocity[0] = 0.0;
  y_velocity[0] = 0.0;
  z_velocity[0] = 0.0;
  
  float my_time = 0.0;
  float freefall_constant = POW(density[0], -0.5);
  float freefall_time_constant = POW(((32 * gravitational_constant) / (3 * pi)), 0.5);
  float dt, density_ratio, timestep_fraction;
  timestep_fraction = 0.1;

  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i = 0;i < 3;i++) {
    grid_dimension[i] = 1;
    grid_start[i] = 0;
    grid_end[i] = 0;
  }

  while (density[0] < 1.e10) {

    dt = timestep_fraction * 
      POW(((3 * pi) / (32 * gravitational_constant * density[0])), 0.5);

    fprintf(stderr, "%"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
            my_time, dt, density[0], (energy[0]*temperature_units));

    density_ratio = POW((freefall_constant - 
                         (0.5 * freefall_time_constant * my_time)), -2.) / density[0];
    
    density[0] *= density_ratio;
    HI_density[0] *= density_ratio;
    HII_density[0] *= density_ratio;
    HM_density[0] *= density_ratio;
    HeI_density[0] *= density_ratio;
    HeII_density[0] *= density_ratio;
    HeIII_density[0] *= density_ratio;
    H2I_density[0] *= density_ratio;
    H2II_density[0] *= density_ratio;
    DI_density[0] *= density_ratio;
    DII_density[0] *= density_ratio;
    HDI_density[0] *= density_ratio;
    e_density[0] *= density_ratio;
    metal_density[0] *= density_ratio;

    energy[0] += (my_chemistry.Gamma - 1) * energy[0] * 
      freefall_time_constant * POW(density[0], 0.5) * dt;

    if (solve_chemistry(my_chemistry, my_units,
                        a_value, dt,
                        grid_rank, grid_dimension,
                        grid_start, grid_end,
                        density, energy,
                        x_velocity, y_velocity, z_velocity,
                        HI_density, HII_density, HM_density,
                        HeI_density, HeII_density, HeIII_density,
                        H2I_density, H2II_density,
                        DI_density, DII_density, HDI_density,
                        e_density, metal_density) == FAIL) {
      fprintf(stderr, "Error in solve_chemistry.\n");
      return FAIL;
    }

    my_time += dt;

  }

  return SUCCESS;
 
 }

 
