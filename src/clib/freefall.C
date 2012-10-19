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
#include "global_data.h"
#include "grackle.h"

//  ENZO Main Program

#ifdef SHARED_LIBRARY
#define MAIN_NAME freefall_main
#else
#define MAIN_NAME main
#endif

Eint32 MAIN_NAME(Eint32 argc, char *argv[])
{

  chemistry_data my_chemistry = set_default_chemistry_parameters();
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

  gr_float energy_units = POW((my_units.length_units / my_units.time_units), 2.0);

  gr_float gravitational_constant = 4.0 * 3.1415926 * 6.6726e-8 * 
    my_units.density_units * POW(my_units.time_units, 2);

  gr_float a_value = 1.0;

  if (initialize_chemistry_data(my_chemistry, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return FAIL;
  }

  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity;
  gr_float *HI_density, *HII_density, *HM_density,
    *HeI_density, *HeII_density, *HeIII_density,
    *H2I_density, *H2II_density,
    *DI_density, *DII_density, *HDI_density,
    *e_density, *metal_density,
    *cooling_time, *temperature, *gamma;
  gr_float tiny_number = 1.e-20;

  gr_int my_size = 1;
  density = new gr_float[my_size];
  energy = new gr_float[my_size];
  x_velocity = new gr_float[my_size];
  y_velocity = new gr_float[my_size];
  z_velocity = new gr_float[my_size];
  HI_density = new gr_float[my_size];
  HII_density = new gr_float[my_size];
  HM_density = new gr_float[my_size];
  HeI_density = new gr_float[my_size];
  HeII_density = new gr_float[my_size];
  HeIII_density = new gr_float[my_size];
  H2I_density = new gr_float[my_size];
  H2II_density = new gr_float[my_size];
  DI_density = new gr_float[my_size];
  DII_density = new gr_float[my_size];
  HDI_density = new gr_float[my_size];
  e_density = new gr_float[my_size];
  metal_density = new gr_float[my_size];
  cooling_time = new gr_float[my_size];
  temperature = new gr_float[my_size];
  gamma = new gr_float[my_size];

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

  gr_float temperature_units = mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;

  energy[0] = 1000. / temperature_units;
  x_velocity[0] = 0.0;
  y_velocity[0] = 0.0;
  z_velocity[0] = 0.0;
  
  gr_float my_time = 0.0;
  gr_float freefall_constant = POW(density[0], -0.5);
  gr_float freefall_time_constant = POW(((32 * gravitational_constant) / (3 * pi)), 0.5);
  gr_float dt, density_ratio, timestep_fraction;
  timestep_fraction = 0.1;

  gr_int grid_rank = 3;
  gr_int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i = 0;i < 3;i++) {
    grid_dimension[i] = 1;
    grid_start[i] = 0;
    grid_end[i] = 0;
  }

  while (density[0] < 1.e10) {

    dt = timestep_fraction * 
      POW(((3 * pi) / (32 * gravitational_constant * density[0])), 0.5);

    if (calculate_cooling_time(my_chemistry, my_units,
        		       a_value, dt,
        		       grid_rank, grid_dimension,
        		       grid_start, grid_end,
        		       density, energy,
        		       x_velocity, y_velocity, z_velocity,
        		       HI_density, HII_density, HM_density,
        		       HeI_density, HeII_density, HeIII_density,
        		       H2I_density, H2II_density,
        		       DI_density, DII_density, HDI_density,
        		       e_density, metal_density, 
        		       cooling_time) == FAIL) {
      fprintf(stderr, "Error in calculate_cooling_time.\n");
      return FAIL;
    }

    if (calculate_temperature(my_chemistry, my_units,
                              grid_rank, grid_dimension,
                              density, energy,
                              HI_density, HII_density, HM_density,
                              HeI_density, HeII_density, HeIII_density,
                              H2I_density, H2II_density,
                              DI_density, DII_density, HDI_density,
                              e_density, metal_density, 
                              temperature) == FAIL) {
      fprintf(stderr, "Error in calculate_temperature.\n");
      return FAIL;
    }

    fprintf(stderr, "%"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
            my_time, dt, density[0], (temperature[0]),
	    (cooling_time[0]*my_units.time_units));

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

 
