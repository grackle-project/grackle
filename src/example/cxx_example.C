/***********************************************************************
/
/ Example executable using libgrackle
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
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

extern "C" {
#include <grackle.h>
}

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

int main(int argc, char *argv[])
{

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.velocity_units = my_units.length_units / my_units.time_units;
  my_units.a_units = 1.0; // units for the expansion factor

  // Second, create a chemistry object for parameters and rate data.
  if (set_default_chemistry_parameters() == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return 0;
  }
  // Set parameter values for chemistry.
  my_chemistry.use_grackle = 1;            // chemistry on
  my_chemistry.with_radiative_cooling = 1; // cooling on
  my_chemistry.primordial_chemistry = 3;   // molecular network with H, He, D
  my_chemistry.metal_cooling = 1;          // metal cooling on
  my_chemistry.UVbackground = 1;           // UV background on
  my_chemistry.grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  gr_float initial_redshift = 100.;
  gr_float a_value = 1. / (1. + initial_redshift);

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units, a_value) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return 0;
  }

  // Allocate field arrays.
  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity,
    *HI_density, *HII_density, *HM_density,
    *HeI_density, *HeII_density, *HeIII_density,
    *H2I_density, *H2II_density,
    *DI_density, *DII_density, *HDI_density,
    *e_density, *metal_density;
  gr_float tiny_number = 1.e-20;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  gr_int field_size = 10;
  gr_int grid_rank = 3;
  // If grid rank is less than 3, set the other dimensions, 
  // start indices, and end indices to 0.
  gr_int grid_dimension[3], grid_start[3], grid_end[3];
  for (int i = 0;i < 3;i++) {
    grid_dimension[i] = 0; // the active dimension not including ghost zones.
    grid_start[i] = 0;
    grid_end[i] = 0;
  }
  grid_dimension[0] = field_size;
  grid_end[0] = field_size - 1;

  density       = new gr_float[field_size];
  energy        = new gr_float[field_size];
  x_velocity    = new gr_float[field_size];
  y_velocity    = new gr_float[field_size];
  z_velocity    = new gr_float[field_size];
  // for primordial_chemistry >= 1
  HI_density    = new gr_float[field_size];
  HII_density   = new gr_float[field_size];
  HeI_density   = new gr_float[field_size];
  HeII_density  = new gr_float[field_size];
  HeIII_density = new gr_float[field_size];
  e_density     = new gr_float[field_size];
  // for primordial_chemistry >= 2
  HM_density    = new gr_float[field_size];
  H2I_density   = new gr_float[field_size];
  H2II_density  = new gr_float[field_size];
  // for primordial_chemistry >= 3
  DI_density    = new gr_float[field_size];
  DII_density   = new gr_float[field_size];
  HDI_density   = new gr_float[field_size];
  // for metal_cooling = 1
  metal_density = new gr_float[field_size];

  // set temperature units
  gr_float temperature_units =  mh * pow(my_units.a_units * 
                                         my_units.length_units /
                                         my_units.time_units, 2) / kboltz;

  int i;
  for (i = 0;i < field_size;i++) {
    density[i] = 1.0;
    HI_density[i] = my_chemistry.HydrogenFractionByMass * density[i];
    HII_density[i] = tiny_number * density[i];
    HM_density[i] = tiny_number * density[i];
    HeI_density[i] = (1.0 - my_chemistry.HydrogenFractionByMass) * density[i];
    HeII_density[i] = tiny_number * density[i];
    HeIII_density[i] = tiny_number * density[i];
    H2I_density[i] = tiny_number * density[i];
    H2II_density[i] = tiny_number * density[i];
    DI_density[i] = 2.0 * 3.4e-5 * density[i];
    DII_density[i] = tiny_number * density[i];
    HDI_density[i] = tiny_number * density[i];
    e_density[i] = tiny_number * density[i];
    // solar metallicity
    metal_density[i] = my_chemistry.SolarMetalFractionByMass * density[i];

    x_velocity[i] = 0.0;
    y_velocity[i] = 0.0;
    z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    energy[i] = 1000. / temperature_units;
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  gr_float dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units,
                      a_value, dt,
                      grid_rank, grid_dimension,
                      grid_start, grid_end,
                      density, energy,
                      x_velocity, y_velocity, z_velocity,
                      HI_density, HII_density, HM_density,
                      HeI_density, HeII_density, HeIII_density,
                      H2I_density, H2II_density,
                      DI_density, DII_density, HDI_density,
                      e_density, metal_density) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return 0;
  }

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(&my_units,
                             a_value,
                             grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, energy,
                             x_velocity, y_velocity, z_velocity,
                             HI_density, HII_density, HM_density,
                             HeI_density, HeII_density, HeIII_density,
                             H2I_density, H2II_density,
                             DI_density, DII_density, HDI_density,
                             e_density, metal_density, 
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return 0;
  }

  // Calculate temperature.
  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(&my_units,
                            grid_rank, grid_dimension,
                            density, energy,
                            HI_density, HII_density, HM_density,
                            HeI_density, HeII_density, HeIII_density,
                            H2I_density, H2II_density,
                            DI_density, DII_density, HDI_density,
                            e_density, metal_density, 
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return 0;
  }

  // Calculate pressure.
  gr_float *pressure;
  pressure = new gr_float[field_size];
  if (calculate_pressure(&my_units,
                         grid_rank, grid_dimension,
                         density, energy,
                         HI_density, HII_density, HM_density,
                         HeI_density, HeII_density, HeIII_density,
                         H2I_density, H2II_density,
                         DI_density, DII_density, HDI_density,
                         e_density, metal_density,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return 0;
  }

  // Calculate gamma.
  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(&my_units,
                      grid_rank, grid_dimension,
                      density, energy,
                      HI_density, HII_density, HM_density,
                      HeI_density, HeII_density, HeIII_density,
                      H2I_density, H2II_density,
                      DI_density, DII_density, HDI_density,
                      e_density, metal_density,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return 0;
  }

  return 1;
}
