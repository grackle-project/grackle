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

#include <grackle.h>

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

int main(int argc, char *argv[])
{

  // Enable output
  grackle_verbose = 1;

  // Set parameters

  // chemistry on
  int use_grackle = 1;
  // cooling on
  int with_radiative_cooling = 1;
  // fully tabular cooling
  int primordial_chemistry = 0;
  // metal cooling on
  int metal_cooling = 1;
  // UV background on
  int UVbackground = 1;
  // cooling data for Haardt & Madau 2012 background
  char *grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5";
  // no dust
  int h2_on_dust = 0;
  // include CMB cooling floor
  int cmb_temperature_floor = 1;
  // set adiabatic index
  double Gamma = 5./3.;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  int comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  double density_units = 1.67e-24;
  double length_units = 1.0;
  double time_units = 1.0e12;
  double velocity_units = length_units / time_units;
  double a_units = 1.0; // units for the expansion factor

  // Set initial expansion factor (for internal units).
  // Set expansion factor to 1 for non-cosmological simulation.
  double initial_redshift = 0.;
  double a_value = 1. / (1. + initial_redshift);

  // Initialize everything.
  if (initialize_grackle(comoving_coordinates,
                         density_units, length_units,
                         time_units, velocity_units,
                         a_units, a_value,
                         use_grackle, with_radiative_cooling,
                         grackle_data_file,
                         primordial_chemistry, metal_cooling,
                         UVbackground, h2_on_dust,
                         cmb_temperature_floor, Gamma) == 0) {
    fprintf(stderr, "Error in initialize_grackle.\n");
    return 0;
  }

  // Allocate field arrays.
  gr_float *density, *energy, *x_velocity, *y_velocity, *z_velocity,
    *metal_density;
  gr_float tiny_number = 1.e-20;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  int grid_rank = 3;
  int grid_dimension[3], grid_start[3], grid_end[3];
  int i;
  for (i = 0;i < 3;i++) {
    grid_dimension[i] = 1; // the active dimension not including ghost zones.
    grid_start[i] = 0;
    grid_end[i] = 0;
  }
  grid_dimension[0] = field_size;
  grid_end[0] = field_size - 1;

  density       = malloc(field_size * sizeof(gr_float));
  energy        = malloc(field_size * sizeof(gr_float));
  x_velocity    = malloc(field_size * sizeof(gr_float));
  y_velocity    = malloc(field_size * sizeof(gr_float));
  z_velocity    = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  metal_density = malloc(field_size * sizeof(gr_float));

  // set temperature units
  double temperature_units = mh * pow(a_units * 
                                      length_units /
                                      time_units, 2) / kboltz;

  double HydrogenFractionByMass = 0.76;
  double SolarMetalFractionByMass = 0.02041;

  for (i = 0;i < field_size;i++) {
    density[i] = 1.0;
    // solar metallicity
    metal_density[i] = SolarMetalFractionByMass * density[i];

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
  double dt = 3.15e7 * 1e6 / time_units;

  if (solve_chemistry_table_(&comoving_coordinates,
                             &density_units, &length_units,
                             &time_units, &velocity_units,
                             &a_units, &a_value, &dt,
                             &grid_rank, grid_dimension,
                             grid_start, grid_end,
                             density, energy,
                             x_velocity, y_velocity, z_velocity,
                             metal_density) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return 0;
  }

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = malloc(field_size * sizeof(gr_float));
  if (calculate_cooling_time_table_(&comoving_coordinates,
                                    &density_units, &length_units,
                                    &time_units, &velocity_units,
                                    &a_units, &a_value,
                                    &grid_rank, grid_dimension,
                                    grid_start, grid_end,
                                    density, energy,
                                    x_velocity, y_velocity, z_velocity,
                                    metal_density,
                                    cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return 0;
  }

  fprintf(stderr, "Cooling time = %le s.\n", cooling_time[0] *
          time_units);

  // Calculate temperature.
  gr_float *temperature;
  temperature = malloc(field_size * sizeof(gr_float));
  if (calculate_temperature_table_(&comoving_coordinates,
                                   &density_units, &length_units,
                                   &time_units, &velocity_units,
                                   &a_units, &a_value,
                                   &grid_rank, grid_dimension,
                                   grid_start, grid_end,
                                   density, energy,
                                   metal_density, 
                                   temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return 0;
  }

  fprintf(stderr, "Temperature = %le K.\n", temperature[0]);

  // Calculate pressure.
  gr_float *pressure;
  pressure = malloc(field_size * sizeof(gr_float));
  if (calculate_pressure_table_(&comoving_coordinates,
                                &density_units, &length_units,
                                &time_units, &velocity_units,
                                &a_units, &a_value,
                                &grid_rank, grid_dimension,
                                grid_start, grid_end,
                                density, energy,
                                pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return 0;
  }

  fprintf(stderr, "Pressure = %le.\n", pressure[0]);

  return 1;
}
