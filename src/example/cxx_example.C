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

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // Set initial redshift (for internal units).
  double initial_redshift = 0.;

  // Check the consistency
  if (gr_check_consistency() != GR_SUCCESS) {
    fprintf(stderr, "Error in gr_check_consistency.\n");
    return EXIT_FAILURE;
  }

  // Enable output
  grackle_verbose = 1;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;
  set_velocity_units(&my_units);

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  chemistry_data *my_grackle_data;
  my_grackle_data = new chemistry_data;
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }
  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  grackle_data->use_grackle = 1;            // chemistry on
  grackle_data->use_isrf_field = 1;
  grackle_data->with_radiative_cooling = 1; // cooling on
  grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
  grackle_data->dust_chemistry = 1;
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 1;           // UV background on
  grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  gr_float tiny_number = 1.e-20;

  // Create struct for storing grackle field data
  grackle_field_data my_fields;
  gr_initialize_field_data(&my_fields);

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start = new int[3];
  my_fields.grid_end = new int[3];
  for (int i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  my_fields.density         = new gr_float[field_size];
  my_fields.internal_energy = new gr_float[field_size];
  my_fields.x_velocity      = new gr_float[field_size];
  my_fields.y_velocity      = new gr_float[field_size];
  my_fields.z_velocity      = new gr_float[field_size];
  // for primordial_chemistry >= 1
  my_fields.HI_density      = new gr_float[field_size];
  my_fields.HII_density     = new gr_float[field_size];
  my_fields.HeI_density     = new gr_float[field_size];
  my_fields.HeII_density    = new gr_float[field_size];
  my_fields.HeIII_density   = new gr_float[field_size];
  my_fields.e_density       = new gr_float[field_size];
  // for primordial_chemistry >= 2
  my_fields.HM_density      = new gr_float[field_size];
  my_fields.H2I_density     = new gr_float[field_size];
  my_fields.H2II_density    = new gr_float[field_size];
  // for primordial_chemistry >= 3
  my_fields.DI_density      = new gr_float[field_size];
  my_fields.DII_density     = new gr_float[field_size];
  my_fields.HDI_density     = new gr_float[field_size];
  // for metal_cooling = 1
  my_fields.metal_density   = new gr_float[field_size];

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = new gr_float[field_size];
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = new gr_float[field_size];

  // radiative transfer ionization / dissociation rate fields (provided in units of [1/s])
  my_fields.RT_HI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeII_ionization_rate = new gr_float[field_size];
  my_fields.RT_H2_dissociation_rate = new gr_float[field_size];
  // radiative transfer heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.RT_heating_rate = new gr_float[field_size];

  // interstellar radiation field strength
  my_fields.isrf_habing = new gr_float[field_size];

  // set temperature units
  double temperature_units = get_temperature_units(&my_units);

  int i;
  for (i = 0;i < field_size;i++) {
    my_fields.density[i] = 1.0;
    my_fields.HI_density[i] = grackle_data->HydrogenFractionByMass *
      my_fields.density[i];
    my_fields.HII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HM_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) *
      my_fields.density[i];
    my_fields.HeII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeIII_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2I_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2II_density[i] = tiny_number * my_fields.density[i];
    my_fields.DI_density[i] = 2.0 * 3.4e-5 * my_fields.density[i];
    my_fields.DII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HDI_density[i] = tiny_number * my_fields.density[i];
    my_fields.e_density[i] = tiny_number * my_fields.density[i];
    // solar metallicity
    my_fields.metal_density[i] = grackle_data->SolarMetalFractionByMass *
      my_fields.density[i];

    my_fields.x_velocity[i] = 0.0;
    my_fields.y_velocity[i] = 0.0;
    my_fields.z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields.internal_energy[i] = 1000. / temperature_units;

    my_fields.volumetric_heating_rate[i] = 0.0;
    my_fields.specific_heating_rate[i] = 0.0;

    my_fields.RT_HI_ionization_rate[i] = 0.0;
    my_fields.RT_HeI_ionization_rate[i] = 0.0;
    my_fields.RT_HeII_ionization_rate[i] = 0.0;
    my_fields.RT_H2_dissociation_rate[i] = 0.0;
    my_fields.RT_heating_rate[i] = 0.0;

    my_fields.isrf_habing[i] = grackle_data->interstellar_radiation_field;
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(&my_units, &my_fields,
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "cooling_time = %24.16g s\n", cooling_time[0] *
          my_units.time_units);

  // Calculate temperature.
  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "temperature = %24.16g K\n", temperature[0]);

  // Calculate pressure.
  gr_float *pressure;
  double pressure_units = my_units.density_units *
    pow(my_units.velocity_units, 2);
  pressure = new gr_float[field_size];
  if (calculate_pressure(&my_units, &my_fields,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "pressure = %24.16g dyne/cm^2\n", pressure[0]*pressure_units);

  // Calculate gamma.
  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(&my_units, &my_fields,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "gamma = %24.16g\n", gamma[0]);

    // Calculate dust temperature.
  gr_float *dust_temperature;
  dust_temperature = new gr_float[field_size];
  if (calculate_dust_temperature(&my_units, &my_fields,
                                 dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "dust_temperature = %24.16g K\n", dust_temperature[0]);

  free_chemistry_data();

  return EXIT_SUCCESS;
}
