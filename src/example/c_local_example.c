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
#define sec_per_Myr 31.5576e12

int main(int argc, char *argv[])
{

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // Set initial redshift (for internal units).
  double initial_redshift = 0.;

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

  // Second, create a chemistry object for parameters.
  chemistry_data *my_grackle_data;
  my_grackle_data = malloc(sizeof(chemistry_data));
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }

  // Set parameter values for chemistry.
  my_grackle_data->use_grackle = 1;            // chemistry on
  my_grackle_data->with_radiative_cooling = 1; // cooling on
  my_grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
  my_grackle_data->metal_cooling = 1;          // metal cooling on
  my_grackle_data->UVbackground = 1;           // UV background on
  my_grackle_data->dust_chemistry = 1;         // dust processes
  my_grackle_data->use_dust_density_field = 1; // follow dust density field
  my_grackle_data->use_isrf_field = 1;         // follow interstellar radiation field
  my_grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012.h5"; // data file

  // Create chemistry data storage object to store rates.
  chemistry_data_storage my_grackle_rates;

  // Finally, initialize the chemistry object.
  if (local_initialize_chemistry_data(my_grackle_data, &my_grackle_rates, &my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  double tiny_number = 1.e-20;

  // Create struct for storing grackle field data
  grackle_field_data my_fields;
  gr_initialize_field_data(&my_fields);

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_start = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_end = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  int i;
  for (i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;

  my_fields.density         = malloc(field_size * sizeof(gr_float));
  my_fields.internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields.x_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields.y_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields.z_velocity      = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields.HI_density      = malloc(field_size * sizeof(gr_float));
  my_fields.HII_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HeI_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HeII_density    = malloc(field_size * sizeof(gr_float));
  my_fields.HeIII_density   = malloc(field_size * sizeof(gr_float));
  my_fields.e_density       = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields.HM_density      = malloc(field_size * sizeof(gr_float));
  my_fields.H2I_density     = malloc(field_size * sizeof(gr_float));
  my_fields.H2II_density    = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields.DI_density      = malloc(field_size * sizeof(gr_float));
  my_fields.DII_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HDI_density     = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields.metal_density   = malloc(field_size * sizeof(gr_float));
  // for use_dust_density_field = 1
  my_fields.dust_density    = malloc(field_size * sizeof(gr_float));

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = malloc(field_size * sizeof(gr_float));

  // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
  my_fields.RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  my_fields.RT_heating_rate = malloc(field_size * sizeof(gr_float));

  // interstellar radiation field strength
  my_fields.isrf_habing = malloc(field_size * sizeof(gr_float));

  // set temperature units
  double temperature_units = get_temperature_units(&my_units);

  for (i = 0;i < field_size;i++) {
    my_fields.density[i] = 1.0;
    my_fields.HI_density[i] = my_grackle_data->HydrogenFractionByMass *
      my_fields.density[i];
    my_fields.HII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HM_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeI_density[i] = (1.0 - my_grackle_data->HydrogenFractionByMass) *
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
    my_fields.metal_density[i] = my_grackle_data->SolarMetalFractionByMass *
      my_fields.density[i];
    // local dust to gas ratio
    my_fields.dust_density[i] = grackle_data->local_dust_to_gas_ratio *
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

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = malloc(field_size * sizeof(gr_float));
  if (local_calculate_cooling_time(my_grackle_data, &my_grackle_rates,
                                   &my_units, &my_fields,
                                   cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Before - cooling_time = %g s.\n", cooling_time[0] *
          my_units.time_units);

  // Calculate temperature.
  gr_float *temperature;
  temperature = malloc(field_size * sizeof(gr_float));
  if (local_calculate_temperature(my_grackle_data, &my_grackle_rates,
                                  &my_units, &my_fields,
                                  temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Before - temperature = %g K.\n", temperature[0]);

  // Calculate pressure.
  gr_float *pressure;
  double pressure_units = my_units.density_units *
    pow(my_units.velocity_units, 2);
  pressure = malloc(field_size * sizeof(gr_float));
  if (local_calculate_pressure(my_grackle_data, &my_grackle_rates,
                               &my_units, &my_fields,
                               pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Before - pressure = %g dyne/cm^2.\n",
          pressure[0]*pressure_units);

  // Calculate gamma.
  gr_float *gamma;
  gamma = malloc(field_size * sizeof(gr_float));
  if (local_calculate_gamma(my_grackle_data, &my_grackle_rates,
                            &my_units, &my_fields,
                            gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Before - gamma = %g.\n", gamma[0]);

  // Calculate dust temperature.
  gr_float *dust_temperature;
  dust_temperature = malloc(field_size * sizeof(gr_float));
  if (local_calculate_dust_temperature(my_grackle_data, &my_grackle_rates,
                                       &my_units, &my_fields,
                                       dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "Before - dust_temperature = %g K.\n", dust_temperature[0]);

  // Evolving the chemistry.
  // some timestep
  double dt = sec_per_Myr / my_units.time_units;
  fprintf(stderr, "Calling solve_chemistry with dt = %g Myr.\n",
          (dt * my_units.time_units / sec_per_Myr));
  if (local_solve_chemistry(my_grackle_data, &my_grackle_rates,
                            &my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }

  // Calculate cooling time.
  if (local_calculate_cooling_time(my_grackle_data, &my_grackle_rates,
                                   &my_units, &my_fields,
                                   cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "After - cooling_time = %g s.\n", cooling_time[0] *
          my_units.time_units);

  // Calculate temperature.
  if (local_calculate_temperature(my_grackle_data, &my_grackle_rates,
                                  &my_units, &my_fields,
                                  temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "After - temperature = %g K.\n", temperature[0]);

  // Calculate pressure.
  if (local_calculate_pressure(my_grackle_data, &my_grackle_rates,
                               &my_units, &my_fields,
                               pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "After - pressure = %g dyne/cm^2.\n",
          pressure[0]*pressure_units);

  // Calculate gamma.
  if (local_calculate_gamma(my_grackle_data, &my_grackle_rates,
                            &my_units, &my_fields,
                            gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "After - gamma = %g.\n", gamma[0]);

  // Calculate dust temperature.
  if (local_calculate_dust_temperature(my_grackle_data, &my_grackle_rates,
                                       &my_units, &my_fields,
                                       dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }
  fprintf(stdout, "After - dust_temperature = %g K.\n", dust_temperature[0]);

  local_free_chemistry_data(my_grackle_data, &my_grackle_rates);

  return EXIT_SUCCESS;
}
