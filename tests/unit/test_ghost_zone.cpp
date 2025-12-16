/***********************************************************************
/
/ This checks that ghost zones are not mutated
/ - this was written a few years before we adopted gtest and was adapted
/   retroactively to work as part of the test suite
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

#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "grtestutils/utils.hpp"
#include "grtestutils/googletest/fixtures.hpp"

#include <grackle.h>

#include <gtest/gtest.h>

#define mh     1.67262171e-24
#define kboltz 1.3806504e-16

typedef int (*property_func)(chemistry_data*, chemistry_data_storage*,
                             code_units*, grackle_field_data*, gr_float*);

typedef std::map<std::string, std::vector<gr_float>> val_vec_map_t;

struct grid_props{
  int dimensions[3];
  int ghost_depth[3];

  int field_size(){
    return dimensions[0] * dimensions[1] * dimensions[2];
  }
};

std::vector<gr_float> init_rand_vals(int field_size,
                                     std::minstd_rand& generator)
{
  std::vector<gr_float> out(field_size);
  for (int i = 0; i < field_size; i++){
    out[i] = (gr_float)grtest::random::uniform_dist_transform(generator);
  }
  return out;
}

class FieldInitHelper{
  // allocates new fields using std::vector and stores the vectors in val_map
public:
  FieldInitHelper(val_vec_map_t& val_map,
                  int field_size,
                  std::minstd_rand& generator)
    : val_map_(val_map), field_size_(field_size), generator_ref_(generator)
  { }

  gr_float* operator()(const std::string& name){
    val_map_[name] = init_rand_vals(field_size_, generator_ref_);
    return val_map_[name].data();
  }

private:
  val_vec_map_t& val_map_;
  int field_size_;
  std::minstd_rand& generator_ref_;
};

// allocates the grackle_field_data struct
void construct_field_data(grackle_field_data& my_fields,
                          grid_props& my_grid_props,
                          const chemistry_data* my_chemistry,
                          const code_units& my_units,
                          val_vec_map_t& val_map,
                          std::minstd_rand& generator){

  gr_float tiny_number = 1.e-20;

  // Create struct for storing grackle field data

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  my_fields.grid_rank = 3;

  for (int i = 0; i < 3; i++){
    int dim = my_grid_props.dimensions[i];
    int ghost_depth = my_grid_props.ghost_depth[i];
    if ( (dim <= 0) | (ghost_depth < 0) | ((2*ghost_depth) >= dim) ){
      abort();
    }
    my_fields.grid_dimension[i] = dim;
    my_fields.grid_start[i] = ghost_depth;
    my_fields.grid_end[i] = dim - ghost_depth - 1;
  }

  int field_size = my_grid_props.field_size();

  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  FieldInitHelper init_helper(val_map, field_size, generator);
  my_fields.density         = init_helper("density");
  my_fields.internal_energy = init_helper("internal_energy");
  my_fields.x_velocity      = init_helper("x_velocity");
  my_fields.y_velocity      = init_helper("y_velocity");
  my_fields.z_velocity      = init_helper("z_velocity");
  // for primordial_chemistry >= 1
  my_fields.HI_density      = init_helper("HI_density");
  my_fields.HII_density     = init_helper("HII_density");
  my_fields.HeI_density     = init_helper("HeI_density");
  my_fields.HeII_density    = init_helper("HeII_density");
  my_fields.HeIII_density   = init_helper("HeIII_density");
  my_fields.e_density       = init_helper("e_density");
  // for primordial_chemistry >= 2
  my_fields.HM_density      = init_helper("HM_density");
  my_fields.H2I_density     = init_helper("H2I_density");
  my_fields.H2II_density    = init_helper("H2II_density");
  // for primordial_chemistry >= 3
  my_fields.DI_density      = init_helper("DI_density");
  my_fields.DII_density     = init_helper("DII_density");
  my_fields.HDI_density     = init_helper("HDI_density");
  // for metal_cooling = 1
  my_fields.metal_density   = init_helper("metal_density");

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = init_helper("volumetric_heating_rate");
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = init_helper("specific_heating_rate");

  // radiative transfer ionization / dissociation rate fields (provided in units of [1/s])
  my_fields.RT_HI_ionization_rate = init_helper("RT_HI_ionization_rate");
  my_fields.RT_HeI_ionization_rate = init_helper("RT_HeI_ionization_rate");
  my_fields.RT_HeII_ionization_rate = init_helper("RT_HeII_ionization_rate");
  my_fields.RT_H2_dissociation_rate = init_helper("RT_H2_dissociation_rate");
  // radiative transfer heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.RT_heating_rate = init_helper("RT_heating_rate");

  // interstellar radiation field strength
  my_fields.isrf_habing = init_helper("isrf_habing");

  // set temperature units
  double temperature_units = get_temperature_units(&my_units);

  int gx = my_grid_props.ghost_depth[0];
  int gy = my_grid_props.ghost_depth[1];
  int gz = my_grid_props.ghost_depth[2];
  int mx = my_grid_props.dimensions[0];
  int my = my_grid_props.dimensions[1];
  int mz = my_grid_props.dimensions[2];

  // initialize values not in ghost zones
  for (int iz = gz; iz < (mz - gz); iz++){
    for (int iy = gy; iy < (my - gy); iy++){
      for (int ix = gx; ix < (mx - gx); ix++){
        int i = ix + mx * (iy + my*iz);
        my_fields.density[i] = 1.0;
        my_fields.HI_density[i] = my_chemistry->HydrogenFractionByMass *
          my_fields.density[i];
        my_fields.HII_density[i] = tiny_number * my_fields.density[i];
        my_fields.HM_density[i] = tiny_number * my_fields.density[i];
        my_fields.HeI_density[i] = (1.0 - my_chemistry->HydrogenFractionByMass)
          * my_fields.density[i];
        my_fields.HeII_density[i] = tiny_number * my_fields.density[i];
        my_fields.HeIII_density[i] = tiny_number * my_fields.density[i];
        my_fields.H2I_density[i] = tiny_number * my_fields.density[i];
        my_fields.H2II_density[i] = tiny_number * my_fields.density[i];
        my_fields.DI_density[i] = 2.0 * 3.4e-5 * my_fields.density[i];
        my_fields.DII_density[i] = tiny_number * my_fields.density[i];
        my_fields.HDI_density[i] = tiny_number * my_fields.density[i];
        my_fields.e_density[i] = tiny_number * my_fields.density[i];
        // solar metallicity
        my_fields.metal_density[i] = my_chemistry->SolarMetalFractionByMass *
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

        my_fields.isrf_habing[i] = my_chemistry->interstellar_radiation_field;
      }
    }
  }

}

// check that all elements in ref and actual that correspond to ghost cells
// are identical.
bool equal_ghost_values(const std::vector<gr_float>& ref,
                        const std::vector<gr_float>& actual,
                        grid_props& my_grid_props)
{
  int gx = my_grid_props.ghost_depth[0];
  int gy = my_grid_props.ghost_depth[1];
  int gz = my_grid_props.ghost_depth[2];
  int mx = my_grid_props.dimensions[0];
  int my = my_grid_props.dimensions[1];
  int mz = my_grid_props.dimensions[2];

  int num_unequal = 0;

  for (int iz = 0; iz < mz; iz++){
    for (int iy = 0; iy < my; iy++){
      for (int ix = 0; ix < mx; ix++){

        int i = ix + mx * (iy + my*iz);
        bool is_ghost = ( (iz < gz) | (iz >= (mz-gz)) |
                          (iy < gy) | (iy >= (my-gy)) |
                          (ix < gx) | (ix >= (mx-gx)) );

        num_unequal += (is_ghost & (ref[i] != actual[i]));
      }
    }
  }
  return num_unequal == 0;
}

bool equal_ghost_values(val_vec_map_t& ref, val_vec_map_t& actual,
                        grid_props& my_grid_props)
{
  for(val_vec_map_t::iterator iter = ref.begin(); iter != ref.end(); ++iter){
    const std::string& key = iter->first;
    if (!equal_ghost_values(ref[key], actual[key], my_grid_props)){
      return false;
    }
  }
  return true;
}

class APIGhostZoneTest: public grtest::ParametrizedConfigPresetFixture
{};

TEST_P(APIGhostZoneTest, GridZoneStartEnd) {

  grid_props my_grid_props = {{5,6,7}, {1,0,2}};

  // alias the pack_ attribute tracked by the fixture
  grtest::GrackleCtxPack& pack = pack_;
  code_units my_units = pack_.initial_units();

  // initialize pseudo random number generator
  std::uint32_t seed = 1379069008;
  std::minstd_rand generator(seed);

  // Create the struct for storing grackle field data.
  grackle_field_data my_fields;
  gr_initialize_field_data(&my_fields);
  int grid_prop_storage[9];
  my_fields.grid_dimension = grid_prop_storage + 0;
  my_fields.grid_start = grid_prop_storage + 3;
  my_fields.grid_end = grid_prop_storage + 6;

  // Initialize the fields in Data in the ghost zone
  // is initialized to random values.

  // Each field data pointer in my_fields actually aliases vector data
  // managed stored in my_field_map.
  // For example: my_fields.density = my_field_map["density"].data()
  val_vec_map_t my_field_map;
  construct_field_data(
    my_fields, my_grid_props, pack.my_chemistry(), my_units, my_field_map, generator
  );

  // orig_field_map_copy is a deepcopy of my_field_map. We will use this as a
  // reference to make sure that no values in the ghost cells are mutated.
  val_vec_map_t orig_field_map_copy = my_field_map;

  int field_size = my_grid_props.field_size();

  /*********************************************************************
  / Now check that ghost zones are not mutated by the chemistry solver or 
  / routines that compute various properties.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (local_solve_chemistry(pack.my_chemistry(), pack.my_rates(),
                            &my_units, &my_fields,
                            dt)!= GR_SUCCESS) {
    FAIL() << "Error running solve_chemistry";
  }

  // check if any ghost values were mutated
  if (!equal_ghost_values(orig_field_map_copy, my_field_map, my_grid_props)) {
    FAIL() << "Some ghost values were modified in solve_chemistry.";
  }

  // Now check what happens when computing various properties
  const char* func_names[5] = {"local_calculate_cooling_time",
                               "local_calculate_temperature",
                               "local_calculate_pressure",
                               "local_calculate_gamma",
                               "local_calculate_dust_temperature"};
  property_func func_ptrs[5] = {&local_calculate_cooling_time,
                                &local_calculate_temperature,
                                &local_calculate_pressure,
                                &local_calculate_gamma,
                                &local_calculate_dust_temperature};
  for (int i = 0; i < 5; i++){
    std::uint32_t seed2 = 1860889605;
    std::minstd_rand generator2(seed2);

    // allocate vector used to hold outputs (initialize with random vals)
    std::vector<gr_float> out_vals = init_rand_vals(field_size, generator2);
    // make a deep copy of out_vals before the calculation
    std::vector<gr_float> pre_calc_copy = out_vals;

    // perform the calculation
    property_func func_ptr = func_ptrs[i];
    if ( (*func_ptr)(pack.my_chemistry(), pack.my_rates(), &my_units, &my_fields, out_vals.data())
         != GR_SUCCESS ) {
      FAIL() << "Error reported by " << func_names;
    }

    // compare the ghost values from before and after the calculation
    if (!equal_ghost_values(out_vals, pre_calc_copy, my_grid_props)) {
      FAIL() << "Some ghost values were modified by " << func_names;
    }
  }

}

using grtest::FullConfPreset;
using grtest::ChemPreset;
using grtest::InitialUnitPreset;

INSTANTIATE_TEST_SUITE_P(
  VaryingPrimordialChem, APIGhostZoneTest,
  ::testing::Values(
    FullConfPreset{ChemPreset::primchem0, InitialUnitPreset::simple_z0},
    FullConfPreset{ChemPreset::primchem1, InitialUnitPreset::simple_z0},
    FullConfPreset{ChemPreset::primchem2, InitialUnitPreset::simple_z0},
    FullConfPreset{ChemPreset::primchem3, InitialUnitPreset::simple_z0})
);
