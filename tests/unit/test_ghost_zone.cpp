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

#include <cstdint>
#include <random>
#include <utility>  // std::move, std::pair
#include <vector>

#include "grtestutils/googletest/fixtures.hpp"
#include "grtestutils/harness/field_container.hpp"
#include "grtestutils/harness/fill_field_vals.hpp"
#include "grtestutils/harness/grackle_ctx_pack.hpp"
#include "grtestutils/harness/status.hpp"
#include "grtestutils/utils.hpp"

#include <grackle.h>

#include <gtest/gtest.h>

typedef int (*property_func)(chemistry_data*, chemistry_data_storage*,
                             code_units*, grackle_field_data*, gr_float*);

void write_random_vals(gr_float* ptr, int size, std::minstd_rand& generator) {
  for (int i = 0; i < size; i++) {
    ptr[i] = (gr_float)grtest::random::uniform_dist_transform(generator);
  }
}

// check that all elements in ref and actual that correspond to ghost cells
// are identical.
testing::AssertionResult equal_ghost_values(
    const gr_float* ref, const gr_float* actual,
    const grtest::GridLayout& my_grid_props) {
  int x_start = my_grid_props.start()[0];
  int y_start = my_grid_props.start()[1];
  int z_start = my_grid_props.start()[2];
  int x_stop = my_grid_props.stop()[0];
  int y_stop = my_grid_props.stop()[1];
  int z_stop = my_grid_props.stop()[2];
  int mx = my_grid_props.dim()[0];
  int my = my_grid_props.dim()[1];
  int mz = my_grid_props.dim()[2];

  auto is_unequal_ghost = [=](int ix, int iy, int iz) -> bool {
    int i = ix + mx * (iy + my * iz);
    bool is_ghost = ((iz < z_start) | (iz >= z_stop) | (iy < y_start) |
                     (iy >= y_stop) | (ix < x_start) | (ix >= x_stop));
    return (is_ghost & (ref[i] != actual[i]));
  };

  int num_unequal = 0;
  for (int iz = 0; iz < mz; iz++) {
    for (int iy = 0; iy < my; iy++) {
      for (int ix = 0; ix < mx; ix++) {
        num_unequal += is_unequal_ghost(ix, iy, iz);
      }
    }
  }

  if (num_unequal == 0) {
    return testing::AssertionSuccess();
  }

  for (int iz = 0; iz < mz; iz++) {
    for (int iy = 0; iy < my; iy++) {
      for (int ix = 0; ix < mx; ix++) {
        if (is_unequal_ghost(ix, iy, iz)) {
          int i = ix + mx * (iy + my * iz);
          return testing::AssertionFailure()
                 << "there are " << num_unequal << " unequal ghost zones.\n"
                 << "the first unequal value is at: " << '{' << ix << ", " << iy
                 << ", " << iz << "} (layout-left)\n"
                 << "  ref =    " << ref[i] << '\n'
                 << "  actual = " << actual[i] << '\n';
        };
      }
    }
  }
  return testing::AssertionFailure() << "something weird happened";
}

testing::AssertionResult equal_ghost_values(
    const grtest::FieldContainer& fc_ref,
    const grtest::FieldContainer& fc_actual,
    const grtest::GridLayout& my_grid_props) {
  for (const auto& name_ptr_pair : fc_ref) {
    const gr_float* ref_ptr = name_ptr_pair.second;
    std::string name = name_ptr_pair.first;
    std::optional<const gr_float*> tmp = fc_actual.find(name);
    if (!tmp.has_value()) {
      return testing::AssertionFailure()
             << "fc_actual has no field named \"" << name << "\"\n"
             << "fc_actual: " << testing::PrintToString(fc_actual);
    }
    const gr_float* actual_ptr = tmp.value();
    testing::AssertionResult rslt =
        equal_ghost_values(ref_ptr, actual_ptr, my_grid_props);
    if (!rslt) {
      return rslt << "\nthis occurred for the \"" << name << "\" field\n"
                  << "fc_ref: " << testing::PrintToString(fc_ref) << '\n'
                  << "fc_actual: " << testing::PrintToString(fc_actual);
    }
  }
  return testing::AssertionSuccess();
}

using APIGhostZoneTest = grtest::ParametrizedConfigPresetFixture;

TEST_P(APIGhostZoneTest, GridZoneStartEnd) {
  int dim[3] = {5, 6, 7};
  int ghostdepth[3] = {1, 0, 2};
  grtest::GridLayout my_grid_props =
      grtest::GridLayout::from_ghostdepth_and_dims(3, ghostdepth, dim);

  // the pack attribute holds grtest::GrackleCtxPack
  std::pair<grtest::FieldContainer, grtest::Status> tmp =
      grtest::FieldContainer::create(pack, my_grid_props);
  if (tmp.second.is_err()) {
    FAIL() << "something went wrong while creating grtest::FieldContainer: "
           << tmp.second.to_string();
  }
  grtest::FieldContainer fc_init = std::move(tmp.first);

  code_units my_units = pack.initial_units();

  // Initialize the fields fc_init (data in the ghost zone is initialized to
  // random values).
  int field_size = fc_init.grid_layout().n_elements();

  uint32_t seed = 1379069008;
  std::minstd_rand generator(seed);
  for (const auto& name_ptr_pair : fc_init) {
    write_random_vals(name_ptr_pair.second, field_size, generator);
  }

  double metallicity = 1.0;
  grtest::SimpleFieldTile field_tile =
      make_simple_tile(pack, my_units, metallicity);
  grtest::Status s = grtest::fill_field_vals(field_tile, fc_init);
  if (s.is_err()) {
    FAIL() << s.to_string();
  }

  // fc is a deepcopy of fc_init. We will only modify fc and use fc_init as the
  // reference to make sure that no values in the ghost cells are mutated.
  grtest::FieldContainer fc = fc_init.clone();

  /*********************************************************************
  / Now check that ghost zones are not mutated by the chemistry solver or
  / routines that compute various properties.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (local_solve_chemistry(pack.my_chemistry(), pack.my_rates(), &my_units,
                            fc.get_ptr(), dt) != GR_SUCCESS) {
    FAIL() << "Error running solve_chemistry";
  }

  // check if any ghost values were mutated
  ASSERT_TRUE(equal_ghost_values(fc, fc_init, my_grid_props))
      << "Some ghost values were modified in solve_chemistry.";

  // Now check what happens when computing various properties
  std::pair<const char*, property_func> name_fn_pairs[] = {
      {"local_calculate_cooling_time", &local_calculate_cooling_time},
      {"local_calculate_temperature", &local_calculate_temperature},
      {"local_calculate_pressure", &local_calculate_pressure},
      {"local_calculate_gamma", &local_calculate_gamma},
      {"local_calculate_dust_temperature", &local_calculate_dust_temperature}};
  for (const auto& [fn_name, fn_ptr] : name_fn_pairs) {
    uint32_t seed2 = 1860889605;
    std::minstd_rand generator2(seed2);

    // allocate vector used to hold outputs (initialize with random vals)
    std::vector<gr_float> out_vals(field_size, 0.0);
    write_random_vals(out_vals.data(), field_size, generator2);
    // make a deep copy of out_vals before the calculation
    std::vector<gr_float> pre_calc_copy = out_vals;

    // perform the calculation
    if ((*fn_ptr)(pack.my_chemistry(), pack.my_rates(), &my_units, fc.get_ptr(),
                  out_vals.data()) != GR_SUCCESS) {
      FAIL() << "Error reported by " << fn_name;
    }

    // compare the ghost values from before and after the calculation
    ASSERT_TRUE(equal_ghost_values(out_vals.data(), pre_calc_copy.data(),
                                   my_grid_props))
        << "Some ghost values were modified by " << fn_name;
  }
}

using grtest::ChemPreset;
using grtest::InitialUnitPreset;
using grtest::ParamConf;

static const ParamConf my_presets_[] = {
    ParamConf::SimplePreset(ChemPreset::primchem0,
                            InitialUnitPreset::simple_z0),
    ParamConf::SimplePreset(ChemPreset::primchem1,
                            InitialUnitPreset::simple_z0),
    ParamConf::SimplePreset(ChemPreset::primchem2,
                            InitialUnitPreset::simple_z0),
    ParamConf::SimplePreset(ChemPreset::primchem3,
                            InitialUnitPreset::simple_z0)};

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, APIGhostZoneTest,
    ::testing::ValuesIn(my_presets_));
