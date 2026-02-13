//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Miscellaneous tests related to (partially) disabling Grackle
///
//===----------------------------------------------------------------------===//

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "grackle.h"
#include "grtestutils/googletest/fixtures.hpp"
#include "grtestutils/harness/grackle_ctx_pack.hpp"
#include "grtestutils/harness/field_container.hpp"
#include "grtestutils/harness/fill_field_vals.hpp"
#include "grtestutils/harness/preset.hpp"
#include "grtestutils/harness/status.hpp"

// this file includes a bunch of assorted API tests
//
// the logic all appears a little repetitive. Usually,
// - we construct a grackle configuration and then check:
//   - changes in the initial and final field (after local_solve_chemistry)
//   - computed derived properties
// - sometimes a single test does this multiple times for slightly different
//   parameter configurations (verify that a single parameter produces the
//   expected result)
//
// I already have some ideas for some slightly better abstractions to further
// simplify the size of our test-cases. But, we should probably wait until
// after we have implemented more test-cases.

// all these tests of the API are forms of property-testing
// - in other words, we are checking invariants that should remain true even
//   when input field values, units, or other configuration parameters are
//   varied (and use reasonable values)
// - we should keep this in mind if we ever adopt a fuzz-testing framework (to
//   hunt for undefined behavior), like https://llvm.org/docs/LibFuzzer.html
//   or https://github.com/google/fuzztest

// this test checks that use_grackle=0 actually disables grackle
// -> personally, I think we should get rid of this parameter, but its useful
//    to check this behavior
// -> todo: we should make separate tests to verify that the
//    `local_calculate_<quan>` functions exit early
TEST(APIMiscTest, UseGrackle0) {
  // setup the context for a grackle solver with metal_cooling
  grtest::ParamConf conf(grtest::ChemPreset::primchem0,
                         grtest::InitialUnitPreset::simple_z0,
                         grtest::make_ParamPair_vec({
                             {"use_grackle", 0},
                         }));
  GRTest_MAKE_CTX_PACK(grtest::GrackleCtxPack ctx_pack, conf);

  // we reuse the initial units as the current units
  code_units my_units = ctx_pack.initial_units();

  // make and fill the field container
  grtest::SimpleFieldTile field_tile =
      grtest::make_simple_tile(ctx_pack, my_units, 1.0);
  std::pair<grtest::FieldContainer, grtest::Status> tmp =
      grtest::create_and_fill_FieldContainer(
          field_tile, ctx_pack, grtest::GridLayout::create_1d<1>());
  if (tmp.second.is_err()) {
    FAIL() << "error while creating & filling grtest::FieldContainer: "
           << tmp.second.to_string();
  }
  grtest::FieldContainer fc = std::move(tmp.first);

  gr_float* eint = fc.find("internal_energy").value_or(nullptr);
  if (eint == nullptr) {
    FAIL() << "Problem accessing the \"internal_energy\" field.";
  }
  gr_float initial_eint = eint[0];

  double dt = 3.15e7 * 1e6 / my_units.time_units;
  if (local_solve_chemistry(ctx_pack.my_chemistry(), ctx_pack.my_rates(),
                            &my_units, fc.get_ptr(), dt) != GR_SUCCESS) {
    FAIL() << "local_solve_chemistry failed";
  }

  EXPECT_EQ(initial_eint, eint[0])
      << "local_solve_chemistry should not change the internal energy when "
      << "grackle is configure with use_grackle=0";
}

// The test suite using this fixture checks behavior when the
// "with_radiative_cooling" parameter is set to 0
class ParametrizedWithRadiativeCooling0
    : public testing::TestWithParam<grtest::ChemPreset> {
protected:
  void SetUp() override {
    // setup this->ctx_pack
    const grtest::ChemPreset chem_preset =
        ParametrizedWithRadiativeCooling0::GetParam();
    grtest::ParamConf conf(chem_preset, grtest::InitialUnitPreset::simple_z0,
                           grtest::make_ParamPair_vec({
                               {"with_radiative_cooling", 0},
                           }));
    GRTest_MAKE_CTX_PACK(ctx_pack, conf);

    // we reuse the initial units as the current units
    my_units = ctx_pack.initial_units();

    // make & fill the field container
    grtest::SimpleFieldTile field_tile =
        grtest::make_simple_tile(ctx_pack, my_units, 1.0);
    grtest::GridLayout layout = grtest::GridLayout::create_1d<1>();
    std::pair<grtest::FieldContainer, grtest::Status> tmp =
        grtest::create_and_fill_FieldContainer(field_tile, ctx_pack, layout);
    if (tmp.second.is_err()) {
      FAIL() << "error while creating & filling grtest::FieldContainer: "
             << tmp.second.to_string();
    }
    fc = std::unique_ptr<grtest::FieldContainer>(
        new grtest::FieldContainer(std::move(tmp.first)));
  }

  code_units my_units;
  grtest::GrackleCtxPack ctx_pack;
  std::unique_ptr<grtest::FieldContainer> fc;
};

TEST_P(ParametrizedWithRadiativeCooling0, SolveChemistry) {
  const grtest::FieldContainer original_fc = fc->clone();

  double dt = 3.15e7 * 1e6 / my_units.time_units;
  if (local_solve_chemistry(ctx_pack.my_chemistry(), ctx_pack.my_rates(),
                            &my_units, fc->get_ptr(), dt) != GR_SUCCESS) {
    FAIL() << "local_solve_chemistry failed";
  }

  // in a cosmological sim, eint may not be exactly the same
  // -> since we convert to and from proper units in-place
  EXPECT_EQ(original_fc.get_or_abort("internal_energy")[0],
            fc->get_or_abort("internal_energy")[0])
      << "local_solve_chemistry should not change the internal energy when "
      << "grackle is configure with with_radiative_cooling=0";

  if (ctx_pack.my_chemistry()->primordial_chemistry > 0) {
    const char* names[] = {"e_density", "HI_density", "HII_density"};
    for (const char* name : names) {
      gr_float before = original_fc.get_or_abort(name)[0];
      gr_float after = fc->get_or_abort(name)[0];
      // we place after in the denominator since before may be very close to 0
      EXPECT_GT(std::fabs((before - after) / after), 1.0e-12)
          << "local_solve_chemistry should have updated the \"" << name
          << "\" species field since the initial conditions were not in "
          << "chemical equilibrium";
    }
  }
}

TEST_P(ParametrizedWithRadiativeCooling0, CalculateCoolingTime) {
  gr_float val = NAN;
  if (local_calculate_cooling_time(ctx_pack.my_chemistry(), ctx_pack.my_rates(),
                                   &my_units, fc->get_ptr(),
                                   &val) != GR_SUCCESS) {
    FAIL() << "local_calculate_cooling_time";
  }
  EXPECT_TRUE(std::isfinite(val))
      << "local_calculate_cooling_time has historically compute values, even "
      << "when grackle is configure with with_radiative_cooling=0";
}

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, ParametrizedWithRadiativeCooling0,
    ::testing::Values(grtest::ChemPreset::primchem0,
                      grtest::ChemPreset::primchem1,
                      grtest::ChemPreset::primchem2,
                      grtest::ChemPreset::primchem3));

// The test suite using this fixture attempt to test that the 3 different
// approaches for disabling metal cooling all produce consistent results
class ParametrizedMetalFreeTest
    : public testing::TestWithParam<grtest::ChemPreset> {
protected:
  // This struct is used to describe a configuration of FieldContainer
  // - the fiducial case enables metal cooling (for the sake of comparison)
  // - the other cases all disable cooling in one way or another
  struct MetalCase {
    const char* descr;
    code_units my_units;
    grtest::GrackleCtxPack& ctx_pack;
    grtest::FieldContainer fc;
  };

  void SetUp() override {
    const grtest::ChemPreset chem_preset =
        ParametrizedMetalFreeTest::GetParam();

    // setup the context for a grackle solver with metal_cooling
    grtest::ParamConf main_conf(chem_preset,
                                grtest::InitialUnitPreset::simple_z0,
                                grtest::make_ParamPair_vec({
                                    {"dust_chemistry", 0},
                                }));
    GRTest_MAKE_CTX_PACK(this->metalcool_ctxpack, main_conf);

    // setup the context for a grackle solver without metal_cooling
    grtest::ParamConf alt_conf(chem_preset,
                               grtest::InitialUnitPreset::simple_z0,
                               grtest::make_ParamPair_vec({
                                   {"metal_cooling", 0},
                                   {"dust_chemistry", 0},
                               }));
    GRTest_MAKE_CTX_PACK(this->nometalcool_ctxpack, alt_conf);

    // all fluid containers have a single element
    const grtest::GridLayout common_layout = grtest::GridLayout::create_1d<1>();

    // initialize the FieldContainer for each MetalCase
    for (int i = 0; i < 4; i++) {
      // get the settings for the current case
      const char* descr;
      grtest::GrackleCtxPack* ctx_pack;
      double local_metallicity = 0.0;
      std::set<std::string> exclude_fields{};

      switch (i) {
        case 0:  // this is the fiducial case where metal cooling is enabled
          descr = "<fiducial case with cooling>";
          ctx_pack = &this->metalcool_ctxpack;
          local_metallicity = 1.0;
          break;
        case 1:
          // in this case, the solver explicitly models metal-cooling, but the
          // local fluid conditions just happen to not have any metals
          // -> this comes up in cosmological simulations at very early cosmic
          //    times (i.e. before there have been any SNe)
          descr = "<chem.metal_cooling=1,fields.metal_density[:]=0.0>";
          ctx_pack = &this->metalcool_ctxpack;
          break;
        case 2:
          // in this case, we explicitly disable metal_cooling
          descr = "<chem.metal_cooling=0>";
          ctx_pack = &this->nometalcool_ctxpack;
          // when chemistry_data::metal_cooling is 0, Grackle expects
          // metal_density to be a nullptr (so, no need to exclude it)
          break;
        default:
          // in this case, metal_cooling is technically enabled, but a nullptr
          // is registered for the metal_density field
          // -> Grackle was explicitly designed to disable metal-cooling (on
          //    the fly) when run in this manner. I suspect that this mode of
          //    operation existed before the chemistry_data::metal_cooling flag
          // -> TODO: (maybe in 4.0), I think we should stop supporting this
          //    mode! If people don't want metal-cooling we should require them
          //    to set the metal_cooling flag to 0. I think that this mode of
          //    operation is somewhat bug-prone
          descr = "<chem.metal_cooling=1,fields.metal_density=nullptr>";
          ctx_pack = &this->metalcool_ctxpack;
          exclude_fields = {"metal_density"};
      }

      // we reuse the initial units as the current units
      code_units my_units = ctx_pack->initial_units();

      // create and fill the field container
      grtest::SimpleFieldTile field_tile =
          make_simple_tile(*ctx_pack, my_units, local_metallicity);
      std::pair<grtest::FieldContainer, grtest::Status> tmp =
          grtest::create_and_fill_FieldContainer(field_tile, *ctx_pack,
                                                 common_layout, exclude_fields);
      if (tmp.second.is_err()) {
        FAIL() << "error while creating & filling grtest::FieldContainer: "
               << tmp.second.to_string();
      }

      // actually create MetalCase instance
      all_metal_cases.emplace_back(
          MetalCase{descr, my_units, *ctx_pack, std::move(tmp.first)});
    }

    fid_metal_cool_case_idx = 0;
    first_nometal_cool = 1;
  }

  // this performs some nonfatal assertions (Ugh... this is clunky)
  void common_checks(const std::vector<double>& vals,
                     const std::string& val_descr) const {
    const char* metal_cool_descr =
        all_metal_cases[fid_metal_cool_case_idx].descr;
    double metal_cool_val = vals[fid_metal_cool_case_idx];
    for (std::size_t i = 0; i < vals.size(); i++) {
      if (i == fid_metal_cool_case_idx) {
        continue;
      }

      EXPECT_GT(std::fabs((vals[i] - metal_cool_val) / metal_cool_val), 1.0e-7)
          << "Comparing " << val_descr << '\n'
          << "    -> " << metal_cool_descr << ": " << metal_cool_val << '\n'
          << "    -> " << all_metal_cases[i].descr << ": " << vals[i] << '\n'
          << "These quantities should be different. Previous Grackle versions "
          << "produced larger relative differences";

      if (i == first_nometal_cool) {
        continue;
      }

      double rel_err = 1.0e-15;
      const char* ref_descr = all_metal_cases[first_nometal_cool].descr;
      double ref_val = vals[first_nometal_cool];

      EXPECT_NEAR(ref_val, vals[i], std::fabs(rel_err * ref_val))
          << "Comparing " << val_descr << '\n'
          << "    -> " << ref_descr << ": " << ref_val << '\n'
          << "    -> " << all_metal_cases[i].descr << ": " << vals[i] << '\n'
          << "These quantities should be the same (within a relative error of "
          << rel_err << ')';
    }
  }

  // holds a fiducial context-pack (where my_chem.metal_cooling == 1)
  grtest::GrackleCtxPack metalcool_ctxpack;
  // holds an alternative context-pack (where my_chem.metal_cooling == 0)
  grtest::GrackleCtxPack nometalcool_ctxpack;

  // holds all MetalCase instances
  std::vector<MetalCase> all_metal_cases;

  // holds the index of the case with metal cooling
  std::size_t fid_metal_cool_case_idx;

  // holds the index of the first case without metal cooling
  std::size_t first_nometal_cool;
};

TEST_P(ParametrizedMetalFreeTest, SolveChemistry) {
  // perform local_solve_chemistry and record the change in internal_energy
  std::vector<double> delta_e(all_metal_cases.size());
  double dt_cgs = 3.15e7 * 1e6;
  for (std::size_t i = 0; i < delta_e.size(); i++) {
    MetalCase& cur_case = all_metal_cases[i];
    gr_float* eint = cur_case.fc.find("internal_energy").value_or(nullptr);
    if (eint == nullptr) {
      FAIL() << "Problem accessing the \"internal_energy\" field for the "
             << cur_case.descr << " case (THIS SHOULD NEVER HAPPEN)";
    }
    gr_float eint0 = eint[0];

    double dt = dt_cgs / cur_case.my_units.time_units;
    if (local_solve_chemistry(cur_case.ctx_pack.my_chemistry(),
                              cur_case.ctx_pack.my_rates(), &cur_case.my_units,
                              cur_case.fc.get_ptr(), dt) != GR_SUCCESS) {
      FAIL() << "local_solve_chemistry failed for the " << cur_case.descr
             << " case";
    }

    delta_e[i] = eint[0] - eint0;
  };

  // now let's compare the change in internal_energy values between the cases
  std::string quan_descr = "change in internal energy from evolution over " +
                           std::to_string(dt_cgs) + " seconds";
  common_checks(delta_e, quan_descr);
}

TEST_P(ParametrizedMetalFreeTest, CalculateTemperature) {
  // perform local_solve_chemistry and record the change in internal_energy
  std::vector<double> temperature(all_metal_cases.size());
  for (std::size_t i = 0; i < temperature.size(); i++) {
    MetalCase& cur_case = all_metal_cases[i];

    gr_float val;
    if (local_calculate_temperature(
            cur_case.ctx_pack.my_chemistry(), cur_case.ctx_pack.my_rates(),
            &cur_case.my_units, cur_case.fc.get_ptr(), &val) != GR_SUCCESS) {
      FAIL() << "local_calculate_temperature failed for the " << cur_case.descr
             << " case";
    }

    temperature[i] = val;
  }

  // now let's compare the computed temperature values between the cases
  common_checks(temperature, "temperature");
}

TEST_P(ParametrizedMetalFreeTest, CalculateCoolingTime) {
  // perform local_solve_chemistry and record the change in internal_energy
  std::vector<double> tcool(all_metal_cases.size());
  for (std::size_t i = 0; i < tcool.size(); i++) {
    MetalCase& cur_case = all_metal_cases[i];

    gr_float val;
    if (local_calculate_cooling_time(
            cur_case.ctx_pack.my_chemistry(), cur_case.ctx_pack.my_rates(),
            &cur_case.my_units, cur_case.fc.get_ptr(), &val) != GR_SUCCESS) {
      FAIL() << "local_calculate_cooling_time failed for the " << cur_case.descr
             << " case";
    }

    tcool[i] = val;
  }

  // now let's compare the computed temperature values between the cases
  common_checks(tcool, "cooling time");
}

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, ParametrizedMetalFreeTest,
    ::testing::Values(grtest::ChemPreset::primchem0,
                      grtest::ChemPreset::primchem1,
                      grtest::ChemPreset::primchem2,
                      grtest::ChemPreset::primchem3));