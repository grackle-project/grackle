//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Miscellaneous tests pertaining to the use of injection pathways
///
//===----------------------------------------------------------------------===//

#include <array>
#include <cstdint>
#include <string>
#include <optional>
#include <utility>

#include "grackle.h"
#include <gtest/gtest.h>
#include "grtestutils/iterator_adaptor.hpp"
#include "grtestutils/preset.hpp"
#include "grtestutils/googletest/assertions.hpp"
#include "grtestutils/googletest/fixtures.hpp"
#include "grtestutils/googletest/matchers.hpp"
#include "grtestutils/ratequery.hpp"
#include "grtestutils/utils.hpp"  // starts_with

using namespace std::string_literals;  // enables s-suffix for std::string

using grtest::ChemPreset;
using grtest::FullConfPreset;
using grtest::InitialUnitPreset;

const std::array<const char*, 12> CANONICAL_INJ_PATH_NAMES = {
    "local_ISM", "ccsn13", "ccsn20", "ccsn25",  "ccsn30",  "fsn13",
    "fsn15",     "fsn50",  "fsn80",  "pisn170", "pisn200", "y19",
};

const std::string INJ_PATH_NAMES_KEY = "inject_model_names"s;
const std::string INJ_PATH_NUCLIDE_GAS_YIELD_PREFIX =
    "inject_path_gas_yield_frac."s;
const std::string INJ_PATH_GRAIN_YIELD_PREFIX =
    "inject_path_grain_yield_frac."s;

using InjPathSuite = grtest::ParametrizedConfigPresetFixture;

static std::optional<std::vector<std::string>> get_field_names(
    const chemistry_data_storage& my_rates, int64_t category) {
  int64_t n = INT64_MAX;
  if (gr_get_field_id(&my_rates, category, nullptr, &n) != GR_SUCCESS) {
    return std::nullopt;
  }
  std::vector<std::string> out;
  out.reserve(n);
  for (int64_t i = 0; i < n; i++) {
    const char* name = nullptr;
    if (gr_get_field_name(&my_rates, category, i, &name) != GR_SUCCESS) {
      return std::nullopt;
    }
    out.emplace_back(name);
  }
  return out;
}

TEST_P(InjPathSuite, FieldNames) {
  std::optional<std::vector<std::string>> str_list_maybe =
      get_field_names(*pack.my_rates(), GR_FIELD_CAT_INJ_PATHWAY);
  ASSERT_NOT_EMPTY(str_list_maybe)
      << "error querying \"" << GR_FIELD_CAT_INJ_PATHWAY << '"';
  std::vector<std::string> str_list =
      std::exchange(str_list_maybe, std::nullopt).value();

  // construct the reference list of field names:
  constexpr std::size_t N_KNOWN_NAMES = CANONICAL_INJ_PATH_NAMES.size();
  std::array<std::string, N_KNOWN_NAMES> known_field_names;
  for (std::size_t i = 0; i < N_KNOWN_NAMES; i++) {
    known_field_names[i] =
        std::string(CANONICAL_INJ_PATH_NAMES[i]) + "_metal_density";
  }

  EXPECT_THAT(str_list, grtest::HoldsUnique());
  EXPECT_THAT(str_list, ::testing::Not(::testing::IsEmpty()));
  EXPECT_THAT(str_list, ::testing::IsSubsetOf(known_field_names));

  // now let's check that `gr_get_field_id` works properly
  // -> in the future, it probably makes some sense to dedicate some tests to
  //    explicitly checking this logic (and decoupling the test from the
  //    inject-pathway logic)
  // -> but, that's not very easy to do right now (i.e. currently, the ONLY
  //    fields accessible through the function are the inject_path fields)
  int64_t list_len = str_list.size();
  for (int64_t i = 0; i < list_len; i++) {
    const char* name = str_list[i].c_str();
    int64_t queried_id = INT64_MAX;
    ASSERT_GR_SUCCESS(gr_get_field_id(pack.my_rates(), GR_FIELD_CAT_INJ_PATHWAY,
                                      name, &queried_id))
        << "problem during call for \"" << name << '"';

    EXPECT_EQ(i, queried_id) << "incorrect id retrieved for\"" << name << '"';
  }
}

TEST_P(InjPathSuite, PathName) {
  // this is a little bit like the previous test, but we use the ratequery
  // machinery rather than the field_id machinery
  // -> the plan is to transition away from this approach!
  std::optional<std::vector<std::string>> str_list_maybe =
      grtest::try_query_str_arr1d(pack.my_rates(), INJ_PATH_NAMES_KEY);
  ASSERT_NOT_EMPTY(str_list_maybe)
      << "error querying \"" << INJ_PATH_NAMES_KEY << '"';
  std::vector<std::string> str_list =
      std::exchange(str_list_maybe, std::nullopt).value();

  EXPECT_THAT(str_list, grtest::HoldsUnique());
  EXPECT_THAT(str_list, ::testing::Not(::testing::IsEmpty()));
  EXPECT_THAT(str_list, ::testing::IsSubsetOf(CANONICAL_INJ_PATH_NAMES));
}

TEST_P(InjPathSuite, Yields) {
  std::optional<std::vector<std::string>> str_list_maybe =
      grtest::try_query_str_arr1d(pack.my_rates(), INJ_PATH_NAMES_KEY);
  ASSERT_NOT_EMPTY(str_list_maybe)
      << "error querying \"" << INJ_PATH_NAMES_KEY << '"';
  std::vector<std::string> str_list =
      std::exchange(str_list_maybe, std::nullopt).value();

  // ideally, we would verify that there is a yield for each kind of relevant
  // nuclide or grain species, but there's no convenient way to do that (yet!)
  // -> for now, lets just settle for checking the most common ones
  std::array<std::string, 2> keys = {INJ_PATH_NUCLIDE_GAS_YIELD_PREFIX + "C",
                                     INJ_PATH_GRAIN_YIELD_PREFIX + "AC_dust"};
  for (const std::string& key : keys) {
    std::optional<std::vector<double>> dbl_arr_maybe =
        grtest::try_query_f64_arr1d(pack.my_rates(), key);
    ASSERT_NOT_EMPTY(dbl_arr_maybe) << "error querying \"" << key << '"';
    std::vector<double> dbl_arr =
        std::exchange(dbl_arr_maybe, std::nullopt).value();
    EXPECT_THAT(dbl_arr, ::testing::SizeIs(str_list.size()))
        << "array associated with \"" << key << "\" should have the same "
        << "length as the array associated with \"" << INJ_PATH_NAMES_KEY
        << '"';
  }
}

// Let's check the sum of all the yields
// - the total yield at an index i should satidfy 0<= total[i] <= 1
// - this a very robust invariant that should be true for all injection pathway
//   information shipped with Grackle
TEST_P(InjPathSuite, YieldSum) {
  std::optional<std::vector<std::string>> str_list_maybe =
      grtest::try_query_str_arr1d(pack.my_rates(), INJ_PATH_NAMES_KEY);
  if (!str_list_maybe.has_value()) {
    GTEST_FAIL() << "error querying \"" << INJ_PATH_NAMES_KEY << '"';
  }
  std::vector<std::string> str_list =
      std::exchange(str_list_maybe, std::nullopt).value();

  // define an inline function that returns a string describing the string
  // associated with index i
  auto describe_index = [&](std::size_t i) -> std::string {
    return "in the current configuration, the index is associated with the \""s +
           str_list[i] + "\" injection pathway"s;
  };
  std::size_t n_pathways = str_list.size();

  // sum up the total yields (for grain species & metal nuclides in gas phase)
  // Reminders: we track yields for each injection pathway
  std::vector<double> totals(str_list.size(), 0.0);

  int n_gas_nuclide_yields = 0;
  int n_grain_yields = 0;
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    const char* descr = nullptr;
    if (grtest::starts_with(pair.name, INJ_PATH_NUCLIDE_GAS_YIELD_PREFIX)) {
      descr = "gas-phase metal nuclide yield fractions";
      n_gas_nuclide_yields++;
    } else if (grtest::starts_with(pair.name, INJ_PATH_GRAIN_YIELD_PREFIX)) {
      descr = "grain-species yield fractions";
      n_grain_yields++;
    } else {
      continue;
    }
    std::optional<std::vector<double>> dbl_arr_maybe =
        grtest::try_query_f64_arr1d(pack.my_rates(), pair.name);
    ASSERT_NOT_EMPTY(dbl_arr_maybe) << "error querying \"" << pair.name << '"';
    const std::vector<double> dbl_arr =
        std::exchange(dbl_arr_maybe, std::nullopt).value();
    for (std::size_t i = 0; i < n_pathways; i++) {
      EXPECT_THAT(dbl_arr[i],
                  ::testing::AllOf(::testing::Ge(0.0), ::testing::Le(1.0)))
          << "the array of " << descr << ", corresponding to the \""
          << pair.name << "\" key, has an invalid value at index " << i << " ("
          << describe_index(i) << ')';
      totals[i] += dbl_arr[i];
    }
  }

  for (std::size_t i = 0; i < n_pathways; i++) {
    EXPECT_THAT(totals[i], ::testing::AnyOf(::testing::Lt(1.0),
                                            ::testing::DoubleNear(1.0, 4e-7)))
        << "the sum of each yield fractions (i.e. fraction of all injected "
        << "non-primordial material) for all relevant grain-species and all "
        << "the gas-phase components of all relevant metal nuclides has an "
        << "invalid value at index " << i << " (" << describe_index(i) << ')';
  }

  EXPECT_GT(n_gas_nuclide_yields, 0)
      << "no keys were encountered for gas-phase metal nuclide yield fractions";
  EXPECT_GT(n_grain_yields, 0)
      << "no keys were encountered for n_gas_nuclide_yields";
}

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, InjPathSuite,
    ::testing::Values(FullConfPreset{ChemPreset::primchem4_dustspecies3,
                                     InitialUnitPreset::simple_z0},
                      FullConfPreset{
                          ChemPreset::primchem4_dustspecies3_allinjectpaths,
                          InitialUnitPreset::simple_z0}));
