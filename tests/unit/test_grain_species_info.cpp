//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Tests grackle::impl::GrainSpeciesInfo
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "LUT.hpp"
#include "dust/grain_species_info.hpp"

namespace {  // stuff in an anonymous namespace is local to this file

/// checks that the string, @p s, ends with the suffix, @p suffix
///
/// @param[in] s the string to check
/// @param[in] suffix the desired behavior
bool str_ends_with(std::string_view s, std::string_view suffix) {
#if __cplusplus >= 202002L
  return s.ends_with(suffix);
#else
  // according to cppreference, this is the behavior of the ends_with method
  // introduced in C++20
  std::size_t s_len = s.size();
  std::size_t suf_len = suffix.size();
  return ((s_len >= suf_len) &&
          (s.compare(s_len - suf_len, std::string_view::npos, suffix) == 0));
#endif
}

/// Represents a species name and its location within the global species list
struct SpNameIndexPair {
  std::string name;
  int index;
};

/// Returns a SpNameIndexPair for every species for which p returns true
///
/// p is a unary function that accepts a string_view and returns a boolean
///
/// @note
/// I think we'll be able to dispose of this in the near future
template <typename UnaryPred>
std::vector<SpNameIndexPair> selected_sp_info_(UnaryPred p) {
  std::vector<SpNameIndexPair> out;
  int count = 0;

#define STRINGIFY_(NAME) #NAME
#define ENTRY(NAME)                                                            \
  if (p(STRINGIFY_(NAME))) {                                                   \
    out.push_back(SpNameIndexPair{STRINGIFY_(NAME), count});                   \
  }                                                                            \
  count++;
#include "field_data_evolved_species.def"
#undef ENTRY
#undef STRINGIFY_

  return out;
}

/// Returns name-index pairs for all grains from the list of Grackle Species
///
/// The index corresponds to the index in global species list
std::vector<SpNameIndexPair> grain_sp_info_from_species_list() {
  auto fn = [](std::string_view s) { return str_ends_with(s, "_dust"); };
  return selected_sp_info_(fn);
}

/// Returns name for every chemical species (in proper order)
std::vector<std::string> chemical_species_list() {
  auto fn = [](std::string_view s) { return !str_ends_with(s, "_dust"); };
  const std::vector<SpNameIndexPair> tmp = selected_sp_info_(fn);
  std::vector<std::string> out;
  for (std::size_t i = 0; i < tmp.size(); i++) {
    out.push_back(tmp[i].name);
  }
  return out;
}

/// returns the number of grain species known to grackle
int number_known_grain_species() {
  // this is a very silly implementation
  return (int)(grain_sp_info_from_species_list().size());
}

/// tracks the max value allowed by Grackle's dust_species runtime parameter
constexpr int MAX_dust_species_VAL = 3;

// down below, we define some machinery to make sure we properly deallocate
// memory when a test fails. This would not be an issue if we were willing to
// make GrainSpeciesInfoDeleter a C++ class with a proper constructor &
// destructor
struct GrainSpeciesInfoDeleter {
  void operator()(grackle::impl::GrainSpeciesInfo* ptr) const {
    grackle::impl::drop_GrainSpeciesInfo(ptr);
    delete ptr;
  }
};

using unique_GrainSpeciesInfo_ptr =
    std::unique_ptr<grackle::impl::GrainSpeciesInfo, GrainSpeciesInfoDeleter>;

/// creates a unique_ptr holding a grackle::impl::GrainSpeciesInfo
///
/// This is useful for preventing memory leaks when tests fail
unique_GrainSpeciesInfo_ptr make_unique_GrainSpeciesInfo(
    int dust_species_param) {
  grackle::impl::GrainSpeciesInfo* ptr = new grackle::impl::GrainSpeciesInfo;
  (*ptr) = grackle::impl::new_GrainSpeciesInfo(dust_species_param);
  return unique_GrainSpeciesInfo_ptr(ptr, GrainSpeciesInfoDeleter());
}

}  // anonymous namespace

// check that OnlyGrainSpLUT contains the entries for every known grain species
TEST(OnlyGrainSpLUTTest, CheckNumEntries) {
  ASSERT_EQ(OnlyGrainSpLUT::NUM_ENTRIES, number_known_grain_species());
}

/// Define a fixture for running parameterized tests of the GrainSpeciesInfo
/// machinery. The tests are parameterized by the dust_species parameter.
class GrainSpeciesInfoTest : public testing::TestWithParam<int> {
protected:
  /// set up the grain_species_info pointer
  ///
  /// @note
  /// We **ONLY** perform setup in this method (rather than in a default
  /// constructor) because we want to perform some basic sanity checks
  void SetUp() override {
    int dust_species = GetParam();
    grain_species_info_ = make_unique_GrainSpeciesInfo(dust_species);

    // perform 3 sanity checks!
    ASSERT_GE(grain_species_info_->n_species, 1)
        << "GrainSpeciesInfo::n_species should be positive";
    ASSERT_LE(grain_species_info_->n_species, number_known_grain_species())
        << "GrainSpeciesInfo::n_species can't exceed the number of known grain "
        << "species";
    ASSERT_NE(grain_species_info_->species_info, nullptr)
        << "GrainSpeciesInfo::species_info can't be a nullptr when "
        << "GrainSpeciesInfo::n_species is positive";
  }

  unique_GrainSpeciesInfo_ptr grain_species_info_;
};

TEST_P(GrainSpeciesInfoTest, CheckOnlyGrainSpeciesLUTConsistency) {
  const int n_species = grain_species_info_->n_species;
  for (int i = 0; i < n_species; i++) {
    // sanity check!
    ASSERT_NE(grain_species_info_->species_info[i].name, nullptr);
    // actual check!
    EXPECT_EQ(i, grain_species_info_->species_info[i].onlygrainsp_idx)
        << "element " << i << " of the GrainSpeciesInfo::species_info array "
        << "doesn't seem to be synchronized with the OnlyGrainSpeciesLUT "
        << "enumeration. At face value (there aren't other related bugs), "
        << "OnlyGrainSpeciesLUT::" << grain_species_info_->species_info[i].name
        << " seems to have a value of "
        << grain_species_info_->species_info[i].onlygrainsp_idx;
  }
}

TEST_P(GrainSpeciesInfoTest, SublimationTemperature) {
  const int n_species = grain_species_info_->n_species;
  for (int i = 0; i < n_species; i++) {
    // sanity check!
    ASSERT_NE(grain_species_info_->species_info[i].name, nullptr);
    // actual check!
    EXPECT_GT(grain_species_info_->species_info[i].sublimation_temperature, 0)
        << "element " << i << " of the GrainSpeciesInfo::species_info array "
        << "holds an invalid sublimation temperature";
  }
}

TEST_P(GrainSpeciesInfoTest, SpeciesLUTCompare) {
  // this is the reference list of all grain species constructed from the same
  // XMacro as that the SpeciesLUT enumeration is built from
  std::vector<SpNameIndexPair> ref_list = grain_sp_info_from_species_list();

  const int n_species = grain_species_info_->n_species;
  for (int i = 0; i < n_species; i++) {
    // sanity check!
    ASSERT_NE(grain_species_info_->species_info[i].name, nullptr);
    // actual check!
    std::string actual_name(grain_species_info_->species_info[i].name);
    EXPECT_EQ(actual_name, ref_list[i].name)
        << "according to the reference list, the grain species at index " << i
        << " should be `" << ref_list[i].name << "` (not `" << actual_name
        << "`)";
    EXPECT_EQ(grain_species_info_->species_info[i].species_idx,
              ref_list[i].index)
        << "the value @ index " << i << " of GrainSpeciesInfo::species_info "
        << "seems to indicate that the SpLUT::" << actual_name << " enumerator "
        << "has a different value compared to the reference list";
  }
}

// check that the grackle::impl::max_ingredients_per_grain_species constant
// indeed holds the maximum number of ingredients that a grain species can
// have.
TEST_P(GrainSpeciesInfoTest, MaxIngredientsPerGrainSpecies) {
  // go through and compute the number of ingredients for each species
  int n_grain_species = grain_species_info_->n_species;
  int max_ingredient_count = 0;
  for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
    max_ingredient_count = std::max(
        max_ingredient_count,
        grain_species_info_->species_info[gsp_idx].n_growth_ingredients);
  }

  if (MAX_dust_species_VAL == GetParam()) {
    EXPECT_EQ(max_ingredient_count,
              grackle::impl::max_ingredients_per_grain_species)
        << "it appears that grackle::impl::max_ingredients_per_grain_species "
        << "should have a value of " << max_ingredient_count;
  } else {
    EXPECT_LE(max_ingredient_count,
              grackle::impl::max_ingredients_per_grain_species)
        << "it appears that grackle::impl::max_ingredients_per_grain_species "
        << "should be at least " << max_ingredient_count;
  }
}

// we choose to use std::pair for this purpose since equality is already
// defined and googletest knows how to compare them
using CoefNamePair = std::pair<int, std::string>;

// we are intentionally not very exhaustive here
std::map<std::string, std::vector<CoefNamePair>> get_ingredients(
    int dust_chemistry_parameter) {
  std::map<std::string, std::vector<CoefNamePair>> out;

  if (dust_chemistry_parameter > 0) {
    out["MgSiO3_dust"] =
        std::vector<CoefNamePair>{{1, "Mg"}, {1, "SiOI"}, {2, "H2O"}};
    out["AC_dust"] = std::vector<CoefNamePair>{{1, "CI"}};
  }
  if (dust_chemistry_parameter > 1) {
    out["SiM_dust"] = std::vector<CoefNamePair>{{1, "SiI"}};
    // skip a few!
    out["Fe3O4_dust"] = std::vector<CoefNamePair>{{3, "Fe"}, {4, "H2O"}};
    // skip a bunch more
  }
  if (dust_chemistry_parameter > 2) {
    std::string names[3] = {"ref_org_dust", "vol_org_dust", "H2O_ice_dust"};
    for (int i = 0; i < 3; i++) {
      out[names[i]] = std::vector<CoefNamePair>();
    }
  };
  return out;
}

// here we are going to some very simple tests that we properly recorded growth
// ingredients.
//
// NOTE: currently, this is extremely crude! We've almost finished machinery
// that will make this a lot easier (and less verbose!)
TEST_P(GrainSpeciesInfoTest, SampledGrainIngredients) {
  int dust_chemistry_parameter = GetParam();

  // get the list of ALL know chemical species names, in the cannonical order
  const std::vector<std::string> chem_species_names = chemical_species_list();
  const std::size_t known_chem_species_count = chem_species_names.size();

  std::map<std::string, std::vector<CoefNamePair>> ref_ingred_map =
      get_ingredients(dust_chemistry_parameter);

  // iterate through each grain species from grain_species_info_
  // - validate that each index of the ingredients is meaningful
  // - if the dust grain is listed in ingred_map, then we will compare the
  //   ingredient list
  int n_grain_species = grain_species_info_->n_species;
  std::size_t n_comparisons = 0;
  for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
    const grackle::impl::GrainSpeciesInfoEntry& species_info =
        grain_species_info_->species_info[gsp_idx];
    std::string name = species_info.name;

    // try to coerce the ingredients to a vector
    std::vector<CoefNamePair> ingred_vec;
    for (int i = 0; i < species_info.n_growth_ingredients; i++) {
      const grackle::impl::GrainGrowthIngredient& cur =
          species_info.growth_ingredients[i];
      ASSERT_TRUE(0 <= cur.species_idx &&
                  cur.species_idx < known_chem_species_count)
          << "for ingredient " << i << " of the \"" << name
          << "\" grain species, "
          << "the species index, " << cur.species_idx
          << ", doesn't correspond to a chemical species";

      ingred_vec.push_back({cur.coef, chem_species_names[cur.species_idx]});
    }

    // now lets check if we can match ingredients
    auto search = ref_ingred_map.find(name);
    if (search == ref_ingred_map.end()) {
      continue;
    }

    const std::vector<CoefNamePair>& ref_vec = search->second;
    EXPECT_EQ(ingred_vec, ref_vec)
        << "there is a disagreement over the number of growth ingredients";
    n_comparisons++;
  }

  EXPECT_EQ(n_comparisons, ref_ingred_map.size())
      << "some of the entries in the reference ingred_map were skipped. This "
      << "is indicative of a problem.";
}

INSTANTIATE_TEST_SUITE_P(
    ,  // <- this comma is meaningful
    GrainSpeciesInfoTest, testing::Range(1, MAX_dust_species_VAL + 1),
    // adjust how the value is formatted in the test name
    [](const testing::TestParamInfo<GrainSpeciesInfoTest::ParamType>& info) {
      int val = info.param;
      std::string name = "DustSpeciesEq" + std::to_string(val);
      return name;
    });

// check the GrainSpeciesInfo object when constructed from dust_species
// parameters that hold extreme values
TEST(GrainSpeciesInfoTestMisc, DustSpeciesExtremeValues) {
  {
    unique_GrainSpeciesInfo_ptr ptr = make_unique_GrainSpeciesInfo(0);
    EXPECT_EQ(ptr->n_species, 0)
        << "GrainSpeciesInfo::n_species should be 0 when the dust_species "
        << "parameter is 0.";
    EXPECT_EQ(ptr->species_info, nullptr)
        << "GrainSpeciesInfo::species_info should be a nullptr when "
        << "dust_species parameter is 0.";
  }

  int invalid_dust_species_values[2] = {-42423, MAX_dust_species_VAL + 1};
  for (auto dust_species_param : invalid_dust_species_values) {
    unique_GrainSpeciesInfo_ptr ptr =
        make_unique_GrainSpeciesInfo(dust_species_param);
    EXPECT_LE(ptr->n_species, -1)
        << "GrainSpeciesInfo::n_species should be negative when the "
        << "dust_species parameter is " << dust_species_param << " (i.e. an "
        << "invalid value).";
    EXPECT_EQ(ptr->species_info, nullptr)
        << "GrainSpeciesInfo::species_info member should be a nullptr when the "
        << "dust_species parameter is " << dust_species_param << " (i.e. an "
        << "invalid value).";
  }

  int n_known_grain_species = number_known_grain_species();
  {
    unique_GrainSpeciesInfo_ptr ptr =
        make_unique_GrainSpeciesInfo(MAX_dust_species_VAL);
    EXPECT_EQ(ptr->n_species, n_known_grain_species)
        << "GrainSpeciesInfo::n_species should specify the number of grain "
        << "species known to grackle, " << n_known_grain_species
        << ", when the "
        << "dust_species parameter holds the max known value, "
        << MAX_dust_species_VAL << '.';
    EXPECT_NE(ptr->species_info, nullptr)
        << "GrainSpeciesInfo::species_info can't be a nullptr when "
        << "GrainSpeciesInfo::n_species is positive";
  }
}
