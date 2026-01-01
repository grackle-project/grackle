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
#include <memory>
#include <vector>
#include <string>
#include <string_view>

#include "LUT.hpp"
#include "dust/grain_species_info.hpp"
#include "utils/FrozenKeyIdxBiMap.hpp"

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

/// Returns name-index pairs for all grains from the list of Grackle Species
///
/// The index corresponds to the index in global species list
std::vector<SpNameIndexPair> grain_sp_info_from_species_list() {
  std::vector<SpNameIndexPair> out;
  int count = 0;

#define STRINGIFY_(NAME) #NAME
#define ENTRY(NAME)                                                            \
  if (str_ends_with(STRINGIFY_(NAME), "_dust")) {                              \
    out.push_back(SpNameIndexPair{STRINGIFY_(NAME), count});                   \
  }                                                                            \
  count++;
#include "field_data_evolved_species.def"
#undef ENTRY
#undef STRINGIFY_

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
  // NOLINTBEGIN(clang-analyzer-cplusplus.NewDeleteLeaks)
  grackle::impl::GrainSpeciesInfo* ptr = new grackle::impl::GrainSpeciesInfo;
  (*ptr) = grackle::impl::new_GrainSpeciesInfo(dust_species_param);
  return unique_GrainSpeciesInfo_ptr(ptr, GrainSpeciesInfoDeleter());
  // NOLINTEND(clang-analyzer-cplusplus.NewDeleteLeaks)
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

// the following test is somewhat redundant with the SpeciesLUTCompare tests
// and is more work to maintain
TEST_P(GrainSpeciesInfoTest, CheckOnlyGrainSpeciesLUTConsistency) {
  // construct a vector with an element for each entry in OnlyGrainSpeciesLUT
  std::vector<SpNameIndexPair> ref_l{
      {"MgSiO3_dust", OnlyGrainSpLUT::MgSiO3_dust},
      {"AC_dust", OnlyGrainSpLUT::AC_dust},
      {"SiM_dust", OnlyGrainSpLUT::SiM_dust},
      {"FeM_dust", OnlyGrainSpLUT::FeM_dust},
      {"Mg2SiO4_dust", OnlyGrainSpLUT::Mg2SiO4_dust},
      {"Fe3O4_dust", OnlyGrainSpLUT::Fe3O4_dust},
      {"SiO2_dust", OnlyGrainSpLUT::SiO2_dust},
      {"MgO_dust", OnlyGrainSpLUT::MgO_dust},
      {"FeS_dust", OnlyGrainSpLUT::FeS_dust},
      {"Al2O3_dust", OnlyGrainSpLUT::Al2O3_dust},
      {"ref_org_dust", OnlyGrainSpLUT::ref_org_dust},
      {"vol_org_dust", OnlyGrainSpLUT::vol_org_dust},
      {"H2O_ice_dust", OnlyGrainSpLUT::H2O_ice_dust},
  };

  const int n_species = grain_species_info_->n_species;
  for (int i = 0; i < n_species; i++) {
    const char* name = grackle::impl::FrozenKeyIdxBiMap_key_from_idx(
        &grain_species_info_->name_map, static_cast<std::uint16_t>(i));

    ASSERT_NE(name, nullptr);      // sanity check!
    EXPECT_EQ(i, ref_l[i].index);  // sanity check!

    // actual check:
    EXPECT_EQ(std::string(name), ref_l[i].name)
        << "the grain species associated with index " << i << " in the "
        << "GrainSpeciesInfo instance doesn't seem to be synchronized with "
        << "the OnlyGrainSpeciesLUT enumeration.";
  }
}

TEST_P(GrainSpeciesInfoTest, SublimationTemperature) {
  const int n_species = grain_species_info_->n_species;
  for (int i = 0; i < n_species; i++) {
    const char* name = grackle::impl::FrozenKeyIdxBiMap_key_from_idx(
        &grain_species_info_->name_map, static_cast<std::uint16_t>(i));
    ASSERT_NE(name, nullptr);  // sanity check!
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
    const char* name_cstr = grackle::impl::FrozenKeyIdxBiMap_key_from_idx(
        &grain_species_info_->name_map, static_cast<std::uint16_t>(i));
    ASSERT_NE(name_cstr, nullptr);  // sanity check!
    // actual check!
    std::string actual_name(name_cstr);
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
  int invalid_dust_species_values[3] = {-42423, 0, MAX_dust_species_VAL + 1};
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

// check that the grackle::impl::max_ingredients_per_grain_species constant
// indeed holds the maximum number of ingredients that a grain species can
// have.
TEST(GrainSpeciesInfoTestMisc, MaxIngredientsPerGrainSpecies) {
  // I really don't understand this linter issue
  // NOLINTBEGIN(clang-analyzer-core.CallAndMessage)

  // construct a GrainSpeciesInfo instance that holds values for every dust
  // grain species
  unique_GrainSpeciesInfo_ptr gsp_info =
      make_unique_GrainSpeciesInfo(MAX_dust_species_VAL);

  // go through and compute the number of ingredients for each species
  int n_grain_species = gsp_info->n_species;
  int max_ingredient_count = 0;
  for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
    max_ingredient_count =
        std::max(max_ingredient_count,
                 gsp_info->species_info[gsp_idx].n_growth_ingredients);
  }

  EXPECT_EQ(max_ingredient_count,
            grackle::impl::max_ingredients_per_grain_species)
      << "it appears that grackle::impl::max_ingredients_per_grain_species "
      << "should have a value of " << max_ingredient_count;
  // NOLINTEND(clang-analyzer-core.CallAndMessage)
}
