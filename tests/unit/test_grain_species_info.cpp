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
#include <ostream>
#include <vector>
#include <string>
#include <string_view>

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
  grackle::impl::GrainSpeciesInfo* ptr = new grackle::impl::GrainSpeciesInfo;
  (*ptr) = grackle::impl::new_GrainSpeciesInfo(dust_species_param);
  return unique_GrainSpeciesInfo_ptr(ptr, GrainSpeciesInfoDeleter());
}

}  // anonymous namespace

// check that OnlyGrainSpLUT contains the entries for every known grain species
TEST(OnlyGrainSpLUTTest, CheckNumEntries) {
  ASSERT_EQ(OnlyGrainSpLUT::NUM_ENTRIES, number_known_grain_species());
}

/// represents the dust_species parameter of chemistry_data
///
/// @note
/// This **ONLY** exists to make the test names more concise
enum class DustSpeciesParam : int {
  DustSpeciesEq1 = 1,
  DustSpeciesEq2 = 2,
  DustSpeciesEq3 = 3,
};

// teach googletest how to print DustSpeciesParam
void PrintTo(const DustSpeciesParam dust_species_param, std::ostream* os) {
  *os << "DustSpeciesEq" << static_cast<int>(dust_species_param);
}

/// Define a fixture for running parameterized tests of the GrainSpeciesInfo
/// machinery. The tests are parameterized by the dust_species parameter.
class GrainSpeciesInfoTest : public testing::TestWithParam<DustSpeciesParam> {
protected:
  /// set up the grain_species_info pointer
  ///
  /// @note
  /// We **ONLY** perform setup in this method (rather than in a default
  /// constructor) because we want to perform some basic sanity checks
  void SetUp() override {
    int dust_species = static_cast<int>(GetParam());
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

TEST_P(GrainSpeciesInfoTest, OnlyGrainSpeciesLUT) {
  // this test seeks to check consistency with the LUT
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

TEST_P(GrainSpeciesInfoTest, SublimationTemp) {
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

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, GrainSpeciesInfoTest,
    testing::Values(DustSpeciesParam::DustSpeciesEq1,
                    DustSpeciesParam::DustSpeciesEq2,
                    DustSpeciesParam::DustSpeciesEq3));

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

// check that the grackle::impl::max_ingredients_per_grain_species constant
// indeed holds the maximum number of ingredients that a grain species can
// have.
TEST(GrainSpeciesInfoTestMisc, MaxIngredientsPerGrainSpecies) {
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
}
