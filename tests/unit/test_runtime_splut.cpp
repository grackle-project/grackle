//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// test logic pertaining to the species lut
///
//===----------------------------------------------------------------------===//

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../grtestutils/RuntimeSpLUT.hpp"
#include "LUT.hpp"

/*
/// Returns a SpNameIndexPair for every species for which p returns true
///
/// p is a unary function that accepts a string_view and returns a boolean
///
/// @note
/// I think we'll be able to dispose of this in the near future
template <typename UnaryPred>
std::vector<std::string> selected_names_(UnaryPred p) {
  std::vector<std::string> out;
#define STRINGIFY_(NAME) #NAME
#define ENTRY(N) if (p(STRINGIFY_(N))) { out.push_back(STRINGIFY_(N)); }

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
*/

TEST(RuntimeSpLUT, Empty) {
  grtest::RuntimeSpLUT lut;
  EXPECT_EQ(lut.size(), 0);
}

TEST(RuntimeSpLUT, PrimordialChemistry1) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(1, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 6);
  EXPECT_THAT(lut.find("e"), ::testing::Optional(SpLUT::e));
  EXPECT_THAT(lut.find("HI"), ::testing::Optional(SpLUT::HI));
  EXPECT_THAT(lut.find("HII"), ::testing::Optional(SpLUT::HII));
  EXPECT_THAT(lut.find("HeI"), ::testing::Optional(SpLUT::HeI));
  EXPECT_THAT(lut.find("HeII"), ::testing::Optional(SpLUT::HeII));
  EXPECT_THAT(lut.find("HeIII"), ::testing::Optional(SpLUT::HeIII));
}

TEST(RuntimeSpLUT, PrimordialChemistry2) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(2, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 9);
  EXPECT_THAT(lut.find("e"), ::testing::Optional(SpLUT::e));
  EXPECT_THAT(lut.find("HI"), ::testing::Optional(SpLUT::HI));
  EXPECT_THAT(lut.find("HII"), ::testing::Optional(SpLUT::HII));
  EXPECT_THAT(lut.find("HeI"), ::testing::Optional(SpLUT::HeI));
  EXPECT_THAT(lut.find("HeII"), ::testing::Optional(SpLUT::HeII));
  EXPECT_THAT(lut.find("HeIII"), ::testing::Optional(SpLUT::HeIII));

  EXPECT_THAT(lut.find("HM"), ::testing::Optional(SpLUT::HM));
  EXPECT_THAT(lut.find("H2I"), ::testing::Optional(SpLUT::H2I));
  EXPECT_THAT(lut.find("H2II"), ::testing::Optional(SpLUT::H2II));
}

TEST(RuntimeSpLUT, PrimordialChemistry3) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(3, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 12);
  EXPECT_THAT(lut.find("e"), ::testing::Optional(SpLUT::e));
  EXPECT_THAT(lut.find("HI"), ::testing::Optional(SpLUT::HI));
  EXPECT_THAT(lut.find("HII"), ::testing::Optional(SpLUT::HII));
  EXPECT_THAT(lut.find("HeI"), ::testing::Optional(SpLUT::HeI));
  EXPECT_THAT(lut.find("HeII"), ::testing::Optional(SpLUT::HeII));
  EXPECT_THAT(lut.find("HeIII"), ::testing::Optional(SpLUT::HeIII));

  EXPECT_THAT(lut.find("HM"), ::testing::Optional(SpLUT::HM));
  EXPECT_THAT(lut.find("H2I"), ::testing::Optional(SpLUT::H2I));
  EXPECT_THAT(lut.find("H2II"), ::testing::Optional(SpLUT::H2II));

  EXPECT_THAT(lut.find("DI"), ::testing::Optional(SpLUT::DI));
  EXPECT_THAT(lut.find("DII"), ::testing::Optional(SpLUT::DII));
  EXPECT_THAT(lut.find("HDI"), ::testing::Optional(SpLUT::HDI));
}

TEST(RuntimeSpLUT, PrimordialChemistry4) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 15);
  EXPECT_THAT(lut.find("e"), ::testing::Optional(SpLUT::e));
  EXPECT_THAT(lut.find("HI"), ::testing::Optional(SpLUT::HI));
  EXPECT_THAT(lut.find("HII"), ::testing::Optional(SpLUT::HII));
  EXPECT_THAT(lut.find("HeI"), ::testing::Optional(SpLUT::HeI));
  EXPECT_THAT(lut.find("HeII"), ::testing::Optional(SpLUT::HeII));
  EXPECT_THAT(lut.find("HeIII"), ::testing::Optional(SpLUT::HeIII));

  EXPECT_THAT(lut.find("HM"), ::testing::Optional(SpLUT::HM));
  EXPECT_THAT(lut.find("H2I"), ::testing::Optional(SpLUT::H2I));
  EXPECT_THAT(lut.find("H2II"), ::testing::Optional(SpLUT::H2II));

  EXPECT_THAT(lut.find("DI"), ::testing::Optional(SpLUT::DI));
  EXPECT_THAT(lut.find("DII"), ::testing::Optional(SpLUT::DII));
  EXPECT_THAT(lut.find("HDI"), ::testing::Optional(SpLUT::HDI));

  EXPECT_THAT(lut.find("DM"), ::testing::Optional(SpLUT::DM));
  EXPECT_THAT(lut.find("HDII"), ::testing::Optional(SpLUT::HDII));
  EXPECT_THAT(lut.find("HeHII"), ::testing::Optional(SpLUT::HeHII));
}
