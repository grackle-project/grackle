//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define some basic tests of the experimental ratequery API
///
//===----------------------------------------------------------------------===//

#include <iterator>
#include <limits>
#include <set>
#include <string>
#include <gtest/gtest.h>
#include "grtestutils/iterator_adaptor.hpp"
#include "grtestutils/googletest/fixtures.hpp"

#include "grackle.h"

/// returns the rateid used to denote invalid rate names
grunstable_rateid_type get_invalid_rateid(const grtest::GrackleCtxPack& pack) {
  // although we don't use pack, yet a forthcoming refactor requires it
  return grunstable_ratequery_id(nullptr);
}

using SimpleRateQueryTest =
    grtest::ConfigPresetFixture<grtest::ChemPreset::primchem4_dustspecies3,
                                grtest::InitialUnitPreset::simple_z0>;

TEST_F(SimpleRateQueryTest, InvalidIthRate) {
  const char* name = grunstable_ith_rate(
      std::numeric_limits<unsigned long long>::max(), nullptr);
  EXPECT_EQ(name, nullptr);
}

TEST_F(SimpleRateQueryTest, EmptyNameQuery) {
  grunstable_rateid_type rateid = grunstable_ratequery_id("");
  EXPECT_EQ(rateid, get_invalid_rateid(pack));
}

TEST_F(SimpleRateQueryTest, InvalidNameQuery) {
  grunstable_rateid_type rateid = grunstable_ratequery_id("NotAValidName");
  EXPECT_EQ(rateid, get_invalid_rateid(pack));
}

TEST_F(SimpleRateQueryTest, PtrInvalidRateId) {
  double* ptr =
      grunstable_ratequery_get_ptr(pack.my_rates(), get_invalid_rateid(pack));
  EXPECT_EQ(ptr, nullptr);
}

// most of the remaining tests (and future planned tests) involve iterating
// through all rates made available via the ratequery interface. To make the
// tests themselves as easy to read as possible, we implate a C++-style
// iterator to wrap part of the interface

using ParametrizedRateQueryTest = grtest::ParametrizedConfigPresetFixture;

TEST_P(ParametrizedRateQueryTest, AllUnique) {
  std::set<std::string> name_set;
  std::set<grunstable_rateid_type> id_set;
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    ASSERT_TRUE(name_set.insert(pair.name).second)
        << "the name, \"" << pair.name << "\" appears more than once";
    ASSERT_TRUE(id_set.insert(pair.id).second)
        << "the id, " << pair.id << " appears more than once";
  }
}

TEST_P(ParametrizedRateQueryTest, ConsistentIDs) {
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    grunstable_rateid_type rateid = grunstable_ratequery_id(pair.name.c_str());
    EXPECT_EQ(rateid, pair.id);
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
                            InitialUnitPreset::simple_z0),
    ParamConf::SimplePreset(ChemPreset::primchem4_dustspecies3,
                            InitialUnitPreset::simple_z0)};

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, ParametrizedRateQueryTest,
    ::testing::ValuesIn(my_presets_));
