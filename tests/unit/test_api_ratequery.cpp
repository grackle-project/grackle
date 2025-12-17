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

#include <limits>
#include <gtest/gtest.h>
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
