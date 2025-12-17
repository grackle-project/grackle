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

// will be implemented (in a much more robust manner) in the near future
unsigned long long grunstable_ratequery_nrates(
    const chemistry_data_storage* my_rates) {
  // current implementation is stupid! (in future, will use my_rates)
  unsigned long long i = 0;
  while (nullptr != grunstable_ith_rate(i, nullptr)) {
    i++;
  }
  return i;
}

// most of the remaining tests (and future planned tests) involve iterating
// through all rates made available via the ratequery interface. To make the
// tests themselves as easy to read as possible, we implate a C++-style
// iterator to wrap part of the interface

struct RateNameIdPair {
  std::string name;
  grunstable_rateid_type rateid;
};

class RQIterator {
  chemistry_data_storage* my_rates_;
  unsigned long long counter_;
  unsigned long long n_rates_;
  RateNameIdPair pair_;

  RQIterator& update_pair_and_ret_(unsigned long long val) {
    if (val < n_rates_) {
      pair_.name = std::string(grunstable_ith_rate(val, &pair_.rateid));
    }
    return *this;
  }

public:
  using iterator_category = std::input_iterator_tag;
  using value_type = RateNameIdPair;
  using difference_type = std::ptrdiff_t;
  using pointer = const RateNameIdPair*;
  using reference = const RateNameIdPair;

  RQIterator(unsigned long long counter, unsigned long long n_rates,
             chemistry_data_storage* my_rates)
      : my_rates_(my_rates), counter_(counter), n_rates_(n_rates) {
    update_pair_and_ret_(counter);
  }

  bool operator==(RQIterator other) const {
    return (counter_ == other.counter_) && (my_rates_ == other.my_rates_);
  }

  bool operator!=(RQIterator other) const { return !(*this == other); }
  reference operator*() const { return pair_; }
  RQIterator& operator++() { return update_pair_and_ret_(++counter_); }

  RQIterator operator++(int) {
    RQIterator ret = *this;
    ++(*this);
    return ret;
  }
};

// used for creating the iterator and within range-based for-loops
class RateQueryRange {
  grtest::GrackleCtxPack& pack_;
  long long n_rates_;

public:
  explicit RateQueryRange(grtest::GrackleCtxPack& pack)
      : pack_(pack), n_rates_(grunstable_ratequery_nrates(pack.my_rates())) {}

  RQIterator begin() { return RQIterator(0, n_rates_, pack_.my_rates()); }
  RQIterator end() { return RQIterator(n_rates_, n_rates_, pack_.my_rates()); }
};

using ParametrizedRateQueryTest = grtest::ParametrizedConfigPresetFixture;

TEST_P(ParametrizedRateQueryTest, AllUnique) {
  std::set<std::string> name_set;
  std::set<grunstable_rateid_type> id_set;
  for (const RateNameIdPair pair : RateQueryRange(pack)) {
    ASSERT_TRUE(name_set.insert(pair.name).second)
        << "the name, \"" << pair.name << "\" appears more than once";
    ASSERT_TRUE(id_set.insert(pair.rateid).second)
        << "the id, " << pair.rateid << " appears more than once";
  }
}

TEST_P(ParametrizedRateQueryTest, ConsistentIDs) {
  for (const RateNameIdPair pair : RateQueryRange(pack)) {
    grunstable_rateid_type rateid = grunstable_ratequery_id(pair.name.c_str());
    EXPECT_EQ(rateid, pair.rateid);
  }
}

using grtest::ChemPreset;
using grtest::FullConfPreset;
using grtest::InitialUnitPreset;

INSTANTIATE_TEST_SUITE_P(
    /* 1st arg is intentionally empty */, ParametrizedRateQueryTest,
    ::testing::Values(
        FullConfPreset{ChemPreset::primchem0, InitialUnitPreset::simple_z0},
        FullConfPreset{ChemPreset::primchem1, InitialUnitPreset::simple_z0},
        FullConfPreset{ChemPreset::primchem2, InitialUnitPreset::simple_z0},
        FullConfPreset{ChemPreset::primchem3, InitialUnitPreset::simple_z0},
        FullConfPreset{ChemPreset::primchem4_dustspecies3,
                       InitialUnitPreset::simple_z0}));
