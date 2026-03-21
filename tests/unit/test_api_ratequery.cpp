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

#include <algorithm>  // std::min, std::max
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "grtestutils/iterator_adaptor.hpp"
#include "grtestutils/googletest/assertions.hpp"
#include "grtestutils/googletest/fixtures.hpp"

#include "grackle.h"
#include "status_reporting.h"

/// returns the rateid used to denote invalid rate names
grunstable_rateid_type get_invalid_rateid(const grtest::GrackleCtxPack& pack) {
  return grunstable_ratequery_id(pack.my_rates(), nullptr);
}

using SimpleRateQueryTest =
    grtest::ConfigPresetFixture<grtest::ChemPreset::primchem4_dustspecies3,
                                grtest::InitialUnitPreset::simple_z0>;

TEST_F(SimpleRateQueryTest, InvalidIthRate) {
  const char* name = grunstable_ith_rate(
      pack.my_rates(), std::numeric_limits<unsigned long long>::max(), nullptr);
  EXPECT_EQ(name, nullptr);
}

TEST_F(SimpleRateQueryTest, IthRateChecks) {
  chemistry_data_storage* my_rates = pack.my_rates();
  unsigned long long n_rates = grunstable_ratequery_nrates(my_rates);
  ASSERT_GT(n_rates, 0);  // <- sanity check!
  EXPECT_NE(nullptr, grunstable_ith_rate(my_rates, n_rates - 1, nullptr));
  EXPECT_EQ(nullptr, grunstable_ith_rate(my_rates, n_rates, nullptr));
  EXPECT_EQ(nullptr, grunstable_ith_rate(my_rates, n_rates + 1, nullptr));
}

TEST_F(SimpleRateQueryTest, EmptyNameQuery) {
  grunstable_rateid_type rateid = grunstable_ratequery_id(pack.my_rates(), "");
  EXPECT_EQ(rateid, get_invalid_rateid(pack));
}

TEST_F(SimpleRateQueryTest, InvalidNameQuery) {
  grunstable_rateid_type rateid =
      grunstable_ratequery_id(pack.my_rates(), "NotAValidName");
  EXPECT_EQ(rateid, get_invalid_rateid(pack));
}

// most of the remaining tests (and future planned tests) involve iterating
// through all rates made available via the ratequery interface. To make the
// tests themselves as easy to read as possible, we implate a C++-style
// iterator to wrap part of the interface

std::string stringify_prop_kind(enum grunstable_ratequery_prop_kind kind) {
  switch (kind) {
    case GRUNSTABLE_QPROP_NDIM:
      return "GRUNSTABLE_QPROP_NDIM";
    case GRUNSTABLE_QPROP_SHAPE:
      return "GRUNSTABLE_QPROP_SHAPE";
    case GRUNSTABLE_QPROP_MAXITEMSIZE:
      return "GRUNSTABLE_QPROP_MAXITEMSIZE";
    case GRUNSTABLE_QPROP_WRITABLE:
      return "GRUNSTABLE_QPROP_WRITABLE";
    case GRUNSTABLE_QPROP_DTYPE:
      return "GRUNSTABLE_QPROP_DTYPE";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

std::string stringify_type(enum grunstable_types kind) {
  switch (kind) {
    case GRUNSTABLE_TYPE_F64:
      return "GRUNSTABLE_TYPE_F64";
    case GRUNSTABLE_TYPE_STR:
      return "GRUNSTABLE_TYPE_STR";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

std::optional<enum grunstable_types> safe_type_enum_cast(long long val) {
  if (static_cast<long long>(GRUNSTABLE_TYPE_F64) == val) {
    return std::make_optional(GRUNSTABLE_TYPE_F64);
  } else if (static_cast<long long>(GRUNSTABLE_TYPE_STR) == val) {
    return std::make_optional(GRUNSTABLE_TYPE_STR);
  } else {
    return std::nullopt;
  }
}

/// return a vector of some invalid rateids
std::vector<grunstable_rateid_type> invalid_rateids(
    grtest::GrackleCtxPack& pack) {
  grunstable_rateid_type max_id =
      std::numeric_limits<grunstable_rateid_type>::lowest();
  grunstable_rateid_type min_id =
      std::numeric_limits<grunstable_rateid_type>::max();

  // iterate over (parameter-name, rate-id) pairs
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    max_id = std::max(max_id, static_cast<grunstable_rateid_type>(pair.id));
    min_id = std::min(min_id, static_cast<grunstable_rateid_type>(pair.id));
  }

  std::vector<grunstable_rateid_type> bad_ids;
  if (min_id > std::numeric_limits<grunstable_rateid_type>::lowest()) {
    bad_ids.push_back(min_id - 1);
  }
  if (max_id < std::numeric_limits<grunstable_rateid_type>::max()) {
    bad_ids.push_back(max_id + 1);
  }
  return bad_ids;
}

TEST_F(SimpleRateQueryTest, PropertyInvalidRateID) {
  std::vector<grunstable_rateid_type> invalid_ids = invalid_rateids(pack);
  if (invalid_ids.empty()) {
    GTEST_SKIP() << "unable to come up with known invalid rate ids to use for "
                 << "the test";
  }

  std::vector<enum grunstable_ratequery_prop_kind> prop_kinds{
      GRUNSTABLE_QPROP_NDIM, GRUNSTABLE_QPROP_SHAPE,
      GRUNSTABLE_QPROP_MAXITEMSIZE};
  std::vector<long long> buf;

  for (grunstable_rateid_type invalid_id : invalid_ids) {
    for (enum grunstable_ratequery_prop_kind kind : prop_kinds) {
      constexpr std::size_t BUF_LEN = 20;  // <- arbitrarily large value
      constexpr long long DEFAULT_VAL = -25634634LL;  // <- arbitrary value
      buf.assign(BUF_LEN, DEFAULT_VAL);
      ASSERT_GR_ERR(grunstable_ratequery_prop(pack.my_rates(), invalid_id, kind,
                                              buf.data()))
          << "executed with invalid_id=" << invalid_id
          << ", kind=" << stringify_prop_kind(kind);
      EXPECT_THAT(buf, testing::Each(testing::Eq(DEFAULT_VAL)))
          << "grunstable_ratequery_prop mutated the ptr even though it "
          << "reported a failure. It was called with an invalid id of "
          << invalid_id << " and a kind of " << stringify_prop_kind(kind);
    }
  }
}

using ParametrizedRateQueryTest = grtest::ParametrizedConfigPresetFixture;

TEST_P(ParametrizedRateQueryTest, InvalidIdCollision) {
  grunstable_rateid_type invalid_id =
      grunstable_ratequery_id(pack.my_rates(), nullptr);
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    EXPECT_NE(invalid_id, pair.id)
        << "there is a collision between the canonical invalid id "
        << "and the id associated with the \"" << pair.name << "\" rate";
  }
}

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
    grunstable_rateid_type rateid =
        grunstable_ratequery_id(pack.my_rates(), pair.name.c_str());
    EXPECT_EQ(rateid, pair.id);
  }
}

TEST_P(ParametrizedRateQueryTest, Property) {
  std::vector<long long> buf;
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    constexpr long long DEFAULT_VAL = -25634634LL;  // <- arbitrary value

    // check ndim
    long long ndim = DEFAULT_VAL;
    EXPECT_GR_SUCCESS(grunstable_ratequery_prop(pack.my_rates(), pair.id,
                                                GRUNSTABLE_QPROP_NDIM, &ndim))
        << "for " << pair;
    ASSERT_GE(ndim, 0LL) << "for " << pair;

    // check shape
    buf.assign((ndim == 0LL) ? 1 : ndim, DEFAULT_VAL);
    EXPECT_GR_SUCCESS(grunstable_ratequery_prop(
        pack.my_rates(), pair.id, GRUNSTABLE_QPROP_SHAPE, buf.data()))
        << "for " << pair;
    if (ndim == 0LL) {
      EXPECT_EQ(buf[0], DEFAULT_VAL)
          << "the buffer passed to grunstable_ratequery_prop was unexpectedly "
          << "modified while querying the shape for the rate " << pair
          << ". It shouldn't be modified since ndim=0.";
    } else {
      EXPECT_THAT(buf, testing::Each(testing::Gt(0)))
          << "buf holds the shape queried for " << pair;
    }

    long long tmp = DEFAULT_VAL;
    EXPECT_GR_SUCCESS(grunstable_ratequery_prop(pack.my_rates(), pair.id,
                                                GRUNSTABLE_QPROP_DTYPE, &tmp))
        << "for " << pair;
    std::optional<enum grunstable_types> dtype_maybe = safe_type_enum_cast(tmp);
    if (!dtype_maybe.has_value()) {
      GTEST_FAIL() << "Error coercing " << tmp << ", the dtype for " << pair
                   << ", to an enum value";
    }
    enum grunstable_types dtype = dtype_maybe.value();

    long long maxitemsize = DEFAULT_VAL;
    EXPECT_GR_SUCCESS(grunstable_ratequery_prop(
        pack.my_rates(), pair.id, GRUNSTABLE_QPROP_MAXITEMSIZE, &maxitemsize))
        << "for " << pair;
    if (dtype == GRUNSTABLE_TYPE_F64) {
      EXPECT_EQ(maxitemsize, sizeof(double)) << "for " << pair;
    } else {
      EXPECT_GT(maxitemsize, 0) << "for " << pair;
    }

    long long writable = DEFAULT_VAL;
    EXPECT_GR_SUCCESS(grunstable_ratequery_prop(
        pack.my_rates(), pair.id, GRUNSTABLE_QPROP_WRITABLE, &writable))
        << "for " << pair;
    EXPECT_THAT(writable, ::testing::AnyOf(0LL, 1LL)) << "for " << pair;
  }
}

/// summarizes details about rate properties
struct RateProperties {
  std::vector<long long> shape;
  std::size_t maxitemsize;
  enum grunstable_types dtype;
  bool writable;

  bool operator==(const RateProperties& other) const {
    return maxitemsize == other.maxitemsize && shape == other.shape &&
           dtype == other.dtype && writable == other.writable;
  }

  long long n_items() const {
    long long n_items = 1LL;
    for (std::size_t i = 0; i < shape.size(); i++) {
      n_items *= shape[i];
    }
    return n_items;  // this is correct even for scalars
  }

  /// teach googletest how to print this type
  friend void PrintTo(const RateProperties& props, std::ostream* os) {
    *os << "{shape={";
    for (std::size_t i = 0; i < props.shape.size(); i++) {
      if (i != 0) {
        *os << ", ";
      }
      *os << props.shape[i];
    }
    *os << "}, maxitemsize=" << props.maxitemsize
        << ", dtype=" << stringify_type(props.dtype)
        << ", writable=" << props.writable << '}';
  }
};

/// construct a RateProperties instance for the specified rate
///
/// returns an empty optional if any there are any issues
///
/// @note
/// After we query each property, we encode an extra sanity-check. We **only**
/// encode this in case there are other underlying issues (so that we don't end
/// up with an extremely crazy set of errors). The Property tests should
/// generally provide more details about these tests. To be clear, user-code
/// should never need to include these sanity check!
std::optional<RateProperties> try_query_RateProperties(
    chemistry_data_storage* my_rates, grunstable_rateid_type rateid) {
  long long ndim = -1LL;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_NDIM,
                                &ndim) != GR_SUCCESS) {
    return std::nullopt;
  }
  if (ndim < 0LL) {
    return std::nullopt;  // sanity-check failed!
  }

  std::vector<long long> shape;
  if (ndim > 0LL) {
    shape.assign(ndim, 0LL);
    if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_SHAPE,
                                  shape.data()) != GR_SUCCESS) {
      return std::nullopt;
    }
  }
  if (std::count_if(shape.begin(), shape.end(),
                    [](long long x) { return x <= 0; }) > 0) {
    return std::nullopt;  // sanity check failed!
  }

  long long maxitemsize = -1LL;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_MAXITEMSIZE,
                                &maxitemsize) != GR_SUCCESS) {
    return std::nullopt;
  }
  if (maxitemsize <= 0LL) {
    return std::nullopt;  // sanity check failed!
  }

  long long dtype_tmp;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_DTYPE,
                                &dtype_tmp) != GR_SUCCESS) {
    return std::nullopt;
  }
  std::optional<enum grunstable_types> dtype = safe_type_enum_cast(dtype_tmp);
  if (!dtype.has_value()) {
    return std::nullopt;  // sanity check failed!
  }

  long long writable;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_WRITABLE,
                                &writable) != GR_SUCCESS) {
    return std::nullopt;
  }
  if ((writable != 0LL) && (writable != 1LL)) {
    return std::nullopt;
  }

  return {RateProperties{shape, static_cast<std::size_t>(maxitemsize),
                         dtype.value(), static_cast<bool>(writable)}};
}

// returns a value that differs from the input
static double remap_value(double in) {
  if (in == 0.0 || !std::isfinite(in)) {
    return 1.0;
  }
  return -in;
}

TEST_P(ParametrizedRateQueryTest, SetAndGet) {
  std::vector<double> initial_buf;
  std::vector<double> post_update_buf;

  // iterate over every known (rate-name, rate-id) pair
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    // get the properties associated with the current rate
    std::optional<RateProperties> maybe_props =
        try_query_RateProperties(pack.my_rates(), pair.id);
    if (!maybe_props.has_value()) {
      GTEST_FAIL()
          << "something went wrong while trying to lookup the properties for "
          << "the " << pair << " rate.";
    }
    RateProperties props = maybe_props.value();
    long long n_items = props.n_items();

    if (!props.writable) {
      continue;
    }

    // load in data associated with the current rate
    initial_buf.assign(n_items, NAN);
    ASSERT_GR_SUCCESS(grunstable_ratequery_get_f64(pack.my_rates(), pair.id,
                                                   initial_buf.data()))
        << "for " << pair;

    // overwrite each entry with a different value
    for (long long i = 0; i < n_items; i++) {
      initial_buf[i] = remap_value(initial_buf[i]);
    }

    // write the new values back to my_rates
    EXPECT_GR_SUCCESS(grunstable_ratequery_set_f64(pack.my_rates(), pair.id,
                                                   initial_buf.data()))
        << "for " << pair;

    // finally, lets check that the values actually got updated
    post_update_buf.assign(n_items, NAN);
    EXPECT_GR_SUCCESS(grunstable_ratequery_get_f64(pack.my_rates(), pair.id,
                                                   post_update_buf.data()))
        << "for " << pair;

    EXPECT_EQ(initial_buf, post_update_buf) << "for " << pair;
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

// now, we are going to check on some well-established rates
// -> this may not be the most maintainable test (we may want to re-evaluate it
//    in the future)

enum RateKind { scalar_f64, simple_1d_rate, k13dd };

static RateProperties RateProperties_from_RateKind(
    const chemistry_data* my_chemistry, RateKind kind) {
  const enum grunstable_types f64dtype = GRUNSTABLE_TYPE_F64;
  switch (kind) {
    case RateKind::scalar_f64: {
      std::vector<long long> shape = {};  // <-- intentionally empty
      return RateProperties{shape, sizeof(double), f64dtype, true};
    }
    case RateKind::simple_1d_rate: {
      std::vector<long long> shape = {my_chemistry->NumberOfTemperatureBins};
      return RateProperties{shape, sizeof(double), f64dtype, true};
    }
    case RateKind::k13dd: {
      std::vector<long long> shape = {my_chemistry->NumberOfTemperatureBins *
                                      14};
      return RateProperties{shape, sizeof(double), f64dtype, true};
    }
  }
  GR_INTERNAL_UNREACHABLE_ERROR()
}

/// returns a map between known rate names and the rate kind
std::map<std::string, RateKind> known_rates() {
  std::map<std::string, RateKind> out;
  // radiative rates:
  for (int i = 24; i < 32; i++) {
    out.insert({"k" + std::to_string(i), RateKind::scalar_f64});
  }

  // standard collisional rates
  for (int i = 1; i < 24; i++) {
    out.insert({"k" + std::to_string(i), RateKind::simple_1d_rate});
  }
  for (int i = 50; i < 58; i++) {
    out.insert({"k" + std::to_string(i), RateKind::simple_1d_rate});
  }
  for (int i = 125; i < 154; i++) {
    out.insert({"k" + std::to_string(i), RateKind::simple_1d_rate});
  }

  // metal chemistry rates
  for (int i = 11; i < 55; i++) {
    out.insert({"kz" + std::to_string(i), RateKind::simple_1d_rate});
  }

  out.insert({"k13dd", RateKind::k13dd});

  return out;
}

using KnownRateQueryTest =
    grtest::ConfigPresetFixture<grtest::ChemPreset::primchem4_dustspecies3,
                                grtest::InitialUnitPreset::simple_z0>;

TEST_F(KnownRateQueryTest, CheckProperties) {
  const std::map<std::string, RateKind> known_rate_map = known_rates();

  // iterate over every known {parameter-name, key-id} pair
  for (const grtest::NameIdPair pair : grtest::RateQueryRange(pack)) {
    // check if the current rateid is in known_rate_map
    std::map<std::string, RateKind>::const_iterator search =
        known_rate_map.find(pair.name);
    if (search == known_rate_map.end()) {
      continue;  // the rateid is **NOT** in known_rate_map
    }
    // construct the expected properties
    RateProperties expected_props =
        RateProperties_from_RateKind(pack.my_chemistry(), search->second);

    // load the actual props
    std::optional<RateProperties> maybe_actual_props =
        try_query_RateProperties(pack.my_rates(), pair.id);
    if (!maybe_actual_props.has_value()) {
      GTEST_FAIL()
          << "something went wrong while trying to lookup the properties for "
          << "the " << pair << " rate.";
    }
    RateProperties actual_props = maybe_actual_props.value();

    EXPECT_EQ(expected_props, actual_props)
        << "this mismatch in the queried properties is for the " << pair
        << " rate.";
  }

  // it might be useful to repeat the tests using another chemistry_data
  // instance with a different NumberOfTemperatureBins value
}
