//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// tests the @ref grackle::impl::FrozenKeyIdxBiMap
///
//===----------------------------------------------------------------------===//

#include <iostream>  // needed to teach googletest how to print
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "support/status_reporting.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "grackle.h"

namespace grackle::impl {

// teach GoogleTest how to print grackle::impl::BiMapMode for shorter test names
void PrintTo(const BiMapMode& mode, std::ostream* os) {
  switch (mode) {
    case BiMapMode::COPIES_KEYDATA:
      *os << "BiMapMode::COPIES_KEYDATA";
      return;
    case BiMapMode::REFS_KEYDATA:
      *os << "BiMapMode::REFS_KEYDATA";
      return;
  }
  GR_INTERNAL_ERROR("should not be reachable");
}

}  // namespace grackle::impl

// allow ourselves to write Optional(val) rather than testing::Optional(val)
using testing::Optional;

// this top test was introduced to provide a more concrete example
// of how we might use FrozenKeyIdxBiMap

TEST(FrozenKeyIdxBiMap, FullExample) {
  // THE SCENARIO: we have a list of unique ordered strings
  //
  // We are going build a FrozenKeyIdxBiMap instance from the following list.
  // The resulting object is a bidirectional map that can both:
  // 1. map a string to its index (at the time of construction) in the list.
  //    - example: "HII" is mapped to 2
  //    - example: "O2II" is mapped to 33
  // 2. perform the reverse mapping (i.e. index -> string)
  //    - example: 2 is mapped to "HII"
  //    - example: 33 is mapped to "O2II"
  //
  // It's worth emphasizing that the mapping is frozen when its constructed &
  // contents can't be changed (even if you reorder the original)
  const char* keys[34] = {
      "e",   "HI",   "HII", "HeI",  "HeII",  "HeIII", "HM",   "H2I",   "H2II",
      "DI",  "DII",  "HDI", "DM",   "HDII",  "HeHII", "CI",   "CII",   "CO",
      "CO2", "OI",   "OH",  "H2O",  "O2",    "SiI",   "SiOI", "SiO2I", "CH",
      "CH2", "COII", "OII", "OHII", "H2OII", "H3OII", "O2II"};

  namespace grimpl = grackle::impl;

  // PART 1: build a FrozenKeyIdxBiMap from this list
  // the 3rd argument tells the string to make copies of each string
  grimpl::FrozenKeyIdxBiMap m = grimpl::FrozenKeyIdxBiMap::create(
      keys, 34, grimpl::BiMapMode::COPIES_KEYDATA);

  // before we use it, we should confirm the constructor succeeded
  if (!grimpl::FrozenKeyIdxBiMap_is_ok(&m)) {
    FAIL() << "creation of the m failed unexpectedly";
  }

  // PART 2: let's show some examples of lookups from names

  // Equivalent Python:  `2 == m["HII"]`
  EXPECT_THAT(m.find("HII"), Optional(2));

  // Equivalent Python/idiomatic C++:  `33 == m["O2II"]`
  EXPECT_THAT(m.find("O2II"), Optional(33));

  // for unknown key, returns AccessRslt{has_value=false, value=<garbage>}
  EXPECT_EQ(m.find("Dummy"), std::nullopt);

  // PART 3: let's show the reverse of the previous lookups
  EXPECT_STREQ("HII", m.inverse_find(2));
  EXPECT_STREQ("O2II", m.inverse_find(33));

  // Behavior is again well-defined when passing an invalid index
  EXPECT_EQ(nullptr, m.inverse_find(131));

  // PART 4: We can also query the length
  EXPECT_EQ(m.size(), 34);
}

// validate basic operations for an empty bimap
TEST(FrozenKeyIdxBiMap, EmptyBasicOps) {
  grackle::impl::FrozenKeyIdxBiMap m = grackle::impl::FrozenKeyIdxBiMap::create(
      nullptr, 0, grackle::impl::BiMapMode::COPIES_KEYDATA);
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&m))
      << "construction of a FrozenKeyIdxBiMap unexpectedly failed";

  EXPECT_EQ(m.size(), 0) << "an empty mapping should have a size of 0";

  EXPECT_EQ(m.find("key"), std::nullopt)
      << "key lookup should always fail for an empty mapping";

  EXPECT_EQ(nullptr, m.inverse_find(0))
      << "index lookup should always fail for an empty mapping";
}

// validate behavior of clone for an empty bimap
TEST(FrozenKeyIdxBiMap, EmptyClone) {
  // make the original
  grackle::impl::FrozenKeyIdxBiMap m = grackle::impl::FrozenKeyIdxBiMap::create(
      nullptr, 0, grackle::impl::BiMapMode::COPIES_KEYDATA);
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&m))
      << "construction of a FrozenKeyIdxBiMap unexpectedly failed";

  // make the clone
  grackle::impl::FrozenKeyIdxBiMap m_clone = m.clone();

  bool success = grackle::impl::FrozenKeyIdxBiMap_is_ok(&m_clone);

  if (!success) {
    FAIL() << "cloning an empty mapping failed!";
  }
}

class BiMapCreate : public testing::TestWithParam<grackle::impl::BiMapMode> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};

TEST_P(BiMapCreate, Simple) {
  const char* keys[] = {"denisty", "internal_energy"};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(keys, 2, GetParam());

  EXPECT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(BiMapCreate, LongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap_detail::KEYLEN_MAX, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(keys, 2, GetParam());

  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(BiMapCreate, TooLongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap_detail::KEYLEN_MAX + 1, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(BiMapCreate, 0LenKey) {
  const char* keys[2] = {"density", ""};
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(BiMapCreate, NullptrWithPosCount) {
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(nullptr, 1, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(BiMapCreate, NotNull0KeyCount) {
  const char* keys[] = {"denisty", "internal_energy"};
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::FrozenKeyIdxBiMap::create(keys, 0, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

INSTANTIATE_TEST_SUITE_P(
    ,  // <- leaving Instantiation name empty
    BiMapCreate,
    testing::Values(grackle::impl::BiMapMode::REFS_KEYDATA,
                    grackle::impl::BiMapMode::COPIES_KEYDATA));

/// helper function to initialize a map from a vector
grackle::impl::FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(
    const std::vector<std::string>& vec_, grackle::impl::BiMapMode mode) {
  std::size_t key_count = vec_.size();

  // create a vector of pointers
  std::vector<const char*> key_ptr_l(key_count, nullptr);
  for (std::size_t i = 0; i < key_count; i++) {
    key_ptr_l[i] = vec_[i].c_str();
  }

  return grackle::impl::FrozenKeyIdxBiMap::create(key_ptr_l.data(), key_count,
                                                  mode);
}

class BiMapGeneral : public testing::TestWithParam<grackle::impl::BiMapMode> {
protected:
  std::vector<std::string> ordered_keys;
  grackle::impl::FrozenKeyIdxBiMap* bimap_p = nullptr;

  // we use SetUp/Teardown instead of constructor and destructor so we can
  // perform some sanity checks with ASSERTIONs
  void SetUp() override {
    ordered_keys =
        std::vector<std::string>{"internal_energy", "density", "metal_density"};

    grackle::impl::BiMapMode mode = GetParam();
    grackle::impl::FrozenKeyIdxBiMap tmp =
        new_FrozenKeyIdxBiMap(ordered_keys, mode);
    ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));

    bimap_p = new grackle::impl::FrozenKeyIdxBiMap;
    (*bimap_p) = std::move(tmp);
  }

  void TearDown() override {
    if (bimap_p != nullptr) {
      delete bimap_p;
    }
  }

  bool ReusesOriginalKeyPtrs(const grackle::impl::FrozenKeyIdxBiMap* p) const {
    for (int i = 0; i < 3; i++) {
      const char* orig_key_ptr = ordered_keys[i].c_str();
      if (p->inverse_find(i) != orig_key_ptr) {
        return false;
      }
    }
    return true;
  }
};

TEST_P(BiMapGeneral, FindContainedKey) {
  EXPECT_THAT(bimap_p->find("density"), Optional(1));
  EXPECT_THAT(bimap_p->find("internal_energy"), Optional(0));
  EXPECT_THAT(bimap_p->find("metal_density"), Optional(2));
}

TEST_P(BiMapGeneral, FindAbsentKey) {
  EXPECT_EQ(bimap_p->find("notAKey"), std::nullopt);
}

TEST_P(BiMapGeneral, FindForbiddenKeys) {
  // let's veryify that trying to find forbidden keys works properly
  // -> the fact that they are forbidden means that they are always absent

  EXPECT_EQ(bimap_p->find(""), std::nullopt);

  std::string key(grackle::impl::bimap_detail::KEYLEN_MAX + 1, 'A');
  EXPECT_EQ(bimap_p->find(key.data()), std::nullopt);
}

TEST_P(BiMapGeneral, KeyFromIdxInvalidIdx) {
  EXPECT_EQ(bimap_p->inverse_find(3), nullptr);
}

TEST_P(BiMapGeneral, KeyFromIdxValidIdx) {
  EXPECT_EQ(std::string(bimap_p->inverse_find(2)),
            std::string("metal_density"));
  EXPECT_EQ(std::string(bimap_p->inverse_find(1)), std::string("density"));
  EXPECT_EQ(std::string(bimap_p->inverse_find(0)),
            std::string("internal_energy"));

  // check whether the bimap is using pointers to the keys used during init
  if (GetParam() == grackle::impl::BiMapMode::REFS_KEYDATA) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(bimap_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(bimap_p));
  }
}

TEST_P(BiMapGeneral, Clone) {
  grackle::impl::FrozenKeyIdxBiMap clone = bimap_p->clone();
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&clone));
  grackle::impl::FrozenKeyIdxBiMap* clone_p = &clone;

  // for the sake of robustly checking everything, we delete bimap_p
  delete bimap_p;
  bimap_p = nullptr;

  EXPECT_THAT(clone_p->find("internal_energy"), Optional(0));
  EXPECT_EQ(clone_p->find("notAKey"), std::nullopt);

  EXPECT_EQ(clone_p->inverse_find(3), nullptr);
  EXPECT_STREQ(clone_p->inverse_find(1), "density");

  // check whether the clone is using pointers to the keys used during init
  if (GetParam() == grackle::impl::BiMapMode::REFS_KEYDATA) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(clone_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(clone_p));
  }
}

INSTANTIATE_TEST_SUITE_P(
    ,  // <- leaving Instantiation name empty
    BiMapGeneral,
    testing::Values(grackle::impl::BiMapMode::REFS_KEYDATA,
                    grackle::impl::BiMapMode::COPIES_KEYDATA));
