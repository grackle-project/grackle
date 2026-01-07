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

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "support/FrozenKeyIdxBiMap.hpp"
#include "grackle.h"

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
  grimpl::FrozenKeyIdxBiMap m = grimpl::new_FrozenKeyIdxBiMap(
      keys, 34, grimpl::BiMapMode::COPIES_KEYDATA);

  // before we use it, we should confirm the constructor succeeded
  if (!grimpl::FrozenKeyIdxBiMap_is_ok(&m)) {
    FAIL() << "creation of the m failed unexpectedly";
  }

  // PART 2: let's show some examples of lookups from names

  // Equivalent Python/idiomatic C++:  `2 == m["HII"]`
  EXPECT_EQ(2, grimpl::FrozenKeyIdxBiMap_idx_from_key(&m, "HII"));
  // Equivalent Python/idiomatic C++:  `33 == m["O2II"]`
  EXPECT_EQ(33, grimpl::FrozenKeyIdxBiMap_idx_from_key(&m, "O2II"));

  // Behavior is well-defined when a key isn't known
  int invalid = grimpl::bimap_detail::INVALID_VAL;
  EXPECT_EQ(invalid, grimpl::FrozenKeyIdxBiMap_idx_from_key(&m, "NotAKey"));

  // PART 3: let's show the reverse of the previous lookups
  EXPECT_STREQ("HII", grimpl::FrozenKeyIdxBiMap_key_from_idx(&m, 2));
  EXPECT_STREQ("O2II", grimpl::FrozenKeyIdxBiMap_key_from_idx(&m, 33));

  // Behavior is again well-defined when passing an invalid index
  EXPECT_EQ(nullptr, grimpl::FrozenKeyIdxBiMap_key_from_idx(&m, 131));

  // PART 4: We can also query the length
  EXPECT_EQ(34, grimpl::FrozenKeyIdxBiMap_size(&m));

  // Finally, to cleanup we will deallocate data tracked internally by `m`
  grimpl::drop_FrozenKeyIdxBiMap(&m);
}

// validate basic operations for an empty bimap
TEST(FrozenKeyIdxBiMap, EmptyBasicOps) {
  grackle::impl::FrozenKeyIdxBiMap m = grackle::impl::new_FrozenKeyIdxBiMap(
      nullptr, 0, grackle::impl::BiMapMode::COPIES_KEYDATA);
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&m))
      << "construction of a FrozenKeyIdxBiMap unexpectedly failed";

  EXPECT_EQ(0, grackle::impl::FrozenKeyIdxBiMap_size(&m))
      << "an empty mapping should have a size of 0";

  EXPECT_EQ(grackle::impl::bimap_detail::INVALID_VAL,
            grackle::impl::FrozenKeyIdxBiMap_idx_from_key(&m, "NotAKey"))
      << "key lookup should always fail for an empty mapping";

  EXPECT_EQ(nullptr, grackle::impl::FrozenKeyIdxBiMap_key_from_idx(&m, 0))
      << "index lookup should always fail for an empty mapping";

  grackle::impl::drop_FrozenKeyIdxBiMap(&m);
}

// validate behavior of clone for an empty bimap
TEST(FrozenKeyIdxBiMap, EmptyClone) {
  // make the original
  grackle::impl::FrozenKeyIdxBiMap m = grackle::impl::new_FrozenKeyIdxBiMap(
      nullptr, 0, grackle::impl::BiMapMode::COPIES_KEYDATA);
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&m))
      << "construction of a FrozenKeyIdxBiMap unexpectedly failed";

  // make the clone
  grackle::impl::FrozenKeyIdxBiMap m_clone =
      grackle::impl::FrozenKeyIdxBiMap_clone(&m);

  bool success = grackle::impl::FrozenKeyIdxBiMap_is_ok(&m_clone);

  grackle::impl::drop_FrozenKeyIdxBiMap(&m);  // drop the original

  if (success) {
    grackle::impl::drop_FrozenKeyIdxBiMap(&m_clone);
  } else {
    FAIL() << "cloning an empty mapping failed!";
  }
}

class FrozenKeyIdxBiMapConstructorSuite
    : public testing::TestWithParam<grackle::impl::BiMapMode> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};

TEST_P(FrozenKeyIdxBiMapConstructorSuite, Simple) {
  const char* keys[] = {"denisty", "internal_energy"};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());

  EXPECT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
  grackle::impl::drop_FrozenKeyIdxBiMap(&tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, LongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap_detail::KEYLEN_MAX, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());

  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
  grackle::impl::drop_FrozenKeyIdxBiMap(&tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, TooLongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap_detail::KEYLEN_MAX + 1, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, 0LenKey) {
  const char* keys[2] = {"density", ""};
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, NullptrWithPosCount) {
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(nullptr, 1, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, NotNull0KeyCount) {
  const char* keys[] = {"denisty", "internal_energy"};
  grackle::impl::FrozenKeyIdxBiMap tmp =
      grackle::impl::new_FrozenKeyIdxBiMap(keys, 0, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

INSTANTIATE_TEST_SUITE_P(
    , /* <- leaving Instantiation name empty */
    FrozenKeyIdxBiMapConstructorSuite,
    testing::Values(grackle::impl::BiMapMode::REFS_KEYDATA,
                    grackle::impl::BiMapMode::COPIES_KEYDATA),
    [](const testing::TestParamInfo<
        FrozenKeyIdxBiMapConstructorSuite::ParamType>& info) {
      if (info.param == grackle::impl::BiMapMode::REFS_KEYDATA) {
        return std::string("BIMAP_REFS_KEYDATA");
      } else {
        return std::string("BIMAP_COPIES_KEYDATA");
      }
    });

/// helper function to initialize a map from a vector
grackle::impl::FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(
    const std::vector<std::string>& vec_, grackle::impl::BiMapMode mode) {
  std::size_t key_count = vec_.size();

  // create a vector of pointers
  std::vector<const char*> key_ptr_l(key_count, nullptr);
  for (std::size_t i = 0; i < key_count; i++) {
    key_ptr_l[i] = vec_[i].c_str();
  }

  return new_FrozenKeyIdxBiMap(key_ptr_l.data(), key_count, mode);
}

class FrozenKeyIdxBiMapGeneralSuite
    : public testing::TestWithParam<grackle::impl::BiMapMode> {
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
    (*bimap_p) = tmp;
  }

  void TearDown() override {
    if (bimap_p != nullptr) {
      grackle::impl::drop_FrozenKeyIdxBiMap(bimap_p);
      delete bimap_p;
    }
  }

  bool ReusesOriginalKeyPtrs(const grackle::impl::FrozenKeyIdxBiMap* p) const {
    for (int i = 0; i < 3; i++) {
      const char* orig_key_ptr = ordered_keys[i].c_str();
      if (grackle::impl::FrozenKeyIdxBiMap_key_from_idx(p, i) != orig_key_ptr) {
        return false;
      }
    }
    return true;
  }
};

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKeyContainedKey) {
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "density"),
            1);
  EXPECT_EQ(
      grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "internal_energy"),
      0);
  EXPECT_EQ(
      grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "metal_density"),
      2);
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKeyAbsentKey) {
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "notAKey"),
            grackle::impl::bimap_detail::INVALID_VAL);
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKeyAbsentIrregularKeys) {
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, ""),
            grackle::impl::bimap_detail::INVALID_VAL);

  std::string key(grackle::impl::bimap_detail::KEYLEN_MAX + 1, 'A');
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, key.data()),
            grackle::impl::bimap_detail::INVALID_VAL);
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdxInvalidIdx) {
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 3), nullptr);
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(
                bimap_p, grackle::impl::bimap_detail::INVALID_VAL),
            nullptr);
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdxValidIdx) {
  EXPECT_EQ(
      std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 2)),
      std::string("metal_density"));
  EXPECT_EQ(
      std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 1)),
      std::string("density"));
  EXPECT_EQ(
      std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 0)),
      std::string("internal_energy"));

  // check whether the bimap is using pointers to the keys used during init
  if (GetParam() == grackle::impl::BiMapMode::REFS_KEYDATA) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(bimap_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(bimap_p));
  }
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, Clone) {
  grackle::impl::FrozenKeyIdxBiMap clone =
      grackle::impl::FrozenKeyIdxBiMap_clone(bimap_p);
  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&clone));
  grackle::impl::FrozenKeyIdxBiMap* clone_p = &clone;

  // for the sake of robustly checking everything, we delete bimap_p
  grackle::impl::drop_FrozenKeyIdxBiMap(bimap_p);
  bimap_p = nullptr;

  EXPECT_EQ(
      grackle::impl::FrozenKeyIdxBiMap_idx_from_key(clone_p, "internal_energy"),
      0);
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_idx_from_key(clone_p, "notAKey"),
            grackle::impl::bimap_detail::INVALID_VAL);

  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(clone_p, 3), nullptr);
  EXPECT_EQ(
      std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(clone_p, 1)),
      std::string("density"));

  // check whether the clone is using pointers to the keys used during init
  if (GetParam() == grackle::impl::BiMapMode::REFS_KEYDATA) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(clone_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(clone_p));
  }

  // finally, cleanup the clone
  grackle::impl::drop_FrozenKeyIdxBiMap(clone_p);
}

INSTANTIATE_TEST_SUITE_P(
    , /* <- leaving Instantiation name empty */
    FrozenKeyIdxBiMapGeneralSuite,
    testing::Values(grackle::impl::BiMapMode::REFS_KEYDATA,
                    grackle::impl::BiMapMode::COPIES_KEYDATA),
    [](const testing::TestParamInfo<FrozenKeyIdxBiMapGeneralSuite::ParamType>&
           info) {
      if (info.param == grackle::impl::BiMapMode::REFS_KEYDATA) {
        return std::string("BIMAP_REFS_KEYDATA");
      } else {
        return std::string("BIMAP_COPIES_KEYDATA");
      }
    });
