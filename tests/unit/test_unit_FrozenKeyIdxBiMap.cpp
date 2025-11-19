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
#include "utils/FrozenKeyIdxBiMap.hpp"
#include "grackle.h"

class FrozenKeyIdxBiMapConstructorSuite :
  public testing::TestWithParam<grackle::impl::BiMapMode> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.

};

TEST_P(FrozenKeyIdxBiMapConstructorSuite, Simple) {
  const char* keys[] = {"denisty", "internal_energy"};

  grackle::impl::FrozenKeyIdxBiMap tmp
    = grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());

  EXPECT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
  grackle::impl::drop_FrozenKeyIdxBiMap(&tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, LongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap::keylen_max, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp
    = grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());

  ASSERT_TRUE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
  grackle::impl::drop_FrozenKeyIdxBiMap(&tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, TooLongKey) {
  const char* first_key = "density";
  std::string long_key(grackle::impl::bimap::keylen_max+1, 'A');
  const char* keys[2] = {first_key, long_key.data()};

  grackle::impl::FrozenKeyIdxBiMap tmp
    = grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, 0LenKey) {
  const char* keys[2] = {"density",""};
  grackle::impl::FrozenKeyIdxBiMap tmp
    = grackle::impl::new_FrozenKeyIdxBiMap(keys, 2, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, NoKeys) {

  grackle::impl::FrozenKeyIdxBiMap tmp
    = grackle::impl::new_FrozenKeyIdxBiMap(nullptr, 0, GetParam());
  ASSERT_FALSE(grackle::impl::FrozenKeyIdxBiMap_is_ok(&tmp));
}

INSTANTIATE_TEST_SUITE_P(
  , /* <- leaving Instantiation name empty */
  FrozenKeyIdxBiMapConstructorSuite,
  testing::Values(grackle::impl::BiMapMode::REFS_KEYDATA, grackle::impl::BiMapMode::COPIES_KEYDATA),
  [](const testing::TestParamInfo<FrozenKeyIdxBiMapConstructorSuite::ParamType>& info){
    if (info.param == grackle::impl::BiMapMode::REFS_KEYDATA) {
      return std::string("BIMAP_REFS_KEYDATA");
    } else {
      return std::string("BIMAP_COPIES_KEYDATA");
    }
  }
);


/// helper function to initialize a map from a vector
grackle::impl::FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(
  const std::vector<std::string>& vec_, grackle::impl::BiMapMode mode)
{
  std::size_t key_count = vec_.size();

  // create a vector of pointers
  std::vector<const char*> key_ptr_l(key_count, nullptr);
  for (std::size_t i = 0; i < key_count; i++) {
    key_ptr_l[i] = vec_[i].c_str();
  }

  return new_FrozenKeyIdxBiMap(key_ptr_l.data(), key_count, mode);
}

class FrozenKeyIdxBiMapGeneralSuite :
  public testing::TestWithParam<grackle::impl::BiMapMode>
{
protected:
  std::vector<std::string> ordered_keys;
  grackle::impl::FrozenKeyIdxBiMap* bimap_p = nullptr;

  // we use SetUp/Teardown instead of constructor and destructor so we can
  // perform some sanity checks with ASSERTIONs
  void SetUp() override {
    ordered_keys = std::vector<std::string>{
      "internal_energy", "density", "metal_density"
    };

    grackle::impl::BiMapMode mode = this->GetParam();
    grackle::impl::FrozenKeyIdxBiMap tmp
      = new_FrozenKeyIdxBiMap(ordered_keys, mode);
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
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "density"),
    1
  );
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "internal_energy"),
    0
  );
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "metal_density"),
    2
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKeyAbsentKey) {
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, "notAKey"),
    grackle::impl::bimap::invalid_val
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKeyAbsentIrregularKeys) {
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, ""),
    grackle::impl::bimap::invalid_val
  );

  std::string key(grackle::impl::bimap::keylen_max+1, 'A');
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(bimap_p, key.data()),
    grackle::impl::bimap::invalid_val
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdxInvalidIdx) {
  EXPECT_EQ(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 3), nullptr);
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_key_from_idx(
      bimap_p, grackle::impl::bimap::invalid_val),
    nullptr
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdxValidIdx) {
  EXPECT_EQ(
    std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 2)),
    std::string("metal_density")
  );
  EXPECT_EQ(
    std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 1)),
    std::string("density")
  );
  EXPECT_EQ(
    std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(bimap_p, 0)),
    std::string("internal_energy")
  );

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
    0
  );
  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_idx_from_key(clone_p, "notAKey"),
    grackle::impl::bimap::invalid_val
  );

  EXPECT_EQ(
    grackle::impl::FrozenKeyIdxBiMap_key_from_idx(clone_p, 3),
    nullptr
  );
  EXPECT_EQ(
    std::string(grackle::impl::FrozenKeyIdxBiMap_key_from_idx(clone_p, 1)),
    std::string("density")
  );

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
  [](const testing::TestParamInfo<FrozenKeyIdxBiMapGeneralSuite::ParamType>& info){
    if (info.param == grackle::impl::BiMapMode::REFS_KEYDATA) {
      return std::string("BIMAP_REFS_KEYDATA");
    } else {
      return std::string("BIMAP_COPIES_KEYDATA");
    }
  }
);
