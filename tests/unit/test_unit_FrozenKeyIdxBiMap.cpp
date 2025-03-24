// to be clear, this file is not directly testing a part of the public API
// -> instead it is testing a well-defined internal building-block

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "FrozenKeyIdxBiMap.h"
#include "grackle.h"

class FrozenKeyIdxBiMapConstructorSuite : 
  public testing::TestWithParam<enum BiMapMode> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.

};

TEST_P(FrozenKeyIdxBiMapConstructorSuite, Simple) {
  FrozenKeyIdxBiMap* tmp = nullptr;
  const char* keys[] = {"denisty", "internal_energy"};

  EXPECT_EQ(new_FrozenKeyIdxBiMap(&tmp, keys, 2, GetParam()), GR_SUCCESS);
  ASSERT_NE(tmp, nullptr);
  drop_FrozenKeyIdxBiMap(tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, LongKey) {
  FrozenKeyIdxBiMap* tmp = nullptr;

  const char* first_key = "density";
  std::string long_key(STRU16MAP_KEYLEN_MAX, 'A');
  const char* keys[2] = {first_key, long_key.data()};
  EXPECT_EQ(new_FrozenKeyIdxBiMap(&tmp, keys, 2, GetParam()), GR_SUCCESS);
  ASSERT_NE(tmp, nullptr);
  drop_FrozenKeyIdxBiMap(tmp);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, TooLongKey) {
  FrozenKeyIdxBiMap* tmp = nullptr;

  const char* first_key = "density";
  std::string long_key(STRU16MAP_KEYLEN_MAX+1, 'A');
  const char* keys[2] = {first_key, long_key.data()};
  EXPECT_NE(new_FrozenKeyIdxBiMap(&tmp, keys, 2, GetParam()), GR_SUCCESS);
  ASSERT_EQ(tmp, nullptr);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, 0LenKey) {
  FrozenKeyIdxBiMap* tmp = nullptr;

  const char* keys[2] = {"density",""};
  EXPECT_NE(new_FrozenKeyIdxBiMap(&tmp, keys, 2, GetParam()), GR_SUCCESS);
  ASSERT_EQ(tmp, nullptr);
}

TEST_P(FrozenKeyIdxBiMapConstructorSuite, NoKeys) {
  FrozenKeyIdxBiMap* tmp = nullptr;

  EXPECT_NE(new_FrozenKeyIdxBiMap(&tmp, nullptr, 0, GetParam()), GR_SUCCESS);
  ASSERT_EQ(tmp, nullptr); // if this fails, we'll leak memory (but not much
                           // can be done)
}

INSTANTIATE_TEST_SUITE_P(
  , /* <- leaving Instantiation name empty */
  FrozenKeyIdxBiMapConstructorSuite,
  testing::Values(BIMAP_REFS_KEYDATA, BIMAP_COPIES_KEYDATA),
  [](const testing::TestParamInfo<FrozenKeyIdxBiMapConstructorSuite::ParamType>& info){
    if (info.param == BIMAP_REFS_KEYDATA) {
      return std::string("BIMAP_REFS_KEYDATA");
    } else {
      return std::string("BIMAP_COPIES_KEYDATA");
    }
  }
);


/// helper function to initialize a map from a vector
int new_FrozenKeyIdxBiMap(FrozenKeyIdxBiMap** out,
                          const std::vector<std::string>& vec_,
                          enum BiMapMode mode)
{
  std::size_t key_count = vec_.size();

  // create a vector of pointers
  std::vector<const char*> key_ptr_l(key_count, nullptr);
  for (std::size_t i = 0; i < key_count; i++) {
    key_ptr_l[i] = vec_[i].c_str();
  }

  return new_FrozenKeyIdxBiMap(out, key_ptr_l.data(), key_count, mode);
}

class FrozenKeyIdxBiMapGeneralSuite :
  public testing::TestWithParam<enum BiMapMode>
{
protected:
  std::vector<std::string> ordered_keys;
  FrozenKeyIdxBiMap* bimap_p = nullptr;

  // I don't think we can use this->GetParam from a constructor... So we use 
  // SetUp/Teardown instead of constructor and destructor instead
  void SetUp() override {
    ordered_keys = std::vector<std::string>{
      "internal_energy", "density", "metal_density"
    };

    enum BiMapMode mode = this->GetParam();

    // since we are in SetUp, (rather than a constructor) we can perform some sanity checks with ASSERTIONs
    ASSERT_EQ(
        new_FrozenKeyIdxBiMap(&this->bimap_p, ordered_keys, mode),
        GR_SUCCESS
    );

  }

  void TearDown() override {
    if (bimap_p != nullptr) { drop_FrozenKeyIdxBiMap(bimap_p); }
  }

  bool ReusesOriginalKeyPtrs(const FrozenKeyIdxBiMap* p) const {
    for (int i = 0; i < 3; i++) {
      const char* orig_key_ptr = ordered_keys[i].c_str();
      if (FrozenKeyIdxBiMap_key_from_idx(p, i) != orig_key_ptr) {return false;}
    }
    return true;
  }

};

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKey_ContainedKey) {
  EXPECT_EQ(FrozenKeyIdxBiMap_idx_from_key(bimap_p, "density"), 1);
  EXPECT_EQ(FrozenKeyIdxBiMap_idx_from_key(bimap_p, "internal_energy"), 0);
  EXPECT_EQ(FrozenKeyIdxBiMap_idx_from_key(bimap_p, "metal_density"), 2);
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKey_AbsentKey) {
  EXPECT_EQ(
    FrozenKeyIdxBiMap_idx_from_key(bimap_p, "notAKey"), STRU16MAP_INVALID_VAL
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, IdxFromKey_AbsentIrregularKeys) {
  EXPECT_EQ(
    FrozenKeyIdxBiMap_idx_from_key(bimap_p, ""), STRU16MAP_INVALID_VAL
  );

  std::string key(STRU16MAP_KEYLEN_MAX+1, 'A');
  EXPECT_EQ(
    FrozenKeyIdxBiMap_idx_from_key(bimap_p, key.data()), STRU16MAP_INVALID_VAL
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdx_InvalidIdx) {
  EXPECT_EQ(FrozenKeyIdxBiMap_key_from_idx(bimap_p, 3), nullptr);
  EXPECT_EQ(
    FrozenKeyIdxBiMap_key_from_idx(bimap_p, STRU16MAP_INVALID_VAL), nullptr
  );
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, KeyFromIdx_ValidIdx) {
  EXPECT_EQ(
    std::string(FrozenKeyIdxBiMap_key_from_idx(bimap_p, 2)),
    std::string("metal_density")
  );
  EXPECT_EQ(
    std::string(FrozenKeyIdxBiMap_key_from_idx(bimap_p, 1)),
    std::string("density")
  );
  EXPECT_EQ(
    std::string(FrozenKeyIdxBiMap_key_from_idx(bimap_p, 0)),
    std::string("internal_energy")
  );

  // check whether the bimap is using pointers to the keys used during init
  if (GetParam() == BIMAP_REFS_KEYDATA) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(bimap_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(bimap_p));
  }
}

TEST_P(FrozenKeyIdxBiMapGeneralSuite, Clone) {
  FrozenKeyIdxBiMap* clone_p = nullptr;
  EXPECT_EQ(FrozenKeyIdxBiMap_clone(&clone_p, bimap_p), GR_SUCCESS);
  ASSERT_NE(clone_p, nullptr);

  // for the sake of robustly checking everything, we delete bimap_p
  drop_FrozenKeyIdxBiMap(bimap_p);
  bimap_p = nullptr;

  EXPECT_EQ(FrozenKeyIdxBiMap_idx_from_key(clone_p, "internal_energy"), 0);
  EXPECT_EQ(
    FrozenKeyIdxBiMap_idx_from_key(clone_p, "notAKey"), STRU16MAP_INVALID_VAL
  );

  EXPECT_EQ(FrozenKeyIdxBiMap_key_from_idx(clone_p, 3), nullptr);
  EXPECT_EQ(
    std::string(FrozenKeyIdxBiMap_key_from_idx(clone_p, 1)),
    std::string("density")
  );

  // check whether the clone is using pointers to the keys used during init
  if (GetParam() == 0) {
    EXPECT_TRUE(ReusesOriginalKeyPtrs(clone_p));
  } else {
    EXPECT_FALSE(ReusesOriginalKeyPtrs(clone_p));
  }

  // finally, cleanup the clone
  drop_FrozenKeyIdxBiMap(clone_p);
}


INSTANTIATE_TEST_SUITE_P(
  , /* <- leaving Instantiation name empty */
  FrozenKeyIdxBiMapGeneralSuite,
  testing::Values(BIMAP_REFS_KEYDATA, BIMAP_COPIES_KEYDATA),
  [](const testing::TestParamInfo<FrozenKeyIdxBiMapGeneralSuite::ParamType>& info){
    if (info.param == BIMAP_REFS_KEYDATA) {
      return std::string("BIMAP_REFS_KEYDATA");
    } else {
      return std::string("BIMAP_COPIES_KEYDATA");
    }
  }
);
