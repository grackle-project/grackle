//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// test the PartMap data type
///
//===----------------------------------------------------------------------===//
#include <iostream>  // needed to teach googletest how to print
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "support/config.hpp"
#include "support/PartMap.hpp"

// teach GoogleTest how to print GRIMPL_NS::partmap::IdxSearch for more
// informative errors (otherwise it just shows the memory's raw byte values)
namespace GRIMPL_NS::partmap {
void PrintTo(const IdxSearch& search, std::ostream* os) {
  bool valid = search.has_val;
  std::string index = (valid) ? std::to_string(search.index) : "<garbage>";
  std::string pd = (valid) ? std::to_string(search.pd) : "<garbage>";
  std::string start_offset =
      (valid) ? std::to_string(search.start_offset) : "<garbage>";
  *os << "{has_val=" << valid << ", index=" << index << ", pd=" << pd
      << ", start_offset=" << start_offset << '}';
}
}  // namespace GRIMPL_NS::partmap

using ::testing::Eq;
using ::testing::Field;
using ::testing::Lt;

// this is a simple case
TEST(PartSeq, Empty) {
  GRIMPL_NS::PartMap m = GRIMPL_NS::new_PartMap(nullptr, nullptr, 0);
  ASSERT_TRUE(grackle::impl::PartMap_is_ok(&m));

  EXPECT_EQ(GRIMPL_NS::PartMap_n_partitions(&m), 0);
  EXPECT_EQ(GRIMPL_NS::PartMap_n_idx(&m), 0);

  EXPECT_THAT(
      GRIMPL_NS::PartMap_part_bounds(&m, 0),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Lt(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Lt(0))));

  using GRIMPL_NS::partmap::IdxSearch;
  EXPECT_THAT(GRIMPL_NS::PartMap_search_idx(&m, 0),
              Field("has_val", &IdxSearch::has_val, Eq(false)));
}

// these act as the names of the partition descriptors that are used in the
// following test-case
//
// Ideally, these would be scoped-enums, but that makes use of PartMap very
// clunky! (The only effective way to use a scoped-enum is to make PartMap
// a class template, where partition_descr
namespace PartitionName {
enum { A, B, C };
}  // namespace PartitionName

// this is the case illustrated in PartMap's docstring
TEST(PartSeq, DocString) {
  const int pds[3] = {PartitionName::A, PartitionName::C, PartitionName::B};
  const int sizes[3] = {4, 2, 3};

  GRIMPL_NS::PartMap m = GRIMPL_NS::new_PartMap(pds, sizes, 3);
  ASSERT_TRUE(grackle::impl::PartMap_is_ok(&m));

  EXPECT_EQ(GRIMPL_NS::PartMap_n_partitions(&m), 3);
  EXPECT_EQ(GRIMPL_NS::PartMap_n_idx(&m), 9);

  EXPECT_THAT(
      GRIMPL_NS::PartMap_part_bounds(&m, PartitionName::A),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(4))));
  EXPECT_THAT(
      GRIMPL_NS::PartMap_part_bounds(&m, PartitionName::C),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(4)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(6))));
  EXPECT_THAT(
      GRIMPL_NS::PartMap_part_bounds(&m, PartitionName::B),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(6)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(9))));

  using GRIMPL_NS::partmap::IdxSearch;
  EXPECT_THAT(
      GRIMPL_NS::PartMap_search_idx(&m, 2),
      ::testing::AllOf(Field("has_val", &IdxSearch::has_val, Eq(true)),
                       Field("index", &IdxSearch::index, Eq(2)),
                       Field("pd", &IdxSearch::pd, Eq(PartitionName::A)),
                       Field("start_offset", &IdxSearch::start_offset, Eq(2))))
      << "index 2 should be in the partition with pd = PartitionName::A & "
         "start_offset "
      << "should be 2";

  EXPECT_THAT(
      GRIMPL_NS::PartMap_search_idx(&m, 5),
      ::testing::AllOf(Field("has_val", &IdxSearch::has_val, Eq(true)),
                       Field("index", &IdxSearch::index, Eq(5)),
                       Field("pd", &IdxSearch::pd, Eq(PartitionName::C)),
                       Field("start_offset", &IdxSearch::start_offset, Eq(1))))
      << "index 5 should be in the partition with pd = PartitionName::C & "
         "start_offset "
      << "should be 1";

  EXPECT_THAT(
      GRIMPL_NS::PartMap_search_idx(&m, 6),
      ::testing::AllOf(Field("has_val", &IdxSearch::has_val, Eq(true)),
                       Field("index", &IdxSearch::index, Eq(6)),
                       Field("pd", &IdxSearch::pd, Eq(PartitionName::B)),
                       Field("start_offset", &IdxSearch::start_offset, Eq(0))))
      << "index 6 should be in the partition with pd = PartitionName::B & "
         "start_offset "
      << "should be 0";

  // extra sanity check!
  EXPECT_THAT(GRIMPL_NS::PartMap_search_idx(&m, 9999),
              Field("is_valid", &IdxSearch::has_val, Eq(false)));
}
