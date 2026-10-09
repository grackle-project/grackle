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

#include "gmock/gmock.h"
#include "support/config.hpp"
#include "support/PartMap.hpp"

using CreateRslt = GRIMPL_NS::Expected<GRIMPL_NS::PartMap, GRIMPL_NS::Error>;

// teach GoogleTest how to print GRIMPL_NS::partmap::IdxSearch for more
// informative errors (otherwise it just shows the memory's raw byte values)
namespace GRIMPL_NS::partmap {
void PrintTo(const IdxSearch& search, std::ostream* os) {
  *os << "{index=" << search.index << ", pd=" << search.pd
      << ", start_offset=" << search.start_offset << '}';
}
}  // namespace GRIMPL_NS::partmap

using ::testing::Eq;
using ::testing::Field;
using ::testing::Lt;
using ::testing::Optional;

// this is a simple case
TEST(PartMap, Empty) {
  CreateRslt rslt = GRIMPL_NS::PartMap::create(nullptr, nullptr, 0);
  ASSERT_TRUE(rslt.has_value());
  GRIMPL_NS::PartMap m = std::move(rslt).value();

  EXPECT_EQ(m.n_partitions(), 0);
  EXPECT_EQ(m.n_idx(), 0);

  EXPECT_THAT(
      m.part_bounds(0),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Lt(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Lt(0))));

  EXPECT_EQ(m.search_idx(0), std::nullopt);
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
TEST(PartMap, DocString) {
  const int pds[3] = {PartitionName::A, PartitionName::C, PartitionName::B};
  const int sizes[3] = {4, 2, 3};

  CreateRslt rslt = GRIMPL_NS::PartMap::create(pds, sizes, 3);
  ASSERT_TRUE(rslt.has_value());
  GRIMPL_NS::PartMap m = std::move(rslt).value();

  EXPECT_EQ(m.n_partitions(), 3);
  EXPECT_EQ(m.n_idx(), 9);

  EXPECT_THAT(
      m.part_bounds(PartitionName::A),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(4))));
  EXPECT_THAT(
      m.part_bounds(PartitionName::C),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(4)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(6))));
  EXPECT_THAT(
      m.part_bounds(PartitionName::B),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(6)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(9))));

  using GRIMPL_NS::partmap::IdxSearch;
  EXPECT_THAT(m.search_idx(2),
              Optional(Eq(IdxSearch{
                  .index = 2, .pd = PartitionName::A, .start_offset = 2})))
      << "index 2 should be in the partition with pd = PartitionName::A & "
         "start_offset "
      << "should be 2";

  EXPECT_THAT(m.search_idx(5),
              Optional(Eq(IdxSearch{
                  .index = 5, .pd = PartitionName::C, .start_offset = 1})))
      << "index 5 should be in the partition with pd = PartitionName::C & "
         "start_offset "
      << "should be 1";

  EXPECT_THAT(m.search_idx(6),
              Optional(Eq(IdxSearch{
                  .index = 6, .pd = PartitionName::B, .start_offset = 0})))
      << "index 6 should be in the partition with pd = PartitionName::B & "
         "start_offset "
      << "should be 0";

  // extra sanity check!
  EXPECT_EQ(m.search_idx(9999), std::nullopt);
}

// test a failure mode of the factory method:
TEST(PartMap, Nullptrs) {
  const int pds[3] = {PartitionName::A, PartitionName::C, PartitionName::B};
  const int sizes[3] = {4, 2, 3};
  {
    CreateRslt rslt = GRIMPL_NS::PartMap::create(nullptr, sizes, 3);
    EXPECT_FALSE(rslt.has_value()) << "pds argument was false";
  }
  {
    CreateRslt rslt = GRIMPL_NS::PartMap::create(pds, nullptr, 3);
    EXPECT_FALSE(rslt.has_value()) << "sizes argument was false";
  }
}