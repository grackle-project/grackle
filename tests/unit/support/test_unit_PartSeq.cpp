//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// test the PartSeq data type
///
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "support/config.hpp"
#include "support/PartSeq.hpp"

using ::testing::Eq;
using ::testing::Field;
using ::testing::Lt;

// this is a simple case
TEST(PartSeq, Empty) {
  GRIMPL_NS::PartSeq s = GRIMPL_NS::PartSeq_from_part_sizes(nullptr, 0);

  EXPECT_EQ(GRIMPL_NS::PartSeq_n_partitions(&s), 0);
  EXPECT_EQ(GRIMPL_NS::PartSeq_n_idx(&s), 0);

  EXPECT_THAT(
      GRIMPL_NS::PartSeq_part_prop(&s, 0),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Lt(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Lt(0))));

  using GRIMPL_NS::RelativePartitionInfo;
  EXPECT_THAT(
      GRIMPL_NS::PartSeq_find_part_info(&s, 0),
      ::testing::AllOf(
          Field("part_id", &RelativePartitionInfo::part_id, Lt(0)),
          Field("start_offset", &RelativePartitionInfo::start_offset, Lt(0))));
}

// this is the case illustrated in PartSeq's docstring
TEST(PartSeq, DocString) {
  const int part_sizes[3] = {4, 2, 3};

  GRIMPL_NS::PartSeq s = GRIMPL_NS::PartSeq_from_part_sizes(part_sizes, 3);

  EXPECT_EQ(GRIMPL_NS::PartSeq_n_partitions(&s), 3);
  EXPECT_EQ(GRIMPL_NS::PartSeq_n_idx(&s), 9);

  EXPECT_THAT(
      GRIMPL_NS::PartSeq_part_prop(&s, 0),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(0)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(4))));
  EXPECT_THAT(
      GRIMPL_NS::PartSeq_part_prop(&s, 1),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(4)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(6))));
  EXPECT_THAT(
      GRIMPL_NS::PartSeq_part_prop(&s, 2),
      ::testing::AllOf(Field("start", &GRIMPL_NS::IdxInterval::start, Eq(6)),
                       Field("stop", &GRIMPL_NS::IdxInterval::stop, Eq(9))));

  using GRIMPL_NS::RelativePartitionInfo;
  EXPECT_THAT(
      GRIMPL_NS::PartSeq_find_part_info(&s, 2),
      ::testing::AllOf(
          Field("part_id", &RelativePartitionInfo::part_id, Eq(0)),
          Field("start_offset", &RelativePartitionInfo::start_offset, Eq(2))))
      << "index 2 should be in the partition with part_id = 0 & start_offset "
      << "should be 2";

  EXPECT_THAT(
      GRIMPL_NS::PartSeq_find_part_info(&s, 6),
      ::testing::AllOf(
          Field("part_id", &RelativePartitionInfo::part_id, Eq(2)),
          Field("start_offset", &RelativePartitionInfo::start_offset, Eq(0))))
      << "index 2 should be in the partition with part_id = 0 & start_offset "
      << "should be 2";
}
