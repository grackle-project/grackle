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
#include "gtest/gtest.h"
#include "support/config.hpp"
#include "support/PartMap.hpp"

using GRIMPL_NS::PartMap;
template <typename PartitionDescrT>
using CreateRslt =
    GRIMPL_NS::Expected<PartMap<PartitionDescrT>, GRIMPL_NS::Error>;

// define a set of enums to act as sample partition descriptors in our tests
enum class PartitionName { A, B, C };
// teach GoogleTest how to print our sample partition descriptors for more
// informative errors
void PrintTo(const PartitionName& pd, std::ostream* os) {
  switch (pd) {
    case PartitionName::A:
      *os << "ParititionName::A";
      return;
    case PartitionName::B:
      *os << "ParititionName::B";
      return;
    case PartitionName::C:
      *os << "ParititionName::C";
      return;
    default:
      *os << "ParititionName::<UNKNOWN>";
      return;
  }
}

// teach GoogleTest how to print GRIMPL_NS::partmap::IdxSearch for more
// informative errors (otherwise it just shows the memory's raw byte values)
namespace GRIMPL_NS::partmap {

template <typename PartitionDescrT>
void PrintTo(const IdxSearch<PartitionDescrT>& search, std::ostream* os) {
  *os << "{index=" << search.index
      << ", pd=" << testing::PrintToString(search.pd)
      << ", start_offset=" << search.start_offset << '}';
}
}  // namespace GRIMPL_NS::partmap

using ::testing::Eq;
using ::testing::Field;
using ::testing::Lt;
using ::testing::Optional;

// this is a simple case
TEST(PartMap, Empty) {
  CreateRslt rslt = PartMap<PartitionName>::create(nullptr, nullptr, 0);
  ASSERT_TRUE(rslt.has_value());
  GRIMPL_NS::PartMap m = std::move(rslt).value();

  EXPECT_EQ(m.n_partitions(), 0);
  EXPECT_EQ(m.n_idx(), 0);

  EXPECT_THAT(m.part_bounds(PartitionName::A),
              ::testing::AllOf(
                  Field("start", &GRIMPL_NS::IndexInterval1D::start, Lt(0)),
                  Field("stop", &GRIMPL_NS::IndexInterval1D::stop, Lt(0))));

  EXPECT_EQ(m.search_idx(0), std::nullopt);
}

// this is the case illustrated in PartMap's docstring
TEST(PartMap, DocString) {
  const PartitionName pds[3] = {PartitionName::A, PartitionName::C,
                                PartitionName::B};
  const int sizes[3] = {4, 2, 3};

  CreateRslt<PartitionName> rslt =
      GRIMPL_NS::PartMap<PartitionName>::create(pds, sizes, 3);
  ASSERT_TRUE(rslt.has_value());
  GRIMPL_NS::PartMap m = std::move(rslt).value();

  EXPECT_EQ(m.n_partitions(), 3);
  EXPECT_EQ(m.n_idx(), 9);

  EXPECT_THAT(m.part_bounds(PartitionName::A),
              ::testing::AllOf(
                  Field("start", &GRIMPL_NS::IndexInterval1D::start, Eq(0)),
                  Field("stop", &GRIMPL_NS::IndexInterval1D::stop, Eq(4))));
  EXPECT_THAT(m.part_bounds(PartitionName::C),
              ::testing::AllOf(
                  Field("start", &GRIMPL_NS::IndexInterval1D::start, Eq(4)),
                  Field("stop", &GRIMPL_NS::IndexInterval1D::stop, Eq(6))));
  EXPECT_THAT(m.part_bounds(PartitionName::B),
              ::testing::AllOf(
                  Field("start", &GRIMPL_NS::IndexInterval1D::start, Eq(6)),
                  Field("stop", &GRIMPL_NS::IndexInterval1D::stop, Eq(9))));

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
  const PartitionName pds[3] = {PartitionName::A, PartitionName::C,
                                PartitionName::B};
  const int sizes[3] = {4, 2, 3};
  {
    CreateRslt<PartitionName> rslt =
        GRIMPL_NS::PartMap<PartitionName>::create(nullptr, sizes, 3);
    EXPECT_FALSE(rslt.has_value()) << "pds argument was false";
  }
  {
    CreateRslt<PartitionName> rslt =
        GRIMPL_NS::PartMap<PartitionName>::create(pds, nullptr, 3);
    EXPECT_FALSE(rslt.has_value()) << "sizes argument was false";
  }
}