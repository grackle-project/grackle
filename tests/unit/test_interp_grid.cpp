//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of the InterpGrid and InterpGridProps logic
///
//===----------------------------------------------------------------------===//

#include <iostream>  // needed to teach googletest how to print
#include <cstdio>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "gtest/gtest.h"
#include "interp_grid.hpp"
#include "support/config.hpp"
#include <gmock/gmock.h>

// teach GoogleTest how to print GRIMPL_NS::InterpDimScale &
// GRIMPL_NS::InterpGridProps for more informative errors (otherwise it just
// shows the memory's raw byte values)
namespace GRIMPL_NAMESPACE_DECL {
void PrintTo(const InterpDimScale& ids, std::ostream* os) {
  *os << "InterpDimScale::Linear(" << ids.count << ", " << ids.start << ", "
      << ids.step << ')';
}

}  // namespace GRIMPL_NAMESPACE_DECL

static std::string matcher_descr_(
    const std::vector<GRIMPL_NS::InterpDimScale>& vec, bool negation) {
  const char* partial = negation ? "doesn't hold" : "holds";
  return partial + testing::PrintToString(vec);
}

// the following defines a custom matcher for checking whether an interpolation
// grid is consistent with a vector of InterpDimScale
//   EXPECT_THAT(arg, HoldsDimScales(vec_of_dim_scales))
MATCHER_P(HoldsDimScales, v, matcher_descr_(v, negation)) {
  int expected_rank = v.size();
  if (!arg) {
    *result_listener << " is invalid";
    return false;
  } else if (arg.rank != expected_rank) {
    *result_listener << " holds " << arg.rank << " dims, rather than "
                     << expected_rank;
    return false;
  } else {
    for (int dim = 0; dim < expected_rank; dim++) {
      printf("comparing dim = %d\n", dim);
      fflush(stdout);
      const GRIMPL_NS::InterpDimScale& dim_scale = v[dim];
      int extent = arg.dimension[dim];
      if (extent != dim_scale.count) {
        *result_listener << " has length " << extent << " along dim " << dim
                         << " (rather than " << dim_scale.count << ')';
        return false;
      }
      double step = arg.parameter_spacing[dim];
      if (step != dim_scale.step) {
        *result_listener << " has spacing " << step << " along dim " << dim
                         << " (rather than " << dim_scale.step << ')';
        return false;
      }
      for (int i = 0; i < extent; i++) {
        printf(".. comparing i = %d\n", i);
        fflush(stdout);
        double expected = dim_scale.access_value(i);
        double actual = arg.parameters[dim][i];
        if (actual != expected) {
          *result_listener << " has a mismatched value at index " << i
                           << ", along dimension " << dim << "(found " << actual
                           << ", not " << expected << ')';
          return false;
        }
      }
    }
    return true;
  }
}

TEST(InterpGridProps, Empty) {
  GRIMPL_NS::InterpGridProps grid_props;
  EXPECT_FALSE(grid_props);
}

TEST(InterpGridProps, Simple1D) {
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(4, 0.0, 1.0)};
  GRIMPL_NS::InterpGridProps grid_props(dim_scales.size(), dim_scales.data());
  EXPECT_TRUE(grid_props);

  EXPECT_EQ(grid_props.data_size, 4);
  EXPECT_THAT(grid_props, HoldsDimScales(dim_scales));

  // the remaining logic is just extra sanity checks
  // (in a sense, we're validating that the HoldsDimScales matcher works
  ASSERT_EQ(grid_props.rank, 1);
  EXPECT_EQ(grid_props.parameter_spacing[0], 1.0);
  EXPECT_EQ(grid_props.dimension[0], 4);
  EXPECT_EQ(grid_props.parameters[0][0], 0.0);
  EXPECT_EQ(grid_props.parameters[0][1], 1.0);
  EXPECT_EQ(grid_props.parameters[0][2], 2.0);
  EXPECT_EQ(grid_props.parameters[0][3], 3.0);
}

TEST(InterpGridProps, Simple2D) {
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(4, 0.0, 1.0),
      GRIMPL_NS::InterpDimScale::Linear(2, 5.0, -10.0),
  };
  GRIMPL_NS::InterpGridProps grid_props(dim_scales.size(), dim_scales.data());
  EXPECT_TRUE(grid_props);

  EXPECT_EQ(grid_props.data_size, 8);
  EXPECT_THAT(grid_props, HoldsDimScales(dim_scales));
}

TEST(InterpGridProps, Simple3D) {
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(4, 0.0, 1.0),
      GRIMPL_NS::InterpDimScale::Linear(2, 5.0, -10.0),
      GRIMPL_NS::InterpDimScale::Linear(5, 30.0, 100.0),
  };
  GRIMPL_NS::InterpGridProps grid_props(dim_scales.size(), dim_scales.data());
  EXPECT_TRUE(grid_props);

  EXPECT_EQ(grid_props.data_size, 40);
  EXPECT_THAT(grid_props, HoldsDimScales(dim_scales));
}

TEST(InterpGridProps, MoveConstruct) {
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(4, 0.0, 1.0),
      GRIMPL_NS::InterpDimScale::Linear(2, 5.0, -10.0),
      GRIMPL_NS::InterpDimScale::Linear(5, 30.0, 100.0),
  };
  GRIMPL_NS::InterpGridProps original(dim_scales.size(), dim_scales.data());
  EXPECT_TRUE(original);

  // grid_props is move-constructed
  GRIMPL_NS::InterpGridProps grid_props(std::move(original));
  EXPECT_TRUE(grid_props);
  EXPECT_FALSE(original) << "the original instance should now be invalid";

  EXPECT_EQ(grid_props.data_size, 40);
  EXPECT_THAT(grid_props, HoldsDimScales(dim_scales));
}

TEST(InterpGridProps, MoveAssign) {
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(4, 0.0, 1.0),
      GRIMPL_NS::InterpDimScale::Linear(2, 5.0, -10.0),
      GRIMPL_NS::InterpDimScale::Linear(5, 30.0, 100.0),
  };
  GRIMPL_NS::InterpGridProps original(dim_scales.size(), dim_scales.data());
  EXPECT_TRUE(original);

  GRIMPL_NS::InterpGridProps grid_props;  // <- default constructed
  grid_props = std::move(original);       // <- move assignment
  EXPECT_TRUE(grid_props);
  EXPECT_FALSE(original) << "the original instance should now be invalid";

  EXPECT_EQ(grid_props.data_size, 40);
  EXPECT_THAT(grid_props, HoldsDimScales(dim_scales));
}

// InterpGrid is less widely used than InterpGridProps

TEST(InterpGrid, Empty) {
  GRIMPL_NS::InterpGridProps interp_grid;
  EXPECT_FALSE(interp_grid);
}

TEST(InterpGrid, NonNullButNoProps) {
  double* data = new double[2];
  data[0] = 1;
  data[1] = 1;
  GRIMPL_NS::InterpGridProps empty_grid_props;
  GRIMPL_NS::InterpGrid interp_grid(std::move(empty_grid_props), data);
  EXPECT_FALSE(interp_grid);
}

TEST(InterpGrid, Simple1D) {
  // create grid_props
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales{
      GRIMPL_NS::InterpDimScale::Linear(2, 0.0, 1.0)};
  GRIMPL_NS::InterpGridProps grid_props(1, dim_scales.data());

  // create pointer with values
  double* data = new double[2];
  data[0] = 1.0;
  data[1] = -1.0;

  // construct the interp_grid
  GRIMPL_NS::InterpGrid interp_grid(std::move(grid_props), data);
  ASSERT_TRUE(interp_grid);
  EXPECT_EQ(interp_grid.data[0], 1.0);
  EXPECT_EQ(interp_grid.data[1], -1.0);
  EXPECT_THAT(interp_grid.props, HoldsDimScales(dim_scales));

  EXPECT_FALSE(grid_props)
      << "the grid_props variable should be in a null-state since it was "
      << "\"consumed\" to create interp_grid";

  // the following should probably be moved to a different test-case
  GRIMPL_NS::InterpGrid interp_grid2(std::move(interp_grid));
  EXPECT_TRUE(interp_grid2);
  EXPECT_FALSE(interp_grid);

  GRIMPL_NS::InterpGrid interp_grid3;
  interp_grid3 = std::move(interp_grid2);
  EXPECT_TRUE(interp_grid3);
  EXPECT_FALSE(interp_grid);
}