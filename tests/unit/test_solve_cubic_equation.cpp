//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Tests grackle::impl::solve_cubic_equation_cpp
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>

#include <array>

#include "solve_cubic_equation.hpp"

namespace {

struct CubicSuccessCase {
  double a;
  double b;
  double c;
  double expected_root;
};

} // anonymous namespace

TEST(SolveCubicEquation, OneRealRootCases) {
  const std::array<CubicSuccessCase, 4> cases{{
    // x^3 + x + 1 = 0
    {0.0, 1.0, 1.0, -0.6823278038280193},
    // (x - 2)(x^2 + 1) = x^3 - 2x^2 + x - 2
    {-2.0, 1.0, -2.0, 2.0},
    // (x + 1)(x^2 - 2x + 2) = x^3 - x^2 + 2
    {-1.0, 0.0, 2.0, -1.0},
    // (x + 2)(x^2 + x + 1) = x^3 + 3x^2 + 3x + 2
    {3.0, 3.0, 2.0, -2.0},
  }};

  for (const auto& one_case : cases) {
    double root = 0.0;
    int rslt = grackle::impl::solve_cubic_equation_cpp(
      one_case.a, one_case.b, one_case.c, root
    );

    EXPECT_EQ(rslt, 0);
    EXPECT_NEAR(root, one_case.expected_root, 1e-12);
  }
}

TEST(SolveCubicEquation, ThreeRealRootsReturnsError) {
  // (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
  double root = 0.0;
  int rslt = grackle::impl::solve_cubic_equation_cpp(-6.0, 11.0, -6.0, root);

  EXPECT_EQ(rslt, 1);
}
