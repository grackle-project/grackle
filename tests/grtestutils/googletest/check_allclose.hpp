//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares assertions for comparing arrays
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_HPP
#define GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_HPP

#include <string>
#include <vector>
#include <gtest/gtest.h>

/// this compares 2 std::vectors
///
/// This draws a lot of inspiration from numpy.testing.assert_allclose
///
/// Parts of this are fairly inefficient, partially because it is adapted from
/// code written from before we adopted googletest
testing::AssertionResult check_allclose(const std::vector<double>& actual,
                                        const std::vector<double>& desired,
                                        double rtol = 0.0, double atol = 0.0,
                                        std::string err_msg = "");

#endif  // GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_HPP
