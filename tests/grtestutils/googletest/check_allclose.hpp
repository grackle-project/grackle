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

#include <vector>

#include <gtest/gtest.h>

#include "../view.hpp"
#include "./check_allclose_detail.hpp"

#define COMPARE_(cmp_fn, actual, desired, selection_mask, idx_mapping)         \
  ::grtest::arraycmp_detail::compare_(::grtest::arraycmp_detail::CmpPack{      \
      {cmp_fn},                                                                \
      {::grtest::arraycmp_detail::PtrPair{actual, desired}},                   \
      selection_mask,                                                          \
      {idx_mapping}})

/// Returns whether 2 pointers are exactly equal
///
/// This draws a lot of inspiration from numpy.testing.assert_array_equal
template <grtest::DataLayout Layout>
testing::AssertionResult check_array_equal(
    const float* actual, const float* desired,
    grtest::IdxMapping<Layout> idx_mapping,
    const bool* selection_mask = nullptr) {
  grtest::arraycmp_detail::FltIsEqual cmp_fn;
  return COMPARE_(cmp_fn, actual, desired, selection_mask, idx_mapping);
}

/// Returns whether 2 pointers are exactly equal
///
/// This draws a lot of inspiration from numpy.testing.assert_array_equal
template <grtest::DataLayout Layout>
testing::AssertionResult check_array_equal(
    const double* actual, const double* desired,
    grtest::IdxMapping<Layout> idx_mapping,
    const bool* selection_mask = nullptr) {
  grtest::arraycmp_detail::FltIsEqual cmp_fn;
  return COMPARE_(cmp_fn, actual, desired, selection_mask, idx_mapping);
}

/// compares 2 pointers
///
/// This draws a lot of inspiration from numpy.testing.assert_allclose
template <grtest::DataLayout Layout>
testing::AssertionResult check_allclose(const float* actual,
                                        const float* desired,
                                        grtest::IdxMapping<Layout> idx_mapping,
                                        double rtol, double atol = 0.0,
                                        const bool* selection_mask = nullptr) {
  grtest::arraycmp_detail::FltIsClose cmp_fn(rtol, atol);
  return COMPARE_(cmp_fn, actual, desired, selection_mask, idx_mapping);
}

/// compares 2 pointers
///
/// This draws a lot of inspiration from numpy.testing.assert_allclose
template <grtest::DataLayout Layout>
testing::AssertionResult check_allclose(const double* actual,
                                        const double* desired,
                                        grtest::IdxMapping<Layout> idx_mapping,
                                        double rtol, double atol = 0.0,
                                        const bool* selection_mask = nullptr) {
  grtest::arraycmp_detail::FltIsClose cmp_fn(rtol, atol);
  return COMPARE_(cmp_fn, actual, desired, selection_mask, idx_mapping);
}

#undef COMPARE_

/// this compares 2 std::vectors
///
/// This draws a lot of inspiration from numpy.testing.assert_allclose
///
/// Parts of this are fairly inefficient, partially because it is adapted from
/// code written from before we adopted googletest
inline testing::AssertionResult check_allclose(
    const std::vector<double>& actual, const std::vector<double>& desired,
    double rtol = 0.0, double atol = 0.0) {
  if (actual.size() != desired.size()) {
    return testing::AssertionFailure()
           << "the compared arrays have different lengths";
  }
  grtest::IdxMapping<grtest::DataLayout::LEFT> idx_mapping(actual.size());
  return check_allclose(actual.data(), desired.data(), idx_mapping, rtol, atol);
}

#endif  // GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_HPP
