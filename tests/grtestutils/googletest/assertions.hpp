//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares some general-purpose assertions for testing grackle
///
//===----------------------------------------------------------------------===//
#ifndef GRTEST_GOOGLETEST_ASSERTIONS_HPP
#define GRTEST_GOOGLETEST_ASSERTIONS_HPP

#include <optional>
#include <string>

#include <gtest/gtest.h>
#include "grackle.h"

namespace grtest {

inline testing::AssertionResult IsGRSUCCESS(const char* expr1, int val1) {
  if (val1 == GR_SUCCESS) {
    return testing::AssertionSuccess();
  }

  std::string descr;
  switch (val1) {
    case GR_SUCCESS: {
      descr = "GR_SUCCESS";
      break;
    }
    case GR_FAIL: {
      descr = "GR_FAIL";
      break;
    }
    default: {
      descr = "not a standard code";
    }
  }
  testing::AssertionResult out = testing::AssertionFailure();
  out << "Evaluated: " << expr1 << '\n'
      << " Expected: GR_SUCCESS (aka " << GR_SUCCESS << ")\n"
      << "   Actual: " << val1 << " (" << descr << ')';
  return out;
}

inline testing::AssertionResult IsGRError(const char* expr1, int val1) {
  if (val1 != GR_SUCCESS) {
    return testing::AssertionSuccess();
  }

  testing::AssertionResult out = testing::AssertionFailure();
  out << "Evaluated: " << expr1 << '\n'
      << " Expected: a value other than GR_SUCCESS\n"
      << "   Actual: " << val1 << " (GR_SUCCESS)";
  return out;
}

}  // namespace grtest

#define EXPECT_GR_SUCCESS(expr) EXPECT_PRED_FORMAT1(::grtest::IsGRSUCCESS, expr)

#define ASSERT_GR_SUCCESS(expr) ASSERT_PRED_FORMAT1(::grtest::IsGRSUCCESS, expr)

#define EXPECT_GR_ERR(expr) EXPECT_PRED_FORMAT1(::grtest::IsGRError, expr)

#define ASSERT_GR_ERR(expr) ASSERT_PRED_FORMAT1(::grtest::IsGRError, expr)

/// applies to a std::optional<T>
#define ASSERT_NOT_EMPTY(expr) ASSERT_NE(expr, std::nullopt)

#endif  // GRTEST_GOOGLETEST_ASSERTIONS_HPP
