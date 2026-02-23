//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/implements a bunch of helper code to assist with implementing
/// logic in check_allclose.hpp
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_DETAIL_HPP
#define GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_DETAIL_HPP

#include <cmath>
#include <string>
#include <type_traits>
#include <variant>
#include <gtest/gtest.h>

#include "../view.hpp"
#include "gtest/gtest.h"

namespace grtest::arraycmp_detail {

template <bool EqualNan, typename T>
[[gnu::always_inline]] bool isclose_(T actual, T desired, double rtol,
                                     double atol) {
  static_assert(std::is_floating_point_v<T>);

  T abs_diff = std::fabs(actual - desired);
  // the following variable is false if actual or desired (or both) is NaN
  bool isclose = abs_diff <= atol + rtol * std::fabs(desired);

  if constexpr (EqualNan) {
    bool actual_isnan = std::isnan(actual);
    bool desired_isnan = std::isnan(desired);
    return (actual_isnan && desired_isnan) || isclose;
  } else {
    return isclose;
  }
}

/// "functor" to check if floating point values are equal within a tolerance
///
/// This effectively implements numpy's isclose function (see
/// https://numpy.org/doc/stable/reference/generated/numpy.isclose.html). As in
/// the original function, the max allowed variations from the relative
/// difference tolerance and the absolute difference tolerance are summed and
/// compared against the absolute difference.
///
/// @note
/// For less experienced C++ developers, `operator()` overloads the "function
/// call operation" (it is analogous to python's `__call__` method)
class FltIsClose {
  double rtol;
  double atol;

public:
  FltIsClose() = delete;
  FltIsClose(double rtol, double atol) : rtol{rtol}, atol{atol} {}

  /// determines whether arguments are equal within the tolerance
  bool operator()(float actual, float desired) const noexcept {
    return isclose_<true>(actual, desired, this->rtol, this->atol);
  }

  /// determines whether arguments are equal within the tolerance
  bool operator()(double actual, double desired) const noexcept {
    return isclose_<true>(actual, desired, this->rtol, this->atol);
  }

  /// describe relationship between values for which this functor returns false
  std::string describe_false() const {
    std::string rtol_str = to_pretty_string(rtol);
    std::string atol_str = to_pretty_string(atol);
    return ("unequal for the tolerance (rtol = " + rtol_str +
            ", atol = " + atol_str + ")");
  }
};

/// "functor" to check if floating point values are exactly equal
///
/// @note
/// For less experienced C++ developers, `operator()` overloads the "function
/// call operation" (it is analogous to python's `__call__` method)
struct FltIsEqual {
  /// determines whether arguments are exactly equal
  bool operator()(float actual, float desired) const noexcept {
    return (actual == desired) || std::isnan(actual) == std::isnan(desired);
  }

  /// determines whether arguments are exactly equal
  bool operator()(double actual, double desired) const noexcept {
    return (actual == desired) || std::isnan(actual) == std::isnan(desired);
  }

  /// describe relationship between values for which this functor returns false
  std::string describe_false() const { return "not exactly equal"; }
};

template <typename T>
using PtrPair = std::pair<const T*, const T*>;

/// Packages up the information for a comparison of 2 pointers
///
/// See the docstring of @ref compare_ for an extended discussion for why
/// this type actually exists.
struct CmpPack {
  std::variant<FltIsEqual, FltIsClose> cmp_fn;
  std::variant<PtrPair<float>, PtrPair<double>> actual_desired_pair;
  const bool* selection_mask;
  std::variant<IdxMapping<DataLayout::LEFT>, IdxMapping<DataLayout::RIGHT>>
      idx_mapping;
};

/// this dispatches the appropriate logic to drive the comparison
///
/// The most pragmatic approach for implementing the underlying comparisons in
/// an extendable manner (without extensive code duplication or sacrificing
/// performance) is to implement them using templates and to make the datatype
/// a template parameter.
///
/// This function was designed in a misguided attempt to shift most of the
/// implementation into source files in order to reduce compile times. In order
/// to hide all calls to a set of templates into a source file, this must be
/// a totally ordinary function that dispatches to the proper templates:
/// - thus, @ref CmpPack as a well-defined type to package up all of the
///   possible type combinations. This is achieved through the use of
///   std::variant (i.e. type-safe unions).
/// - the idea is that callers package up `CmpPack`, call this function, and
///   then function unpacks the values from `CmpPack` and dispatches to the
///   appropriate template function.
///
/// With the benefit of hindsight, this was probably all a mistake... Reducing
/// the compilation cost was probably **NOT** worth the added complexity (if
/// nothing else, we probably should have measured it first...)
testing::AssertionResult compare_(CmpPack pack);

}  // namespace grtest::arraycmp_detail

#endif  // GRTESTUTILS_GOOGLETEST_CHECK_ALLCLOSE_DETAIL_HPP