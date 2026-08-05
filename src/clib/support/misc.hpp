//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define miscellaneous internal logic
///
//===----------------------------------------------------------------------===//

#ifndef SUPPORT_MISC_HPP
#define SUPPORT_MISC_HPP

#include "./config.hpp"

#include <type_traits>

namespace GRIMPL_NAMESPACE_DECL {

/// @brief a branchless choice between 2 values
///
/// The returned value is analogous to the expression `(cond) ? a : b`, except
/// that:
/// - @p a and @p b are selected based on 2 multiplications and an addition,
///   rather than a branch instruction
/// - don't forget that there is a big difference between writing
///   `branchless_choice(cond, expr_a, expr_b)` and  `(cond) ? expr_a : expr_b`,
///   where `expr_a` and `expr_b` are expressions (e.g. `x + 5`). In the former
///   case, both expressions are evaluated. In the latter case, only one
///   expression is considered
/// - this function breaks down if a or b isn't finite
///
/// @note
/// This is primarily useful in loops where a branch could prevent the use of
/// SIMD instructions (or to avoid divergence of threads in a GPU kernel)
template <typename T>
GRIMPL_FORCE_INLINE T branchless_choice(bool cond, T a, T b) {
  // some sanity checks
  using sanitized_T = std::remove_reference_t<std::remove_cv_t<T>>;
  static_assert(std::is_same_v<sanitized_T, T>);
  static_assert(std::is_arithmetic_v<T>);

  return cond * a + (!cond) * b;
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_MISC_HPP