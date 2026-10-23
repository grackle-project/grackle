//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines some miscellaneous math functions
///
//===----------------------------------------------------------------------===//
#ifndef MATH_MISC_HPP
#define MATH_MISC_HPP

#include "support/config.hpp"

#include <cmath>

namespace GRIMPL_NAMESPACE_DECL {

///@{
/// Implements alternatives to fmax that take more args
///
/// @note
/// These functions **only** exists to support transcription
inline float fmax(float a, float b, float c) {
  return std::fmax(a, std::fmax(b, c));
}
inline double fmax(double a, double b, double c) {
  return std::fmax(a, std::fmax(b, c));
}
///@}

///@{
/// Implements alternatives to fmin that take more args
///
/// @note
/// These functions **only** exists to support transcription
inline float fmin(float a, float b, float c) {
  return std::fmin(a, std::fmin(b, c));
}
inline float fmin(float a, float b, float c, float d) {
  return std::fmin(std::fmin(a, b), std::fmin(c, d));
}
inline double fmin(double a, double b, double c) {
  return std::fmin(a, std::fmin(b, c));
}
inline double fmin(double a, double b, double c, double d) {
  return std::fmin(std::fmin(a, b), std::fmin(c, d));
}
///@}

///@{
/// @brief Returns the closest value to the value ``v``, that lies within the
///     bounded interval [``lo``, ``hi``].
///
/// @todo Consider replacing these functions with ``std::clamp`` in the future.
///     The main concern with doing so is that the ``<algorithm>`` header
///     (where ``std::clamp`` is declared) is fairly template-heavy and could
///     slightly worsen compilation times
inline int clamp(int v, int lo, int hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
inline long clamp(long v, long lo, long hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
inline long long clamp(long long v, long long lo, long long hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
///@}

///@{
/// @brief Returns the closest value to the floating point value ``v``, that
///     lies within the bounded interval [``lo``, ``hi``].
///
/// In contrast to ``std::clamp``, this functions provide consistent behvaior
/// when a NaN is passed into the function.
/// - For context, passing a NaN to ``std::clamp`` triggers undefined behavior
/// - in the future, we might want to revist the precise behavior when a NaN
///   is passed
inline float fclamp(float v, float lo, float hi) {
  return std::fmax(lo, std::fmin(v, hi));
}
inline double fclamp(double v, double lo, double hi) {
  return std::fmax(lo, std::fmin(v, hi));
}
///@}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MATH_MISC_HPP