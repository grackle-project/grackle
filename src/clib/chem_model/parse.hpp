//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares parsing logic
///
//===----------------------------------------------------------------------===//
#ifndef CHEM_MODEL_PARSE_HPP
#define CHEM_MODEL_PARSE_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>  // std::tie
#include <utility>

#include "support/config.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace parse {

/// @brief represents a parsing error
class Err {
  std::size_t pos_;    // <- string position associated with errors
  const char* descr_;  // <- nullptr cooresponds to invalid_char_err

  constexpr Err(std::size_t pos, const char* descr) noexcept
      : pos_(pos), descr_(descr) {}

public:
  /// coerce to a string representation
  ///
  /// @param parsed_s Is the string being parsed where the error arose
  std::string to_string(std::string_view parsed_s) const;

  constexpr bool is_InvalidChar() const { return descr_ == nullptr; }

  /// @brief factory method for an invalid character error
  ///
  /// @param pos the position where the error occurred
  static constexpr Err InvalidChar(std::size_t pos) noexcept {
    return {pos, nullptr};
  }

  /// @brief factory method for a custom error message
  ///
  /// @param pos the position where the error occurred
  /// @param descr a non-null literal string specifying an error message
  static constexpr Err Custom(std::size_t pos, const char* descr) noexcept {
    return {pos, descr};
  }
};

// these functions assumes that c is an ascii character
constexpr bool islower_(char c) { return 'a' <= c && c <= 'z'; }
constexpr bool isupper_(char c) { return 'A' <= c && c <= 'Z'; }
constexpr bool isdigit_(char c) { return '0' <= c && c <= '9'; }

/// @brief try to parse a right-justified roman numeral from the string @p s
///
/// This function tries to parse the largest, valid roman numeral. At the time
/// of writing, the numeral with the maximum value that can be parsed is 99.
///
/// @param s string to be parsed
/// @param min_start lowest index in @p s that can be a digit of the numeral
/// @param s_stop Index in @p s where the substring ends.
/// @returns On success, produces a pair holding index of @p s where the numeral
///     starts and the value of the numeral
inline constexpr Expected<std::pair<std::size_t, int>, Err> roman_from_right(
    std::string_view s, std::size_t min_start = 0,
    std::size_t s_stop = std::string_view::npos) {
  // summarizes roman numeral table for a given decimal place
  struct DecimalInfo {
    char tally;          // <- acts like tally marks
    char tally_cluster;  // <- acts like the 5-tally mark
    int pow_of_10;
  };
  constexpr unsigned N = 3;
  DecimalInfo di_table[N] = {{'I', 'V', 1}, {'X', 'L', 10}, {'C', '\0', 100}};

  std::size_t pos = std::min(s_stop, s.size());
  int total = 0;

  for (unsigned row = 0; row < N - 1 && pos > min_start; row++) {
    const DecimalInfo& di = di_table[row];
    bool maybe_subtractive = (pos >= min_start + 1) && (s[pos - 2] == di.tally);

    if (maybe_subtractive && s[pos - 1] == di_table[row + 1].tally) {
      total += 9 * di.pow_of_10;
      pos -= 2;
    } else if (maybe_subtractive && s[pos - 1] == di.tally_cluster) {
      total += 4 * di.pow_of_10;
      pos -= 2;
    } else {
      for (int c = 0; c < 3 && pos > min_start && di.tally == s[pos - 1]; c++) {
        total += di.pow_of_10;
        pos--;
      }
      if (pos > min_start && s[pos - 1] == di.tally_cluster) {
        total += 5 * di.pow_of_10;
        pos--;
      }
    }
  }

  if (total == 0) {
    return Unexpected(Err::InvalidChar(s.size() - 1));
  }
  return std::pair<std::size_t, int>(pos, total);
}

/// @brief try to parse a positive integer from the string @p s
///
/// This assumes decimal notation and forbids leading zeros. For the sake of
/// clarity, this is a very deliberate choice.
///
/// @param s string to be parsed
/// @param pos Index in @p s of where the integer starts.
/// @param s_stop Index in @p s where the substring ends
/// @returns On success, produces a pair holding the first position after the
///     parsed integer and the parsed integer.
template <typename IntegerT>
constexpr Expected<std::pair<std::size_t, IntegerT>, Err> positive_int(
    std::string_view s, std::size_t pos = 0,
    std::size_t s_stop = std::string_view::npos) {
  constexpr IntegerT max_val = std::numeric_limits<IntegerT>::max();
  constexpr IntegerT max_div_10 = max_val / 10;
  constexpr IntegerT digit_sentinal = max_val - max_div_10;

  std::size_t stop = std::min(s_stop, s.size());

  if (pos >= s_stop) {
    return Unexpected(Err::Custom(s_stop - 1, "s_stop must exceed pos"));
  } else if (s[pos] == '0') {
    return Unexpected(Err::Custom(pos, "leading digit of a number can't be 0"));
  }
  IntegerT accum = 0;
  std::size_t orig = pos;
  for (pos = orig; pos < stop && isdigit_(s[pos]); pos++) {
    IntegerT digit = static_cast<IntegerT>(s[pos] - '0');
    if (accum > max_div_10 || (accum == max_div_10 && digit > digit_sentinal)) {
      return Unexpected(Err::Custom(orig, "number exceeds max allowed value"));
    }
    accum = accum * static_cast<IntegerT>(10) + digit;
  }
  return std::pair<std::size_t, IntegerT>(pos, accum);
}

}  // namespace parse
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_PARSE_HPP
