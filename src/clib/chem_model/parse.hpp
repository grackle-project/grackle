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
#include <limits>  // std::numeric_limits
#include <string>
#include <string_view>
#include <tuple>  // std::tie
#include <utility>

#include "support/config.hpp"
#include "support/expected.hpp"
#include "support/inplace_vector.hpp"

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
  // pos tracks the most recently parsed character in the numeral
  std::size_t pos = std::min(s_stop, s.size());

  // lambda function that checks if `c` is the last character before `pos`
  auto is_last = [&](char c) {
    return (pos >= (min_start + 1) && s[pos - 1] == c);
  };

  // lambda function that reads all symbols (i.e. updates `pos` accordingly)
  // that contribute to a single decimal place digit & returns the digit
  auto read_decimal_digit = [&](char sym_1, char sym_5, char sym_10) -> int {
    bool maybe_subtractive = (pos >= (min_start + 2)) && (s[pos - 2] == sym_1);
    if (maybe_subtractive && is_last(sym_10)) {
      pos -= 2;
      return 9;
    } else if (maybe_subtractive && is_last(sym_5)) {
      pos -= 2;
      return 4;
    } else {
      int sum = 0;
      for (int i = 0; i < 3 && is_last(sym_1); i++) {
        pos--;
        sum++;
      }
      if (is_last(sym_5)) {
        pos--;
        sum += 5;
      }
      return sum;
    }
  };

  int total = 0;
  total += read_decimal_digit('I', 'V', 'X');
  total += 10 * read_decimal_digit('X', 'L', 'C');
  // total += 100 * read_decimal_digit('C', 'D', 'M');  // <- not needed!

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

/// @brief try to extract characters in the format of a nuclide symbol
///
/// This function deliberately doesn't check whether extracted symbols map to
/// to a known nuclide symbol
///
/// @param s string to be parsed
/// @param pos Index in @p s of where the symbol starts.
/// @param s_stop Index in @p s where the substring ends
/// @returns On success, produces a pair holding the first position after the
///     symbol-candidate and a view of the symbol-candidate.
constexpr Expected<std::pair<std::size_t, std::string_view>, Err>
potential_nuclide_sym(std::string_view s, std::size_t pos = 0,
                      std::size_t s_stop = std::string_view::npos) {
  std::size_t stop = std::min(s_stop, s.size());
  if (!isupper_(s[pos])) {
    return Unexpected(Err::InvalidChar(pos));
  } else if (pos + 1 == stop || !islower_(s[pos + 1])) {
    return std::pair<std::size_t, std::string_view>(pos + 1, s.substr(pos, 1));
  } else {
    return std::pair<std::size_t, std::string_view>(pos + 2, s.substr(pos, 2));
  }
}

// encodes a coefficient-string_view pair
using CoefSViewPair = std::pair<uint8_t, std::string_view>;

inline constexpr std::size_t MaxSymsPerFormula = 4;

/// Represents the parsed components of a chemical formula
struct FormulaComps {
  /// tracks the value of the roman numeral
  uint8_t roman_numeral_val;

  /// tracks subscripts of a nuclide symbol and the symbol itself
  InplaceVector<CoefSViewPair, MaxSymsPerFormula> subscript_symbol_pairs;
};

/// @brief try to parse components of a chemical formula from the string @p s
///
/// The function fails if all characters in the range
///   [`s.data() + pos`, `s.data() + stop`)
/// aren't all part of a single chemical formula
///
/// @param[in]  s string to be parsed.
/// @param[in]  pos Index in @p s of where the formula starts.
/// @param[in]  stop Index in @p s where the substring ends
constexpr inline Expected<FormulaComps, Err> formula_comps(
    std::string_view s, std::size_t pos = 0,
    std::size_t stop = std::string_view::npos) {
  // try to parse a roman numeral starting from the right edge of s
  Expected<std::pair<std::size_t, int>, Err> retval =
      roman_from_right(s, 1, stop);
  if (!retval.has_value()) {
    return Unexpected(
        Err::Custom(pos, "must end in a roman numeral (with value under 100)"));
  }

  FormulaComps comps;
  // update `stop` to point to the leftmost character in the roman numeral (we
  // took care to ensure this is always at least pos+1) and store numeral value
  std::tie(stop, comps.roman_numeral_val) = retval.value();
  while (pos < stop) {
    // try to parse a nuclide symbol
    Expected<std::pair<std::size_t, std::string_view>, Err> rv_a =
        potential_nuclide_sym(s, pos, stop);
    if (!rv_a.has_value()) {
      return Unexpected(rv_a.error());
    }
    std::string_view sym;
    std::tie(pos, sym) = rv_a.value();

    // determine subscript associated with the symbol
    Expected<std::pair<std::size_t, uint8_t>, Err> rv_b =
        positive_int<uint8_t>(s, pos, stop);
    if (pos != stop && !rv_b.has_value() && !rv_b.error().is_InvalidChar()) {
      return Unexpected(rv_b.error());
    }
    uint8_t subscript;
    std::tie(pos, subscript) = rv_b.value_or({pos, 1});

    // record the subscript-symbol pair
    if (!comps.subscript_symbol_pairs.try_emplace_back(subscript, sym)) {
      return Unexpected(Err::Custom(pos, "too many nuclide symbols"));
    }
  }

  return comps;
}

}  // namespace parse
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_PARSE_HPP
