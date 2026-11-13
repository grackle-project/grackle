//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of parsing roman numerals
///
//===----------------------------------------------------------------------===//

#include <array>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "chem_model/parse.hpp"
#include "support/expected.hpp"

using RomanRslt =
    GRIMPL_NS::Expected<std::pair<std::size_t, int>, GRIMPL_NS::parse::Err>;

/// a roman numeral string specifying an integer ``i`` is stored at index ``i``
/// of this array.
///
/// We only store values for the first 99 roman numerals. Index 0 holds
/// a dummy value
static constexpr std::array<const char*, 100> roman_numeral_strs = {
    "",       "I",       "II",     "III",     "IV",       "V",      "VI",
    "VII",    "VIII",    "IX",     "X",       "XI",       "XII",    "XIII",
    "XIV",    "XV",      "XVI",    "XVII",    "XVIII",    "XIX",    "XX",
    "XXI",    "XXII",    "XXIII",  "XXIV",    "XXV",      "XXVI",   "XXVII",
    "XXVIII", "XXIX",    "XXX",    "XXXI",    "XXXII",    "XXXIII", "XXXIV",
    "XXXV",   "XXXVI",   "XXXVII", "XXXVIII", "XXXIX",    "XL",     "XLI",
    "XLII",   "XLIII",   "XLIV",   "XLV",     "XLVI",     "XLVII",  "XLVIII",
    "XLIX",   "L",       "LI",     "LII",     "LIII",     "LIV",    "LV",
    "LVI",    "LVII",    "LVIII",  "LIX",     "LX",       "LXI",    "LXII",
    "LXIII",  "LXIV",    "LXV",    "LXVI",    "LXVII",    "LXVIII", "LXIX",
    "LXX",    "LXXI",    "LXXII",  "LXXIII",  "LXXIV",    "LXXV",   "LXXVI",
    "LXXVII", "LXXVIII", "LXXIX",  "LXXX",    "LXXXI",    "LXXXII", "LXXXIII",
    "LXXXIV", "LXXXV",   "LXXXVI", "LXXXVII", "LXXXVIII", "LXXXIX", "XC",
    "XCI",    "XCII",    "XCIII",  "XCIV",    "XCV",      "XCVI",   "XCVII",
    "XCVIII", "XCIX"};

TEST(Roman, EmptyString) {
  EXPECT_FALSE(GRIMPL_NS::parse::roman_from_right("").has_value());
}

TEST(Roman, InvalidMinStart) {
  std::string s = "I";
  for (std::size_t min_start : {s.size(), s.size() + 2}) {
    EXPECT_FALSE(GRIMPL_NS::parse::roman_from_right(s, min_start).has_value())
        << "for min_start = " << min_start;
  }
}

using testing::Pair;

TEST(Roman, SimpleVal) {
  for (int i = 1; i <= 99; i++) {
    std::string_view s = roman_numeral_strs[i];
    RomanRslt retval = GRIMPL_NS::parse::roman_from_right(s);
    ASSERT_TRUE(retval.has_value())
        << "error parsing the roman numeral: \"" << s << '"';
    EXPECT_THAT(retval.value(), Pair(0, i))
        << "wrong answer when parsing the roman numeral: \"" << s << '"';
  }
}

TEST(Roman, CheckIgnoreExtraI) {
  // we want to make sure that an extra I prepended to the start of the roman
  // numeral gets  properly ignored
  for (int i = 1; i <= 99; i++) {
    if (i <= 2 || i == 5 || i == 10) {
      continue;
    }
    std::string s = std::string("I") + roman_numeral_strs[i];
    RomanRslt retval = GRIMPL_NS::parse::roman_from_right(s);
    ASSERT_TRUE(retval.has_value())
        << "error parsing roman numeral @ end of the string: \"" << s << '"';
    EXPECT_THAT(retval.value(), Pair(1, i))
        << "parsing numeral @ the end of the the string, \"" << s
        << "\" gives the wrong answer";
  }
}

TEST(Roman, CheckIgnoreExtraV) {
  // we want to make sure that an extra V prepended to the start of the roman
  // numeral gets  properly ignored
  for (int i = 1; i <= 99; i++) {
    if (i <= 3) {
      continue;
    }
    std::string s = std::string("V") + roman_numeral_strs[i];
    RomanRslt retval = GRIMPL_NS::parse::roman_from_right(s);
    ASSERT_TRUE(retval.has_value())
        << "error parsing roman numeral @ end of the string: \"" << s << '"';
    EXPECT_THAT(retval.value(), Pair(1, i))
        << "parsing numeral @ the end of the the string, \"" << s
        << "\" gives the wrong answer";
  }
}

TEST(Roman, CheckIgnoreExtraX) {
  // we want to make sure that an extra X prepended to the start of the roman
  // numeral gets  properly ignored
  for (int i = 1; i <= 99; i++) {
    if (i <= 29 || (50 <= i && i <= 59) || i == 100) {
      continue;
    }
    std::string s = std::string("X") + roman_numeral_strs[i];
    RomanRslt retval = GRIMPL_NS::parse::roman_from_right(s);
    ASSERT_TRUE(retval.has_value())
        << "error parsing roman numeral @ end of the string: \"" << s << '"';
    EXPECT_THAT(retval.value(), Pair(1, i))
        << "parsing numeral @ the end of the the string, \"" << s
        << "\" gives the wrong answer";
  }
}

TEST(Roman, CheckIgnoreExtraL) {
  // we want to make sure that an extra L prepended to the start of the roman
  // numeral gets properly ignored
  for (int i = 1; i <= 99; i++) {
    if (i <= 49) {
      continue;
    }
    std::string s = std::string("L") + roman_numeral_strs[i];
    RomanRslt retval = GRIMPL_NS::parse::roman_from_right(s);
    ASSERT_TRUE(retval.has_value())
        << "error parsing roman numeral @ end of the string: \"" << s << '"';
    EXPECT_THAT(retval.value(), Pair(1, i))
        << "parsing numeral @ the end of the the string, \"" << s
        << "\" gives the wrong answer";
  }
}