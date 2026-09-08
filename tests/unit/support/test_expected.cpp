//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of Expected
///
//===----------------------------------------------------------------------===//

#include <cerrno>   // errno
#include <cstdlib>  // strtol
#include <string>

#include <gtest/gtest.h>
#include "support/expected.hpp"

using GRIMPL_NS::Expected;    // <- this wouldn't be necessary in real code
using GRIMPL_NS::Unexpected;  // <- this wouldn't be necessary in real code

// here is a sample usecase
// ========================
enum class ParseErr { OutOfRange, EmptyInput, InvalidInput };

Expected<unsigned long, ParseErr> parse_ulong(const std::string& s) {
  if (s.size() == 0) {
    return Unexpected(ParseErr::EmptyInput);
  }
  char* end;
  errno = 0;
  unsigned long v = std::strtoul(s.data(), &end, 10);
  if (errno == ERANGE) {
    return Unexpected(ParseErr::OutOfRange);
  } else if (s.data() + s.size() != end) {
    return Unexpected(ParseErr::InvalidInput);
  } else {
    return v;
  }
}

TEST(ExpectedSample, Successful) {
  std::string s = "1234";
  Expected<unsigned long, ParseErr> rslt = parse_ulong(s);
  ASSERT_TRUE(rslt.has_value());
  ASSERT_TRUE((bool)rslt);
  EXPECT_EQ(rslt.value(), 1234UL);
  EXPECT_EQ(*rslt, 1234UL);
}

TEST(ExpectedSample, ErrEmptyInput) {
  Expected<unsigned long, ParseErr> rslt = parse_ulong("");
  ASSERT_FALSE(rslt.has_value());
  ASSERT_FALSE((bool)rslt);
  EXPECT_EQ(rslt.error(), ParseErr::EmptyInput);
}

TEST(ExpectedSample, ErrInvalidInput) {
  Expected<unsigned long, ParseErr> rslt = parse_ulong("NotANumber");
  ASSERT_FALSE(rslt.has_value());
  ASSERT_FALSE((bool)rslt);
  EXPECT_EQ(rslt.error(), ParseErr::InvalidInput);
}