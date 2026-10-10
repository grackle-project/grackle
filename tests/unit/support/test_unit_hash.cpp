//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of the hash function
///
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>
#include "support/fnv1a_hash.hpp"
#include <iomanip>
#include <ostream>
#include <string_view>

namespace GRIMPL_NAMESPACE_DECL {
/// Teach GTest how to print HashRsltPack
/// @note it's important this is in the same namespace as HashRsltPack
void PrintTo(const HashRsltPack& pack, std::ostream* os) {
  *os << "{success: " << pack.success << ", keylen: " << pack.keylen
      << ", hash: 0x" << std::setfill('0')
      << std::setw(8)  // u32 has 8 hex digits
      << std::hex << pack.hash << "}";
}

bool operator==(const HashRsltPack& a, const HashRsltPack& b) {
  return a.success == b.success && a.keylen == b.keylen && a.hash == b.hash;
}

}  // namespace GRIMPL_NAMESPACE_DECL

// the test answers primarily came from Appendix C of
// https://datatracker.ietf.org/doc/html/draft-eastlake-fnv-17

using GRIMPL_NS::FNV1aHasher;

TEST(FNV1a, EmptyString) {
  GRIMPL_NS::HashRsltPack expected{true, 0, 0x811c9dc5ULL};
  ASSERT_EQ(FNV1aHasher<>::calc(""), expected);
  ASSERT_EQ(FNV1aHasher<>::calc(std::string_view("")), expected);
}

TEST(FNV1a, aString) {
  GRIMPL_NS::HashRsltPack expected{true, 1, 0xe40c292cULL};
  ASSERT_EQ(FNV1aHasher<>::calc("a"), expected);
  ASSERT_EQ(FNV1aHasher<>::calc(std::string_view("a")), expected);
}

TEST(FNV1a, foobarString) {
  GRIMPL_NS::HashRsltPack expected{true, 6, 0xbf9cf968ULL};
  ASSERT_EQ(FNV1aHasher<>::calc("foobar"), expected);
  ASSERT_EQ(FNV1aHasher<>::calc(std::string_view("foobar")), expected);
}

TEST(FNV1a, MaxSizeString) {
  constexpr int MaxKeyLen = 6;  // <- exactly matches the key's length
  GRIMPL_NS::HashRsltPack expected{true, MaxKeyLen, 0xbf9cf968ULL};
  ASSERT_EQ(FNV1aHasher<>::calc("foobar"), expected);
  ASSERT_EQ(FNV1aHasher<>::calc(std::string_view("foobar")), expected);
}

TEST(FNV1a, TooLongString) {
  constexpr int MaxKeyLen = 5;  // <- shorter than the queried key
  ASSERT_FALSE(FNV1aHasher<MaxKeyLen>::calc("foobar").success);
  ASSERT_FALSE(
      FNV1aHasher<MaxKeyLen>::calc(std::string_view("foobar")).success);
}
