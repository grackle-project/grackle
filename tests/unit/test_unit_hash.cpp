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

namespace grackle::impl {
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

}  // namespace grackle::impl

// the test answers primarily came from Appendix C of
// https://datatracker.ietf.org/doc/html/draft-eastlake-fnv-17

TEST(FNV1a, EmptyString) {
  grackle::impl::HashRsltPack expected{true, 0, 0x811c9dc5ULL};
  ASSERT_EQ(grackle::impl::fnv1a_hash(""), expected);
}

TEST(FNV1a, aString) {
  grackle::impl::HashRsltPack expected{true, 1, 0xe40c292cULL};
  ASSERT_EQ(grackle::impl::fnv1a_hash("a"), expected);
}

TEST(FNV1a, foobarString) {
  grackle::impl::HashRsltPack expected{true, 6, 0xbf9cf968ULL};
  ASSERT_EQ(grackle::impl::fnv1a_hash("foobar"), expected);
}

TEST(FNV1a, MaxSizeString) {
  constexpr int MaxKeyLen = 6;  // <- exactly matches the key's length
  grackle::impl::HashRsltPack expected{true, 6, 0xbf9cf968ULL};
  ASSERT_EQ(grackle::impl::fnv1a_hash("foobar"), expected);
}

TEST(FNV1a, TooLongString) {
  constexpr int MaxKeyLen = 5;  // <- shorter than the queried key
  ASSERT_FALSE(grackle::impl::fnv1a_hash<MaxKeyLen>("foobar").success);
}
