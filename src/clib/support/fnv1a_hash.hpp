//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// implements the 32-bit fnv-1a hash function
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_FNV1A_HASH_HPP
#define SUPPORT_FNV1A_HASH_HPP

#include <cstdint>
#include <limits>

namespace grackle::impl {

/// Holds the result of a call to fnv1a_hash
struct HashRsltPack {
  bool success;
  std::uint16_t keylen;
  std::uint32_t hash;
};

/// calculate 32-bit FNV-1a hash of key and measures the key's length.
///
/// @tparam MaxKeyLen the max number of characters in key (excluding '\0'). By
///     default, it's the largest value HashRsltPack::keylen holds. A smaller
///     value can be specified as an optimization.
/// @param key the null-terminated string. Behavior is deliberately undefined
///     when passed a `nullptr`
///
/// @note
/// The current implementation prioritizes convenience. We may want to evaluate
/// whether alternatives (e.g. fxhash) are faster or have fewer collisions with
/// our typical keys.
///
/// @warning
/// Obviously this is @b NOT cryptographically secure
template <int MaxKeyLen = std::numeric_limits<std::uint16_t>::max()>
HashRsltPack fnv1a_hash(const char* key) {
  static_assert(
      0 <= MaxKeyLen && MaxKeyLen <= std::numeric_limits<std::uint16_t>::max(),
      "MaxKeyLen can't be encoded by HashRsltPack");

  constexpr std::uint32_t prime = 16777619;
  constexpr std::uint32_t offset = 2166136261;

  std::uint32_t hash = offset;
  for (int i = 0; i <= MaxKeyLen; i++) {  // the `<=` is intentional
    if (key[i] == '\0') {
      return HashRsltPack{true, static_cast<std::uint16_t>(i), hash};
    }
    hash = (hash ^ key[i]) * prime;
  }
  return HashRsltPack{false, 0, 0};
}

}  // namespace grackle::impl

#endif  // SUPPORT_FNV1A_HASH_HPP
