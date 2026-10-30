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
#include <string_view>

#include "support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Holds the result of a call to fnv1a_hash
struct HashRsltPack {
  bool success;
  std::uint16_t keylen;
  std::uint32_t hash;
};

/// @brief collects methods for computing a key's length and 32-bit FNV-1a hash
///
/// @tparam MaxKeyLen the max number of characters in key (excluding '\0'). By
///     default, it's the largest value HashRsltPack::keylen holds. A smaller
///     value can be specified as an optimization.
///
/// @note
/// This hash function prioritizes convenience. We may want to
/// evaluate whether alternatives (e.g. fxhash) are faster or have fewer
/// collisions with our typical keys.
///
/// @warning
/// Obviously this is @b NOT cryptographically secure
template <int MaxKeyLen = std::numeric_limits<std::uint16_t>::max()>
struct FNV1aHasher {
  static_assert(0 <= MaxKeyLen &&
                    MaxKeyLen <= std::numeric_limits<std::uint16_t>::max(),
                "MaxKeyLen can't be encoded by HashRsltPack");

  static constexpr std::uint32_t FNV1A_PRIME = 16777619;
  static constexpr std::uint32_t FNV1A_OFFSET = 2166136261;

  /// @brief calculate the hash value
  ///
  /// @param key the null-terminated string. Behavior is deliberately undefined
  ///     when passed a `nullptr`
  static constexpr HashRsltPack calc(const char* key) {
    std::uint32_t hash = FNV1A_OFFSET;
    for (int i = 0; i <= MaxKeyLen; i++) {  // the `<=` is intentional
      if (key[i] == '\0') {
        return {.success = true,
                .keylen = static_cast<std::uint16_t>(i),
                .hash = hash};
      }
      hash = (hash ^ static_cast<std::uint8_t>(key[i])) * FNV1A_PRIME;
    }
    return {.success = false, .keylen = 0, .hash = 0};
  }

  static constexpr HashRsltPack calc(std::string_view key) {
    int len = key.size();
    if (len > MaxKeyLen) {
      return {.success = false, .keylen = 0, .hash = 0};
    }
    std::uint32_t hash = FNV1A_OFFSET;
    for (int i = 0; i < len; i++) {
      hash = (hash ^ static_cast<std::uint8_t>(key[i])) * FNV1A_PRIME;
    }
    return {
        .success = true, .keylen = static_cast<uint16_t>(len), .hash = hash};
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_FNV1A_HASH_HPP
