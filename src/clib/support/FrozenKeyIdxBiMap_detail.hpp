//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines the hash table machinery used to implement the internals of the
/// FrozenKeyIdxBiMap type
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_FROZENKEYIDXBIMAP_DETAIL_HPP
#define SUPPORT_FROZENKEYIDXBIMAP_DETAIL_HPP

#include <algorithm>  // std:::min
#include <cstdint>
#include <cstring>  // std::memcmp
#include <limits>
#include "status_reporting.h"

#include "fnv1a_hash.hpp"

namespace grackle::impl {

/// This namespace holds a generic constants and typedefs
namespace bimap_detail {
/// specifies an invalid value of the map
///
/// In other words, a map must have fewer entries than the maximum possible u16
/// value
inline constexpr uint16_t INVALID_VAL = std::numeric_limits<uint16_t>::max();

/// specifies maximum allowed length of a key (excluding the null terminator).
///
/// @note
/// While the value may seem low, it's probably large enough. Restricting
/// strings to 30 elements (including the terminator) allows for a hypothetical
/// optimization (see the @ref FrozenKeyIdxBiMap docstring for info).
inline constexpr uint16_t KEYLEN_MAX = 29;

/// the type used for indexing the rows of the BiMap's internal hash table
///
/// @note
/// This does not necessarily need to be a uint16_t (even if the internal
/// values of the hash table are uint16_t). It just needs to be able to hold
/// the maximum possible value
typedef uint16_t rowidx_type;

}  // namespace bimap_detail

// -----------------------------------------------------------------

/// Encloses code and logic pertaining to the capacity of a FrozenKeyIdxBiMap
namespace bimap_cap_detail {

/// the load factor specifies the fraction of the capacity of the Hash
/// table that is filled. This should be an integer.
///
/// Generally, the larger this is, the fewer collisions there are, but the more
/// memory is required. Lookups probably get slower if it's too big
inline constexpr int INVERSE_LOAD_FACTOR = 2;
static_assert(INVERSE_LOAD_FACTOR > 1);

/// list of allowed capacities. These are all prime numbers that
/// have nearly constant linear spacing (it may make more sense to have
/// logarithmic spacing)
inline constexpr uint32_t CAPACITIES_LIST[] = {
    // each number increases by ~10 in this first batch
    7, 19, 31, 41, 53, 61, 71, 83, 89, 101, 113, 127, 139, 149, 163, 173, 181,
    191, 199, 211, 223, 233, 241, 251,
    // probably won't use these (so we increase spacing:
    293, 401, 503, 601, 701, 797, 907, 997};

inline constexpr int N_CAPACITIES = sizeof(CAPACITIES_LIST) / sizeof(uint32_t);

/// compute the maximum number of keys
inline bimap_detail::rowidx_type max_key_count() {
  return std::min<bimap_detail::rowidx_type>(
      std::numeric_limits<bimap_detail::rowidx_type>::max(),
      CAPACITIES_LIST[N_CAPACITIES - 1] / INVERSE_LOAD_FACTOR);
}

/// compute the capacity of the map (a value of 0 indicates that a large enough
/// capacity can't be found)
///
/// @param key_count the desired number of keys (should be positive)
inline uint16_t calc_map_capacity(int key_count) {
  uint64_t c = INVERSE_LOAD_FACTOR * static_cast<uint64_t>(key_count);
  for (int i = 0; i < N_CAPACITIES; i++) {  // binary search may be faster
    if (c < CAPACITIES_LIST[i]) {
      return static_cast<uint16_t>(CAPACITIES_LIST[i]);
    }
  }
  return 0;
}

}  // namespace bimap_cap_detail

// -----------------------------------------------------------------

/// Holds machinery for hash tablea that implement FrozenKeyIdxBiMap
///
/// This machinery is specialized for (string, u16) key-value pairs
namespace bimap_StrU16_detail {

/// entry in a hash table
///
/// This acts as a (key,value) pair with a little extra metadata. A hash table
/// is fundamentally an array of these instances
///
/// @note members are ordered to minimize the struct size (i.e. smallest members
/// listed first) to pack as many entries into a cacheline as possible
struct Row {
  /// specifies the value associated with the current key
  uint16_t value;

  /// specifies the length of the key (not including the '\0')
  ///
  /// @note Tracked for short-circuiting comparisons (while probing collisions)
  uint16_t keylen;
  /// identifies the address of this entry's key
  const char* key;
};

static void overwrite_row(Row* row, const char* key, uint16_t keylen,
                          uint16_t value, bool copy_key_data) {
  GR_INTERNAL_REQUIRE(row->keylen == 0, "Sanity check failed!");
  row->value = value;
  row->keylen = keylen;
  const char* key_ptr = key;
  if (copy_key_data) {
    std::size_t total_len = keylen + 1;  // <- add 1 to account for '\0'
    char* ptr = new char[total_len];
    std::memcpy(ptr, key, total_len);
    key_ptr = ptr;
  }
  *row = Row{value, keylen, key_ptr};
}

/// represents the result of an internal search for a key
///
/// @note
/// As a rule of thumb, it's generally better (for compiler optimization) to
/// return a struct of integers than rely on modifying pointer arguments
struct SearchRslt {
  /// specifies value found by the search (or @ref bimap_detail::INVALID_VAL)
  uint16_t val;
  /// specified the number of probes before the search returned
  int probe_count;
  /// specify the index of the "row" corresponding to the search result
  int rowidx;
};

/// Search for the row matching key. The search ends when a match is found, an
/// an empty row is found, or the function has probed `max_probe` entries
///
/// @param rows an array of rows to search to be compared
/// @param key the key to be compared
/// @param capacity the length of the rows array
/// @param max_probe the maximum number of rows to check before giving up
///
/// @important
/// The behavior is undefined if @p key is a @c nullptr, @p keylen is 0, or
/// @p keylen exceeds @p bimap::keylen
///
/// @note
/// This is declared as `static inline` to facilitate inlining within
/// FrozenKeyIdxBiMap's interface API.
inline SearchRslt search(const Row* rows, const char* key, int capacity,
                         int max_probe) {
  GR_INTERNAL_REQUIRE(key != nullptr, "Major programming oversight");
  max_probe = (max_probe <= 0 || max_probe > capacity) ? capacity : max_probe;

  HashRsltPack h = fnv1a_hash<bimap_detail::KEYLEN_MAX>(key);
  int i = -1;  // <- set to a dummy value
  int launched_probes = 0;
  if (h.keylen > 0 && h.success && max_probe > 0) {
    int guess_i = static_cast<int>(h.hash % capacity);  // <- initial guess

    do {  // circularly loop over rows to search for key (start at guess_i)
      i = (guess_i + launched_probes) % capacity;
      launched_probes++;  // <- about to perform a new probe
      const Row& r = rows[i];

      if (r.keylen == h.keylen && std::memcmp(r.key, key, h.keylen) == 0) {
        return SearchRslt{r.value, launched_probes, i};  // match found!
      }

      // check if rows[i] is empty or if we have hit the limit on searches
    } while (rows[i].keylen != 0 && launched_probes < max_probe);
  }

  return SearchRslt{bimap_detail::INVALID_VAL, launched_probes, i};
}

}  // namespace bimap_StrU16_detail
}  // namespace grackle::impl
#endif  // SUPPORT_FROZENKEYIDXBIMAP_DETAIL_HPP
