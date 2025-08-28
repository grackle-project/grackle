// See LICENSE file for license and copyright information

/// @file FrozenKeyIdxBiMap_details.h
/// @brief Defines the hashtable machinery used to help implement the internals
///        of the FrozenKeyIdxBiMap type
///
/// In more detail, the type internally uses a hash table to associate strings
/// with a u16 (a 16-bit unsigned integer or uint16_t).

#ifndef FROZEN_KEY_IDX_BIMAP_DETAILS_H
#define FROZEN_KEY_IDX_BIMAP_DETAILS_H

#include <stdint.h> // uint16_t, UINT16_MAX
#include <string.h> // memcmp
#include "status_reporting.h"

#ifdef __cplusplus
extern "C" {
#endif

// this lists an invalid value of the map (we state that you can't store the
// max u16 value)
#define STRU16MAP_INVALID_VAL UINT16_MAX

// this is the max number of characters allowed in a key (excluding the null
// character). While this may seem a little low, it will enable some really
// cool optimizations
#define STRU16MAP_KEYLEN_MAX 29


// keylen is the length of key (excluding the null character)
struct hash_rslt_pair_{ uint16_t keylen; uint32_t hash; };

/// calculate hash of key (and measure the length of the key). A key length of
/// 0 indicates an error
///
/// @param key a null-terminated string
///
/// We use 32bit FNV-1 hash (this should be checked!). We may want to consider
/// a faster hash (maybe fxhash) or check hash functions against the known keys
static inline struct hash_rslt_pair_ hash_from_str_(const char* key) {

  const uint32_t prime = 16777619;
  const uint32_t offset = 2166136261;

  // initialize to a value denoting an error
  struct hash_rslt_pair_ out = { 0, 0 };

  // here we need to short-circuit
  if (key[0] == '\0') { return out; }

  uint32_t hash = offset;
  for (int i = 0; i <= STRU16MAP_KEYLEN_MAX; i++) { // the `<=` is intentional
    char cur = key[i];
    if (cur == '\0') {
      out.keylen = i; // since i is 0-indexed, we don't subtract by one
      out.hash = hash;
      return out;
    }

    // update hash
    hash = hash * prime;
    hash = hash ^ cur;
  }

  return out; // remember, this was initialized to denote a problem
}

/// represents the result of an internal searching for a key
///
/// @note
/// As a rule of thumb, it's generally better (for compiler optimization) to
/// return a struct of integers than rely modifying pointer arguments 
struct StrU16Search_{
  /// specifies the value found by the search (or STRU16MAP_INVALID_VAL)
  uint16_t val;
  /// specified the number of probes before returning
  int probe_count;
  /// specify the index of the "row" corresponding to the search result
  int rowidx;
};

/// entry in the hash table.
///
/// the members are ordered to minimize the struct size (i.e. smallest members
/// are listed first) try to pack as many entries into a cacheline as possible
struct StrU16Row_{

  /// specifies the value associated with the current key 
  uint16_t value;

  /// specifies the length of the key (not including a null character)
  ///
  /// Included to short-circuit comparisons (to try to speed up probing when
  /// collisions occur)
  uint16_t keylen;

  /// identifies the address of this entry's key
  const char* key;  
};

/// Search for the row matching key. The search ends when a match is found, the
/// an empty row is found, or the function has probed `max_probe` entries
///
/// @param rows an array of rows to search to be compared
/// @param key the key to be compared
/// @param capacity the length of the rows array
/// @param max_probe the maximum number of rows to check before giving up
///
/// @important
/// The behavior is undefined if key is NULL, keylen is 0, keylen exceeds
/// STRU16MAP_KEYLEN_MAX or strlen(key) != keylen
///
/// @note
/// This is declared as `static inline` to facillitate inlining within
/// FrozenKeyIdxBiMap's interface API.
static inline struct StrU16Search_ StrU16Map_find_match_or_empty_(
  const struct StrU16Row_* rows, const char* key, int capacity, int max_probe
) {
  if ((max_probe < 0) || (max_probe > capacity)) { max_probe = capacity; }

  GR_INTERNAL_REQUIRE(key != NULL, "A nullptr key is forbidden");
  struct hash_rslt_pair_ hash_rslt = hash_from_str_(key);
  const uint16_t keylen = hash_rslt.keylen;
  int i = (int)(hash_rslt.hash % capacity); // <- 1st index to check

  for (int probe_count = 1; probe_count <= max_probe; probe_count++) {
    const struct StrU16Row_ row = rows[i];

    if (rows[i].keylen == 0) { // order matters (to handle hash_rslt.keylen==0)
      struct StrU16Search_ out = {STRU16MAP_INVALID_VAL, probe_count, i};
      return out;
    } else if ((row.keylen == keylen) && (memcmp(row.key, key, keylen) == 0)) {
      struct StrU16Search_ out = {rows[i].value, probe_count, i};
      return out;
    }

    i = (i != 0) ? i - 1 : capacity - 1; // prep for next pass thru loop
  }

  struct StrU16Search_ out = {STRU16MAP_INVALID_VAL, max_probe, i};
  return out;
}


#ifdef __cplusplus
} // extern "C"
#endif

#endif /* FROZEN_KEY_IDX_BIMAP_DETAILS_H */
