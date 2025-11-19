//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the internal FrozenKeyIdxBiMap type
///
/// All insertions occur during initialization. This simplifies a lot of
/// bookkeeping.
///
/// The underlying implementation of this type is *highly* suboptimal (it's
/// simplistic at the cost of speed). PR #270 introduce a drop-in replacement
/// that maintains the exact same API
///
//===----------------------------------------------------------------------===//
#ifndef UTILS_FROZENKEYIDXBIMAP_HPP
#define UTILS_FROZENKEYIDXBIMAP_HPP

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "grackle.h"
#include "status_reporting.h"

// the motivation for these constants are provided in PR #270 (they are related
// to some optimizations in the FrozenKeyIdxBiMap implementation)
namespace grackle::impl::bimap {
/// specifies an invalid value of the map (we state that you can't store the
/// maximum u16 value)
inline constexpr std::uint16_t invalid_val =
    std::numeric_limits<std::uint16_t>::max();

/// specifies maximum allowed length of a key (excluding the null character).
inline constexpr std::uint16_t keylen_max = 29;
}

// these are just here for to make it easier for us to adopt changes from PR
// #270 (then we can delete these macros)
#define STRU16MAP_INVALID_VAL grackle::impl::bimap::invalid_val
#define STRU16MAP_KEYLEN_MAX grackle::impl::bimap::keylen_max

namespace grackle::impl {

enum class BiMapMode {
  REFS_KEYDATA = 0,
  COPIES_KEYDATA = 1
};

/// @brief This is a bidirectional map (bimap). It is specialized to map `n`
/// unique string keys to unique indexes with values of `0` through `n-1` and
/// vice versa. The ordering of keys is set at initialization and frozen.
///
/// This is primarily intended to be used in the implementation of Maps of
/// arrays (where the values could be part of a single contiguous array or are
/// individual arrays), but this but may be broadly useful for other
/// applications.
///
/// This operates in 2 modes:
/// 1. @ref BiMapMode::REFS_KEYDATA This is the default, where we operate
///    under the assumption that the allocations holding the string characters
///    outlive  the bimap. In this mode the bimap is intended to hold
///    string-literals. (which are live for the entirety of a program). This
///    minimizes memory usage.
/// 2. @ref BiMapMode::COPIES_KEYDATA Under this mode, the bimap copies the
///    data of all keys. This is useful for testing purposes. In the long-term,
///    if we allow dynamic extension of chemistry networks, it will also be
///    useful. If we are implement the optimizations described down below
///    (where we directly embed the string in the hash-table-rows), this will
///    probably be a quite a bit faster
///
/// Replacement in PR #270
/// ======================
/// The current implementation is extremely oversimplified and inefficient! It
/// doesn't even use a hash table. The purpose is to create a simple abstract
/// data structure for which the implementation will be dramatically improved
/// by PR #270 (but the interface won't be touched at all).
///
/// The PR with the improved version, also updates this docstring with a
/// detailed explanation of design decisions (like why the contents
/// are "frozen") and highlights a number of potential improvements.
///
/// > [!note]
/// > The contents of this struct should be considered an implementation
/// > detail! Always prefer the associated functions (they are defined in such
/// > a way that they should be inlined
struct FrozenKeyIdxBiMap{
  // don't forget to update FrozenKeyIdxBiMap_clone when changing members

  /// the number of contained strings
  int length;
  /// array of keys
  const char** keys;
  /// indicates whether the map "owns" the memory holding the characters in
  /// each key or just references it
  BiMapMode mode;
};

/// Constructs a new FrozenKeyIdxBiMap
///
/// @param[out] out Pointer where the allocated type is stored
/// @param[in]  keys Sequence of 1 or more unique strings. Each string must
///     include at least 1 non-null character and be null-terminated
/// @param[in]  key_count The length of keys
/// @param[in]  mode specifies handling of keys. This will be passed on to any
///     clones that are made.
inline int new_FrozenKeyIdxBiMap(
  FrozenKeyIdxBiMap** out, const char* keys[], int key_count,
  BiMapMode mode
) {

  // check the specified keys
  long long max_keys = static_cast<long long>(bimap::invalid_val) - 1LL;
  if (key_count < 1 || static_cast<long long>(key_count) > max_keys) {
    return GrPrintAndReturnErr(
      "key_count must be positive and cannot exceed %lld", max_keys
    );
  } else if (keys == nullptr) {
    return GrPrintAndReturnErr("keys must not be a nullptr");
  }
  for (int i = 0; i < key_count; i++) {
    GR_INTERNAL_REQUIRE(keys[i] != nullptr, "Can't specify a nullptr key");
    std::size_t n_chrs_without_nul = std::strlen(keys[i]);
    if (n_chrs_without_nul == 0 || n_chrs_without_nul > bimap::keylen_max) {
      return GrPrintAndReturnErr(
        "calling strlen on \"%s\", the key @ index %d, yields 0 or a length "
        "exceeding %d",
        keys[i], i, bimap::keylen_max
      );
    }
    // check uniqueness
    for (int j = 0; j < i; j++) {
      if (strcmp(keys[i], keys[j]) == 0) {
        return GrPrintAndReturnErr("\"%s\" key repeats", keys[i]);
      }
    }
  }

  // now, actually construct the result
  const char** out_keys = nullptr;
  switch (mode) {
    case BiMapMode::REFS_KEYDATA: {
      out_keys = new const char*[key_count];
      for (int i = 0; i < key_count; i++) {
        out_keys[i] = keys[i];
      }
      break;
    }
    case BiMapMode::COPIES_KEYDATA: {
      char** tmp_keys = new char*[key_count];
      for (int i = 0; i < key_count; i++) {
        std::size_t n_chrs_without_nul = std::strlen(keys[i]);
        tmp_keys[i] = new char[n_chrs_without_nul + 1];
        std::memcpy(tmp_keys[i], keys[i], n_chrs_without_nul + 1);
      }
      out_keys = (const char**)tmp_keys;
      break;
    }
    default: return GrPrintAndReturnErr("unknown mode");
  }

  (*out) = new FrozenKeyIdxBiMap;
  (*out)->length = key_count;
  (*out)->keys = out_keys;

  return GR_SUCCESS;
}

/// Destroys the specified FrozenKeyIdxBiMap
inline void drop_FrozenKeyIdxBiMap(FrozenKeyIdxBiMap* ptr) {
  if (ptr->mode == BiMapMode::COPIES_KEYDATA) {
    for (int i = 0; i < ptr->length; i++) {
      delete[] ptr->keys[i];
    }
  }
  delete[] ptr->keys;
  delete ptr;
}


/// Makes a clone of the specified FrozenKeyIdxBiMap (the clone inherites the
/// original BiMapMode).
int FrozenKeyIdxBiMap_clone(FrozenKeyIdxBiMap** out,
                            const FrozenKeyIdxBiMap* ptr) {
  return new_FrozenKeyIdxBiMap(out, ptr->keys, ptr->length, ptr->mode);
};

/// returns the value associated with the key or (if the key can't be found)
/// @ref grackle::impl::bimap::invalid_val
inline std::uint16_t FrozenKeyIdxBiMap_idx_from_key(
  const FrozenKeyIdxBiMap* map, const char* key
)
{
  GR_INTERNAL_REQUIRE(key != nullptr, "A nullptr key is forbidden");
  for (int i = 0; i < map->length; i++) {
    if (std::strcmp(map->keys[i], key) == 0) {
      return static_cast<std::uint16_t>(i);
    }
  }
  return bimap::invalid_val;
}

/// checks if the map contains a key
inline int FrozenKeyIdxBiMap_contains(const FrozenKeyIdxBiMap* map,
                                      const char* key) {
  return FrozenKeyIdxBiMap_idx_from_key(map, key) != bimap::invalid_val;
}

inline int FrozenKeyIdxBiMap_size(const FrozenKeyIdxBiMap* map)
{ return map->length; }

/// Return the ith key (this is effectively a reverse lookup)
inline const char* FrozenKeyIdxBiMap_key_from_idx(
  const FrozenKeyIdxBiMap* map, std::uint16_t i
) {
  if (i >= map->length) { return nullptr; }
  return map->keys[i]; // this can't be a nullptr
}

}  // namespace grackle::impl

#endif  // UTILS_FROZENKEYIDXBIMAP_HPP
