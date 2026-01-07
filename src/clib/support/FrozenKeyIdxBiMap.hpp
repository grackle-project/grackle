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
/// simplistic at the cost of speed). PR #484 introduce a drop-in replacement
/// the API won't change that is much faster
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_FROZENKEYIDXBIMAP_HPP
#define SUPPORT_FROZENKEYIDXBIMAP_HPP

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "grackle.h"
#include "status_reporting.h"

namespace grackle::impl {

// the motivation for these constants are provided in PR #484 (they are related
// to some optimizations in the FrozenKeyIdxBiMap implementation)
namespace bimap_detail {
/// specifies an invalid value of the map (we state that you can't store the
/// maximum u16 value)
inline constexpr uint16_t INVALID_VAL = std::numeric_limits<uint16_t>::max();

/// specifies maximum allowed length of a key (excluding the null character).
inline constexpr uint16_t KEYLEN_MAX = 29;
}  // namespace bimap_detail

// the following doxygen comment block logically groups every all parts of
// the (internal) API for Grackle's (internal) FrozenKeyIdxBiMap. It's useful
// when generating a doxygen webpage

/// @defgroup bimap-grp FrozenKeyIdxBiMap Data Type
///
/// FrozenKeyIdxBiMap provides specialized mapping functionality for internal
/// use within Grackle. The functionality is useful as a building-block for
/// runtime lookup-tables and other data types with map-like interface.
///
/// The data type was implemented in a C-style. The FrozenKeyIdxBiMap struct
/// should be treated as an opaque type that is operated upon by a set of
/// associated functions. More idiomatic C++ (or languages like Rust & Swift),
/// the associated functions would be attached to the struct as methods
/** @{*/

/// describes the operating modes of @ref FrozenKeyIdxBiMap
enum class BiMapMode {
  /// The preferred default mode, where the creation of a BiMap involves making
  /// copies of each key (cleaning up a BiMap will deallocate the copies)
  ///
  /// In general, this is much safer, and it will be @b very useful in the
  /// longer-term if we allow dynamic extension of chemistry networks. If we
  /// adopt the embedded-key optimization (discussed in the FrozenKeyIdxBiMap),
  /// this mode will probably be significantly faster.
  COPIES_KEYDATA = 1,
  /// This mode aims to reduce memory usage by having the BiMap reference
  /// external keys. In other words, the BiMap won't attempt to manage
  /// allocations holding each character in a string.
  ///
  /// @warning
  /// For safety this should @b ONLY be used when all keys are immutable
  /// string-literals (i.e. when the strings are valid for program's duration)
  REFS_KEYDATA = 0,
};

/// @brief A bidirectional map (bimap), specialized to map @c n unique string
/// keys to unique indexes with values of @c 0 through @c (n-1) and
/// vice versa. The ordering & values of keys are set at creation and frozen.
///
/// This type is useful in a number of scenarios. For example, it can be used
/// to implement a type representing a Map of arrays (where the values could
/// be part of a single contiguous array or are individual arrays).
///
/// This type operates in 2 modes: @ref BiMapMode::COPIES_KEYDATA and
/// @ref BiMapMode::REFS_KEYDATA. Their docstrings provide more context. When
/// in doubt, prefer the former mode.
///
/// @par Replacement in PR #484
/// The current implementation is extremely oversimplified and inefficient! It
/// doesn't even use a hash table. The purpose is to create a simple abstract
/// data structure for which the implementation will be dramatically improved
/// by PR #484 (but the interface won't be touched at all).
///
/// @par
/// The PR with the improved version, also updates this docstring with a
/// detailed explanation of design decisions (like why the contents
/// are "frozen") and highlights a number of potential improvements.
///
/// @note
/// The contents of this struct should be considered an implementation
/// detail! Always prefer the associated functions (they are defined in such
/// a way that they should be inlined)
struct FrozenKeyIdxBiMap {
  // don't forget to update FrozenKeyIdxBiMap_clone when changing members

  /// the number of contained strings
  int length;
  /// array of keys
  const char** keys;
  /// specifies ownership of keys, @see BiMapMode
  BiMapMode mode;
};

/// Constructs a new FrozenKeyIdxBiMap
///
/// @param[in]  keys Sequence of 1 or more unique strings. Each string must
///     include at least 1 non-null character and be null-terminated
/// @param[in]  key_count The length of keys
/// @param[in]  mode specifies handling of keys. This will be passed on to any
///     clones that are made.
///
/// @note
/// Callers should pass the returned value to @ref FrozenKeyIdxBiMap_is_ok
/// to check whether there was an error during creation. This is pretty
/// ugly/clunky, but it's the only practical way to achieve comparable behavior
/// to other internal data types. The best alternatives involve things like
/// std::optional or converting this type to a simple C++ class.
inline FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(const char* keys[],
                                               int key_count, BiMapMode mode) {
  // this will be returned if there is an error
  FrozenKeyIdxBiMap erroneous_obj{-1, nullptr, BiMapMode::REFS_KEYDATA};

  if (keys == nullptr && key_count == 0) {
    return FrozenKeyIdxBiMap{0, nullptr, mode};
  }

  // check the specified keys
  long long max_keys = static_cast<long long>(bimap_detail::INVALID_VAL) - 1LL;
  if (key_count < 1 || static_cast<long long>(key_count) > max_keys) {
    GrPrintErrMsg("key_count must be positive and cannot exceed %lld",
                  max_keys);
    return erroneous_obj;
  } else if (keys == nullptr) {
    GrPrintErrMsg("keys must not be a nullptr");
    return erroneous_obj;
  }
  for (int i = 0; i < key_count; i++) {
    GR_INTERNAL_REQUIRE(keys[i] != nullptr, "Can't specify a nullptr key");
    std::size_t n_chrs_without_nul = std::strlen(keys[i]);
    if (n_chrs_without_nul == 0 ||
        n_chrs_without_nul > bimap_detail::KEYLEN_MAX) {
      GrPrintErrMsg(
          "calling strlen on \"%s\", the key @ index %d, yields 0 or a length "
          "exceeding %d",
          keys[i], i, bimap_detail::KEYLEN_MAX);
      return erroneous_obj;
    }
    // check uniqueness
    for (int j = 0; j < i; j++) {
      if (strcmp(keys[i], keys[j]) == 0) {
        GrPrintErrMsg("\"%s\" key repeats", keys[i]);
        return erroneous_obj;
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
    default: {
      GrPrintErrMsg("unknown mode");
      return erroneous_obj;
    }
  }

  return FrozenKeyIdxBiMap{/* length = */ key_count,
                           /* keys = */ out_keys,
                           /* mode = */ mode};
}

/// checks whether a creational function produced a valid bimap
///
/// @param[in] ptr Points to the object being checked
/// @return true if the value is ok or false if the value is invalid
///
/// @important
/// The interface of @ref FrozenKeyIdxBiMap sets values in a very particular
/// way to signal that FrozenKeyIdxBiMap is in an invalid state. This function
/// @b ONLY checks for that particular signature.
inline bool FrozenKeyIdxBiMap_is_ok(const FrozenKeyIdxBiMap* ptr) {
  return (ptr->length != -1);
}

/// Destroys the internal data tracked by an instance
///
/// @param[in] ptr A non-null pointer to a valid bimap instance
///
/// @warning
/// As with any C datatype, care is required to avoid issues with internal
/// dangling pointers. YOU SHOULD ONLY CALL THIS ONCE for a given instance
/// (and only if the instance was properly by the interface)
/// - while some efforts are made to reduce the possiblity of issues, some
///   things just can't be avoided (especially when it comes to shallow copies)
/// - here's a problematic example:
///   @code{.cpp}
///   FrozenKeyIdxBiMap bimap = new_FrozenKeyIdxBiMap( /*<args...>*/ );
///   // (the FrozenKeyIdxBiMap_is_ok check is elided for brevity)
///
///   // you should generally avoid shallow copies (if possible)
///   FrozenKeyIdxBiMap shallow_cpy = bimap;
///
///   // problems arise below (if we swap order, the 2nd call is still bad)
///   drop_FrozenKeyIdxBiMap(&shallow_cpy);  // <- this is OK
///   drop_FrozenKeyIdxBiMap(&bimap);        // <- this is BAD
///   @endcode
inline void drop_FrozenKeyIdxBiMap(FrozenKeyIdxBiMap* ptr) {
  if (ptr->mode == BiMapMode::COPIES_KEYDATA) {
    for (int i = 0; i < ptr->length; i++) {
      delete[] ptr->keys[i];
    }
  }
  if (ptr->keys != nullptr) {
    delete[] ptr->keys;
  }
}

/// Makes a clone of the specified FrozenKeyIdxBiMap
///
/// The clone inherits the original's BiMapMode value. If it held
/// BiMapMode::COPIES_KEYDATA, then fresh copies of the strings are made
///
/// @note
/// Callers should pass the returned value to @ref FrozenKeyIdxBiMap_is_ok
/// to check whether there was an error during creation. This is pretty
/// ugly/clunky, but it's the only practical way to achieve comparable behavior
/// to other internal data types. The best alternatives involve things like
/// std::optional or converting this type to a simple C++ class.
inline FrozenKeyIdxBiMap FrozenKeyIdxBiMap_clone(const FrozenKeyIdxBiMap* ptr) {
  return new_FrozenKeyIdxBiMap(ptr->keys, ptr->length, ptr->mode);
};

namespace bimap {

/// holds the result of a call to @ref FrozenKeyIdxBiMap_get
///
/// @note This C-style approximation of std::optional<uint16_t>
struct AccessRslt {
  /// Indicates whether the value member is valid
  bool has_value;
  /// the loaded value (if has_value is false then this holds garbage)
  uint16_t value;
};

}

/// lookup the value associated with the key
///
/// This is the analog to calling `map[key]` in python.
///
/// @param[in] map A pointer to a valid bimap
/// @param[in] key A null-terminated string
///
/// @return An instance of @ref bimap::AccessRslt that encodes the value (if
///     the key is present)
inline bimap::AccessRslt FrozenKeyIdxBiMap_get(
    const FrozenKeyIdxBiMap* map, const char* key) {
  GR_INTERNAL_REQUIRE(key != nullptr, "A nullptr key is forbidden");
  for (int i = 0; i < map->length; i++) {
    if (std::strcmp(map->keys[i], key) == 0) {
      return bimap::AccessRslt{true, static_cast<uint16_t>(i)};
    }
  }
  return bimap::AccessRslt{false, 0};
}

/// returns whether the map contains the key
///
/// @param[in] map A pointer to a valid bimap
/// @param[in] key A null-terminated string
inline bool FrozenKeyIdxBiMap_contains(const FrozenKeyIdxBiMap* map,
                                       const char* key) {
  return FrozenKeyIdxBiMap_get(map, key).has_value;
}

/// return the number of keys in the map
///
/// @param[in] map A pointer to a valid bimap
inline int FrozenKeyIdxBiMap_size(const FrozenKeyIdxBiMap* map) {
  return map->length;
}

/// Return the key associated with the specified value
///
/// For some context, if this function returns a string `s` for some index `i`,
/// then a call to @ref FrozenKeyIdxBiMap_get that passes `s` will
/// return `i`
///
/// This is intended for use in situations where you briefly need the string
/// (i.e. and you plan to stop using the pointer before or at the same time as
/// the @p map is destroyed). In more detail:
/// - If the @p map was constructed in @ref BiMapMode::COPIES_KEYDATA mode,
///   returned strings have the same lifetime as @p map (i.e. they are
///   deallocated when the contents of @p map are deallocated).
/// - Otherwise, the returned string's allocation is externally managed. But,
///   any scenario where the allocation doesn't live at least as long as @p map,
///   is ill-formed
///
/// @param[in] map A pointer to a valid bimap
/// @param[in] idx The index to check
/// @return The pointer to the appropriate key
inline const char* FrozenKeyIdxBiMap_inverse_get(const FrozenKeyIdxBiMap* map,
                                                 uint16_t idx) {
  if (idx >= map->length) {
    return nullptr;
  }
  return map->keys[idx];  // this can't be a nullptr
}

/** @}*/  // end of group

}  // namespace grackle::impl

#endif  // SUPPORT_FROZENKEYIDXBIMAP_HPP
