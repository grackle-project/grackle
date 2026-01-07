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
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_FROZENKEYIDXBIMAP_HPP
#define SUPPORT_FROZENKEYIDXBIMAP_HPP

#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "status_reporting.h"
#include "FrozenKeyIdxBiMap_detail.hpp"

namespace grackle::impl {

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
/// @par Implementation Notes
/// At the time of writing, the type is primarily implemented in terms of
/// a hash table that uses open-addressing with linear probing to resolve
/// collisions. The implementation heavily draws from logic I wrote for Enzo-E:
///   https://github.com/enzo-project/enzo-e/blob/main/src/Cello/view_StringIndRdOnlyMap.hpp
/// (More details are provided below under C++ considerations)
///
/// @par Why Frozen?
/// The contents are "frozen" for 3 primary reasons:
/// 1. It drastically simplifies the implementation (we don't have to worry
///    about deletion -- which can be quite messy)
/// 2. Linear-probing generally provides better data locality than other hash
///    collision resolution techniques, but generally has other drawbacks.
///    Freezing the contents lets us mitigate many drawbacks (mostly related to
///    the deletion operation)
/// 3. It could let us make copy operations cheaper. If we know the map won't
///    change, we could just use reference counting.
///
/// @par Consideration: Reference Counting
/// The original C++ leverages @c std::shared_ptr to achieve reference counting
/// (and reduce the cost of copying). Theoretically, I would like to see us use
/// some kind of reference-counting too. But this is tricky in library code,
/// given the diversity of threading libraries that are not formally
/// interoperable. I think the only way to properly do this would be to come up
/// with a system for allowing registration of locks/atomics with Grackle as a
/// whole.
///
/// @par C++ Considerations
/// It would definitely be worth evaluating whether we should embrace C++
/// in order to convert this to a full-blown class and adopt characteristics
/// present in the original Enzo-E version:
/// - most importantly, it would greatly reduce the chance of memory leaks
/// - (much less importantly) it would be a lot more ergonomic (& less clunky)
/// - But, for reasons expressed above, I am concerned about using
///   @c std::shared_ptr for reference counting.
///
/// @par
/// I would be stunned if <tt> std::map<std::string, uint16_t> </tt> or
/// <tt> std::map<const char*, uint16_t> </tt> is faster than the internal
/// hash table since @c std::map is usually implemented as a tree.
///
/// @par Potential Improvements
/// Simple Ideas:
/// - We could be smarter about the order that we insert keys into the table
///   (in the constructor) to minimize the search time.
/// - We might be able to come up with a better hash function
///
/// @par
/// A more ambitious idea is to embed string allocations within the rows for
/// @ref BiMapMode::COPIES_KEYDATA mode. This is possible thanks to the fact
/// that we use @ref bimap_detail::KEYLEN_MAX to limit the size of keys.
/// - Essentially, we would replace @ref bimap_StrU16_detail::Row with
///   something like the following:
///   @code{.cpp}
///   struct alignas(32) PackedRow { char data[32]; };
///
///   bool is_empty(PackedRow& r) { return data[0] == '\0' }
///   const char* get_key(PackedRow r) { return r.data; }
///   uint16_t get_val(Packed r) {
///     stdd::uint16_t o;
///     std::memcpy(&o, r.data+30, 2);
///     return o;
///   }
///   @endcode
/// - additional context about the preceding snippet:
///   - when empty, a @c PackedRow::data is filled with '\0'
///   - otherwise, @c PackedRow::data encodes the key-value pair:
///     - data[0:30] is the null-terminated key string ('\0' fills unused space)
///     - data[30:32] encodes the 16-bit value
///   - @c alignas(32) primarily ensures better cacheline alignment.
/// - Benefits of this change:
///   1. better locality (if @c PackedRow is in the cache, so is the key-string)
///   2. probing can use memcmp without a checking whether a row is empty
///   3. with a little extra care, we could use the forced alignment of
///      @c PackedRow::data to compare strings with SIMD instructions
///
/// @note
/// The contents of this struct should be considered an implementation
/// detail! Always prefer the associated functions (they are defined in such
/// a way that they should be inlined)
struct FrozenKeyIdxBiMap {
  // don't forget to update FrozenKeyIdxBiMap_clone when changing members

  /// the number of contained strings
  bimap_detail::rowidx_type length;
  /// the number of elements in table_rows
  bimap_detail::rowidx_type capacity;
  /// max number of rows that must be probed to determine if a key is contained
  bimap_detail::rowidx_type max_probe;
  /// specifies ownership of keys, @see BiMapMode
  BiMapMode mode;

  /// actual hash table data
  bimap_StrU16_detail::Row* table_rows;
  /// tracks the row indices to make iteration faster
  bimap_detail::rowidx_type* ordered_row_indices;
};

/// Create an invalid FrozenKeyIdxBiMap
///
/// @note
/// ugh, it's unfortunate that we need to make this... but for now it's useful.
/// Ideally, we would refactor so that we can get rid of this function. A
/// useful compromise might simply put it within the bimap_detail namespace
inline FrozenKeyIdxBiMap mk_invalid_FrozenKeyIdxBiMap() {
  return FrozenKeyIdxBiMap{bimap_detail::INVALID_VAL,
                           bimap_detail::INVALID_VAL,
                           0,
                           BiMapMode::REFS_KEYDATA,
                           nullptr,
                           nullptr};
}

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
FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(const char* const keys[], int key_count,
                                        BiMapMode mode);

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
  return ptr->length != bimap_detail::INVALID_VAL;
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
  if (FrozenKeyIdxBiMap_is_ok(ptr)) {
    if (ptr->length > 0) {
      if (ptr->mode == BiMapMode::COPIES_KEYDATA) {
        for (bimap_detail::rowidx_type i = 0; i < ptr->capacity; i++) {
          bimap_StrU16_detail::Row* row = ptr->table_rows + i;
          // casting from (const char*) to (char*) should be legal (as long as
          // there were no bugs modifying the value of ptr->mode)
          if (row->keylen > 0) {
            delete[] row->key;
          }
        }
      }
      delete[] ptr->table_rows;
      delete[] ptr->ordered_row_indices;
    }  // ptr->length > 0
    (*ptr) = mk_invalid_FrozenKeyIdxBiMap();
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
FrozenKeyIdxBiMap FrozenKeyIdxBiMap_clone(const FrozenKeyIdxBiMap* ptr);

namespace bimap {

/// holds the result of a call to @ref FrozenKeyIdxBiMap_find
///
/// @note
/// This is a C-style approximation of std::optional<uint16_t>. Additionally,
/// the choice to make value a uint16_t is motivated by PR #484
struct AccessRslt {
  /// Indicates whether the value member is valid
  bool has_value;
  /// the loaded value (if has_value is false then this holds garbage)
  uint16_t value;
};

}  // namespace bimap

/// lookup the value associated with the key
///
/// This is the analog to calling `map[key]` in python.
///
/// @param[in] map A pointer to a valid bimap
/// @param[in] key A null-terminated string
///
/// @return An instance of @ref bimap::AccessRslt that encodes the value (if
///     the key is present)
inline bimap::AccessRslt FrozenKeyIdxBiMap_find(const FrozenKeyIdxBiMap* map,
                                                const char* key) {
  uint16_t tmp = bimap_StrU16_detail::search(map->table_rows, key,
                                             map->capacity, map->max_probe)
                     .val;
  return bimap::AccessRslt{tmp != bimap_detail::INVALID_VAL, tmp};
}

/// returns whether the map contains the key
///
/// @param[in] map A pointer to a valid bimap
/// @param[in] key A null-terminated string
inline bool FrozenKeyIdxBiMap_contains(const FrozenKeyIdxBiMap* map,
                                       const char* key) {
  return FrozenKeyIdxBiMap_find(map, key).has_value;
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
/// then a call to @ref FrozenKeyIdxBiMap_find that passes `s` will
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
inline const char* FrozenKeyIdxBiMap_inverse_find(const FrozenKeyIdxBiMap* map,
                                                  uint16_t idx) {
  if (idx >= map->length) {
    return nullptr;
  }
  const char* out = map->table_rows[map->ordered_row_indices[idx]].key;
  GR_INTERNAL_REQUIRE(out != nullptr, "logical error: string can't be nullptr");
  return out;
}

/** @}*/  // end of group

namespace bimap_detail {

/// a helper function used to actually allocate memory for FrozenKeyIdxBiMap
inline FrozenKeyIdxBiMap alloc(uint16_t length, uint16_t capacity,
                               BiMapMode mode) {
  // it would be nice to handle allocate all pointers as a single block of
  // memory, but that gets tricky. Essentially, we would allocate uninitialized
  // memory and manually use placement-new (and the corresponding `delete`)
  using bimap_detail::rowidx_type;
  using bimap_StrU16_detail::Row;
  FrozenKeyIdxBiMap out = {
      /*length=*/length,
      /*capacity=*/capacity,
      /*max_probe=*/capacity,
      /*mode=*/mode,
      /*table_rows=*/(capacity > 0) ? new Row[capacity] : nullptr,
      /*ordered_row_indices=*/(length > 0) ? new rowidx_type[length] : nullptr};
  for (uint16_t i = 0; i < capacity; i++) {
    out.table_rows[i].keylen = 0;
  }
  return out;
}

}  // namespace bimap_detail

inline FrozenKeyIdxBiMap new_FrozenKeyIdxBiMap(const char* const keys[],
                                               int key_count, BiMapMode mode) {
  int64_t max_len = static_cast<int64_t>(bimap_cap_detail::max_key_count());
  if (keys == nullptr && key_count == 0) {
    return bimap_detail::alloc(0, 0, mode);
  } else if (keys == nullptr) {
    GrPrintErrMsg("keys must not be a nullptr");
    return mk_invalid_FrozenKeyIdxBiMap();
  } else if (key_count < 1 || static_cast<int64_t>(key_count) > max_len) {
    GrPrintErrMsg("key_count must be positive & can't exceed %lld", max_len);
    return mk_invalid_FrozenKeyIdxBiMap();
  }

  // based on the preceding check, this shouldn't be able to fail
  bimap_detail::rowidx_type capacity =
      bimap_cap_detail::calc_map_capacity(key_count);
  GR_INTERNAL_REQUIRE(capacity > 0, "something went wrong");

  // let's validate the keys
  for (int i = 0; i < key_count; i++) {
    GR_INTERNAL_REQUIRE(keys[i] != nullptr, "Can't specify a nullptr key");
    std::size_t n_chrs_without_nul = std::strlen(keys[i]);
    if (n_chrs_without_nul == 0 ||
        n_chrs_without_nul > bimap_detail::KEYLEN_MAX) {
      GrPrintErrMsg(
          "calling strlen on \"%s\", the key @ index %d, yields 0 or a length "
          "exceeding %d",
          keys[i], i, bimap_detail::KEYLEN_MAX);
      return mk_invalid_FrozenKeyIdxBiMap();
    }
    // check uniqueness
    for (int j = 0; j < i; j++) {
      if (strcmp(keys[i], keys[j]) == 0) {
        GrPrintErrMsg("\"%s\" key repeats", keys[i]);
        return mk_invalid_FrozenKeyIdxBiMap();
      }
    }
  }

  // now, that we know we will succeed, lets construct the bimap
  FrozenKeyIdxBiMap out = bimap_detail::alloc(key_count, capacity, mode);

  // now it's time to fill in the array
  int max_probe_count = 1;
  for (int i = 0; i < key_count; i++) {
    // search for the first empty row
    bimap_StrU16_detail::SearchRslt search_rslt = bimap_StrU16_detail::search(
        out.table_rows, keys[i], capacity, capacity);
    // this should be infallible (especially after we already did some checks)
    GR_INTERNAL_REQUIRE(search_rslt.probe_count != 0, "sanity check failed");

    // now we overwrite the row
    bimap_StrU16_detail::overwrite_row(out.table_rows + search_rslt.rowidx,
                                       keys[i], std::strlen(keys[i]), i,
                                       mode == BiMapMode::COPIES_KEYDATA);
    out.ordered_row_indices[i] = search_rslt.rowidx;

    max_probe_count = std::max(max_probe_count, search_rslt.probe_count);
  }
  out.max_probe = max_probe_count;

  return out;
}

inline FrozenKeyIdxBiMap FrozenKeyIdxBiMap_clone(const FrozenKeyIdxBiMap* ptr) {
  FrozenKeyIdxBiMap out =
      bimap_detail::alloc(ptr->length, ptr->capacity, ptr->mode);
  out.max_probe = ptr->max_probe;

  if (ptr->length == 0 || !FrozenKeyIdxBiMap_is_ok(ptr)) {
    return out;
  }

  // give the compiler/linter a hint that out.table_rows is not a nullptr
  // (this is guaranteed by the preceding early exit)
  GR_INTERNAL_REQUIRE(
      (out.table_rows != nullptr) && (out.ordered_row_indices != nullptr),
      "something is very wrong!");

  bool copy_key_data = out.mode == BiMapMode::COPIES_KEYDATA;
  for (bimap_detail::rowidx_type i = 0; i < ptr->capacity; i++) {
    const bimap_StrU16_detail::Row& ref_row = ptr->table_rows[i];
    if (ref_row.keylen > 0) {
      bimap_StrU16_detail::overwrite_row(out.table_rows + i, ref_row.key,
                                         ref_row.keylen, ref_row.value,
                                         copy_key_data);
    }
  }

  for (bimap_detail::rowidx_type i = 0; i < ptr->length; i++) {
    out.ordered_row_indices[i] = ptr->ordered_row_indices[i];
  }
  return out;
};

}  // namespace grackle::impl

#endif  // SUPPORT_FROZENKEYIDXBIMAP_HPP
