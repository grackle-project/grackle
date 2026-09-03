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
#include <optional>
#include <utility>  // std::swap

#include "config.hpp"
#include "support/status_reporting.hpp"
#include "FrozenKeyIdxBiMap_detail.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// the following doxygen comment block logically groups every all parts of
// the (internal) API for Grackle's (internal) FrozenKeyIdxBiMap. It's useful
// when generating a doxygen webpage

/// @defgroup bimap-grp FrozenKeyIdxBiMap Data Type
///
/// FrozenKeyIdxBiMap provides specialized mapping functionality for internal
/// use within Grackle. The functionality is useful as a building-block for
/// runtime lookup-tables and other data types with map-like interface.
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
class FrozenKeyIdxBiMap {
  // don't forget to update the clone method when changing data members

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

  /// a helper function used to actually allocate memory for FrozenKeyIdxBiMap
  void alloc_(uint16_t target_length, uint16_t target_capacity,
              BiMapMode target_mode) {
    // it would be nice to handle allocate all pointers as a single block of
    // memory, but that gets tricky. Essentially, we would allocate
    // uninitialized memory and manually use placement-new (and the
    // corresponding `delete`)
    using bimap_detail::rowidx_type;
    using bimap_StrU16_detail::Row;
    length = target_length;
    capacity = target_capacity;
    max_probe = target_capacity;
    mode = target_mode;
    table_rows = (target_capacity > 0) ? new Row[target_capacity] : nullptr;
    ordered_row_indices =
        (target_length > 0) ? new rowidx_type[target_length] : nullptr;
    for (uint16_t i = 0; i < target_capacity; i++) {
      table_rows[i].keylen = 0;
    }
  }

  /// @brief Create an invalid FrozenKeyIdxBiMap
  ///
  /// @note
  /// ugh, it's unfortunate that we need to make this... but for now it's
  /// useful. Ideally, we will get rid of this function.
  static FrozenKeyIdxBiMap make_invalid_() {
    FrozenKeyIdxBiMap out;
    out.length = bimap_detail::INVALID_VAL;
    out.capacity = bimap_detail::INVALID_VAL;
    out.max_probe = 0;
    out.mode = BiMapMode::REFS_KEYDATA;
    out.table_rows = nullptr;
    out.ordered_row_indices = nullptr;
    return out;
  }

public:  // interface methods
  /// @brief Factory Method (constructs a new FrozenKeyIdxBiMap)
  ///
  /// @param[in]  keys Sequence of 1 or more unique strings. Each string must
  ///     include at least 1 non-null character and be null-terminated
  /// @param[in]  key_count The length of keys
  /// @param[in]  mode specifies handling of keys. This will be passed on to any
  ///     clones that are made.
  ///
  /// @note
  /// Callers should pass the returned value to @ref FrozenKeyIdxBiMap::is_ok
  /// to check whether there was an error during creation. This is pretty
  /// ugly/clunky, but it's the only practical way to achieve comparable
  /// behavior to other internal data types. The best alternatives involve the
  /// use of std::optional or C++23's std::expected.
  static FrozenKeyIdxBiMap create(const char* const keys[], int key_count,
                                  BiMapMode mode);

  /// @brief Default Constructor
  ///
  /// Returns an instance holding 0 elements.
  FrozenKeyIdxBiMap() { alloc_(0, 0, BiMapMode::REFS_KEYDATA); }

  // for now, lets disble copy construction and assignment (its usually a
  // mistake when that happens and this is important to transitioning towards
  // acting more like a class)
  FrozenKeyIdxBiMap(const FrozenKeyIdxBiMap&) = delete;
  FrozenKeyIdxBiMap& operator=(const FrozenKeyIdxBiMap&) = delete;

  /// @brief Move Constructor
  ///
  /// Constructs a new instance and transfers the contents from @p other into
  /// the new instance. @p other is left in an undefined state.
  ///
  /// @param other The source of contents for the newly constructed instance
  FrozenKeyIdxBiMap(FrozenKeyIdxBiMap&& other) noexcept : FrozenKeyIdxBiMap() {
    swap(other);
  }

  /// @brief Move Assignment
  ///
  /// Transfers the contents from @p other into `this`. @p other is left in an
  /// undefined state.
  ///
  /// @param other The source of contents for the newly constructed instance
  FrozenKeyIdxBiMap& operator=(FrozenKeyIdxBiMap&& other) {
    if (this != &other) {
      swap(other);
    }
    return *this;
  }

  /// @brief Destuctor
  inline ~FrozenKeyIdxBiMap() noexcept;

  /// @brief Makes a clone of this
  ///
  /// The clone inherits the original's BiMapMode value. If it held
  /// BiMapMode::COPIES_KEYDATA, then fresh copies of the strings are made
  ///
  /// @warning
  /// Callers should pass the returned value to @ref FrozenKeyIdxBiMap::is_ok
  /// to check whether there was an error during creation. This is pretty
  /// ugly/clunky, but it's the only practical way to achieve comparable
  /// behavior to other internal data types. The best alternatives involve
  /// things like std::optional or C++23's std::expected.
  ///
  /// @note
  /// If we wanted slightly more idiomatic C++, we would fold this into the
  /// copy constructor and copy assignment methods.
  FrozenKeyIdxBiMap clone() const;

  /// @brief swaps contents
  void swap(FrozenKeyIdxBiMap& other) noexcept {
    std::swap(length, other.length);
    std::swap(capacity, other.capacity);
    std::swap(max_probe, other.max_probe);
    std::swap(mode, other.mode);
    std::swap(table_rows, other.table_rows);
    std::swap(ordered_row_indices, other.ordered_row_indices);
  }

  /// @brief checks whether a creational function produced a valid container
  ///
  /// @return true if the container is valid (otherwise, it returns false)
  bool is_ok() const noexcept { return length != bimap_detail::INVALID_VAL; }

  /// @brief lookup the value associated with the specified key
  ///
  /// This is the analog to calling `map[key]` in python.
  ///
  /// @param[in] key A null-terminated string
  /// @return An optional that contains the value if the key can be found
  std::optional<uint16_t> find(const char* key) const noexcept {
    uint16_t tmp =
        bimap_StrU16_detail::search(table_rows, key, capacity, max_probe).val;
    bool success = tmp != bimap_detail::INVALID_VAL;
    return (success) ? std::make_optional(tmp) : std::nullopt;
  }

  /// @brief Return the key associated with the specified value
  ///
  /// For some context, if this function returns a string `s` for some index
  /// `i`, then a call to @ref FrozenKeyIdxBiMap::find that passes `s` will
  /// return `i`
  ///
  /// This is intended for use in situations where you briefly need the string
  /// (i.e. and you plan to stop using the pointer before or at the same time as
  /// `this` is destroyed). In more detail:
  /// - If `this` was constructed in @ref BiMapMode::COPIES_KEYDATA mode,
  ///   returned strings have the same lifetime as `this` (i.e. they are
  ///   deallocated when the contents of `this` are deallocated).
  /// - Otherwise, the returned string's allocation is externally managed. But,
  ///   any scenario where the allocation doesn't live at least as long as
  ///   `this`, is ill-formed
  ///
  /// @param[in] idx The index to check
  /// @return The pointer to the appropriate key
  const char* inverse_find(int idx) const noexcept {
    if (idx >= length || idx < 0) {
      return nullptr;
    }
    const char* out = table_rows[ordered_row_indices[idx]].key;
    GR_INTERNAL_REQUIRE(out != nullptr,
                        "logical error: string can't be nullptr");
    return out;
  }

  /// @brief return the number of keys in the map
  int size() const noexcept { return length; }
};

/** @}*/  // end of group

inline FrozenKeyIdxBiMap::~FrozenKeyIdxBiMap() noexcept {
  if (is_ok()) {
    if (length > 0) {
      if (mode == BiMapMode::COPIES_KEYDATA) {
        for (bimap_detail::rowidx_type i = 0; i < capacity; i++) {
          bimap_StrU16_detail::Row* row = table_rows + i;
          // casting from (const char*) to (char*) should be legal (as long as
          // there were no bugs modifying the value of ptr->mode)
          if (row->keylen > 0) {
            delete[] row->key;
          }
        }
      }
      delete[] table_rows;
      delete[] ordered_row_indices;
    }  // ptr->length > 0
  }
}

inline FrozenKeyIdxBiMap FrozenKeyIdxBiMap::create(const char* const keys[],
                                                   int key_count,
                                                   BiMapMode mode) {
  int64_t max_len = static_cast<int64_t>(bimap_cap_detail::max_key_count());
  if (keys == nullptr && key_count == 0) {
    FrozenKeyIdxBiMap out;
    out.alloc_(0, 0, mode);
    return out;
  } else if (keys == nullptr) {
    GrPrintErrMsg("keys must not be a nullptr");
    return FrozenKeyIdxBiMap::make_invalid_();
  } else if (key_count < 1 || static_cast<int64_t>(key_count) > max_len) {
    GrPrintErrMsg("key_count must be positive & can't exceed %lld",
                  static_cast<long long int>(max_len));
    return FrozenKeyIdxBiMap::make_invalid_();
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
      return FrozenKeyIdxBiMap::make_invalid_();
    }
    // check uniqueness
    for (int j = 0; j < i; j++) {
      if (strcmp(keys[i], keys[j]) == 0) {
        GrPrintErrMsg("\"%s\" key repeats", keys[i]);
        return FrozenKeyIdxBiMap::make_invalid_();
      }
    }
  }

  // now, that we know we will succeed, lets construct the bimap
  FrozenKeyIdxBiMap out;
  out.alloc_(key_count, capacity, mode);

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

inline FrozenKeyIdxBiMap FrozenKeyIdxBiMap::clone() const {
  if (!is_ok()) {
    return FrozenKeyIdxBiMap::make_invalid_();
  }

  FrozenKeyIdxBiMap out;
  out.alloc_(length, capacity, mode);
  out.max_probe = max_probe;

  if (length == 0) {
    return out;
  }

  // give the compiler/linter a hint that out.table_rows is not a nullptr
  // (this is guaranteed by the preceding early exit)
  GR_INTERNAL_REQUIRE(
      (out.table_rows != nullptr) && (out.ordered_row_indices != nullptr),
      "something is very wrong!");

  bool copy_key_data = out.mode == BiMapMode::COPIES_KEYDATA;
  for (bimap_detail::rowidx_type i = 0; i < capacity; i++) {
    const bimap_StrU16_detail::Row& ref_row = table_rows[i];
    if (ref_row.keylen > 0) {
      bimap_StrU16_detail::overwrite_row(out.table_rows + i, ref_row.key,
                                         ref_row.keylen, ref_row.value,
                                         copy_key_data);
    }
  }

  for (bimap_detail::rowidx_type i = 0; i < length; i++) {
    out.ordered_row_indices[i] = ordered_row_indices[i];
  }
  return out;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_FROZENKEYIDXBIMAP_HPP
