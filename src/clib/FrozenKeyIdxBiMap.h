// See LICENSE file for license and copyright information

/// @file FrozenKeyIdxBiMap.h
/// @brief Declares the internal FrozenKeyIdxBiMap type
///
///
/// All insertions occur during initialization. This simplifies a lot of
/// bookkeeping.
///
/// If necessary, a number of optimizations could be made to the implementation
/// that might make key lookups faster. These optimizations could take
/// advantage of the following factors:
///    - Assumptions about the max key size and the max capacity of the map.
///      For example, if the max key size never exceeds ~22 characters and
///      there are never more than ~128 entries, it would probably be optimal
///      to store the strings in-place (improving cache locallity).
///    - Alternatively, making assumptions 

#ifndef FROZEN_KEY_IDX_BIMAP_HPP
#define FROZEN_KEY_IDX_BIMAP_HPP

#include <stddef.h>
#include <stdint.h>

#include "FrozenKeyIdxBiMap_details.h"

#ifdef __cplusplus
extern "C" {
#endif

/// the type used for indexing the rows of the BiMap's internal hash table
///
/// @note
/// This does not necessarily need to be a uint16_t (even if the internal
/// values of the hash table are uint16_t).
typedef uint16_t BiMap_rowidx_type;

/// @def BIMAP_ROWIDX_TYPE_MAX
/// @brief the maximum value representable by BiMap_rowidx_type
#define BIMAP_ROWIDX_TYPE_MAX UINT16_MAX


enum BiMapMode {
  BIMAP_REFS_KEYDATA = 0,
  BIMAP_COPIES_KEYDATA = 1
};

/// @brief This is a bidirectional map (bimap). It is specialized to map `n`
///     unique string keys to unique indexes with values of `0` thru `n-1` and 
///     vice-versa. The ordering of keys is set at initialization and frozen.
///
/// This is primarily intended to be used in the implementation of Maps of 
/// arrays (where the values could be part of a single contiguous array or are
/// individual arrays), but this but may be broadly useful for other
/// applications.
///
/// This operates in 2 modes:
/// 1. BIMAP_REFS_KEYDATA: This is the default, where we operate under the
///    assumption that the allocations holding the string characters outlive
///    the bimap. In this mode the bimap  is intended to hold string-literals
///    (which are live for the entirety of a program). This minimizes memory
///    usage.
/// 2. BIMAP_COPIES_KEYDATA: Under this mode, the bimap copies the data of all
///    of the keys. This is useful for testing purposes. In the long-term, if
///    we allow dynamic extension of chemistry networks, it will also be
///    useful. If we are implement the optimizations described down below
///    (where we directly embed the string in the hash-table-rows), this will
///    probably be a quite a bit faster
///
/// This is based on the logic I wrote for Enzo-E.
///   https://github.com/enzo-project/enzo-e/blob/main/src/Cello/view_StringIndRdOnlyMap.hpp
/// (More details are provided below under C++ considerations)
///
/// Why Frozen?
/// ===========
/// The contents are "frozen" for a 3 primary reasons:
/// - It drastically simplifies the implementation (we don't have to worry
///   about deletion -- which can be quite messy)
/// - Linear-probing generally provides better data locallity than other hash
///   collision resolution techniques, but generally has other drawbacks.
///   Freezing the contents let's us mitigate many drawbacks (mostly related to
///   the deletion operation)
/// - It could let us make copy operations cheaper. If we know the map won't
///   change, we could just use reference counting.
///
/// Improvements: Reference Counting
/// ================================
/// The original C++ leverages std::shared_ptr to achieve reference counting
/// (and reduce the cost of copying). Theoretically, I would like to see us use
/// some kind of reference-counting too. But this is tricky in library code,
/// given the diversity of threading libraries that are not formally
/// interoperable. I think the only way to properly do this would be to come up
/// with a system for allowing registration of locks/atomics with Grackle as a
/// whole.
/// 
/// C++ Considerations
/// ==================
/// As we starting using a C++ compiler, I would like us to eventually embrace
/// the class structure present in the original Enzo-E version (it would
/// greatly reduce the chance of memory leaks). But, for reasons expressed
/// above, I am concerned about using std::shared_ptr for reference counting.
///
/// I would be stunned if std::map<std::string, uint16_t> or
/// std::map<const char*, uint16_t> is faster than the internal hash table
/// since std::map is usually implemented as a tree.
///
/// Potential Improvements
/// ======================
/// There is definitely room for optimizing this implementation:
/// - We could be smarter about the order that we insert keys into the table
///   (in the constructor) to minimize the search time.
/// - We might be able to come up with a better hash function
/// - We can achieve even better locality, in BIMAP_COPIES_KEYDATA mode, thanks
///   to our use of STRU16MAP_KEYLEN_MAX (which is currently 29). Once we
///   shift everything to C++17 (it is possible now, but requires care), we
///   could define:
///   ```{.cpp}
///   struct alignas(32) packedrow_ { char data[32]; };
///   bool is_empty(packed_row_ r) { return data[0] == '\0' }
///   const char* get_key(packed_row_ r) { return r.data; }
///   uint16_t get_val(packed_row_ r)
///   { uint16_t o; memcpy(&o, r.data+30, 2); return o; }
///   ```
///   This is useful since it improves locality of the string. `alignas(32)` is
///   present to help ensure better cacheline alignment and with a little extra
///   care, it lets us use SIMD operations for faster probing
///
/// > [!note]
/// > The contents of this struct should be considered an implementation
/// > detail! Always prefer the associated functions (they are defined in such
/// > a way that they should be inlined
struct FrozenKeyIdxBiMap{
  // don't forget to update FrozenKeyIdxBiMap_clone when changing members

  /// the number of contained strings
  BiMap_rowidx_type length;
  /// the number of elements in table_rows
  BiMap_rowidx_type capacity;
  /// max number of rows that must be probed to determine if a key is contained
  BiMap_rowidx_type max_probe;
  /// indicates whether the map "owns" the memory holding the characters in
  /// each key or just references it
  enum BiMapMode mode;

  /// actual hash table data
  struct StrU16Row_* table_rows;
  /// tracks the row indices to make iteration faster
  BiMap_rowidx_type* ordered_row_indices;
};

typedef struct FrozenKeyIdxBiMap FrozenKeyIdxBiMap;

/// Constructs a new FrozenKeyIdxBiMap
///
/// @param[out] out Pointer where the allocated type is stored
/// @param[in]  keys Sequence of 1 or more unique strings. Each string must
///     include at least 1 non-null character and be null-terminated
/// @param[in]  key_count The length of keys
/// @param[in]  mode specifies handling of keys. This will be passed on to any
///     clones that are made.
int new_FrozenKeyIdxBiMap(
  FrozenKeyIdxBiMap** out, const char* keys[], int key_count,
  enum BiMapMode mode
);

/// Destroys the specified FrozenKeyIdxBiMap
void drop_FrozenKeyIdxBiMap(FrozenKeyIdxBiMap*);

/// Makes a clone of the specified FrozenKeyIdxBiMap (the clone inherites the
/// original BiMapMode).
int FrozenKeyIdxBiMap_clone(FrozenKeyIdxBiMap** out,
                            const FrozenKeyIdxBiMap* ptr);

/// returns the value associated with the key (or STRU16MAP_INVALID_VAL)
///
/// @note
/// denoting this as static is the only effective way to define inline in C
static inline uint16_t FrozenKeyIdxBiMap_idx_from_key(
  const FrozenKeyIdxBiMap* map, const char* key
)
{
  return StrU16Map_find_match_or_empty_(
    map->table_rows, key, map->capacity, map->max_probe
  ).val;
}

/// checks if the map contains a key
static inline int FrozenKeyIdxBiMap_contains(const FrozenKeyIdxBiMap* map,
                                             const char* key) {
  return FrozenKeyIdxBiMap_idx_from_key(map, key) != STRU16MAP_INVALID_VAL;
}

static inline int FrozenKeyIdxBiMap_size(const FrozenKeyIdxBiMap* map)
{ return map->length; }

/// Return the ith key (this is effectively a reverse lookup)
static inline const char* FrozenKeyIdxBiMap_key_from_idx(
  const FrozenKeyIdxBiMap* map, uint16_t i
) {
  if (i >= map->length) { return NULL; }
  const char* out = map->table_rows[map->ordered_row_indices[i]].key;
  GR_INTERNAL_REQUIRE(out != NULL, "The string can't be NULL - logical error");
  return out;
}



#ifdef __cplusplus
} // extern "C"
#endif

#endif /* FROZEN_KEY_IDX_BIMAP_HPP */
