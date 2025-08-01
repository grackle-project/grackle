// See LICENSE file for license and copyright information

/// @file FrozenKeyIdxBiMap.c
/// @brief Implements parts of FrozenKeyIdxBiMap type

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "FrozenKeyIdxBiMap.h"
#include "status_reporting.h"
#include "grackle.h" // GR_SUCCESS

/// @def    INVERSE_LOAD_FACTOR
/// @brief  the load factor specifies the fraction of the capcity of the Hash
///         table that is filled. This should be an integer.
///
/// This should not be 1. Generally, the larger this is, the fewer collisions
/// there are, but the more memory is required. If this is too large, it may
/// actually slow down lookups. 
#define INVERSE_LOAD_FACTOR 2

// this is the list of allowed capacities. These are all prime numbers that
// have nearly constant linear spacing (it may make more sense to have
// logarithmic spacing)
static const uint32_t allowed_capacities_l_[] = {
  7,  19,  31,  41, 53, 61,  71,  83, 89, 101, 113, 127, 139, 149, 163, 173,
  181, 191, 199, 211, 223, 233, 241, 251
};

static const int allowed_capacities_len_ = (
  sizeof(allowed_capacities_l_) / sizeof(uint32_t)
);

/// compute the maximum number of keys
static inline BiMap_rowidx_type get_max_key_count_(void) {
  uint16_t max_from_allowed_capacity = (
    allowed_capacities_l_[allowed_capacities_len_ - 1] / INVERSE_LOAD_FACTOR
  );

  // while the following may seem silly, it is important to keep it here
  if (BIMAP_ROWIDX_TYPE_MAX < max_from_allowed_capacity) {
    return BIMAP_ROWIDX_TYPE_MAX;
  } else {
    return max_from_allowed_capacity;
  }
}

/// compute the capacity of the map
///
/// @param num_keys the desired number of keys (should be positive)
///
/// @returns the capacity. A value of 0 means that there is not enough space
static inline uint16_t calc_map_capacity_(int key_count) {
  uint64_t target_capacity = INVERSE_LOAD_FACTOR * (uint64_t)key_count;
  // binary search may be faster
  for (int i = 0; i < allowed_capacities_len_; i++) {
    if (target_capacity < allowed_capacities_l_[i]) {
      return (uint16_t)allowed_capacities_l_[i];
    }
  }
  return (uint16_t)0;
}

static int allocate_FrozenKeyIdxBiMap_(FrozenKeyIdxBiMap** out,
                                       uint16_t length,
                                       uint16_t capacity) {
  // when we shift to C++, it would be nice to perform one call to malloc (so
  // that the memory is all contiguous) and then use placement new
  *out = (FrozenKeyIdxBiMap*)malloc(sizeof(FrozenKeyIdxBiMap));
  (*out)->length = length;
  (*out)->capacity = capacity;
  // we use calloc so keylen is always 0
  (*out)->table_rows = (struct StrU16Row_*)calloc(
    capacity, sizeof(struct StrU16Row_)
  );
  (*out)->ordered_row_indices = (uint16_t*)malloc(sizeof(uint16_t) * length);
  // most platforms don't let us gracefully handle failure when it comes to
  // allocations, so we don't bother...
  return GR_SUCCESS;
}

static void StrU16Row_overwrite_(
  struct StrU16Row_* row, const char* key, size_t keylen, uint16_t value,
  enum BiMapMode mode
) {
  GR_INTERNAL_REQUIRE(row->keylen == 0, "Sanity check failed!");
  row->value = value;
  row->keylen = keylen;
  if (mode == BIMAP_COPIES_KEYDATA) {
    size_t total_len = keylen + 1; // <- add 1 to account for null-terminator
    char* ptr = (char*)malloc(sizeof(char) * total_len);
    memcpy(ptr, key, total_len);
    row->key = ptr;
  } else {
    row->key = key;
  }
}

// do some checks on keys and get keylen
int check_key_get_len_(const char* keys[], int i, size_t* keylen) {
  GR_INTERNAL_REQUIRE(keys[i] != NULL, "Can't specify a nullptr key");
  *keylen = strlen(keys[i]); // <- excludes the null terminator from len
  if ((*keylen == 0) || (*keylen > STRU16MAP_KEYLEN_MAX)) {
    return GrPrintAndReturnErr(
      "calling strlen on \"%s\", the key @ index %d, yields 0 or a length "
      "exceeding %d",
      keys[i], i, STRU16MAP_KEYLEN_MAX
    );
  }
  return GR_SUCCESS;
}

#ifdef __cplusplus
extern "C" {
#endif

int new_FrozenKeyIdxBiMap(
  FrozenKeyIdxBiMap** out, const char* keys[], int key_count,
  enum BiMapMode mode
) {
  if ((mode != BIMAP_REFS_KEYDATA) && (mode != BIMAP_COPIES_KEYDATA)) {
    return GrPrintAndReturnErr(
      "mode must be BIMAP_REFS_KEYDATA or BIMAP_COPIES_KEYDATA"
    );
  } else if ( (key_count < 1) || (key_count > get_max_key_count_()) ) {
    return GrPrintAndReturnErr(
      "key_count must be positive and cannot exceed %lld",
      (long long)get_max_key_count_()
    );
  }

  uint16_t capacity = calc_map_capacity_(key_count);
  // the following shouldn't be possible based on the prior check
  GR_INTERNAL_REQUIRE(capacity > 0, "something went wrong");

  int ec = allocate_FrozenKeyIdxBiMap_(out, key_count, capacity);
  if ((*out) == NULL) { return GrPrintAndReturnErr("something went wrong"); }

  (*out)->mode = mode;
  
  int max_probe_count = 1;

  for (int i = 0; i < key_count; i++) {
    size_t keylen;
    int tmp = check_key_get_len_(keys, i, &keylen);
    if (tmp != GR_SUCCESS) {
      drop_FrozenKeyIdxBiMap(*out);
      *out = NULL;
      return tmp;
    }

    // insert the row
    struct StrU16Search_ search_rslt = StrU16Map_find_match_or_empty_(
      (*out)->table_rows, keys[i], capacity, capacity
    );
    StrU16Row_overwrite_((*out)->table_rows + search_rslt.rowidx, keys[i],
                         keylen, i, mode);
    (*out)->ordered_row_indices[i] = search_rslt.rowidx;

    if (max_probe_count < search_rslt.probe_count) {
      max_probe_count = search_rslt.probe_count;
    }
  }
  (*out)->max_probe = max_probe_count;

  return GR_SUCCESS;
}

int FrozenKeyIdxBiMap_clone(FrozenKeyIdxBiMap** out,
                            const FrozenKeyIdxBiMap* ptr) {
  int ec = allocate_FrozenKeyIdxBiMap_(out, ptr->length, ptr->capacity);
  if (ec != GR_SUCCESS) { return ec; }

  (*out)->length = ptr->length; // (technically this is also redundant)
  (*out)->capacity = ptr->capacity; // (technically this is also redundant)
  (*out)->max_probe = ptr->max_probe;
  (*out)->mode = ptr->mode;
  
  if (!ptr->mode) {
    memcpy((*out)->table_rows, ptr->table_rows,
           sizeof(struct StrU16Row_)*ptr->capacity);
  } else {
    const int capacity = ptr->capacity; // this may help compiler
    for (int i = 0; i < capacity; i++) {
      struct StrU16Row_ ref_row = ptr->table_rows[i];
      if (ref_row.keylen > 0) {
        StrU16Row_overwrite_((*out)->table_rows + i, ref_row.key,
                             ref_row.keylen, ref_row.value, 1);
      } else {
        (*out)->table_rows[i] = ref_row; // technically unnecessary
      }
    }
  }

  memcpy((*out)->ordered_row_indices, ptr->ordered_row_indices,
         sizeof(BiMap_rowidx_type) * ptr->capacity);
  return GR_SUCCESS;
}

void drop_FrozenKeyIdxBiMap(FrozenKeyIdxBiMap* ptr) {
  if (ptr->mode) {
    BiMap_rowidx_type capacity = ptr->capacity;
    for (BiMap_rowidx_type i = 0; i < capacity; i++) {
      struct StrU16Row_* row = ptr->table_rows+i;
      // casting from (const char*) to (char*) should be legal (as long as
      // there were no bugs modifying the value of ptr->mode)
      if (row->keylen > 0) { free((char*)row->key); }
    }
  }
  free(ptr->table_rows);
  free(ptr->ordered_row_indices);
  free(ptr);
}

#ifdef __cplusplus
} // extern "C"
#endif
