// See LICENSE file for license and copyright information

/// @file gr_fields_private.h
/// @brief Declares private details related to the gr_fields type

#ifndef GR_FIELDS_PRIVATE_H
#define GR_FIELDS_PRIVATE_H

#include "grackle.h"
#include "fields/experimental_api.h"
#include "fields/internal_field_data.h"
#include "FrozenKeyIdxBiMap.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_GRID_RANK 3

/// Tracks the underlying arrays used for grid properties (so that we don't need
/// to create heap allocations for each of them)
struct GridPropData{
  /// specifies the length of each dimension
  int dimension[MAX_GRID_RANK];

  /// specifies the min valid (0-based) index along each dimension
  int start[MAX_GRID_RANK];

  /// specifies the max valid (0-based) index along each dimension.
  ///
  /// In the future, it would be more idiomatic if this tracked the minimum
  /// invalid, non-negative index along each dimension (in this case, we should
  /// rename this `stop`)
  int end[MAX_GRID_RANK];
};
typedef struct GridPropData GridPropData;


/// aggregate properties about the managed field data allocations (if gr_fields
/// was allocated in such a manner) 
struct FieldDataAllocProps {
  gr_float* ptr; // <-- NULL if allocations are unmanaged
  uint64_t elem_per_field; // <- 0 if allocations are unmanaged
};
typedef struct FieldDataAllocProps FieldDataAllocProps;


/// This struct is used to track all of the necessary information to map between
/// the grFields_<fn> API and the internal_field_data type
struct gr_fields_interface_ {
  /// associates field-names with indexes (the indexes are essentially
  FrozenKeyIdxBiMap* field_names;

  /// tracks the field slots that can be filled in. Each entry corresponds to an
  /// entry in field_data
  gr_float*** field_slots;

  /// when gr_fields is configured to manage
  FieldDataAllocProps field_data_alloc_props;

  /// tracks the storage of grid_data.
  ///
  /// In the future, when the `internal_field_data` is defined separately from
  /// `grackle_field_data`, this member should be part of internal_field_data
  GridPropData grid_props;

  /// the underlying field-data representation used within core grackle routines
  internal_field_data field_data;
};

/// internal method to initialize gr_fields (it is wrapped by the public API)
///
/// @param[out] out holds the created pointer
/// @param[in]  field_names specifies the field names that are to be used. If
///     this function is successful, this entity is "consumed" (see below)
/// @param[in]  flags Bitwise ORed creation flags. If no flags are used, this
///     should be 0
/// @param[in]  nelements_per_field When the GR_FCREATE_MANAGED_FIELDDATA flag
///     is provided, this argument specifies the number of elements allocated
///     per field. Otherwise, this should be zero
///
/// On success, ``field_names`` is "consumed" by the function and the resulting
/// gr_fields instance takes ownership of the pointer (i.e. the caller should
/// never touch it again & deallocation is managed by the resulting instance).
/// On failure, the caller is responsible for cleaning up ``field_names``
///
/// @note
/// The interface could be improved. While we currently just support a single
/// flag, we could imagine a future where we use flags to specify the
/// memory-space or a specialized implementation. But maybe we should just
/// simplify the implementation for now?
int grFields_from_field_names_(gr_fields** out, FrozenKeyIdxBiMap* field_names,
                               uint64_t flags, uint64_t nelements_per_field);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* GR_FIELDS_PRIVATE_H */
