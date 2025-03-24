// See LICENSE file for license and copyright information

/// @file internal_field_data.h
/// @brief Declares the internal_field_data type

#ifndef INTERNAL_FIELD_DATA_H
#define INTERNAL_FIELD_DATA_H

#include "grackle.h"

#ifdef __cplusplus
extern "C" {
#endif

/// defines the the internal_field_data type.
/// - The intention is to use this type internally  the transcribed Grackle
///   routines.
/// - This type is **only** directly accessible through the grField api
/// - This type can also be internally constructed from grackle_field_data
/// - In a future major release, we may remove grackle_field_data
///
/// @important
/// The contents of this type will **NEVER** be exposed as part of the public
/// API. This provides flexibility for us to refactor this type without breaking
/// the API.
///
/// @note
/// Accordingly, the fact that we choose to define internal_field_data as an
/// alias of the grackle_field_data type is just a temporary convience,
/// reflecting the fact that grackle_field_data and internal_field_data are
/// presently identical. This won't always be the case.
///
/// In the future, internal_field_data will be an independently defined struct
/// (if we wanted, it could even be a C++ class).
/// - internal_field_data will provides access to all of the information
///   currently tracked by grackle_field_data.
/// - initially, the way data is organized within internal_field_data will
///   remain similar to that of grackle_field_data, but that will change with
///   time.
/// - as we add new fields to grackle_field_data, they will always be added to
///   the internal_field_data type (as a general rule, they won't be added to
///   grackle_field_data, but exceptions can be made...)
typedef grackle_field_data internal_field_data;

/// construct an instance of internal_field_data from grackle_field_data
///
/// @note
/// Since these 2 types are currently aliases, the implementation is currently
/// extremely trivial. This will change in the future
static inline internal_field_data make_internal_from_grackle_field_data_(
  grackle_field_data* my_fields
)
{
  internal_field_data out = *my_fields;
  return out;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // INTERNAL_FIELD_DATA_H
