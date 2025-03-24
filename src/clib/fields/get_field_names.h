// See LICENSE file for license and copyright information

/// @file get_field_names.h
/// @brief Declares functions to get the list of field names.

#ifndef GET_FIELD_NAMES_H
#define GET_FIELD_NAMES_H

#include "FrozenKeyIdxBiMap.h"
#include "grackle.h"

#ifdef __cplusplus
extern "C" {
#endif

/// an internal function to determine the list of required fields
///
/// @param[out] out Holds the field names
/// @param[in]  my_chemistry Specifies the configuration
/// @param[in]  my_rates Specifies additional configuration information
int get_field_names_from_local_(
  FrozenKeyIdxBiMap** out, const chemistry_data* my_chemistry
);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* GET_FIELD_NAMES_H */
