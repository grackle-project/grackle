/***********************************************************************
/
/ Declare utility functions used internally by Grackle to encapsulate
/ logic for determining data files
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef DATA_FILE_UTILS_H
#define DATA_FILE_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/// temporarily declared for testing purposes
///
/// (in the future, we probably don't want tests to directly use this)
char* calc_checksum_str_(const char* fname);

/// used as the return type when determining the self-shielding location
///
/// if ``path`` is ``NULL``, then there is an error. This struct should NEVER
/// be exposed as part of the public API
struct generic_file_props {
  const char* path;
  int path_requires_dealloc;
  const char* checksum;
  int checksum_requires_dealloc;
};


/// Determines the path to the data file.
///
/// @param[in]  grackle_data_file specified grackle data file
/// @param[in]  grackle_data_file_options specifies how to interpret the first
///     argument
/// @param[in]  registry When testing, this can be used to specify a file
///     registry. In production, this should always be NULL.
///
/// @note
/// If this functionality ever gets exposed as part of the public API, we
/// should:
/// - stop using generic_data_file_props as a return type.
/// - make it possible for the caller to pre-allocate any buffers to hold the
///   file path and the computed checksum (in that case, we should consider
///   adopting an interface sorta like snprintf)
/// - want to remove registry as an argument
struct generic_file_props determine_data_file_(const char* grackle_data_file,
                                               int grackle_data_file_options,
                                               const char** registry);

/// Deallocates the memory held within a given ``struct generic_file_props``
void free_generic_file_props_(struct generic_file_props* ptr);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif /* DATA_FILE_UTILS_H */
