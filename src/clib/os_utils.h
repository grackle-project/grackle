/***********************************************************************
/
/ Declare utility functions used internally by Grackle to related to
/ path manipulation and OS-specific functionality
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#ifndef OS_UTILS_H
#define OS_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/// join together fragments of a string into 1 newly allocated string
///
/// @note
/// in principle, we could give sep == '\0' special significance
char* join_parts_(char sep, const char** parts, int nparts);

/// represents the known platform types (that produce different results)
enum platform_kind { platform_kind_generic_unix, platform_kind_unknown };

/// function that returns the appropriate platform enum
enum platform_kind get_platform_(void);

/// get the Grackle data directory
char* get_data_dir_(enum platform_kind kind);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif /* OS_UTILS_H */
