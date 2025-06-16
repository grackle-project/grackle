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

/// a portable version of strdup, which is provided on posix and in C23
char* my_strdup_(const char* src);

/// For a string ``s`` that starts with prefix ``prefix``, this returns
/// the first character in ``s`` after the prefix. Otherwise, it returns NULL.
///
/// If the returned non-NULL ptr points to a '\0' character, then both strings
/// are identical.
///
/// @param s the full string that may begin with the prefix
/// @param prefix the prefix that the full string may begin with
///
/// @return ``NULL`` either argument was ``NULL`` or if ``path`` does not
///   start with ``prefix``. Otherwise, this returns ``path + strlen(prefix)``
const char* post_prefix_ptr_(const char* s, const char* prefix);

/// join together fragments of a string into 1 newly allocated string
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
