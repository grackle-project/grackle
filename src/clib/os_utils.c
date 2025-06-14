/***********************************************************************
/
/ Implement utility functions used internally by Grackle to related to
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

#include <errno.h>  // ERANGE
#include <stdio.h>  // fprintf, stderr
#include <stdlib.h>  // malloc, realloc, free, abort
#include <string.h>  // memcpy, strlen

#include "os_utils.h"
#include "grackle.h" // grackle_verbose

/// Just like getenv, except it returns NULL in place of strings of length 0.
static const char* getenv_nonempty_(const char* name) {
  const char* out = getenv(name);
  return ((out == NULL) || (out[0] == '\0')) ? NULL : out;
}

char* my_strdup_(const char* src) {
  size_t len_with_nul = strlen(src) + 1;
  char* out = malloc(sizeof(char) * len_with_nul);
  memcpy(out, src, len_with_nul);
  return out;
}

const char* post_prefix_ptr_(const char* s, const char* prefix) {
  if ((s == NULL) || (prefix == NULL)) return NULL;

  // these lengths don't include the null terminator
  size_t len_s = strlen(s);
  size_t len_prefix = strlen(prefix);
  if ((len_s < len_prefix) || (len_prefix == 0)) return NULL;

  if (memcmp(s, prefix, len_prefix) != 0) return NULL;

  return s + len_prefix;
}

char* join_parts_(char sep, const char** parts, int nparts) {
  if (nparts < 2) return NULL;

  // in principle, we could give sep == '\0' special significance

  size_t total_len = 0;
  for (int i = 0; i < nparts; i++) {
    if (parts[i] == NULL) return NULL;
    total_len += strlen(parts[i]); // we don't include the nul-terminator
  }
  total_len += (nparts - 1); // account for the size of sep and
  total_len++; // account for trailing nul-terminator

  char* out = malloc(total_len);
  size_t cur_offset = 0;
  for (int i = 0; i < nparts; i++) {
    if (i > 0) {
      out[cur_offset] = sep;
      cur_offset++;
    }
    size_t cur_part_len = strlen(parts[i]);
    memcpy(out + cur_offset, parts[i], cur_part_len);
    cur_offset += cur_part_len;
  }
  out[cur_offset] = '\0';
  if ((cur_offset+1) != total_len) abort();
  return out;
}


// Platform-Specific Stuff
// -----------------------

enum platform_kind get_platform_(void) {
#if defined(PLATFORM_GENERIC_UNIX) && defined(PLATFORM_MACOS)
  #error "more than 1 platform macro was defined"
#elif defined(PLATFORM_GENERIC_UNIX)
  return platform_kind_generic_unix;
#elif defined(PLATFORM_MACOS)
  return platform_kind_macos;
#else
  return platform_kind_unknown;
#endif
}

// define a function to get the home directory

/// returns the user's home directory
///
/// If it is defined with a non-empty value, the function honors the value in
/// the ``HOME`` environment variable. Otherwise, the function falls back to
/// fetching the value using platform specific apis.
///
/// @return a string pointing to the current user's home directory. ``NULL`` is
///     returned if there was an error. The caller is always responsible for
///     deallocating this string.
static char* get_home_dir(void);

#if defined(PLATFORM_GENERIC_UNIX) || defined(PLATFORM_MACOS)

// assume a posix-platform, the following headers are all standard 


#include <sys/types.h> // uid_t
#include <unistd.h> // getuid, sysconf
#include <pwd.h> // getpwuid, struct passwd

static char* get_home_dir(void)
{
  // first, try to get the value set in the environment
  const char* env_str = getenv_nonempty_("HOME");
  if (env_str != NULL) return my_strdup_(env_str);

  // fall back to checking the user database (standard on posix systems)

  // ask the system for an upper limit on the buffersize to hold the results
  const long initial_bufsize_guess = sysconf(_SC_GETPW_R_SIZE_MAX);

  // If the system can't give a firm answer, we guess.
  long bufsize = (initial_bufsize_guess == -1) ? 2048 : initial_bufsize_guess;
  char* buffer = NULL;

  struct passwd pwd, *result;
  int return_code;

  do {
    if (buffer == NULL) { // our 1st attempt
      buffer = malloc(sizeof(char)*bufsize);
    } else { // our next attempt
      bufsize *= 2;
      char* tmp = realloc(buffer, sizeof(char)*bufsize);
      if (tmp == NULL) break;
      buffer = tmp;
    }
    return_code = getpwuid_r(getuid(), &pwd, buffer, bufsize, &result);
  } while ((return_code == ERANGE) && (bufsize < 1000000000));

  if (return_code != 0) {
    free(buffer);
    fprintf(stderr, "ERROR while determining the HOME directory\n");
    return NULL;
  }

  char* out = my_strdup_(pwd.pw_dir);
  free(buffer);
  return out;
}
#else
static char* get_home_dir(void) {
  fprintf(stderr,
          "Don't know how to determine HOME directory on current platform\n");
  return NULL;
}
#endif

/// Returns a string specifying the default data directory
///
/// All of these choices are inspired by the API description of the
/// platformdirs python package
/// * we only looked at online documentation:
///   https://platformdirs.readthedocs.io/en/latest/
/// * we have NOT read any source code
static char* default_data_dir_(enum platform_kind kind) {
  const char* appname = "grackle";
  switch(kind) {

    case platform_kind_unknown: {
      fprintf(
        stderr,
        ("ERROR: can't infer default data dir on unknown platform.\n"
         " -> can only infer data directories on macOS and unix systems\n")
      );
      return NULL;
    }

    case platform_kind_macos: {
      // https://developer.apple.com/library/archive/documentation/FileManagement/Conceptual/FileSystemProgrammingGuide/MacOSXDirectories/MacOSXDirectories.html
      char* home_dir = get_home_dir();
      const char * parts[3] = {
        home_dir, "Library/Application Support", appname
      };
      char* out = join_parts_('/', parts, 3);
      free(home_dir);
      return out;
    }

    case platform_kind_generic_unix: {
      // https://specifications.freedesktop.org/basedir-spec/latest/
      const char* env_str = getenv_nonempty_("XDG_DATA_HOME");

      // check if we need to fall back to the default
      const char* dflt = "~/.local/share";
      if (env_str == NULL) {
        env_str = dflt;
      } else if ((env_str[0] != '~') && (env_str[0] != '/')) {
        // this is what the specification tells us to do
        fprintf(stderr,
                "WARNING: ignoring XDG_DATA_HOME because it doesn't hold an "
                "absolute path\n");
        env_str = dflt;
      }

      // now actually infer the absolute path
      if (env_str[0] == '~') {
        if (post_prefix_ptr_(env_str, "~/") == NULL) {
          fprintf(stderr,
                  "ERROR: can't expand env-variable, XDG_DATA_HOME when it "
                  "starts with `~user/` or just contains `~`\n");
          return NULL;
        }

        char* home_dir = get_home_dir();
        const char* parts[3] = {home_dir, env_str + 1, appname};
        char* out = join_parts_('/', parts, 3);
        free(home_dir);
        return out;

      } else {
        const char* parts[2] = {env_str, appname};
        char* out = join_parts_('/', parts, 2);
        return out;

      }
    }
  
  }

  fprintf(stderr,
          "ERROR: This part of the function should be unreachable! Did you add "
          "a new platform_kind and forget to update the function?\n");
  abort();
}

char* get_data_dir_(enum platform_kind kind) {
  const char* env_str = getenv_nonempty_("GRACKLE_DATA_DIR");
  char* out;
  const char* description;
  if (env_str != NULL) {
    out = my_strdup_(env_str);
    description = "from the `GRACKLE_DATA_DIR` environment variable";
  } else {
    if (grackle_verbose) {
      fprintf(stdout,
              ("INFO: looking up system-default for the data directory since "
               "`GRACKLE_DATA_DIR` env variable is empty\n"));
      fflush(stdout); // flush in case we run into an error in the next call
    }
    out = default_data_dir_(kind);
    description = "inferred from the system defaults";
  }

  // confirm we are providing an absolute path
  if (out[0] != '/') {
    fprintf(stderr,
            "ERROR: the data-directory %s, `%s` is not an absolute path\n",
            description, out);
    free(out);
    return out;
  }
  if (grackle_verbose) {
    fprintf(stdout, "INFO: the data-directory (%s) is: `%s`\n",
            description, out);
  }
  return out;
}
