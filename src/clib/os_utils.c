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

#include <stdbool.h>
#include <stddef.h>  // size_t
#include <stdio.h>   // fprintf, stdout, fflush
#include <stdlib.h>  // malloc, realloc, free, abort
#include <string.h>  // memcpy, strlen

#include "os_utils.h"
#include "grackle.h"           // grackle_verbose
#include "status_reporting.h"  // GrPrintErrMsg

/// Just like getenv, except it returns NULL in place of strings of length 0.
static const char* getenv_nonempty_(const char* name) {
  const char* out = getenv(name);
  return ((out == NULL) || (out[0] == '\0')) ? NULL : out;
}

char* join_parts_(char sep, const char** parts, int nparts) {
  if (nparts < 1) return NULL;

  size_t total_len = 0;
  for (int i = 0; i < nparts; i++) {
    if (parts[i] == NULL) return NULL;
    total_len += strlen(parts[i]);  // we don't include the nul-terminator
  }
  total_len += (nparts - 1);  // account for the size of sep and
  total_len++;                // account for trailing nul-terminator

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
  GR_INTERNAL_REQUIRE((cur_offset + 1) == total_len, "sanity-check failed!");
  return out;
}

enum platform_kind get_platform_(void) {
#if defined(PLATFORM_GENERIC_UNIX)
  return platform_kind_generic_unix;
#else
  return platform_kind_unknown;
#endif
}

struct DataDirSpec {
  const char* envvar;
  const char* suffix_path;
  bool fatal_if_relative_path;
};

char* get_data_dir_(enum platform_kind kind) {
  const char* msg_common = "can't infer the data directory";

  if (kind == platform_kind_unknown) {
    GrPrintErrMsg("%s on an unknown platform.\n", msg_common);
    return NULL;
  }

  // following list specifies the order of directories in order of precedence
  // -> https://specifications.freedesktop.org/basedir-spec/latest/ instructs
  //    us to simply skip over XDG_DATA_HOME if it specifies a relative path
  // if you look back through the commit-history, you'll see that earlier drafts
  // of this functionality:
  // -> used "$HOME/Library/Application Support/grackle" on macOS. For
  //    simplicity, we now treat macOS like any other Unix
  // -> used standard posix functionality to infer $HOME if the envvar wasn't
  //    defined
  static const struct DataDirSpec data_dir_chain[] = {
      {"GRACKLE_DATA_DIR", NULL, true},
      {"XDG_DATA_HOME", "grackle", false},
      {"HOME", ".local/share/grackle", true},
  };

  const size_t n_entries = sizeof(data_dir_chain) / sizeof(struct DataDirSpec);
  for (size_t i = 0; i < n_entries; i++) {
    const char* env_str = getenv_nonempty_(data_dir_chain[i].envvar);

    if (env_str == NULL) {
      if (grackle_verbose) {
        fprintf(stdout, "INFO: %s using \"%s\" envvar (it isn't defined)\n",
                msg_common, data_dir_chain[i].envvar);
        fflush(stdout);  // flush in case we run into an error in the next call
      }
      continue;

    } else if (env_str[0] == '~') {
      GrPrintErrMsg("%s when \"%s\" envvar path starts with '~'\n", msg_common,
                    data_dir_chain[i].envvar);
      return NULL;

    } else if (env_str[0] != '/' && data_dir_chain[i].fatal_if_relative_path) {
      GrPrintErrMsg("%s when \"%s\" envvar specifies a relative path\n",
                    msg_common, data_dir_chain[i].envvar);
      return NULL;

    } else if (env_str[0] != '/') {
      if (grackle_verbose) {
        fprintf(stdout, "INFO: %s using \"%s\" envvar (it's a relative path)\n",
                msg_common, data_dir_chain[i].envvar);
        fflush(stdout);  // flush in case we run into an error in the next call
      }
      continue;

    } else {  // we are done!
      const char* parts[2] = {env_str, data_dir_chain[i].suffix_path};
      int nparts = 1 + (parts[1] != NULL);
      char* out = join_parts_('/', parts, nparts);

      if (grackle_verbose) {
        fprintf(stdout, "INFO: the data-directory is: `%s`\n", out);
      }
      return out;
    }
  }
  GrPrintErrMsg("%s: this probably means the \"HOME\" isn't an env variable\n",
                msg_common);
  return NULL;
}
