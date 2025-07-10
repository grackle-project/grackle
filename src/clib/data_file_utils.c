/***********************************************************************
/
/ Implement logic used internally by Grackle to determine data files.
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <limits.h> // CHAR_BIT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data_file_utils.h"
#include "file_registry.h"
#include "sha256.h"
#include "status_reporting.h"
#include "os_utils.h"

#include "grackle.h" // get_grackle_version


#define CKSUM_ALGORITHM "sha256"
#define CKSUM_STR_PREFIX CKSUM_ALGORITHM ":"
#define CKSUM_DIGEST_N_BYTES 32
#define CKSUM_DIGEST_N_HEXDIGITS (2*CKSUM_DIGEST_N_BYTES)

// confirm a byte is 8 bits
// -> the C standard technically allows it to be larger. But on any modern
//    POSIX system (or even Windows) it must be 8 bytes
// -> this scenario only comes up on highly specialize DSP hardware
#if CHAR_BIT != 8
  #error "our assumption that a byte is 8 bits is violated"
#endif

/// returns whether 2 null-terminated checksum strings are equal
///
/// A checksum string is composed of lowercase letters & has 2 parts:
/// - a prefix that includes the name of a hash algorithm used to compute the
///   checksum followed by a colon (e.g. `md5:`, `sha1:`, `sha256:`)
/// - the suffix that specifies the actual values of the checksum as a string
///   of hexadecimal digits.
static int cksum_str_eq_(const char* lhs, const char*rhs) {
  if ((lhs == NULL) || (rhs == NULL)) return 0;
  return strcmp(lhs, rhs) == 0;
}

/// Converts a checksum digest into a hexadecimal string
///
/// @param[in] digest is an array of bytes where each byte has an
///     arbitrary value from 0 to 255
/// @param[in] digest_len is the length of digest
/// @param[out] str is an emtpy array of length `2*digest_len + 1`
///     entries. At the conclusion of this operation,
///     `str[i*2:i*2+2]` specifies the value of `digest[i]` in
///     hexadecimal notation. `str[digest_len*2]` will be assigned
///     the null terminator.
static void convert_to_hex_(const unsigned char* digest, int digest_len, char* str) {

  // some important context: the standard does not specify whether `char` is 
  //   signed or unsigned and the call to snprintf will only work if we
  //   pass the values of each byte as an unsigned char.
  //
  // Thus: we need to explicitly reinterpret the value of each element digest
  //       as an unsigned char.
  // - there are rules in the C & C++ standard for this topic. Consider
  //   an object of T1. We want to access it through a value of type T2.
  //   For arbitrary types, this "type punning" is undefined behavior
  //   (i.e. standards compliant compilers are free to do whatever they
  //   want when they encounter undefined behavior without any consistency)
  // - C++ is generally stricter about this topic (e.g. they forbid using
  //   unions to reinterpret values)
  // - there are exceptions when it comes to `unsigned char` & `char`
  // - discussions of these rules for C & C++ are found at
  //   https://en.cppreference.com/w/c/language/object#Strict_aliasing
  //   https://en.cppreference.com/w/cpp/language/reinterpret_cast#Type_aliasing

  for (int i = 0; i < digest_len; i++){
    // while it may seem like there are faster ways to do this, please don't
    // change this without including a reference or argument explaining why
    // your approach won't invoke undefined behavior.

    char elem = digest[i];
#ifdef __cplusplus
    unsigned char *uchar_ptr = reinterpret_cast<unsigned char*>(&elem);
#else
    unsigned char *uchar_ptr = (unsigned char*)(&elem);
#endif
    snprintf(str + 2*i, 3, "%02hhx", *uchar_ptr);
  }
}

char* calc_checksum_str_(const char* fname) {

  FILE* fp = fopen(fname, "rb");
  if (!fp) {
    GrPrintErrMsg("can't open `%s` to calculate checksum. Does it exist?\n",
                  fname);
    return NULL;
  }

  hash_state ctx;
  sha256_init(&ctx);

  const size_t CHUNKSIZE = 4096;
  char* buffer = malloc(CHUNKSIZE);

  int any_data_read = 0;
  size_t cur_chunk_len;
  do {
    cur_chunk_len = fread(buffer, 1, CHUNKSIZE, fp);
    if (cur_chunk_len != 0) {
      sha256_process(&ctx, (unsigned char*)buffer,
                     (unsigned long)cur_chunk_len);
      any_data_read = 1;
    }
  } while(cur_chunk_len == CHUNKSIZE);
  free(buffer);
  fclose(fp);

  if (!any_data_read) {
    GrPrintErrMsg("`%s` specifies a path to an empty file\n", fname);
    return NULL;
  }

  unsigned char digest[CKSUM_DIGEST_N_BYTES];
  sha256_done(&ctx, digest);

  // now we just need to convert all the bytes to a string of hexadecimal
  // digits for the sake of comparison
  const char prefix[] = CKSUM_STR_PREFIX;
  size_t prefix_len = strlen(prefix); // excludes nul character
  size_t out_length = prefix_len + (2 * sizeof(digest)) + 1; // add 1 for nul

  char* out = malloc(out_length);
  memcpy(out, prefix, prefix_len);
  convert_to_hex_(digest, sizeof(digest), out+prefix_len);
  return out;
}


static struct generic_file_props file_from_data_dir_(
  const char* grackle_data_file, int grackle_data_file_options,
  int calculate_checksum, const char** registry
)
{
  // initialize out in a format that denotes an error (if its never modified)
  struct generic_file_props out = {NULL, 0, NULL, 0};

  // first, let's check if the specified file name is known to Grackle
  const char* expected_cksum_str = expected_file_cksum_(grackle_data_file,
                                                        registry);
  if (expected_cksum_str == NULL) {
    // in the future, depending on the value of grackle_data_file_options,
    // we may want to provide special handling for the case where
    // grackle_data_file starts with `"user-data/..."

    GrPrintErrMsg(
      "can't load %s from data directory (it isn't in the file registry)\n",
      grackle_data_file);
    return out;
  }

  // now it's time to construct the full path to the file (if it exists)
  grackle_version version_info = get_grackle_version();

  // get the data_directory
  char* data_dir_path = get_data_dir_(get_platform_());
  if (data_dir_path == NULL) return out;

  const char* path_parts[4] = {
    data_dir_path, "data-store-v1", version_info.version, grackle_data_file
  };
  char* full_path = join_parts_('/', path_parts, 4);
  free(data_dir_path);

  if (calculate_checksum == 0) { // skip the checksum calculation
    out.path = full_path;
    out.path_requires_dealloc = 1;
    return out;
  }

  char* measured_cksum_str = calc_checksum_str_(full_path);

  if (measured_cksum_str == NULL) {
    return out;
  } else if (cksum_str_eq_(measured_cksum_str, expected_cksum_str) == 0) {
    GrPrintErrMsg(
        "The measured checksums doesn't match expectations\n"
        "    -> measured: \"%s\"\n"
        "    -> expected: \"%s\"\n"
        "    -> path: `%s`\n"
        "  This error is indicative of 1 of 3 possible scenarios:\n"
        "     1. There is a bug in the core Grackle library.\n"
        "     2. There is a bug in the Grackle's data-file management tool.\n"
        "     3. It isn't Grackle's fault. Either the datafile was\n"
        "        corrupted or its the fault of the user/some other tool\n"
        "        that tried to modify the file.\n",
        measured_cksum_str, expected_cksum_str, full_path);
    free(measured_cksum_str);
    free(full_path);
  } else {
    out.path = full_path;
    out.path_requires_dealloc = 1;
    out.checksum = measured_cksum_str;
    out.checksum_requires_dealloc = 1;
  }
  return out;
}


struct generic_file_props determine_data_file_(const char* grackle_data_file,
                                               int grackle_data_file_options,
                                               const char** registry)
{
  // initialize output struct in a format that will denote an error (if is
  // never modified)
  struct generic_file_props out = {NULL, 0, NULL, 0};

  if (grackle_data_file == NULL) {
    GrPrintErrMsg("grackle_data_file must not be NULL\n");
    return out;
  }

  if (grackle_data_file_options == -1) {
    grackle_data_file_options = GR_DFOPT_FULLPATH_NO_CKSUM;  // the legacy case
  }

  switch (grackle_data_file_options) {
    case GR_DFOPT_FULLPATH_NO_CKSUM: {
      out.path = grackle_data_file;
      return out;
    }
    case GR_DFOPT_MANAGED: {
      return file_from_data_dir_(grackle_data_file, grackle_data_file_options,
                                 1, registry);
    }
    case GR_DFOPT_MANAGED_NO_CKSUM: {
      return file_from_data_dir_(grackle_data_file, grackle_data_file_options,
                                 0, registry);
    }
    default: {
      GrPrintErrMsg("grackle_data_file_options has an unexpected value: %d\n",
                    grackle_data_file_options);
      return out;
    }
  }
}

void free_generic_file_props_(struct generic_file_props* ptr) {
  if (ptr != NULL) {
    if (ptr->path_requires_dealloc){
      free((char*)ptr->path);
    }
    ptr->path = NULL;
    if (ptr->checksum_requires_dealloc){
      free((char*)ptr->checksum);
    }
    ptr->checksum = NULL;
  }
}
