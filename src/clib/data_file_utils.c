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

#include <ctype.h> // tolower
#include <limits.h> // CHAR_BIT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "picohash.h"

#include "data_file_utils.h"
#include "os_utils.h"

#include "grackle.h" // get_grackle_version


#define CKSUM_ALGORITHM "sha1"
#define CKSUM_STR_PREFIX CKSUM_ALGORITHM ":"
#define CKSUM_DIGEST_N_BYTES PICOHASH_SHA1_DIGEST_LENGTH
//#define CKSUM_DIGEST_N_BYTES 20
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
/// A checksum string consists of 2 parts:
/// - a prefix that includes the name of a hash algorthim used to compute the
///   checksum followed by a colon (e.g. `md5:`, `sha1:`, `sha256:`)
/// - the suffix that specifies the actual values of the checksum as a string
///   of hexadecimal digits.
///
/// @note
/// We could make this faster by encoding the checksum as an array of bytes
/// (rather than a string of hexadecimal digits).
///  - This would involve half the memory and we wouldn't need to worry about
///    case-insensitivity.
///  - But it's not worth the effort to do this to perform just a single
///    checksum comparison. (we need to compute the string-representation
///    anyway in order to effectively communicate issues with library users)
static int cksum_str_eq_(const char* lhs, const char*rhs)
{
  // locales could theoretically be an issue here... (but we should be fine)
  // as long as the strings only contain latin letters (without modifiers)
  // and arabic numerals
  if ((lhs == NULL) || (rhs == NULL)) return 0;

  size_t len = strlen(lhs); // excludes trailing '\0'
  if ((len == 0) || (len != strlen(rhs))) return 0;

  int neq = 0;
  for (size_t i = 0; i < len; i++){
    neq += (tolower(lhs[i]) == tolower(rhs[i]));
  }
  return (len == (size_t)neq);
}

/// abort the program with an error message if the checksum string
/// isn't valid
///
/// we abort, rather than return NULL because there is a programming error
/// (and people can simply avoid this error by running their program without
/// any cksum calculations)
///
/// behavior is undefined when cksum_str is NULL
void assert_valid_cksum_str_(const char* cksum_str,
                             const char* cksum_origin_descr,
                             const char* extra_fmt_arg) {
  char* err = NULL;

  const char* hexstr_start = post_prefix_ptr_(cksum_str, CKSUM_STR_PREFIX);
  const char* colon_pos = strchr(cksum_str, ':');

  // ignore '\0' in length calculation
  size_t hexstr_len = (hexstr_start == NULL) ? 0 : strlen(hexstr_start);

  if ((hexstr_start == NULL) && (colon_pos == NULL)){
    err = my_strdup_(
      "no prefix specifying an algorithm name (e.g. \"" CKSUM_STR_PREFIX "\")"
    );
  } else if (hexstr_start == NULL) {
    err = my_strdup_(
      "the algorithm name (i.e. characters before the colon), doesn't match"
      "  \"" CKSUM_ALGORITHM "\""
    );
  } else if (hexstr_len != CKSUM_DIGEST_N_HEXDIGITS) {
    const char fmt[] = "it should have %d characters after the prefix, not %d";
    int sz = snprintf(err, 0, fmt, CKSUM_DIGEST_N_HEXDIGITS, (int)hexstr_len);
    err = malloc(sz+1);
    snprintf(err, sz, fmt, CKSUM_DIGEST_N_HEXDIGITS, (int)hexstr_len);

  } else {
    const char hexdigits[] = "0123456789abcdefABCDEF";
    int bad_digits = 0;
    for (int i = 0; i < CKSUM_DIGEST_N_HEXDIGITS; i++) {
      bad_digits += (strchr(hexdigits, hexstr_start[i]) == NULL);
    }
    if (bad_digits) {
      err = strdup(
        "the characters after the prefix include non-hexadecimal digit(s)"
      );
    }
  }

  // let's perform some sanity checks on the contents of this string!
  if (err != NULL) {
    const char* extra_fmt = (extra_fmt_arg == NULL) ? "" : extra_fmt_arg;
    fprintf(
      stderr,
      ("INTERNAL ERROR: There is a problem with a checksum string\n"
       "  string value: \"%s\"\n"
       "  origin: %s %s\n"
       "  issue: %s\n"),
      cksum_str, cksum_origin_descr, extra_fmt_arg, err);
    free(err);
    abort();
  }
}

// ===========================================
// LOGIC TO BE RELOCATED BEGIN
// ===========================================
// the logic in the section will be relocated to a separate header file

typedef struct { const char* fname; const char* cksum; } registry_entry;

static registry_entry file_registry[] = {
  // for now, this includes a subset (will be autogenerated in the future)
  {"CloudyData_UVB=FG2011.h5", "sha1:5b3423fb5cb96d6f8fae65655e204f1f82a276fa"},
  {"CloudyData_UVB=HM2012.h5", "sha1:3ae95f71926aa9543964fbd41c5e53a42345c19c"},
};

/// return the full checksum string of the file if it is in the registry (or
/// NULL if there isn'a match)
static inline const char* expected_file_cksum_(const char* fname) {
  if (fname == NULL) return NULL;

  const size_t n_entries = sizeof(file_registry) / sizeof(registry_entry);
  const char* cksum_str = NULL;
  for (size_t i = 0; i < n_entries; i++) {
    if (strcmp(fname, file_registry[i].fname) == 0) {
      return file_registry[i].cksum;
    }
  }
  return NULL;
}
// ===========================================
// LOGIC TO BE RELOCATED END
// ===========================================

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
static void convert_to_hex_(char* digest, int digest_len, char* str) {

  // some important context: the standard does not specify whether `char` is 
  //   signed or unsigned and the call to snprintf will only only work if we
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
    snprintf(str + 2*i, 3, "%02hhx\n", *uchar_ptr);
  }
}

/// calculate the checksum for the specified file
static char* calc_checksum_str_(const char* fname) {

  FILE* fp = fopen(fname, "rb");
  if (!fp) {
    fprintf(stderr,
            ("ERROR: unable to open `%s` to calculate checksum. Does the file "
             "actually exist?"),
            fname);
    return NULL;
  }

  picohash_ctx_t ctx;
  picohash_init_sha1(&ctx);

  const size_t CHUNKSIZE = 4096;
  char* buffer = malloc(CHUNKSIZE);

  int any_data_read = 0;
  int cur_len;
  do {
    cur_len = fread(buffer, 1, CHUNKSIZE, fp);
    if (cur_len != 0) {
      picohash_update(&ctx, buffer, cur_len);
      any_data_read = 1;
    }
  } while(cur_len == CHUNKSIZE);
  free(buffer);
  fclose(fp);

  if (!any_data_read) {
    fprintf(stderr, "ERROR: `%s` either specifies a path to an empty file\n",
            fname);
  }

  char digest[PICOHASH_SHA1_DIGEST_LENGTH];
  picohash_final(&ctx, digest);

  // now we just need to convert all of the bytes to a string of hexadecimal
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
  const char* grackle_data_file, int grackle_data_file_options
)
{
  // initialize output struct in a format that will denote an error (if is
  // never modified)
  struct generic_file_props out = {NULL, 0, NULL, 0};

  // first, let's check if the specified file name is known to Grackle
  const char* expected_cksum_str = expected_file_cksum_(grackle_data_file);
  if (expected_cksum_str == NULL) {
    // in the future, depending on the value of grackle_data_file_options,
    // we may want to provide special handling for the case where
    // grackle_data_file starts with `"user-data/..."

    fprintf(stderr,
            "ERROR: can't load %s from data directory, no such file is in "
            "the file registry\n",
            grackle_data_file);
    return out;
  }

  // sanity check that checksum from the file registry was properly formatted
  assert_valid_cksum_str_(expected_cksum_str,
                          "from the file-registry for the file named",
                          grackle_data_file);

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

  char* measured_cksum_str = calc_checksum_str_(full_path);

  if (measured_cksum_str == NULL) {
    return out;
  } else if (cksum_str_eq_(measured_cksum_str, expected_cksum_str) == 0) {
    fprintf(stderr,
            "ERROR: the measured checksums doesn't match expectations\n"
            "     -> measured: \"%s\"\n"
            "     -> expected: \"%s\"\n"
            "     -> path: `%s`\n"
            "  This error is indicative of 1 of 3 scenarios:\n"
            "     1. There is a bug in the core Grackle library for locating\n"
            "        the file or computing the checksum\n"
            "     2. There is a bug in the Grackle's data-file management\n"
            "        tool.\n"
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
                                               int grackle_data_file_options)
{
  // initialize output struct in a format that will denote an error (if is
  // never modified)
  struct generic_file_props out = {NULL, 0, NULL, 0};

  if (grackle_data_file == NULL) {
    fprintf(stderr, "grackle_data_file must not be NULL\n");
    return out;

  } else if (grackle_data_file_options == -1) { // the legacy case!
    out.path = grackle_data_file;
    return out;
  } else if (grackle_data_file_options == 1) {
    return file_from_data_dir_(grackle_data_file, grackle_data_file_options);
  } else {
    fprintf(stderr, "grackle_data_file_options has an unexpected value: %d\n",
            grackle_data_file_options);
    return out;

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
