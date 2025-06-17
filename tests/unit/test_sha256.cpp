// testing sha256

#include "gtest/gtest.h"
#include "sha256.h"
#include "grtest_os.hpp"

#include <filesystem>
#include <string>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdio>
#include <filesystem>

#define CKSUM_ALGORITHM "sha256"
#define CKSUM_STR_PREFIX CKSUM_ALGORITHM ":"
#define CKSUM_DIGEST_N_BYTES 32
#define CKSUM_DIGEST_N_HEXDIGITS (2*CKSUM_DIGEST_N_BYTES)

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
static void convert_to_hex_(unsigned char* digest, int digest_len, char* str) {

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

  Hash_state ctx;
  sha256_init(&ctx);

  const size_t CHUNKSIZE = 4096;
  char* buffer = (char*)malloc(CHUNKSIZE);

  int any_data_read = 0;
  int cur_len;
  do {
    cur_len = fread(buffer, 1, CHUNKSIZE, fp);
    if (cur_len != 0) {
      sha256_process(&ctx, (unsigned char*)buffer,
                     (unsigned long)cur_len);
      any_data_read = 1;
    }
  } while(cur_len == CHUNKSIZE);
  free(buffer);
  fclose(fp);

  if (!any_data_read) {
    fprintf(stderr, "ERROR: `%s` either specifies a path to an empty file\n",
            fname);
  }

  unsigned char digest[CKSUM_DIGEST_N_BYTES];
  sha256_done(&ctx, digest);

  // now we just need to convert all of the bytes to a string of hexadecimal
  // digits for the sake of comparison
  const char prefix[] = CKSUM_STR_PREFIX;
  size_t prefix_len = strlen(prefix); // excludes nul character
  size_t out_length = prefix_len + (2 * sizeof(digest)) + 1; // add 1 for nul

  char* out = (char*)malloc(out_length);
  memcpy(out, prefix, prefix_len);
  convert_to_hex_(digest, sizeof(digest), out+prefix_len);
  return out;
}

static std::string get_sha256_cksum(const std::string& str) {
  // we are going to work out of a temporary dir
  grtest::TempDir tmpdir = grtest::TempDir::create();
  std::string path = (tmpdir.get_path() / "my-file.txt").string();

  std::FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) { return ""; }
  std::fprintf(fp, "%s", str.c_str());
  std::fclose(fp);

  char* tmp = calc_checksum_str_(path.c_str());
  std::string out;
  if (tmp != nullptr) {
    out = tmp;
    std::free(tmp);
  }
  return out;
}

TEST(sha256, abc) {
  std::string message = "abc";
  std::string answer =
      "sha256:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
  std::string result = get_sha256_cksum(message);
  ASSERT_EQ(result, answer);
}

TEST(sha256, other) {
  std::string message =
      "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  std::string answer =
      "sha256:248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1";
  std::string result = get_sha256_cksum(message);
  ASSERT_EQ(result, answer);
}

TEST(sha256, other2) {
  std::string message =
      "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmnoijklmnop"
      "jklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
  std::string answer =
      "sha256:cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac45037afee9d1";
  std::string result = get_sha256_cksum(message);
  ASSERT_EQ(result, answer);
}
