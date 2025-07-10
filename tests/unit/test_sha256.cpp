// testing sha256

#include "gtest/gtest.h"
#include "data_file_utils.h"
#include "grtest_os.hpp"

#include <filesystem>
#include <string>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cstdio>

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
