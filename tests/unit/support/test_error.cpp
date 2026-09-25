//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of @ref GRIMPL_NS::Error
///
//===----------------------------------------------------------------------===//

#include "support/error.hpp"
#include <gtest/gtest.h>

#include <string>
#include <optional>

using Error = GRIMPL_NS::Error;

static Error erroneous_open(std::string path) {
  return Error::msg_literal("File not found");
}

TEST(ErrorFmt, ErroneousOpen) {
  Error obj = erroneous_open("path/to/file");
  std::string err_msg = std::format("{}", obj);
  ASSERT_EQ(err_msg, "File not found");
}

static Error erroneous_read_data(std::string path) {
  Error err = erroneous_open(path);
  err.context("problem loading data from {}", path);
  return err;
}

TEST(ErrorFmt, ErroneousReadData) {
  Error obj = erroneous_read_data("path/to/file");
  std::string err_msg = std::format("\n{}", obj);
  const char* expected = R"""(
problem loading data from path/to/file

Caused By:
     File not found)""";
  ASSERT_EQ(err_msg, expected);
}

static Error erroneous_read_InterpTable(std::string path) {
  Error err = erroneous_read_data(path);
  err.context_literal("unable to create InterpTable");
  return err;
}

TEST(ErrorFmt, ErroneousReadInterpTable) {
  Error obj = erroneous_read_InterpTable("path/to/file");
  std::string err_msg = std::format("\n{}", obj);
  const char* expected = R"""(
unable to create InterpTable

Caused By:
  1: problem loading data from path/to/file
  2: File not found)""";
  ASSERT_EQ(err_msg, expected);
}

TEST(Error, EquivMsg) {
  // let's verify the consistency of results for 3 different approaches that
  // should encode equivalent error messages

  constexpr int N_APPROACH = 4;
  Error errs[N_APPROACH] = {
      Error::msg("{} is a {}", 1, "number"),
      Error::msg("{}", std::string("1 is a number")),
      Error::msg("1 is a number"),
      Error::msg_literal("1 is a number"),
  };

  for (int i = 0; i < N_APPROACH; i++) {
    std::string msg = std::format("{}", Error::msg("{} is a {}", 1, "number"));
    EXPECT_EQ(msg, "1 is a number") << "issue with approach " << i;
  }
}

// there's an assumption that f is open in binary mode
static std::optional<std::string> read_full_file_(std::FILE& f) {
  // get number of bytes in the file and jump back to start!
  if (std::fseek(&f, 0, SEEK_END) != 0) {
    return std::nullopt;
  }
  long n_chars = std::ftell(&f);
  if (n_chars < 0) {  // denotes an error
    return std::nullopt;
  }
  if (std::fseek(&f, 0, SEEK_SET) != 0) {
    return std::nullopt;
  }

  std::string out;
  out.resize(n_chars);
  for (long i = 0; i < n_chars; i++) {
    out[i] = std::fgetc(&f);  // <- this is inefficient
  }

  return {out};
}

TEST(Error, Write) {
  Error err = Error::msg("{} is a {}", 1, "number");

  // write error message to a temporary file
  std::FILE* fp = std::tmpfile();
  ASSERT_TRUE(fp != nullptr) << "can't make temporary file";
  err.write(fp);  // <- by default, a newline is appended to the message

  // read the error message from the temporary file
  std::optional<std::string> maybe_content = read_full_file_(*fp);
  std::fclose(fp);
  ASSERT_TRUE(maybe_content.has_value()) << "issue reading from tmp file";
  std::string msg = maybe_content.value();

  EXPECT_EQ(msg, "1 is a number\n");
}