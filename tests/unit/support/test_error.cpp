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

#include <format>
#include <string>

using Error = GRIMPL_NS::Error;

static Error erroneous_open(std::string path) {
  return Error::msg_literal("File not found");
}

TEST(Error, ErroneousOpen) {
  Error obj = erroneous_open("path/to/file");
  std::string err_msg = std::format("{}", obj);
  ASSERT_EQ(err_msg, "File not found");
}

static Error erroneous_read_data(std::string path) {
  Error err = erroneous_open(path);
  err.context("problem loading data from " + path);
  return err;
}

TEST(Error, ErroneousReadData) {
  Error obj = erroneous_read_data("path/to/file");
  std::string err_msg = std::format("{}", obj);
  const char* expected = R"""(problem loading data from path/to/file

Caused By:
     File not found)""";
  ASSERT_EQ(err_msg, expected);
}

static Error erroneous_read_InterpTable(std::string path) {
  Error err = erroneous_read_data(path);
  err.context("unable to create InterpTable");
  return err;
}

TEST(Error, ErroneousReadInterpTable) {
  Error obj = erroneous_read_InterpTable("path/to/file");
  std::string err_msg = std::format("{}", obj);
  const char* expected = R"""(unable to create InterpTable

Caused By:
  1: problem loading data from path/to/file
  2: File not found)""";
  ASSERT_EQ(err_msg, expected);
}
