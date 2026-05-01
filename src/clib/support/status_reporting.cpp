
//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// @brief Implements status reporting functionality
///
//===----------------------------------------------------------------------===//

#include <cstdarg>  // std::va_list, va_copy, va_end, va_start
#include <cstdio>   // std::fflush, std::fprintf, stderr, std::vsnprintf
#include <cstdlib>  // std::abort
#include <vector>

#include "status_reporting.hpp"
#include "grackle.h" // GR_FAIL

// this is the internal routine that everything else dispatches to
static void vprint_err_(
  int internal_error, const struct grimpl_source_location_ locinfo,
  const char* msg, std::va_list vlist
) {
  const char* santized_func_name = (locinfo.fn_name == nullptr)
    ? "{unspecified}" : locinfo.fn_name;

  const char* fallback_msg_ = "{nullptr encountered instead of error message}";
  std::vector<char> dynamic_msg_buf;
  const char* msg_buf;
  if (msg == nullptr) {
    msg_buf = fallback_msg_;
  } else {
    // make a copy of the variadic function arguments
    std::va_list vlist_copy;
    va_copy(vlist_copy, vlist);

    // get the total size of the formatted message
    std::size_t msg_len = std::vsnprintf(nullptr, 0, msg, vlist_copy) + 1;
    va_end(vlist_copy);

    // allocate the buffer to hold the message
    dynamic_msg_buf.resize(msg_len);

    // actually format the message
    std::vsnprintf(dynamic_msg_buf.data(), msg_len, msg, vlist);
    va_end(vlist);

    msg_buf = dynamic_msg_buf.data();
  }

  const char* descr = (internal_error == 1) ? "FATAL" : "ERROR";

  std::fprintf(
    stderr, "Grackle-%s %s:%d in %s] %s\n", descr, locinfo.file,
    locinfo.lineno, santized_func_name, msg_buf
  );
}

void grimpl_abort_with_internal_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  std::va_list args;
  va_start(args, msg);
  vprint_err_(1, locinfo, msg, args);
  // while stderr should flush by default, people may overwrite this behavior.
  // We want to force flushing here since we are aborting the program
  std::fflush(stderr);
  std::abort();
}

int grimpl_print_and_return_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  std::va_list args;
  va_start(args, msg);
  vprint_err_(0, locinfo, msg, args);
  va_end(args);
  return GR_FAIL;
}

void grimpl_print_err_msg_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  va_list args;
  va_start(args, msg);
  vprint_err_(0, locinfo, msg, args);
  va_end(args);
}

