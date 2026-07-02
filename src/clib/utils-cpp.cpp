//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements C++ utilities used in the routines transcribed from Fortran
///
//===----------------------------------------------------------------------===//

#include <cstdarg>
#include <cstdio>

#include "utils-cpp.hpp"

int eprintf(const char* format, ...) {
  std::va_list args;
  va_start(args, format);
  int out = vfprintf(stderr, format, args);
  va_end(args);
  return out;
}
