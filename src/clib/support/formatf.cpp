//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement machinery related to printf-style string formatting
///
//===----------------------------------------------------------------------===//

#include <cstdarg>  // va_start, va_copy, va_end, std::va_list
#include <cstdio>   // std::vsnprintf

#include "./formatf.hpp"

namespace GRIMPL_NAMESPACE_DECL {

std::string str_formatf(const char* s, ...) {
  std::va_list vlist;
  va_start(vlist, s);
  std::string out = vstr_formatf(s, vlist);
  va_end(vlist);
  return out;
}

std::string vstr_formatf(const char* s, std::va_list vlist) {
  std::va_list vlist_copy;
  va_copy(vlist_copy, vlist);

  // call vsnprintf to get the size of the output buffer
  std::size_t sz_without_terminator = std::vsnprintf(nullptr, 0, s, vlist);
  va_end(vlist);

  // initialize the std::string with `sz_without_terminator` characters (it's
  // filled without ' ' chars). In practice, this allocates a buffer with
  // `sz_without_terminator + 1` characters and the final character is '\0'
  std::string out(sz_without_terminator, ' ');

  // call std::vsnprintf a 2nd time to write the formatted string to `out`
  // - `out.data()` points to the underlying buffer that will be mutated
  //   (note: mutating `out.data()` was illegal before C++17)
  // - the 2nd argument of std::vsnprintf is the size of the output buffer. We
  //   **NEED** to specify the size including the terminator character in this
  //   argument, otherwise the formatted message will be truncated
  // - In this function call, std::vsnprintf may try to overwrite the value at
  //   `out.data()[sz_without_terminator]` so that it stores '\0' (i.e. it
  //   doesn't know that a '\0' is already stored there). This is explicitly
  //   allowed by the C++ standard.
  std::vsnprintf(out.data(), sz_without_terminator + 1, s, vlist_copy);
  va_end(vlist_copy);
  return out;
}

}  // namespace GRIMPL_NAMESPACE_DECL