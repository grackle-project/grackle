//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement logic of the @ref Error type
///
//===----------------------------------------------------------------------===//

#include "./error.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void Error::write(std::FILE* stream, bool append_newline) const {
  // C++23 could just use std::print/std::println for this purpose
  // -> the current implementation could also be much more efficient!
  const char* suffix = append_newline ? "\n" : "";
  std::string tmp = std::format("{}{}", *this, suffix);
  std::fprintf(stream, "%s", tmp.c_str());
}

}  // namespace GRIMPL_NAMESPACE_DECL