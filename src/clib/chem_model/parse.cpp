//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements some parsing logic
///
//===----------------------------------------------------------------------===//

#include <format>

#include "parse.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace parse {

std::string Err::to_string(std::string_view parsed_s) const {
  const char* descr;
  if (descr_ == nullptr) {
    descr = "unexpected character";
  } else {
    descr = descr_;
  }

  // -> first line specifies the parsing error
  // -> second line shows the line being parsed
  // -> final line uses a carrot, '^', to point to point to the position where
  //    the parsing error actually occured
  std::size_t padding_width = pos_ + 1;
  return std::format("parsing error: {}\n   |{}\n   |{:>{}c}",
                     descr,     // <- describes the error
                     parsed_s,  // <- shows the line being parsed
                     '^', padding_width);
}

}  // namespace parse
}  // namespace GRIMPL_NAMESPACE_DECL