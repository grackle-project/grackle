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

#include <iomanip>  // std::setw
#include <sstream>  // std::ostringstream
#include <utility>  // std::move

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

  std::size_t padding_width = pos_ + 1;
  std::ostringstream s;
  // -> first line specifies the parsing error
  s << "parsing error: " << descr
    << '\n'
    // -> second line shows the line being parsed
    << "   |" << parsed_s
    << '\n'
    // -> final line uses a carrot, '^', to point to the position where the
    //    parsing error actually occured
    << "    " << std::setw(padding_width) << '^';

  // we try to move string buffer out of the stringstream rather than copy it
  return std::move(s).str();
}

}  // namespace parse
}  // namespace GRIMPL_NAMESPACE_DECL