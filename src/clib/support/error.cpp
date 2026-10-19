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

#include <sstream>

#include "./error.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void Error::write(std::FILE* stream, bool append_newline) const {
  const char* suffix = append_newline ? "\n" : "";
  std::string tmp = to_string();
  std::fprintf(stream, "%s%s", tmp.c_str(), suffix);
}

std::string Error::to_string() const {
  if (impl_.get() == nullptr) {  // <- possible after move operation
    return GRIMPL_NS::Error::DFLT_MSG_;
  }
  // I suspect that this implementation is quite inefficient
  // -> I THINK std::ostringstream is just wrapping a std::string and every time
  //    we add more to the output we are reallocating the buffer underpinning
  //    std::string
  // -> ideally, we would determine the full length of the output string,
  //    allocate that and then fill the buffer. Unfortunately there is no
  //    concise way to this. I believe C++ 20's std::format machinery implements
  //    this exact strategy (we can't currently use that machinery for reasons
  //    explained in comments in the header file)
  // -> in practice, I think that it's probably ok to punt on dealing with this
  //    inefficiency... We are only executing this logic when we want to
  //    present an error message to an end user (in which case the program is
  //    realistically going to stop executing). And if we wait long enough,
  //    we can start using std::format
  std::ostringstream s;
  s << impl_->get_string_view();

  // print out the chain of causes (if any)
  ErrImpl_* c = impl_->err_cause_.get();
  if (c != nullptr) {
    s << "\n\nCaused By:";
    for (int count = 1; c != nullptr; c = c->err_cause_.get(), count++) {
      if (count > 1 || c->err_cause_.get() != nullptr) {
        s << "\n  " << count << ": " << c->get_string_view();
      } else {
        s << "\n     " << c->get_string_view();
      }
    }
  }
  // we move the final string out of std::ostringstream, rather than copy it,
  // to try to avoid a large heap allocation
  return std::move(s).str();
}

}  // namespace GRIMPL_NAMESPACE_DECL