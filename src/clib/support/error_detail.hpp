//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement @ref ErrorImpl Type
///
//===----------------------------------------------------------------------===//

#include <memory>
#include <string>
#include <string_view>

#include "./config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Internal state of an error
struct ErrImpl_ {
  // these attributes are used to encode the error information
  // - only one of these is used at a time. msg_literal exists as a minor
  //   optimization (to reduce heap allocated memory)
  // - in the future, we could add additional representations for parameterizing
  //   common classes of errors (perhaps an out-of-range error) that requires
  //   less memory than a full string
  // - if we do that, we may want to track the error info inside a std::variant
  //   (i.e. a type-safe union)
  std::string_view msg_literal;
  std::string msg;

  /// @brief may point to the cause of this error
  std::shared_ptr<ErrImpl_> err_cause_;

  // in the future, we could also consider tracking information like the
  // location (file, line number, function name?) where the error occurred or
  // perhaps even a stack-trace

  /// @brief return a view of the error message
  ///
  /// @warning
  /// Undefined behavior occurs if the returned string_view outlives ``this``
  ///
  /// @note
  /// If we make it possible to encode errors in alternative ways (without
  /// internally tracking a string message), we may need to replace this
  /// function
  std::string_view get_string_view() const {
    if (msg_literal.size() != 0) {
      return msg_literal;
    } else {
      return msg;
    }
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL