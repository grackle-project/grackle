//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare the @ref Error type
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_ERROR_HPP
#define SUPPORT_ERROR_HPP

#include <cstdio>
#include <format>
#include <memory>  // std::shared_ptr;
#include <string>

#include "./config.hpp"
#include "./error_detail.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Represents a generic Error type
///
/// @note
/// The output format and api takes a little inspiration from the anyhow
/// rust crate. The idea of wrapping the implementation in a shared pointer
/// was inspired by the Error type in the jiff rust crate
class Error {
  friend struct std::formatter<GRIMPL_NS::Error>;

  // This is wrapped by a shared_ptr (i.e. its heap-allocated with
  // atomic reference counting) in order to make it easier to
  std::shared_ptr<ErrImpl_> impl_;  ///< internal representation of Error

  Error() = default;  // <- this is private to force use of factory methods
public:
  Error(Error&&) = default;
  Error(const Error&) = default;
  Error& operator=(Error&&) = default;
  Error& operator=(const Error&) = default;
  ~Error() = default;

  /// @brief wrap the existing err info in additional context
  ///
  /// @note returns a reference to this for convenience
  Error& context(std::string msg) {
    std::shared_ptr<ErrImpl_> tmp =
        std::make_shared<ErrImpl_>(nullptr, std::move(msg), nullptr);
    tmp->err_cause_ = this->impl_;
    this->impl_ = tmp;
    return *this;
  }

  // factory methods (we may add more in the future!)
  /// @brief Construct an error from an arbitrary message
  static Error msg(std::string msg) {
    Error out;
    out.impl_ = std::make_shared<ErrImpl_>(nullptr, std::move(msg), nullptr);
    return out;
  }

  /// @brief Construct an error from a string-literal message
  static Error msg_literal(const char* msg_literal) {
    Error out;
    out.impl_ = std::make_shared<ErrImpl_>(msg_literal, "", nullptr);
    return out;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

// by specializing std::formatter for GRIMPL_NS::Error, you can use std::format
// to get a string representation of GRIMPL_NS::Error.
template <>
struct std::formatter<GRIMPL_NS::Error> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) {
    return ctx.begin();
  }

  template <class FmtContext>
  auto format(const GRIMPL_NS::Error& error, FmtContext& ctx) const {
    using OutT = typename FmtContext::iterator;
    OutT out = std::format_to(ctx.out(), "{}", *error.impl_);

    // print out the chain of causes (if any)
    GRIMPL_NS::ErrImpl_* c = error.impl_->err_cause_.get();
    if (c != nullptr) {
      ctx.advance_to(out);
      out = std::format_to(ctx.out(), "\n\nCaused By:");
      for (int count = 1; c != nullptr; c = c->err_cause_.get(), count++) {
        ctx.advance_to(out);
        if (count > 1 || c->err_cause_.get() != nullptr) {
          out = std::format_to(out, "\n  {}: {}", count, *c);
        } else {
          out = std::format_to(out, "\n     {}", *c);
        }
      }
    }
    return out;
  }
};

#endif  // SUPPORT_ERROR_HPP