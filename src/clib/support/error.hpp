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
#include <type_traits>
#include <string>
#include <string_view>

#include "./config.hpp"
#include "./error_detail.hpp"

namespace GRIMPL_NAMESPACE_DECL {

namespace error_detail {

// std::format_string was retroactively exposted in C++ 20.
// -> this facillitates compile-time for wrapped calls to std::format and
//    std::vformat that explicitly check (at compile-time) whether the format
//    string is consistent with the number of specified values to be formatted
//    and with their types (compiler manually perform comparable checks for
//    printf)
// -> this ifdef statement directly recommended by the report introducing this
//    retroactive change
// -> when we adopt C++ 23 as a minimum version, we can assume that
//    std::format_string is always provided
#if __cpp_lib_format >= 202207L
template <typename... Args>
using my_format_string = std::format_string<std::type_identity_t<Args>...>;
#else
template <typename... Args>
using my_format_string = std::string_view;
#endif
}  // namespace error_detail

/// @brief Represents a generic Error type
///
/// @note
/// The output format and api takes a little inspiration from the anyhow
/// rust crate. The idea of wrapping the implementation in a shared pointer
/// was inspired by the Error type in the jiff rust crate
class Error {
  friend struct std::formatter<GRIMPL_NS::Error>;

  static constexpr const char* DFLT_MSG_ = "<UNKNOWN ERROR>";

  // This is wrapped by a shared_ptr (i.e. its heap-allocated with
  // atomic reference counting) in order to make it easier to
  std::shared_ptr<ErrImpl_> impl_;  ///< internal representation of Error

  Error() = default;  // <- this is private to force use of factory methods

  Error& context_helper_(ErrImpl_ new_ctx_err) {
    if (impl_.get() == nullptr) {  // <- possible after move-operation
      *this = Error::msg_literal(Error::DFLT_MSG_);
    }
    std::shared_ptr<ErrImpl_> tmp = std::make_shared<ErrImpl_>(new_ctx_err);
    tmp->err_cause_ = this->impl_;
    this->impl_ = tmp;
    return *this;
  }

public:
  Error(Error&&) = default;
  Error(const Error&) = default;
  Error& operator=(Error&&) = default;
  Error& operator=(const Error&) = default;
  ~Error() = default;

  /// @brief wraps the existing err information in additional context
  ///
  /// This uses C++'s modern string formatting syntax (equivalent to python's
  /// formatting mini-language). Given a error object `err`, invoking
  ///    ``err.context("{} is a {}", 1, "number");``
  /// introduces context comparable to invoking
  ///    ``err.context_literal("1 is a number");``
  /// (under the hood, the way
  ///
  /// @param fmt The format-string. This **MUST** be a string literal.
  /// @param args optional arguments to be formatted
  ///
  /// @warning
  /// Passing a non-literal string as @p fmt introduces undefined behavior.
  /// (While older C++ compilers may compile the code, newer compilers will
  /// explicitly refuse to compile the program).
  ///
  /// @note
  /// The proper way to wrap an error object `err` in a context message encoded
  /// in a std::string object called `s` (this object may have been dynamically
  /// constructed) is to call `err.context("{}", s);`
  template <typename... Args>
  Error& context(error_detail::my_format_string<Args...> fmt, Args&&... args) {
#if __cpp_lib_format >= 202207L
    std::string_view fmt_sv = fmt.get();
#else
    std::string_view& fmt_sv = fmt;
#endif
    // in the future (with a little refactoring):
    //   if (sizeof...(Types) == 0) -> we can skip allocating a std::string
    std::string msg = std::vformat(fmt_sv, std::make_format_args(args...));
    return context_helper_(ErrImpl_(nullptr, std::move(msg), nullptr));
  }

  /// @brief wrap the existing err info in additional context
  ///
  /// returns a reference to `this` for convenience
  ///
  /// @note
  /// This exists because the vast majority of Grackle's error messages are
  /// string literals (if we don't coerce to a std::string that reduces heap
  /// usage)
  Error& context_literal(const char* msg) {
    return context_helper_(ErrImpl_{msg, "", nullptr});
  }

  /// @brief convenience method to make it easier to write to stderr
  void write(std::FILE* stream, bool append_newline = true) const;

  // factory methods (we may add more in the future!)
  // ================================================

  /// @brief Construct an error object by formatting an error message
  ///
  /// This uses C++'s modern string formatting syntax (equivalent to python's
  /// formatting mini-language). For example, invoking
  ///    ``Error::msg("{} is a {}", 1, "number");``
  /// represent a error-message analogous to
  ///    ``Error::msg_literal("1 is a number");``
  /// (under the hood, the internal representation is different)
  ///
  /// @param fmt The format-string. This **MUST** be a string literal.
  /// @param args optional arguments to be formatted
  /// @returns An error object
  ///
  /// @warning
  /// Passing a non-literal string as @p fmt introduces undefined behavior.
  /// (While older C++ compilers may compile the code, newer compilers will
  /// explicitly refuse to compile the program).
  ///
  /// @note
  /// The proper way to create an error object encoding a message copied from
  /// a string `s` (this object may have been dynamically constructed) is to
  /// call `Error::msg("{}", s);`
  template <typename... Args>
  static Error msg(error_detail::my_format_string<Args...> fmt,
                   Args&&... args) {
    // in future: if sizeof...(Types) == 0, we can skip allocating a std::string
#if __cpp_lib_format >= 202207L
    std::string_view fmt_sv = fmt.get();
#else
    std::string_view& fmt_sv = fmt;
#endif
    // in the future (with a little refactoring):
    //   if (sizeof...(Types) == 0) -> we can skip allocating a std::string
    Error out;
    std::string msg = std::vformat(fmt_sv, std::make_format_args(args...));
    out.impl_ = std::make_shared<ErrImpl_>(nullptr, std::move(msg), nullptr);
    return out;
  }

  /// @brief Construct an error from a string-literal message
  ///
  /// @note
  /// This exists because the vast majority of Grackle's error messages are
  /// string literals (if we don't coerce to a std::string that reduces heap
  /// usage)
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
    if (error.impl_.get() == nullptr) {  // <- possible after move operation
      return std::format_to(ctx.out(), "{}", GRIMPL_NS::Error::DFLT_MSG_);
    }

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