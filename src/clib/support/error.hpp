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

#include <cstdarg>
#include <cstdio>
#include <format>
#include <memory>  // std::shared_ptr;
#include <type_traits>
#include <string>
#include <string_view>

// it would be a lot simpler (and idiomatic) to implement string-conversion
// logic in terms of c++ 20's std::format machinery, but this was one of the
// last C++ 20 features that compilers implemented
//
// There's an added wrinkle that if we use std::format we would really like
// to use machinery retroactively introduced to the C++ 20 standard in a "defect
// report" called P2508R1 (A link to this report can be found here:
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p2508r1.html)
// I'm pretty sure this isn't a big deal -- I think this report was written
// before std::format was implemented in most cases (so if an implementation
// supports std::format, it probably includes this retroactive change)
//
// -> to support this we need to require minimum compiler versions based on when
//    the associated standard library added support
//    -> gcc 13: libstdc++ first added support for std::format in that release
//       https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#status.iso.2020
//    -> clang 17: libcxx made std::format non-experimental in that release
//       https://releases.llvm.org/17.0.1/projects/libcxx/docs/ReleaseNotes.html
//    -> (I think clang15 actually experimentally supported it)
//    -> I think apple-clang supports it starting with xcode 15.3
//       (https://developer.apple.com/documentation/Xcode-Release-Notes/xcode-15_3-release-notes)
// -> a more informative table can be found here: (see the row for P2508R1):
//    https://en.cppreference.com/cpp/compiler_support/23
//    -> don't be surprised that the table technically describes C++ 23 features
//       -- P2508r1 is DEFINITELY applicable for C++20
//    -> interestingly that table claims that xcode 14.0.3 (if you hover over
//       the entry it says 14.3) supports the feature even though it's not in
//       the release notes
//
// ASIDE: from reading the standard, it seems like this *SHOULD* all be as easy
//        as checking whether the __cpp_lib_format macro is defined and if it
//        has a value >=202207L, but clang's libcxx runtime library doesn't
//        define this macro at all in certain versions (e.g. clang++ 17 and 18)
//        because some edge cases aren't fully implemented (this probably also
//        applies to apple-clang)

#include "./config.hpp"
#include "./error_detail.hpp"
#include "./formatf.hpp"

namespace GRIMPL_NAMESPACE_DECL {

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
  /// The additional context is specified as a printf-style message
  ///
  /// @param format The format-string. This **MUST** be a string literal.
  /// @param ... optional arguments to be formatted
  /// @return Returns a reference to `this` for convenience (e.g. to facillitate
  ///     chaining of operations)
  ///
  /// @warning
  /// Passing a non-literal string as @p format introduces undefined behavior.
  /// Compilers should generally warn about this
  ///
  /// @note
  /// The proper way to wrap an error object `err` in a context message encoded
  /// in a `std::string` object called `s` (this object may have been
  /// dynamically constructed) is to call `err.contextf("%s", s.data());`
  ///
  /// Implementaton Note
  /// ------------------
  /// For less experienced c++ developers: ``[[gnu::format(printf, 2, 3)]]``
  /// is a compiler-specific attribute that instructs gcc (or clang++) to check
  /// at compile-time that the arguments are consistent with the printf style
  /// format string `format`. It **SHOULD** also check that `format` is a
  /// string-literal. If a compiler doesn't recognize the attribute, it simply
  /// ignores (that's mandated by the C++ standard)
  [[gnu::format(printf, 2, 3)]] Error& contextf(const char* format, ...) {
    // the gnu::format attribute is told that s is argument 2 (rather than arg
    // 1) and that the first variadic argument is argument 3 (rather than arg 2)
    // because `this` is an implicit 1st argument for non-static member
    // functions like this
    std::va_list vlist;
    va_start(vlist, format);
    std::string msg = vstr_formatf(format, vlist);
    va_end(vlist);
    Error out;
    return context_helper_(ErrImpl_("", std::move(msg), nullptr));
  }

  /// @brief wrap the existing err info in additional context
  ///
  /// This method exists to allow
  ///
  /// @return Returns a reference to `this` for convenience (e.g. to facillitate
  ///     chaining of operations)
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

  /// @brief get a string representation of the error
  std::string to_string() const;

  // factory methods (we may add more in the future!)
  // ================================================

  /// @brief Construct an error from a printf-style message
  ///
  /// @param format The format-string. This **MUST** be a string literal.
  /// @param ... optional arguments to be formatted
  /// @returns An error object
  ///
  /// @warning
  /// Passing a non-literal string as @p format introduces undefined behavior.
  /// Compilers should generally warn about this
  ///
  /// @note
  /// The proper way to create an error object encoding a message copied from
  /// a std::string `s` (this object may have been dynamically constructed) is
  /// to call `Error::msgf("%s", s.c_str());`
  ///
  /// Implementaton Note
  /// ------------------
  /// For less experienced c++ developers: ``[[gnu::format(printf, 1, 2)]]``
  /// is a compiler-specific attribute that instructs gcc (or clang++) to check
  /// at compile-time that the arguments are consistent with the printf style
  /// format string `format`. It **SHOULD** also check that `format` is a
  /// string-literal. If a compiler doesn't recognize the attribute, it simply
  /// ignores (that's mandated by the C++ standard)
  [[gnu::format(printf, 1, 2)]] static Error msgf(const char* format, ...) {
    // theoretically, we could iterate through characters in s. If we encounter
    // any occurrence of `%` other than `%%`, then we continue to implement
    // the current behavior. Otherwise, we could just treat it as a literal
    std::va_list vlist;
    va_start(vlist, format);
    std::string msg = vstr_formatf(format, vlist);
    va_end(vlist);
    Error out;
    out.impl_ = std::make_shared<ErrImpl_>("", std::move(msg), nullptr);
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