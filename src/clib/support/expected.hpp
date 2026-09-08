//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// A lightweight backport of std::expected from c++ 20
///
//===----------------------------------------------------------------------===//

#ifndef SUPPORT_EXPECTED_HPP
#define SUPPORT_EXPECTED_HPP

#include <concepts>
#include <type_traits>
#include <variant>

#include "config.hpp"

// this file implements a minimal subset of std::expected (& std::unexpected)
// from C++ 23 (when we eventually transition to C++23, the plan is to start
// using the machinery provided by the standard library)
//
// This is intended to help us with internal management of errors during
// grackle's initialization:
// -> Under our current error handling strategy:
//    -> a function reports presence of an error to caller by either:
//       1. directly return an error code from the function (in this case, the
//          output computed by a function is often reported by modifying an
//          argument or pointer)
//       2. returned objects have an explicit error-state
//       3. returned objects are wrapped in std::optional
//    -> we provide informative error messages, we generally print out the
//       error message immediately as they occurs
// -> Issues with this strategy:
//    -> approaches 1 & 2 are somewhat undesirable for reporting errors to a
//       caller of a function because it means that various objects need to
//       have additional error-states
//    -> The practice of eagerly printing error-messages is actually quite
//       annoying in cases where you want to call a function and check for a
//       particular failure-mode. This comes up more than you might think,
//       especially in parsing code. It also makes it much harder to write unit
//       tests that explicitly check that certain error cases are handled
// -> Having a function return an Expected object is just like returning a
//    std::optional, except you provide error information when there is a
//    problem, instead of returning an empty container. When coupled with a
//    well-defined Error class, this mechanism allows us to defer the act of
//    reporting detailed errors
//
// For more information about the std::expected API, see:
// https://en.cppreference.com/cpp/utility/expected

namespace GRIMPL_NAMESPACE_DECL {

// This facillitate creation of an Expected object holding an error
template <class E>
  requires std::copyable<E>
class Unexpected {
  E err_;

public:
  constexpr explicit Unexpected(const E& err) : err_{err} {};
  constexpr explicit Unexpected(E&& err) : err_{err} {};
  constexpr E& error() { return err_; }
  constexpr const E& error() const { return err_; }
};

// the standard library declares the following deduction guide (but I don't
// think we actually need it)
//     template<class E> Unexpected(E) -> Unexpected<E>;

template <class T, class E>
  requires std::copyable<E>
class Expected {
  static_assert(!std::is_void_v<T>, "T isn't currently allowed to be void");

  std::variant<T, E> u_;  // a type-safe union holding either a value or error

public:
  // the following allow implicit casts from T and Unexpected<E> objects
  // NOLINTBEGIN(google-explicit-constructor)
  constexpr Expected(T v) : u_(std::in_place_index<0>, v) {}
  constexpr Expected(Unexpected<E> e)
      : u_(std::in_place_index<1>, std::move(e.error())) {}
  // NOLINTEND(google-explicit-constructor)

  // define basic operations (if supported by the type T)
  constexpr Expected()
    requires std::default_initializable<T>
  = default;
  constexpr Expected(const Expected&)
    requires std::copy_constructible<T>
  = default;
  constexpr Expected& operator=(const Expected&)
    requires std::assignable_from<T&, const T&>
  = default;

  // the following pair of operations are skipped for now (they're tricky!)
  Expected(Expected&&) = delete;
  Expected& operator=(Expected&&) = delete;

  constexpr bool has_value() const { return u_.index() == 0; }
  constexpr explicit operator bool() const { return u_.index() == 0; }

  constexpr T& value() { return std::get<0>(u_); }
  constexpr const T& value() const { return std::get<0>(u_); }
  constexpr E& error() { return std::get<1>(u_); }
  constexpr const E& error() const { return std::get<1>(u_); }

  constexpr T value_or(T dflt_val) const
    requires std::copy_constructible<T>
  {
    return has_value() ? value() : dflt_val;
  }

  // usage of the following functions trigger undefined behavior if `this`
  // doesn't contain a value (the fact that the current implementation aborts
  // the program is an implementation detail that can/will change)
  constexpr T& operator*() noexcept { return std::get<0>(u_); }
  constexpr const T& operator*() const noexcept { return std::get<0>(u_); }
  constexpr T* operator->() noexcept { return &std::get<0>(u_); }
  constexpr const T* operator->() const noexcept { return &std::get<0>(u_); }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_EXPECTED_HPP