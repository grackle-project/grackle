//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare and implement InplaceVector
///
//===----------------------------------------------------------------------===//

#ifndef SUPPORT_INPLACE_VECTOR
#define SUPPORT_INPLACE_VECTOR

#include <cstddef>  // std::size_t
#include <utility>  // std::forward

#include "support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief A crude, HIGHLY simplified backport of C++26's `std::inplace_vector`
///
/// Unlike a `std::vector` (which always embeds a pointer to heap-allocated
/// storage) the storage is always allocated in place (embeded within an
/// an `InplaceVector`) with a maximum capacity of `N` elements.
///
/// For our purposes, we simply implement this as a C-style array with an
/// associated length (we always default-initialize ``T``)
///
/// For more details about each method see
/// https://en.cppreference.com/cpp/container/inplace_vector
template <class T, std::size_t N>
class InplaceVector {
public:
  // these are type definitions that all standard library containers possess
  // (we are defining these to make use of googletest machinery that depends
  // on these typedefs)
  using value_type = T;
  using reference = T&;
  using const_reference = T&;
  using iterator = T*;
  using const_iterator = const T*;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;

private:
  std::size_t len_ = 0;
  T data_[N];

public:
  // the constructor and assignment operations are implicitly defined
  // (the default constructor builds an initially empty array)

  /// @brief Overloads array subscript operation
  ///
  /// @note For less-experienced C++ developers, this is like python's
  /// `__getitem__` magic method
  constexpr const T& operator[](std::size_t idx) const { return data_[idx]; }
  constexpr T& operator[](std::size_t idx) { return data_[idx]; }

  constexpr std::size_t size() const { return len_; }
  static constexpr std::size_t capacity() noexcept { return N; }
  constexpr bool empty() const { return len_ == 0; }

  // these are used for iterating
  constexpr iterator begin() { return data_; }
  constexpr const_iterator begin() const { return data_; }
  constexpr iterator end() { return begin() + len_; }
  constexpr const_iterator end() const { return begin() + len_; }

  /// @brief try to construct an element in-place at the back
  ///
  /// @note std::inplace_vector actually returns a std::optional<T&>, but we
  /// can't do that before C++26
  template <class... Args>
  constexpr bool try_emplace_back(Args&&... args) {
    bool has_space = len_ < N;
    if (has_space) {
      data_[len_++] = T(std::forward<Args>(args)...);
    }
    return has_space;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_INPLACE_VECTOR
