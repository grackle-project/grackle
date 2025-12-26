//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the SimpleVec type
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_SIMPLEVEC_HPP
#define SUPPORT_SIMPLEVEC_HPP

#include <limits>
// #include <type_traits>

#include "status_reporting.h"

namespace grackle::impl {

/// This is a **VERY** simplified version of std::vector that is intended to be
/// very C-like.
///
/// It is C-like in the sense that:
/// - it essentially acts like a struct with a bunch of associated functions
///   (the associated functions act like methods).
/// - **IMPORTANTLY:** the caller is responsible for explicitly calling the
///   destructor-function.
///
/// @par Motivation
/// This type is motivated by the fact that to dynamically build up objects,
/// you commonly need to construct arrays where the instances are not
/// well-known ahead of time. This comes up with enough frequency that it is
/// useful to define an abstraction for this data-structure
///
/// @todo
/// We should **STRONGLY** consider replacing this type with std::vector. In my
/// opinion, the primary "cost" is that std::vector is more "contagious."
/// Unlike the status quo where we can extract the underlying pointer (and take
/// ownership of it), storage can't be taken from a std::vector; you need to
/// either allocate a new pointer or continue carrying around the underlying
/// vector)
///
/// @par Justification for making this a class template
/// The alternatives to making this a template are *MUCH* less desirable:
/// 1. We could use a macro to define a version of this type and all associated
///    functions for each contained type. This is messy, and we would need to
///    come up with unique names for each version of the primary struct. In
///    practice, the template is doing this under the hood
/// 2. We could convert data from `T*` to `void**`. While this is doable, it
///    will require extra memory. Every time we push back a value, we would
///    need to copy that value into newly allocated memory and cast a pointer
///    to that memory to `void*`. For cleanup, you would need to cast each
///    `void*` back to the original type before calling `delete`
template <typename T>
struct SimpleVec {
  // declare the struct-members
  // - default-member initialization is used to ensure that this struct is in
  //   a valid state without calling an explicit constructor function
  int capacity = 0;
  int len = 0;
  T* data = nullptr;

  // in practice, the following logic prevents other structs from directly
  // embedding an instance of this type as a data member (a pointer to this
  // struct must be stored)
  // - There is always a risk of mistakes with a C-like API when it comes to
  //   dangling pointers. But, the SimpleVec_extract_ptr_and_make_empty
  //   function definitely amplifies the risk.
  // - By forcing the use of pointers to this type, we are mitigating this risk
  //   to an extent at the cost of an extra pointer indirection.
  // - (if we're willing to fully embrace C++, we could do a **LOT** better)
  SimpleVec() = default;
  SimpleVec(const SimpleVec&) = delete;
  SimpleVec(SimpleVec&&) = delete;
  SimpleVec& operator=(const SimpleVec&) = delete;
  SimpleVec& operator=(SimpleVec&&) = delete;
};

/// Deletes internal data within a vec (acts like a destructor)
///
/// As per usual, this does not try to directly deallocate the pointer itself
template <typename T>
void drop_SimpleVec(SimpleVec<T>* vec) {
  if (vec->data != nullptr) {
    delete[] vec->data;
  }
  vec->data = nullptr;
  vec->capacity = 0;
  vec->len = 0;
}

/// returns the internal data pointer (the caller takes ownership of it) and
/// modifies the provided vec so that it's equivalent to an empty vector
///
/// After this function executes, passing the argument to drop_SimpleVec or to
/// SimpleVec_push_back will **NOT** affect the extracted pointer.
template <typename T>
T* SimpleVec_extract_ptr_and_make_empty(SimpleVec<T>* vec) {
  T* ptr = vec->data;
  vec->data = nullptr;
  vec->capacity = 0;
  vec->len = 0;
  return ptr;
}

/// Deletes internal data within a vec
template <typename T>
void SimpleVec_push_back(SimpleVec<T>* vec, T value) {
  if (vec->capacity == 0) {
    vec->capacity = 5;
    vec->data = new T[vec->capacity];
  } else if (vec->capacity == vec->len) {
    int max_cap = std::numeric_limits<int>::max();
    GR_INTERNAL_REQUIRE(vec->capacity < max_cap, "should never happen!");
    int new_cap =
        ((max_cap / 2) <= vec->capacity) ? max_cap : vec->capacity * 2;
    T* new_data = new T[new_cap];
    for (int i = 0; i < vec->len; i++) {
      new_data[i] = vec->data[i];
    }
    delete[] vec->data;
    vec->data = new_data;
    vec->capacity = new_cap;
  }

  vec->data[vec->len] = value;
  ++(vec->len);
}

/// Deletes internal data within a vec
template <typename T>
int SimpleVec_len(const SimpleVec<T>* vec) {
  return vec->len;
}

}  // namespace grackle::impl

#endif  // SUPPORT_SIMPLEVEC_HPP
