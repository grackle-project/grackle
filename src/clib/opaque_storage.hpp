//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define the struct for holding the opaque storage
///
//===----------------------------------------------------------------------===//

#ifndef OPAQUE_STORAGE_HPP
#define OPAQUE_STORAGE_HPP

#include "internal_types.hpp"

/// a struct that used to wrap some private storage details
///
/// The short-term goal is to have @ref chemistry_data_storage struct hold an
/// opaque pointer to an instance of this struct. In more detail:
/// - the public headers will have a forward declaration of this struct
///   (allowing us to declare pointers to `struct gr_opaque_storage`)
/// - the struct is **ONLY** fully defined in a private header
/// - an opaque pointer can only be allocated/destroyed/dereferenced by
///   routines defined by the library that defines the underlying struct
///
/// Opaque pointers are a common tool for defining C APIs. A famous example is
/// the `FILE*` pointer used in the standard C library. Opaque pointers are
/// primarily used as a mechanism for hiding internal implementation details.
/// This provides 2 relevant concrete benefits:
/// 1. It facilitates internal refactoring (and adding/removing features) while
///    maintaining a stable API and ABI
/// 2. It allows us to expose a C API, even if the internals of the struct are
///    implemented in a different language like C++, Fortran, Cython, Rust, ...
///
/// In an ideal world, the Grackle API would treat the
/// @ref chemistry_data_storage struct as an opaque pointer. However, that
/// change can only be made in a Major Release (as it would break any existing
/// code that allocates the struct on the stack). Furthermore, if we'll need to
/// "accessor functions" for every detail we want to expose about the internal
/// state. The gr_opaque_storage struct can be used to help us gradually
/// transition towards this case
struct gr_opaque_storage {
  // in the future, we may want refactor the following set of members into
  // a separate datatype that takes full responsibility for "normal"
  // collisional rates

  /// holds the collision rate tables
  grackle::impl::CollisionalRxnRateCollection* kcol_rate_tables;

  /// a list of the indices that are actualy used from kcol_rate_tables in the
  /// current calculation
  int* used_kcol_rate_indices;
  /// length of used_kcol_rate_indices
  int n_kcol_rate_indices;
};

#endif /* OPAQUE_STORAGE_HPP */
