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

#include "grackle.h"
#include "dust/grain_species_info.hpp"
#include "grain_metal_inject_pathways.hpp"
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

  /// @name rxn-rate-tables-group
  /// These members are all used to describe the tables of "ordinary"
  /// collisional reaction rates. In this context, "ordinary" means that the
  /// that are each interpolated with respect to the ln(T_gas) grid specified
  /// by the input parameters (here "ln" stands for natural log)
  ///@{
  /// holds the collision rate tables
  grackle::impl::CollisionalRxnRateCollection* kcol_rate_tables;
  /// a list of the indices that are actualy used from kcol_rate_tables in the
  /// current calculation
  int* used_kcol_rate_indices;
  /// length of used_kcol_rate_indices
  int n_kcol_rate_indices;
  ///@}

  /// tracks the grid of values used for interpolating grain-species specific
  /// coefficients for computing rates of H2 formation on dust grains. (At the
  /// time of writing, chemistry_data_storage::h2dustS &
  /// chemistry_data_storage::h2dustC hold the interpolated values)
  ///
  /// In more detail, this represents a 2D grid, where ln(Tdust) varies along
  /// axis 0 and ln(Tgas) varies along axis 1. Below we highlight a few notes
  /// about the properties of these grid:
  /// - Importantly, the set of ln(Tgas) values is very same grid of ln(Tgas)
  ///   values that is used for 1D interpolation of the majority of other rates
  /// - For context, all of other uses of this ln(Tgas) grid have historically
  ///   avoided the explicit construction of the ln(Tgas) grid
  /// - Furthermore, interpolation of the h2dust table (which holds
  ///   coefficients for computing H2 formation on dust grains) has been
  ///   interpolated on the same 2D [ln(Tdust), ln(Tgas)] grid. We have also
  ///   historically avoided the explicit construction of the grid while
  ///   interpolating h2dust).
  /// - In the future we should give some thought to whether we really want to
  ///   track this grid. An argument could be made for breaking this up and
  ///   explicitly tracking a grid of ln(Tgas) values and a grid of ln(Tdust)
  ///   values (for debugging purposes).
  gr_interp_grid_props h2dust_grain_interp_props;

  /// Tracks basic information about each relevant grain species
  ///
  /// > [!note]
  /// > A case could be made that we may not want to directly store a
  /// > GrainSpeciesInfo instances directly inside of this struct (since it
  /// > contains some extra information that is unnecessary during the
  /// > calculations). An alternative would be to briefly initialize an
  /// > instance during setup and then repack the data.
  grackle::impl::GrainSpeciesInfo* grain_species_info;

  /// Tracks metal and grain yields for each modeled injection pathway as well
  /// as other grain properties
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props;
};

#endif /* OPAQUE_STORAGE_HPP */
