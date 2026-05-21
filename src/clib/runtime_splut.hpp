//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines logic pertaining to tracking information about the species lut at
/// runtime.
///
//===----------------------------------------------------------------------===//

#ifndef GRACKLE_RUNTIME_SPLUT_HPP
#define GRACKLE_RUNTIME_SPLUT_HPP

#include "dust/grain_species_info.hpp"
#include "support/config.hpp"
#include "support/PartMap.hpp"
#include "support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Describes the kinds of dynamically evolved species
///
/// This is to be used with an instance of @ref PartMap.
///
/// @note
/// FILTER is a temporary kind right now that corresponds to the gap in species
/// when dust_chemistry == 1. After we finish implementing the use of this type
/// throughout the codebase, the immediate goal is to get rid of this category.
///
/// Motivation
/// ----------
/// Before creating this construct, we have been tracking the mapping between
/// all dynamically evolved species names and corresponding indices within a
/// single compile-time LUT (lookup table). This approach basically just works
/// if Grackle's configurations can be placed on a 1D slider, where each
/// setting just adds more species to the network of the previous configuration
/// (they are kinda like nesting dolls).
///
/// The short/medium term plan is to move to a system where we track separate
/// compile-time LUTs for each of these kinds of species and use an instance of
/// @ref PartMap to tell us how these kind-specific LUTs are stitched togeth.
///
/// Longer term, the goal is to move to a more dynamical system for tracking
/// chemical networks (while retaining the hardcoded networks for the simplest
/// primordial chemistry networks to avoid losing speed). This is a stepping
/// stone towards that goal (it will help both systems coexist as we develop
/// the more dynamic system).
struct SpKind {
  enum {
    CHEMICAL = 0,  ///< primordial or metal species (includes electrons)
    FILLER = 1,    ///< right now, this is temporary
    DUST = 2,      ///< single grain species (or perhaps grain ensemble)
  };
};

namespace species_detail {

inline constexpr int primordial_count(int primordial_chemistry) {
  switch (primordial_chemistry) {
    case 0:
      return 0;
    case 1:
      return 6;
    case 2:
      return 9;
    case 3:
      return 12;
    case 4:
      return 15;
    default: {
      // earlier error checking should prevent this
      GR_INTERNAL_UNREACHABLE_ERROR();
    }
  }
}

}  // namespace species_detail

/// construct a @ref PartMap that describes the kinds of species that are
/// dynamically evolved. See @ref SpKind for more information about the species
/// kinds.
///
/// @par Future Thoughts
/// It's a little weird to pass a @ref GrainSpeciesInfo directly into this
/// function. We may want to revisit that choice in the future, (especially if
/// other dust models are implemented that affect the dynamical species fields)
inline PartMap make_species_kind_map(int primordial_chemistry,
                                     int metal_chemistry,
                                     const GrainSpeciesInfo* grain_info) {
  // species counts in various categories
  int n_primordial = species_detail::primordial_count(primordial_chemistry);
  int n_metal_normal = (metal_chemistry == 1) ? 19 : 0;
  int n_metal_accounting;  // these include Mg, Al, S, and Fe
  int n_filler;  // empty filler slots b/t chemical and dust grain species

  if (grain_info == nullptr) {
    n_metal_accounting = 0;
    n_filler = 0;
  } else if (grain_info->dust_species_parameter == 1) {
    n_metal_accounting = 1;  // we explicitly track Mg
    n_filler = 3;            // we track empty filler slots for Al, S, Fe
  } else {
    n_metal_accounting = 4;  // we explicitly track Mg, Al, S, Fe
    n_filler = 0;
  }

  // now, let's construct the actual partition maps
  int n_chemical_sp = n_primordial + n_metal_normal + n_metal_accounting;
  int n_grain_sp = (grain_info == nullptr) ? 0 : grain_info->n_species;

  constexpr int N_KINDS = 3;
  int species_kinds[N_KINDS] = {SpKind::CHEMICAL, SpKind::FILLER, SpKind::DUST};
  int sizes[N_KINDS] = {n_chemical_sp, n_filler, n_grain_sp};
  return new_PartMap(species_kinds, sizes, N_KINDS);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // GRACKLE_RUNTIME_SPLUT_HPP