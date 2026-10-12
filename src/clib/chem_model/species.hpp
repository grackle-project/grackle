//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// declare information pertaining to dynamically evolved species
///
//===----------------------------------------------------------------------===//

#ifndef CHEM_MODEL_SPECIES_HPP
#define CHEM_MODEL_SPECIES_HPP

#include "../dust/grain_species_info.hpp"
#include "./chemsp_props.hpp"
#include "../support/config.hpp"
#include "../support/error.hpp"
#include "../support/expected.hpp"
#include "../support/FrozenKeyIdxBiMap.hpp"
#include "../support/PartMap.hpp"
#include "chem_model/nuclide_model.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Describes the kinds of dynamically evolved species
///
/// This is to be used with an instance of @ref PartMap.
///
/// @note
/// FILLER_A and FILLER_B are temporary kinds right now that corresponds to the
/// gap in species when dust_chemistry == 1. After we finish phasing out SpLUT
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
/// compile-time LUTs for each of these kinds of species (see
/// @ref PrimordialSpLUT, @ref MetalSpLUT, and @ref DustSpLUT) and use an
/// instance of @ref PartMap to tell us how these kind-specific LUTs are
/// stitched together.
///
/// Longer term, the goal is to move to a more dynamical system for tracking
/// chemical networks (while retaining the hardcoded networks for the simplest
/// primordial chemistry networks to avoid losing speed). This is a stepping
/// stone towards that goal (it will help both systems coexist as we develop
/// the more dynamic system). When we get to that point, we may want to
/// consider replacing @ref SpKind::PRIMORDIAL and @ref SpKind::METAL with a
/// single generic kind that includes both categories (perhaps ``CHEMICAL``?)
enum class SpKind {
  PRIMORDIAL,  ///< a primordial species (includes electrons)
  METAL,       ///< a metal species
  DUST,        ///< a single grain species (or perhaps grain ensemble)
  FILLER_A,    ///< right now, this is temporary
  FILLER_B     ///< right now, this is temporary
};

/// @brief aggregates queryable information about the species model
///
/// The species model describes dynamically evolved species. Importantly, each
/// species index is associated with a species name.
///
/// Whereas @ref kind_map has applications during the core calculations, the
/// other data members are primarily intended to help set up a Grackle
/// configuration
struct SpInitializeInfo {
  /// @brief describes the kinds of species that are dynamically evolved
  ///
  /// Specifically, this maps between species index and the species kind. See
  /// @ref SpKind for more information about the species kinds.
  PartMap<SpKind> kind_map;

  /// @brief maps between species names and indices
  ///
  /// @note
  /// As long the codebase continues to make use of @ref SpLUT there may be
  /// ranges of indices (denoted by the @ref SpKind::FILLER_A and
  /// @ref SpKind::FILLER_B partitions in @ref kind_map) that hold garbage
  /// values. This is necessary to ensure consistency of index-name mappings
  /// with @ref SpLUT in certain cases when grackle is configured to not use
  /// all possible primordial species or not use all possible metal species
  FrozenKeyIdxBiMap name_map;

  /// @brief holds a ChemSpeciesProps for chemical species
  ///
  /// The value at index `i` describes the species associated with that index
  /// by @ref name_map. Importantly, accessing an index that isn't associated
  /// with a @ref SpKind::PRIMORDIAL or @ref SpKind::METAL (see @ref kind_map)
  /// will trigger undefined behavior.
  ///
  /// @note
  /// At the time of writing this does **NOT** specify entries for dust species
  ChemSpeciesPropsVec chem_species_props;

  /// @brief associates properties with each "nuclide symbol"
  NuclideModel nuclide_model;
};

/// @brief infer the species properties
///
/// Under the hood this makes use of @ref SpeciesBuilder. However, this function
/// takes care to ensure that the mapping between a species name and its index
/// remain consistent with compile-time lookup tables.
///
/// @par Future Thoughts
/// It's a little weird to pass a @ref GrainSpeciesInfo directly into this
/// function. We may want to revisit that choice in the future, (especially if
/// other dust models are implemented that affect the dynamical species fields)
Expected<SpInitializeInfo, Error> infer_species_maps(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEMISTRY_SPECIES_HPP