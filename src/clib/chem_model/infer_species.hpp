//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// <add a short description>
///
//===----------------------------------------------------------------------===//

// TODO: maybe rename this file after the primary function declared here?

#ifndef CHEM_MODEL_SPECIES_HPP
#define CHEM_MODEL_SPECIES_HPP

#include "../dust/grain_species_info.hpp"
#include "../support/config.hpp"
#include "../support/error.hpp"
#include "../support/expected.hpp"
#include "../support/FrozenKeyIdxBiMap.hpp"
#include "../support/PartMap.hpp"

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

/// @brief aggregates objects build by @ref infer_species_maps
struct SpInitializeInfo {
  PartMap<SpKind> kind_map;
  FrozenKeyIdxBiMap name_map;
};

/// @brief infer the species properties
///
/// In more detail, this constructs:
/// 1. a @ref PartMap that describes the kinds of species that are dynamically
///    evolved dynamically evolved. See @ref SpKind for more information about
///    the species kinds.
/// 2. a @ref FrozenKeyIdxBiMap mapping species names and indices. This is
///    primarily intended to help with setting things up. Note: for as long the
///    codebase continues to make use of @ref SpLUT there may be ranges of
///    indices (denoted by the @ref SpKind::FILLER_A and @ref SpKind::FILLER_B
///    partitions in the PartMap) that hold garbage values.
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