//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of @ref make_consistent
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#ifndef MAKE_CONSISTENT_HPP
#define MAKE_CONSISTENT_HPP

#include "grackle.h"
#include "field_adaptor.hpp"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Enforces consistency of species fields with expected total abundances
/// for indices in the specified index-range
///
/// The function produces undefined behavior if called for a grackle
/// configuration that doesn't dynamically evolve any species densities
///
/// @note
/// This logic encapsulated by this function needs to be **HEAVILY** reworked
/// (e.g. the way we handle dust growth doesn't make much sense)
///
/// @param[in] idx_range index range for which the operation is performed
/// @param[in] dom  equivalent to `nH_cgs / rhoH_proper_codeU`.  Multiplying
///     this quantity by Hydrogen mass density, in (proper) density code units,
///     will produce the Hydrogen number density with units of cm^-3
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] inject_pathway_props holds data about the modelled injection
///     pathways for all metals and grain species
/// @param[in] my_fields specifies all field data
/// @param[inout] sp_densities Specifies densities of the various species that
///     Grackle evolves in a format that allows the values to be accessed with
///     the index lookup table. Wherever possible, data should be be accessed
///     through this argument, rather than with @p my_fields
void make_consistent(const IndexRange& idx_range, double dom,
                     const chemistry_data* my_chemistry,
                     const GrainMetalInjectPathways* inject_pathway_props,
                     grackle_field_data* my_fields,
                     SpeciesMultiView<gr_float> sp_densities);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* MAKE_CONSISTENT_HPP */
