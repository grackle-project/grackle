//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of make_consistent_g
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#ifndef MAKE_CONSISTENT_HPP
#define MAKE_CONSISTENT_HPP

#include "grackle.h"
#include "inject_model/grain_metal_inject_pathways.hpp"

namespace grackle::impl {

/// Enforces consistency of species fields with expected total abundances
///
/// @param[in] imetal specifies whether or not the caller provided a metal
///     field (0 = no, 1 = yes)
/// @param[in] dom
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] inject_pathway_props holds data about the modelled injection
///     pathways for all metals and grain species
/// @param[in] my_fields specifies the field data
void make_consistent(
    int imetal, double dom, chemistry_data* my_chemistry,
    const grackle::impl::GrainMetalInjectPathways* inject_pathway_props,
    grackle_field_data* my_fields);

}  // namespace grackle::impl

#endif /* MAKE_CONSISTENT_HPP */
