//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines miscellaneous machinery. In the long term, we can hopefully remove
/// this machinery
///
//===----------------------------------------------------------------------===//
#ifndef INJECT_MODEL_MISC_HPP
#define INJECT_MODEL_MISC_HPP

#include "opaque_storage.hpp"
#include "raw_data.hpp"  // grackle::impl::inj_model_input::N_Injection_Pathways
#include "status_reporting.h"

namespace grackle::impl {

/// this is a hacky short-term helper function
///
/// @note
/// Because (at the time of writing), when `metal_chemistry > 0` and
/// `multi_metals == 0`, Grackle simply uses the `metal_density` field rather
/// than a pointer from `my_fields->inject_pathway_metal_density`, this function
/// returns 0 in that case.
inline int get_n_inject_pathway_density_ptrs(
    const chemistry_data_storage* my_rates) {
  const gr_opaque_storage* tmp = my_rates->opaque_storage;
  GR_INTERNAL_REQUIRE(tmp != nullptr, "sanity-check failed!");

  if (tmp->inject_pathway_props == nullptr) {
    return 0;
  }
  int n_pathways = tmp->inject_pathway_props->n_pathways;
  // at the time of writing, the only way `n_pathways` can have a value of 1 is
  // when `metal_chemistry > 0 && multi_metals == 0`
  // -> as we note in the docstring, we return 0 in this case
  return (n_pathways == 1) ? 0 : n_pathways;
}

/// Retrieves a pointer that corresponds to the array of user-specified fields
/// holding the cumulative metal density injected by each relevant injection
/// pathway
///
/// @note
/// At the time of writing, when `metal_chemistry > 0 && multi_metals == 0`,
/// Grackle assumes that 100% of the density within `my_fields->metal_density`
/// is injected by the single pathway. In this scenario, we return a pointer
/// `my_fields->metal_density`.
/// - theoretically, this could introduce slightly unexpected behavior if a
///   routine tried to modify `my_fields->metal_density` and expected the
///   injection pathway metal density to remain unchanged. In practice, this
///   never comes up
/// - in the longer term we plan to remove this specialized behavior from
///   Grackle altogether (i.e. users will need to specify separate injection
///   pathway metal densities in all cases). It will lead to less surprising
///   behavior and will become important once we start reading injection
///   pathway details from user-specified HDF5 files
inline const gr_float* const* get_inject_pathway_metal_density(
    const chemistry_data* my_chem, grackle_field_data* my_fields) {
  if (my_chem->metal_chemistry <= 0) {
    return nullptr;
  } else if (my_chem->multi_metals == 0) {
    return &my_fields->metal_density;
  } else {
    return my_fields->inject_pathway_metal_density;
  }
}

}  // namespace grackle::impl

#endif  // INJECT_MODEL_MISC_HPP
