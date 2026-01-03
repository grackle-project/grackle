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

#include "raw_data.hpp"  // grackle::impl::inj_model_input::N_Injection_Pathways

namespace grackle::impl {

/// this is a hacky short-term helper function
///
/// @todo
/// Get rid of this function. This will be necessary in the future in order to
/// start reading injection pathway details from user-specified HDF5 files.
inline int get_n_inject_pathway_density_ptrs(const chemistry_data* my_chem) {
  bool nonzero = (my_chem->metal_chemistry > 0) && (my_chem->multi_metals > 0);
  return (nonzero) ? grackle::impl::inj_model_input::N_Injection_Pathways : 0;
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
