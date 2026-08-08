//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define the methods associated with gr_opaque_storage
///
//===----------------------------------------------------------------------===//

#include "opaque_storage.hpp"
#include "support/config.hpp"

gr_opaque_storage::~gr_opaque_storage() {
  if (kcol_rate_tables != nullptr) {
    // delete contents of kcol_rate_tables
    drop_CollisionalRxnRateCollection(kcol_rate_tables);
    // delete kcol_rate_tables, itself
    delete kcol_rate_tables;
  }

  if (used_kcol_rate_indices != nullptr) {
    // since used_kcol_rate_indices are just integers, we can directly
    // deallocate them
    delete[] used_kcol_rate_indices;
  }

  // all of the InterpGrid & InterpGridProps data members have destructors that
  // handle their deallocation

  if (grain_species_info != nullptr) {
    // delete contents of grain_species_info
    grackle::impl::drop_GrainSpeciesInfo(grain_species_info);
    // delete grain_species_info, itself
    delete grain_species_info;
  }

  if (inject_pathway_props != nullptr) {
    // delete contents of inject_pathway_props
    grackle::impl::drop_GrainMetalInjectPathways(inject_pathway_props);
    // delete inject_pathway_props, itself
    delete inject_pathway_props;
  }

  if (registry != nullptr) {
    // delete will trigger the Registry's destructor
    delete registry;
  }
}
