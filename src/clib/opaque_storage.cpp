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

  // since h2dust_grain_interp_props isn't a pointer, there is nothing more to
  // allocate right here

  GRIMPL_NS::free_interp_grid_(&LH2);
  GRIMPL_NS::free_interp_grid_(&LHD);

  // we deal with freeing other interp grids inside of
  // free_misc_species_cool_rates

  GRIMPL_NS::free_interp_grid_(&alphap);

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
    // delete contents of registry
    grackle::impl::ratequery::drop_Registry(registry);
    // delete registry, itself
    delete registry;
  }
}