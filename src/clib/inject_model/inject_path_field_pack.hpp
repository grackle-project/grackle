//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the InjectPathFieldPack type
///
//===----------------------------------------------------------------------===//
#ifndef INJECT_MODEL_INJECT_PATH_FIELD_PACK_HPP
#define INJECT_MODEL_INJECT_PATH_FIELD_PACK_HPP

#include "raw_data.hpp"  // grackle::impl::inj_model_input::N_Injection_Pathways

namespace grackle::impl {

/// A temporary type that is used to collect the injection pathway metal
/// density fields in the order that is assumed throughout Grackle
///
/// I have some ideas that will let us dispose of this type in the near future.
struct InjectPathFieldPack {
  /// Specifies the number of injection pathways that should have data in
  /// the fields member.
  int n_fields;

  /// holds pointers to the various injection model density fields
  const gr_float* const* fields;
};

/// Setup an InjectPathFieldPack instance
inline InjectPathFieldPack setup_InjectPathFieldPack(
    const chemistry_data* my_chem, const grackle_field_data* my_fields) {
  if ((my_chem->metal_chemistry > 0) && (my_chem->multi_metals == 1)) {
    grackle_field_data* my_fields_ = const_cast<grackle_field_data*>(my_fields);
    my_fields_->inject_pathway_metal_density[0] =
        my_fields_->local_ISM_metal_density;
    my_fields_->inject_pathway_metal_density[1] =
        my_fields_->ccsn13_metal_density;
    my_fields_->inject_pathway_metal_density[2] =
        my_fields_->ccsn20_metal_density;
    my_fields_->inject_pathway_metal_density[3] =
        my_fields_->ccsn25_metal_density;
    my_fields_->inject_pathway_metal_density[4] =
        my_fields_->ccsn30_metal_density;
    my_fields_->inject_pathway_metal_density[5] =
        my_fields_->fsn13_metal_density;
    my_fields_->inject_pathway_metal_density[6] =
        my_fields_->fsn15_metal_density;
    my_fields_->inject_pathway_metal_density[7] =
        my_fields_->fsn50_metal_density;
    my_fields_->inject_pathway_metal_density[8] =
        my_fields_->fsn80_metal_density;
    my_fields_->inject_pathway_metal_density[9] =
        my_fields_->pisn170_metal_density;
    my_fields_->inject_pathway_metal_density[10] =
        my_fields_->pisn200_metal_density;
    my_fields_->inject_pathway_metal_density[11] =
        my_fields_->y19_metal_density;

    return InjectPathFieldPack{inj_model_input::N_Injection_Pathways,
                               my_fields->inject_pathway_metal_density};
  }

  InjectPathFieldPack out{0, nullptr};
  if ((my_chem->metal_chemistry > 0) && (my_chem->multi_metals == 0)) {
    out.n_fields = 1;
    out.fields = &my_fields->metal_density;
  }

  return out;
}

}  // namespace grackle::impl

#endif /* INJECT_MODEL_INJECT_PATH_FIELD_PACK_HPP */
