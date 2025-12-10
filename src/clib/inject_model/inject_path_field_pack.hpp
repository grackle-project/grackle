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
  const gr_float* fields[inj_model_input::N_Injection_Pathways];
};

/// Setup an InjectPathFieldPack instance
inline InjectPathFieldPack setup_InjectPathFieldPack(
    const chemistry_data* my_chem, const grackle_field_data* my_fields) {
  if ((my_chem->metal_chemistry > 0) && (my_chem->multi_metals == 1)) {
    return InjectPathFieldPack{{
        my_fields->local_ISM_metal_density,
        my_fields->ccsn13_metal_density,
        my_fields->ccsn20_metal_density,
        my_fields->ccsn25_metal_density,
        my_fields->ccsn30_metal_density,
        my_fields->fsn13_metal_density,
        my_fields->fsn15_metal_density,
        my_fields->fsn50_metal_density,
        my_fields->fsn80_metal_density,
        my_fields->pisn170_metal_density,
        my_fields->pisn200_metal_density,
        my_fields->y19_metal_density,
    }};
  }

  InjectPathFieldPack out;
  for (int i = 0; i < inj_model_input::N_Injection_Pathways; i++) {
    out.fields[i] = nullptr;
  }

  if ((my_chem->metal_chemistry > 0) && (my_chem->multi_metals == 0)) {
    out.fields[my_chem->metal_abundances] = my_fields->metal_density;
  }
  return out;
}

}  // namespace grackle::impl

#endif /* INJECT_MODEL_INJECT_PATH_FIELD_PACK_HPP */
