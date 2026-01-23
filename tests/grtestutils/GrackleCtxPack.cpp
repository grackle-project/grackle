//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Help implement GrackleCtxPack
///
//===----------------------------------------------------------------------===//

#include "GrackleCtxPack.hpp"

namespace grtest {

std::pair<GrackleCtxPack, Status> GrackleCtxPack::create(
    const SimpleConfPreset& preset) {
  // allocate chemistry_data and set the defaults
  std::unique_ptr<chemistry_data> my_chem(new chemistry_data);
  if (local_initialize_chemistry_parameters(my_chem.get()) != GR_SUCCESS) {
    return {GrackleCtxPack(),
            error::Adhoc("initialize_chemistry_parameters failed")};
  }

  // lookup the parameters associated with the preset
  std::pair<std::vector<ParamPair>, Status> chem_preset_rslt =
      get_chem_preset_vals_(preset.chemistry);
  if (!chem_preset_rslt.second.is_ok()) {
    return {GrackleCtxPack(), chem_preset_rslt.second};
  }
  const std::vector<ParamPair>& params = chem_preset_rslt.first;

  // update my_chem with values from the preset
  param_detail::StrAllocTracker str_allocs;
  Status status =
      set_params(params.begin(), params.end(), *my_chem, &str_allocs);
  if (!status.is_ok()) {
    return {GrackleCtxPack(), status};
  }

  // set up chemistry_data_storage
  code_units initial_unit = setup_initial_unit(preset.unit);
  std::unique_ptr<chemistry_data_storage> my_rates(new chemistry_data_storage);
  if (local_initialize_chemistry_data(my_chem.get(), my_rates.get(),
                                      &initial_unit) != GR_SUCCESS) {
    return {GrackleCtxPack(), error::Adhoc("initialize_chemistry_data failed")};
  }

  std::pair<GrackleCtxPack, Status> out{GrackleCtxPack(), OkStatus()};
  out.first.initial_units_ = initial_unit;
  out.first.str_allocs_ = std::move(str_allocs);
  out.first.my_chemistry_ = std::move(my_chem);
  out.first.my_rates_ = std::move(my_rates);
  return out;
}

}  // namespace grtest
