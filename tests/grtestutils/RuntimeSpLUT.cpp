//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements RuntimeSpLUT
///
//===----------------------------------------------------------------------===//

#include "RuntimeSpLUT.hpp"
namespace grtest {
std::optional<RuntimeSpLUT> RuntimeSpLUT::create(
    int primordial_chemistry, int metal_chemistry,
    const GRIMPL_NS::GrainSpeciesInfo* grain_info) {
  GRIMPL_NS::PartitionedKeyIdxBiMapRslt intermediate =
      grackle::impl::make_runtime_splut(primordial_chemistry, metal_chemistry,
                                        grain_info);

  if (intermediate.success) {
    RuntimeSpLUT tmp;
    tmp.partition_map_ = intermediate.part_map;
    tmp.bimap_ = intermediate.bimap;
    return {tmp};
  }
  return std::nullopt;
}

// this is a little silly (but it's useful)
std::optional<RuntimeSpLUT> RuntimeSpLUT::create(
    const chemistry_data* my_chemistry) {
  int dust_species_param = my_chemistry->dust_species;
  int n_grain_species = GRIMPL_NS::get_n_grain_species(dust_species_param);
  if (n_grain_species < 0) {
    return std::nullopt;
  } else if (n_grain_species == 0) {
    return RuntimeSpLUT::create(my_chemistry->primordial_chemistry,
                                my_chemistry->metal_chemistry, nullptr);
  } else {
    GRIMPL_NS::GrainSpeciesInfo tmp =
        GRIMPL_NS::new_GrainSpeciesInfo(dust_species_param);
    std::optional<RuntimeSpLUT> out =
        RuntimeSpLUT::create(my_chemistry->primordial_chemistry,
                             my_chemistry->metal_chemistry, &tmp);
    GRIMPL_NS::drop_GrainSpeciesInfo(&tmp);
    return out;
  }
}
}  // namespace grtest