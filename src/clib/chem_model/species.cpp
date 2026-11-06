//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// define logic pertaining to dynamically evolved species
///
//===----------------------------------------------------------------------===//

#include "chem_model/species.hpp"
#include "chem_model/species_builder.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<SpInitializeInfo, Error> infer_species_maps(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info) {
  if ((primordial_chemistry == 0) &&
      (metal_chemistry > 0 || grain_info != nullptr)) {
    return Unexpected(Error::msgf(
        "can't have metal or grain species if primordial_chemistry==0"));
  }

  // configure primordial and metal species from the preset
  Expected<SpeciesBuilder, Error> tmp =
      SpeciesBuilder::from_preset(primordial_chemistry, metal_chemistry);
  if (!tmp.has_value()) {
    return Unexpected(tmp.error());
  }
  SpeciesBuilder builder = std::move(tmp).value();

  // perform configurations based dust chemistry
  // ===========================================
  // in the future, we may consider giving GrainSpeciesInfo a flag to perform
  // this logic itself... This might be useful as the number of dust models
  // continues to grow. However, we this would introduce a wrinkle:
  // -> how do we robustly check that other metal fields is already specified?
  //    - a simple solution: add a method to SpeciesBuilder that returns whether
  //      any metal species are present
  // -> More generally: how do we handle the case when there are multiple sets
  //    of metal fields?
  //    - I'm very hesitant to add the ability to start introspecting lists of
  //      fields in the builder (that's going to complicate the design)
  //    - perhaps we could make it possible to register a field more than once?
  //      If you think about all of the consequences of this approach, you'll
  //      realize that it's only tractable in a scenario where the modelled
  //      metal species (and their order) is modelled dynamically (this is a
  //      bit of a ways off, but not entirely crazy)
  if (grain_info != nullptr) {
    int dust_species_parameter = grain_info->dust_species_parameter();
    if (metal_chemistry == 0) {
      return Unexpected(
          Error::msg_literal("dust chemistry requires metal fields & "
                             "metal_chemistry is disabled"));
    }

    // append any extra metal species
    if (dust_species_parameter >= 1) {
      builder.chemical<SpKind::METAL>(SpecialCaseSpecies::Mg);
    }
    if (dust_species_parameter >= 2) {
      builder.chemical<SpKind::METAL>(SpecialCaseSpecies::Al)
          .chemical<SpKind::METAL>(SpecialCaseSpecies::S)
          .chemical<SpKind::METAL>(SpecialCaseSpecies::Fe);
    }

    // now register the dust species
    const FrozenKeyIdxBiMap& grainsp_name_map = grain_info->name_map();
    int n_grain_species = grain_info->n_species();
    for (int i = 0; i < n_grain_species; i++) {
      builder.dust(grainsp_name_map.inverse_find(i));
    }
  }

  return builder.build();
}

}  // namespace GRIMPL_NAMESPACE_DECL