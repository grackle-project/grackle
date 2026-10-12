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

namespace GRIMPL_NAMESPACE_DECL {

Expected<SpInitializeInfo, Error> infer_species_maps(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info, const NuclideModel& nuclide_model) {
  SpeciesBuilder builder;

  // configure the primordial species
  // ================================
  if ((primordial_chemistry == 0) &&
      (metal_chemistry > 0 || grain_info != nullptr)) {
    return Unexpected(Error::msgf(
        "can't have metal or grain species if primordial_chemistry==0"));
  }

  // a case could be made for introducing a method for introducing a method
  // (perhaps called `SpeciesBuilder::primordial_from_primordial_chem`) that
  // implements this logic for adding the primordial species based on the value
  // of the primordial_chemistry parameter
  if (primordial_chemistry >= 1) {
    builder.chemical<SpKind::PRIMORDIAL>(SpecialCaseSpecies::e)
        .chemical<SpKind::PRIMORDIAL>("HI")
        .chemical<SpKind::PRIMORDIAL>("HII")
        .chemical<SpKind::PRIMORDIAL>("HeI")
        .chemical<SpKind::PRIMORDIAL>("HeII")
        .chemical<SpKind::PRIMORDIAL>("HeIII");
  }
  if (primordial_chemistry >= 2) {
    builder.chemical<SpKind::PRIMORDIAL>(SpecialCaseSpecies::HM)
        .chemical<SpKind::PRIMORDIAL>("H2I")
        .chemical<SpKind::PRIMORDIAL>("H2II");
  }
  if (primordial_chemistry >= 3) {
    builder.chemical<SpKind::PRIMORDIAL>("DI")
        .chemical<SpKind::PRIMORDIAL>("DII")
        .chemical<SpKind::PRIMORDIAL>("HDI");
  }
  if (primordial_chemistry >= 4) {
    builder.chemical<SpKind::PRIMORDIAL>(SpecialCaseSpecies::DM)
        .chemical<SpKind::PRIMORDIAL>("HDII")
        .chemical<SpKind::PRIMORDIAL>("HeHII");
  }
  if (primordial_chemistry > 4 || primordial_chemistry < 0) {
    return Unexpected(Error::msgf("invalid primordial_chemsitry value: %d",
                                  primordial_chemistry));
  }

  // configure metal species pertaining to the metal_chemistry flag
  // ==============================================================
  // a case could be made for introducing a method for introducing a method
  // (perhaps called `SpeciesBuilder::metal_from_metal_chem`) that
  // implements this logic for adding the primordial species based on the value
  // of the metal_chemistry parameter
  if (metal_chemistry == 0) {
    // do nothing!
  } else if (metal_chemistry == 1) {
    builder.chemical<SpKind::METAL>("CI")
        .chemical<SpKind::METAL>("CII")
        .chemical<SpKind::METAL>("COI")
        .chemical<SpKind::METAL>("CO2I")
        .chemical<SpKind::METAL>("OI")
        .chemical<SpKind::METAL>("OHI")
        .chemical<SpKind::METAL>("H2OI")
        .chemical<SpKind::METAL>("O2I")
        .chemical<SpKind::METAL>("SiI")
        .chemical<SpKind::METAL>("SiOI")
        .chemical<SpKind::METAL>("SiO2I")
        .chemical<SpKind::METAL>("CHI")
        .chemical<SpKind::METAL>("CH2I")
        .chemical<SpKind::METAL>("COII")
        .chemical<SpKind::METAL>("OII")
        .chemical<SpKind::METAL>("OHII")
        .chemical<SpKind::METAL>("H2OII")
        .chemical<SpKind::METAL>("H3OII")
        .chemical<SpKind::METAL>("O2II");
  } else {
    return Unexpected(
        Error::msgf("invalid metal_chemistry value: %d", metal_chemistry));
  }

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

  return builder.build(nuclide_model);
}

}  // namespace GRIMPL_NAMESPACE_DECL