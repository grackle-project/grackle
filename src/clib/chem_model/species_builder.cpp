
//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement logic for SpeciesBuilder
///
//===----------------------------------------------------------------------===//

#include <array>
#include <type_traits>
#include <vector>

#include "LUT.hpp"
#include "chem_model/infer_species.hpp"
#include "chem_model/species_builder.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// these are dummy names that won't collide with known species names
// -> NaSp stands for "Not a Species"
// -> this is a kludge until we remove all occurrences of SpLUT
static constexpr std::array<const char*, 40> dummy_names_ = {
    "_NaSp00", "_NaSp01", "_NaSp02", "_NaSp03", "_NaSp04", "_NaSp05", "_NaSp06",
    "_NaSp07", "_NaSp08", "_NaSp09", "_NaSp10", "_NaSp11", "_NaSp12", "_NaSp13",
    "_NaSp14", "_NaSp15", "_NaSp16", "_NaSp17", "_NaSp18", "_NaSp19", "_NaSp20",
    "_NaSp21", "_NaSp22", "_NaSp23", "_NaSp24", "_NaSp25", "_NaSp26", "_NaSp27",
    "_NaSp28", "_NaSp29", "_NaSp30", "_NaSp31", "_NaSp32", "_NaSp33", "_NaSp34",
    "_NaSp35", "_NaSp36", "_NaSp37", "_NaSp38", "_NaSp39"};

Expected<SpInitializeInfo, Error> SpeciesBuilder::build() const {
  if (primordial_vec_.empty() && metal_vec_.empty() && dust_vec_.empty()) {
    return SpInitializeInfo();  // <- return species maps that are empty
  }

  std::vector<SpKind> sp_kinds;
  std::vector<int> sp_kind_counts;
  std::vector<const char*> sp_names;

  // define a lambda function that infers the name of a species from either of
  // the types contained by a ChemEntryVariant_
  auto get_species_name = [](const auto& entry) -> const char* {
    using T = std::decay_t<decltype(entry)>;
    if constexpr (std::is_same_v<T, std::string>) {
      return entry.c_str();
    } else if constexpr (std::is_same_v<T, SpecialCaseSpecies>) {
      switch (entry) {
        case SpecialCaseSpecies::e:
          return "e";
        case SpecialCaseSpecies::HM:
          return "HM";
        case SpecialCaseSpecies::DM:
          return "DM";
        case SpecialCaseSpecies::Mg:
          return "Mg";
        case SpecialCaseSpecies::Al:
          return "Al";
        case SpecialCaseSpecies::S:
          return "S";
        case SpecialCaseSpecies::Fe:
          return "Fe";
        default:
          GR_INTERNAL_ERROR("should be unreachable");
      }
    } else {
      GR_INTERNAL_ERROR("should be unreachable");
    }
  };

  SpKind chemical_kinds[2] = {SpKind::PRIMORDIAL, SpKind::METAL};
  for (SpKind kind : chemical_kinds) {
    const std::vector<ChemEntryVariant_>& v =
        (kind == SpKind::PRIMORDIAL) ? primordial_vec_ : metal_vec_;
    sp_kinds.push_back(kind);
    int cur_spkind_count = v.size();
    sp_kind_counts.push_back(cur_spkind_count);
    for (const ChemEntryVariant_& entry : v) {
      // it's important that we pass by reference to an entry inside of the
      // vector (if we make a copy of a std::string, there will be issues with
      // dangling pointers)
      sp_names.emplace_back(std::visit(get_species_name, entry));
    }

    // check if we need to append any filler fields
    // -> this is needed to acheive alignment with indices of SpLUT when
    //    - not using all possible primordial species & at least 1 metal species
    //    - not using all metal species & at least 1 dust species (historically
    //      this only happens if using "Mg", but not "Al", "S", and "Fe")
    // -> we can stop doing this after we remove SpLUT from the codebase
    int expected_spkind_count = cur_spkind_count;
    SpKind filler_kind;
    if (kind == SpKind::PRIMORDIAL && !metal_vec_.empty()) {
      expected_spkind_count = PrimordialSpLUT::NUM_ENTRIES;
      filler_kind = SpKind::FILLER_A;
    } else if (kind == SpKind::METAL && !dust_vec_.empty()) {
      expected_spkind_count = MetalSpLUT::NUM_ENTRIES;
      filler_kind = SpKind::FILLER_B;
    }

    if (cur_spkind_count > expected_spkind_count) {
      GR_INTERNAL_ERROR(
          "this situation shouldn't occur "
          "(cur_spkind_count=%d, expected spkind_count=%d)",
          cur_spkind_count, expected_spkind_count);
    } else if (cur_spkind_count < expected_spkind_count) {
      int n_filler = expected_spkind_count - cur_spkind_count;
      sp_kinds.push_back(filler_kind);
      sp_kind_counts.push_back(n_filler);
      for (int i = 0; i < n_filler; i++) {
        sp_names.push_back(dummy_names_[sp_names.size()]);
      }
    }
  }

  if (!dust_vec_.empty()) {
    sp_kinds.push_back(SpKind::DUST);
    sp_kind_counts.push_back(dust_vec_.size());
    for (const std::string& name : dust_vec_) {
      sp_names.push_back(name.c_str());
    }
  }

  // Now, let's construct the actual objects that we return
  // ======================================================
  Expected<PartMap<SpKind>, Error> kind_map_rslt = PartMap<SpKind>::create(
      sp_kinds.data(), sp_kind_counts.data(), sp_kinds.size());
  if (!kind_map_rslt.has_value()) {
    return Unexpected(
        kind_map_rslt.error().context_literal("error building kind_map"));
  }
  Expected<FrozenKeyIdxBiMap, Error> name_map_rslt = FrozenKeyIdxBiMap::create(
      sp_names.data(), sp_names.size(), BiMapMode::COPIES_KEYDATA);
  if (!name_map_rslt.has_value()) {
    return Unexpected(
        name_map_rslt.error().context_literal("error building name_map"));
  }

  return SpInitializeInfo{.kind_map = std::move(kind_map_rslt).value(),
                          .name_map = std::move(name_map_rslt).value()};
}

Expected<SpInitializeInfo, Error> infer_species_maps(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info) {
  SpeciesBuilder builder;

  // configure the primordial species
  // ================================
  if ((primordial_chemistry == 0) &&
      (metal_chemistry > 0 || grain_info != nullptr)) {
    return Unexpected(Error::msgf(
        "can't have metal or grain species if primordial_chemistry==0"));
  }

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

  return builder.build();
}

}  // namespace GRIMPL_NAMESPACE_DECL