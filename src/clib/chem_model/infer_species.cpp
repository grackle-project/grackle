
//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// <add a short description>
///
//===----------------------------------------------------------------------===//

#include <array>
#include <vector>

#include "chem_model/infer_species.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {
/// Encodes the choice for the species that are only used to track the mass
/// in a subset of metal nuclides
///
/// @note
/// At the time of writing, the choice is tied to the dust model
enum class MetalNuclideMassSpChoice { NONE, ONLY_Mg, ALL };

static MetalNuclideMassSpChoice from_dust_species_param(int dust_species) {
  switch (dust_species) {
    case 1:
      return MetalNuclideMassSpChoice::ONLY_Mg;
    case 2:
      [[fallthrough]];
    case 3:
      return MetalNuclideMassSpChoice::ALL;
    default:
      return MetalNuclideMassSpChoice::NONE;
  }
}

struct CanonicalSpList {
  const char* const* entries;
  int actual_len;
  int max_len;
};

static constexpr std::array<const char*, 15> canonical_primsp_ = {
    "e",    "HI", "HII", "HeI", "HeII", "HeIII", "HM",   "H2I",
    "H2II", "DI", "DII", "HDI", "DM",   "HDII",  "HeHII"};

/// @param[in] primordial_chemistry The primordial_chemistry configuration
///     parameter
Expected<CanonicalSpList, Error> canonical_prim_sp(int primordial_chemistry) {
  CanonicalSpList out;
  out.entries = canonical_primsp_.data();
  out.max_len = canonical_primsp_.size();
  if (primordial_chemistry == 1) {
    out.actual_len = 6;
  } else if (primordial_chemistry == 2) {
    out.actual_len = 9;
  } else if (primordial_chemistry == 3) {
    out.actual_len = 12;
  } else if (primordial_chemistry == 4) {
    out.actual_len = canonical_primsp_.size();
  } else {
    return Unexpected(Error::msgf("invalid primordial chemsitry value: %d",
                                  primordial_chemistry));
  }
  return out;
}

static constexpr std::array<const char*, 23> canonical_metalsp_ = {
    "CI",    "CII",   "COI",   "CO2I", "OI",   "OHI",  "H2OI", "O2I",
    "SiI",   "SiOI",  "SiO2I", "CHI",  "CH2I", "COII", "OII",  "OHII",
    "H2OII", "H3OII", "O2II",  "Mg",   "Al",   "S",    "Fe"};

/// @param[in] metal_chemistry The metal_chemistry configuration parameter
/// @param[in] choice The choice pertaining the metal nuclide mass-tracking
///     species
static Expected<CanonicalSpList, Error> canonical_metal_sp(
    int metal_chemistry, MetalNuclideMassSpChoice choice) {
  CanonicalSpList out;
  out.entries = canonical_metalsp_.data();
  out.max_len = canonical_metalsp_.size();
  if (metal_chemistry == 0 && choice != MetalNuclideMassSpChoice::NONE) {
    return Unexpected(Error::msg_literal(
        "dust chemistry requires metal fields & metal_chemistry is disabled"));
  } else if (metal_chemistry == 0) {
    out.actual_len = 0;
  } else if (metal_chemistry == 1) {
    int min_length = canonical_metalsp_.size() - 4;
    switch (choice) {
      case MetalNuclideMassSpChoice::NONE:
        out.actual_len = min_length;
        break;
      case MetalNuclideMassSpChoice::ONLY_Mg:
        out.actual_len = min_length + 1;
        break;
      case MetalNuclideMassSpChoice::ALL:
        out.actual_len = min_length + 4;
        break;
      default:
        GR_INTERNAL_UNREACHABLE_ERROR();
    }
  } else {
    return Unexpected(Error::msg_literal("not a valid metal_chemistry value"));
  }
  return out;
}

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

Expected<SpInitializeInfo, Error> infer_species_maps(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info) {
  // ugh, this function is so ugly...
  // -> in the future, we'll probably be better off adopting a builder pattern
  std::vector<SpKind> sp_kinds;
  std::vector<int> sp_kind_counts;
  std::vector<const char*> sp_names;

  // primordial species
  // ==================
  // lookup the list of primordial species
  Expected<CanonicalSpList, Error> prim_list_rslt =
      canonical_prim_sp(primordial_chemistry);
  if (!prim_list_rslt.has_value()) {
    return Unexpected(prim_list_rslt.error());
  }
  const CanonicalSpList& prim_list = prim_list_rslt.value();
  if (prim_list.actual_len == 0) {
    return Unexpected(Error::msg_literal("no primordial species"));
  }

  // append all of the primordial species names
  sp_kinds.push_back(SpKind::PRIMORDIAL);
  sp_kind_counts.push_back(prim_list.actual_len);
  for (int i = 0; i < prim_list.actual_len; i++) {
    sp_names.push_back(prim_list.entries[i]);
  }

  // metal species
  // =============
  // lookup the list of metal species names
  MetalNuclideMassSpChoice choice =
      (grain_info == nullptr)
          ? MetalNuclideMassSpChoice::NONE
          : from_dust_species_param(grain_info->dust_species_parameter());
  Expected<CanonicalSpList, Error> metal_list_rslt =
      canonical_metal_sp(metal_chemistry, choice);
  if (!metal_list_rslt.has_value()) {
    return Unexpected(metal_list_rslt.error());
  }
  const CanonicalSpList& metal_list = metal_list_rslt.value();

  // check if we need to append any filler fields
  // -> reminder: this is necessary to acheive alignment with indices of SpLUT
  //    if we aren't using all known primordial species
  // -> we can stop doing this after we remove SpLUT
  if (metal_list.actual_len > 0 && prim_list.actual_len != prim_list.max_len) {
    int n_filler = prim_list.max_len - prim_list.actual_len;
    sp_kinds.push_back(SpKind::FILLER_A);
    sp_kind_counts.push_back(n_filler);
    for (int i = 0; i < n_filler; i++) {
      sp_names.push_back(dummy_names_[sp_names.size()]);
    }
  }

  // append all of the metal species names
  sp_kinds.push_back(SpKind::METAL);
  sp_kind_counts.push_back(metal_list.actual_len);
  for (int i = 0; i < metal_list.actual_len; i++) {
    sp_names.push_back(metal_list.entries[i]);
  }

  // dust species
  // ============
  int n_grain_species = (grain_info == nullptr) ? 0 : grain_info->n_species();

  // check if we need to append any filler fields
  // -> reminder: this is necessary to acheive alignment with indices of SpLUT
  //    if we aren't using all known metal species (this happens with
  //    choice == MetalNuclideMassSpChoice::ONLY_Mg)
  // -> we can stop doing this after we remove SpLUT
  if (n_grain_species > 0 && metal_list.actual_len != metal_list.max_len) {
    int n_filler = metal_list.max_len - metal_list.actual_len;
    sp_kinds.push_back(SpKind::FILLER_B);
    sp_kind_counts.push_back(n_filler);
    for (int i = 0; i < n_filler; i++) {
      sp_names.push_back(dummy_names_[sp_names.size()]);
    }
  }

  // append all of the dust species names
  sp_kinds.push_back(SpKind::DUST);
  sp_kind_counts.push_back(n_grain_species);
  if (grain_info != nullptr) {
    const FrozenKeyIdxBiMap& grainsp_name_map = grain_info->name_map();
    for (int i = 0; i < n_grain_species; i++) {
      sp_names.push_back(grainsp_name_map.inverse_find(i));
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

}  // namespace GRIMPL_NAMESPACE_DECL