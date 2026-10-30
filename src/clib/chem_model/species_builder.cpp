
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
#include <utility>  // std::move, std::tuple
#include <vector>

#include "LUT.hpp"
#include "chem_model/chemsp_props.hpp"
#include "chem_model/species_builder.hpp"
#include "chem_model/nuclide_model.hpp"
#include "chem_model/species.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<SpeciesBuilder, Error> SpeciesBuilder::from_preset(
    int primordial_chemistry, int metal_chemistry) {
  SpeciesBuilder builder;
  // in this function, we're going to hold ourselves to just use
  // SpeciesBuilder's public interface functions
  // -> the one exception is that we will reserve space in the underlying
  //    vectors to reduce reallocations
  builder.primordial_vec_.reserve(PrimordialSpLUT::NUM_ENTRIES);
  builder.metal_vec_.reserve(MetalSpLUT::NUM_ENTRIES);

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

  return builder;
}

namespace {  // stuff inside an anonymous namespace is local to this file

// these are dummy names that won't collide with known species names
// -> NaSp stands for "Not a Species"
// -> this is a kludge until we remove all occurrences of SpLUT
constexpr std::array<const char*, 40> dummy_names_ = {
    "_NaSp00", "_NaSp01", "_NaSp02", "_NaSp03", "_NaSp04", "_NaSp05", "_NaSp06",
    "_NaSp07", "_NaSp08", "_NaSp09", "_NaSp10", "_NaSp11", "_NaSp12", "_NaSp13",
    "_NaSp14", "_NaSp15", "_NaSp16", "_NaSp17", "_NaSp18", "_NaSp19", "_NaSp20",
    "_NaSp21", "_NaSp22", "_NaSp23", "_NaSp24", "_NaSp25", "_NaSp26", "_NaSp27",
    "_NaSp28", "_NaSp29", "_NaSp30", "_NaSp31", "_NaSp32", "_NaSp33", "_NaSp34",
    "_NaSp35", "_NaSp36", "_NaSp37", "_NaSp38", "_NaSp39"};

struct SpecialCaseProps {
  const char* name;
  std::string_view formula_with_equiv_nuclides;
  int charge;
};

SpecialCaseProps special_case_props_(SpecialCaseSpecies s) {
  switch (s) {
    case SpecialCaseSpecies::e:
      return {"e", "not-a-formula", -1};
    case SpecialCaseSpecies::HM:
      return {"HM", "HI", -1};
    case SpecialCaseSpecies::DM:
      return {"DM", "DI", -1};
    case SpecialCaseSpecies::Mg:
      return {"Mg", "MgI", 0};
    case SpecialCaseSpecies::Al:
      return {"Al", "AlI", 0};
    case SpecialCaseSpecies::S:
      return {"S", "SI", 0};
    case SpecialCaseSpecies::Fe:
      return {"Fe", "FeI", 0};
    default:
      GR_INTERNAL_ERROR("should be unreachable");
  }
}

// here we define the infer_ function that uses the value within a
// SpeciesBuilder::ChemEntryVariant_ to infer the ChemSpeciesProps & get the
// species_name
//
// we need to define overloads where the type of the first argument varies
// -> we need an overload for each type a SpeciesBuilder::ChemEntryVariant_ can
//    hold

using InferRslt = Expected<ChemSpeciesProps, Error>;
static InferRslt infer_(const std::string& formula, const char*& species_name,
                        const NuclideModel& n_model) {
  species_name = formula.c_str();
  return ChemSpeciesProps::from_formula(formula, n_model);
}

static InferRslt infer_(SpecialCaseSpecies s, const char*& species_name,
                        const NuclideModel& n_model) {
  SpecialCaseProps case_props = special_case_props_(s);
  species_name = case_props.name;
  if (s == SpecialCaseSpecies::e) {
    return ChemSpeciesProps::free_electron();
  }
  Expected<ChemSpeciesProps, Error> out = ChemSpeciesProps::from_formula(
      case_props.formula_with_equiv_nuclides, n_model);
  if (out.has_value()) {
    out->elementary_charge = case_props.charge;
  }
  return out;
}

}  // anonymous namespace

Expected<SpInitializeInfo, Error> SpeciesBuilder::build() const {
  if (primordial_vec_.empty() && metal_vec_.empty() && dust_vec_.empty()) {
    return SpInitializeInfo();  // <- return species maps that are empty
  }

  // construct the nuclide model (in the future, the builder could conceivably
  // support customizing the model)
  NuclideModel nuclide_model;

  std::vector<SpKind> sp_kinds;
  std::vector<int> sp_kind_counts;
  std::map<std::size_t, ChemSpeciesProps> chemsp_props;
  std::vector<const char*> sp_names;

  // define a lambda function that records the species name and inferrred
  // ChemSpeciesProps from either of the types contained by a ChemEntryVariant_
  auto record_ = [&chemsp_props, &sp_names,
                  &nuclide_model](const auto& entry) -> Expected<void, Error> {
    // get the species name and infer the ChemSpeciesProps
    const char* sp_name = nullptr;
    Expected<ChemSpeciesProps, Error> props_rslt =
        infer_(entry, sp_name, nuclide_model);
    if (!props_rslt.has_value()) {
      return Unexpected(props_rslt.error());
    }

    // now actually record the information
    std::size_t sp_index = sp_names.size();
    bool emplace_success =
        chemsp_props.emplace(sp_index, std::move(props_rslt).value()).second;
    if (!emplace_success) {
      return Unexpected(
          Error::msg_literal("unable to insert chem species properties"));
    }
    sp_names.push_back(sp_name);

    return {};  // <- indicates success
  };

  SpKind chemical_kinds[2] = {SpKind::PRIMORDIAL, SpKind::METAL};
  for (SpKind kind : chemical_kinds) {
    const std::vector<ChemEntryVariant_>& v =
        (kind == SpKind::PRIMORDIAL) ? primordial_vec_ : metal_vec_;

    sp_kinds.push_back(kind);
    int cur_spkind_count = v.size();
    sp_kind_counts.push_back(cur_spkind_count);

    for (int i = 0; i < cur_spkind_count; i++) {
      Expected<void, Error> rslt = std::visit(record_, v[i]);
      if (!rslt.has_value()) {
        const char* kind_name =
            (kind == SpKind::PRIMORDIAL) ? "primordial" : "metal";
        return Unexpected(rslt.error().contextf(
            "problem recording %s species number %d", kind_name, i));
      }
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
  ChemSpeciesPropsVec chem_species_props;
  chem_species_props.m_ = std::move(chemsp_props);

  return SpInitializeInfo{.kind_map = std::move(kind_map_rslt).value(),
                          .name_map = std::move(name_map_rslt).value(),
                          .chem_species_props = std::move(chem_species_props),
                          .nuclide_model = std::move(nuclide_model)};
}

}  // namespace GRIMPL_NAMESPACE_DECL