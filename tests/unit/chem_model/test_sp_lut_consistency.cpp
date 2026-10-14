//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// test that the order of species are consistent with the LUTs
///
//===----------------------------------------------------------------------===//

#include <array>
#include <cstdint>
#include <iostream>  // needed to teach googletest how to print
#include <vector>
#include <set>
#include <utility>  // std::tuple

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "chem_model/chemsp_props.hpp"
#include "chem_model/species.hpp"
#include "dust/grain_species_info.hpp"
#include "LUT.hpp"
#include "gtest/gtest.h"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"
#include "support/expected.hpp"

using GRIMPL_NS::Error;
using GRIMPL_NS::Expected;

struct ParameterSpec {
  int primordial_chemistry;
  int metal_chemistry;
  int dust_species;
};

// teach GoogleTest how to print our sample partition descriptors for more
// informative errors
void PrintTo(ParameterSpec spec, std::ostream* os) {
  *os << "{prim=" << spec.primordial_chemistry
      << ",metal=" << spec.metal_chemistry << ",dustsp=" << spec.dust_species
      << '}';
}

namespace GRIMPL_NAMESPACE_DECL {
// teach GoogleTest how to print our SpKind enum for more informative errors
void PrintTo(const SpKind& kind, std::ostream* os) {
  switch (kind) {
    case SpKind::PRIMORDIAL:
      *os << "SpKind::PRIMORDIAL";
      return;
    case SpKind::METAL:
      *os << "SpKind::METAL";
      return;
    case SpKind::DUST:
      *os << "SpKind::DUST";
      return;
    case SpKind::FILLER_A:
      *os << "SpKind::FILLER_A";
      return;
    case SpKind::FILLER_B:
      *os << "SpKind::FILLER_B";
      return;
    default:
      *os << "SpKind::<UNKNOWN>";
      return;
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

struct SpTriple {
  std::string name;
  int splut_idx;
  int kindlut_idx;
};

struct ExpectedSpData {
  std::vector<SpTriple> primordial;
  std::vector<SpTriple> metal;
  std::vector<SpTriple> dust;
  // none of the strings in the following set should be included
  std::set<std::string> omit_set;
};

ExpectedSpData get_expected_sp_data(ParameterSpec spec) {
  ExpectedSpData out;
  std::vector<SpTriple> omit;

  {
    std::vector<SpTriple>& v =
        (spec.primordial_chemistry >= 1) ? out.primordial : omit;
    v.emplace_back("e", SpLUT::e, PrimordialSpLUT::e);
    v.emplace_back("HI", SpLUT::HI, PrimordialSpLUT::HI);
    v.emplace_back("HII", SpLUT::HII, PrimordialSpLUT::HII);
    v.emplace_back("HeI", SpLUT::HeI, PrimordialSpLUT::HeI);
    v.emplace_back("HeII", SpLUT::HeII, PrimordialSpLUT::HeII);
    v.emplace_back("HeIII", SpLUT::HeIII, PrimordialSpLUT::HeIII);
  }

  {
    std::vector<SpTriple>& v =
        (spec.primordial_chemistry >= 2) ? out.primordial : omit;
    v.emplace_back("HM", SpLUT::HM, PrimordialSpLUT::HM);
    v.emplace_back("H2I", SpLUT::H2I, PrimordialSpLUT::H2I);
    v.emplace_back("H2II", SpLUT::H2II, PrimordialSpLUT::H2II);
  }

  {
    std::vector<SpTriple>& v =
        (spec.primordial_chemistry >= 3) ? out.primordial : omit;
    v.emplace_back("DI", SpLUT::DI, PrimordialSpLUT::DI);
    v.emplace_back("DII", SpLUT::DII, PrimordialSpLUT::DII);
    v.emplace_back("HDI", SpLUT::HDI, PrimordialSpLUT::HDI);
  }

  {
    std::vector<SpTriple>& v =
        (spec.primordial_chemistry >= 4) ? out.primordial : omit;
    v.emplace_back("DM", SpLUT::DM, PrimordialSpLUT::DM);
    v.emplace_back("HDII", SpLUT::HDII, PrimordialSpLUT::HDII);
    v.emplace_back("HeHII", SpLUT::HeHII, PrimordialSpLUT::HeHII);
  }

  {
    std::vector<SpTriple>& v = (spec.metal_chemistry == 1) ? out.metal : omit;
    v.emplace_back("CI", SpLUT::CI, MetalSpLUT::CI);
    v.emplace_back("CII", SpLUT::CII, MetalSpLUT::CII);
    v.emplace_back("COI", SpLUT::COI, MetalSpLUT::COI);
    v.emplace_back("CO2I", SpLUT::CO2I, MetalSpLUT::CO2I);
    v.emplace_back("OI", SpLUT::OI, MetalSpLUT::OI);
    v.emplace_back("OHI", SpLUT::OHI, MetalSpLUT::OHI);
    v.emplace_back("H2OI", SpLUT::H2OI, MetalSpLUT::H2OI);
    v.emplace_back("O2I", SpLUT::O2I, MetalSpLUT::O2I);
    v.emplace_back("SiI", SpLUT::SiI, MetalSpLUT::SiI);
    v.emplace_back("SiOI", SpLUT::SiOI, MetalSpLUT::SiOI);
    v.emplace_back("SiO2I", SpLUT::SiO2I, MetalSpLUT::SiO2I);
    v.emplace_back("CHI", SpLUT::CHI, MetalSpLUT::CHI);
    v.emplace_back("CH2I", SpLUT::CH2I, MetalSpLUT::CH2I);
    v.emplace_back("COII", SpLUT::COII, MetalSpLUT::COII);
    v.emplace_back("OII", SpLUT::OII, MetalSpLUT::OII);
    v.emplace_back("OHII", SpLUT::OHII, MetalSpLUT::OHII);
    v.emplace_back("H2OII", SpLUT::H2OII, MetalSpLUT::H2OII);
    v.emplace_back("H3OII", SpLUT::H3OII, MetalSpLUT::H3OII);
    v.emplace_back("O2II", SpLUT::O2II, MetalSpLUT::O2II);
  }

  {
    std::vector<SpTriple>& v = (spec.dust_species >= 1) ? out.metal : omit;
    v.emplace_back("Mg", SpLUT::Mg, MetalSpLUT::Mg);
  }

  {
    std::vector<SpTriple>& v = (spec.dust_species >= 2) ? out.metal : omit;
    v.emplace_back("Al", SpLUT::Al, MetalSpLUT::Al);
    v.emplace_back("S", SpLUT::S, MetalSpLUT::S);
    v.emplace_back("Fe", SpLUT::Fe, MetalSpLUT::Fe);
  }

  {
    std::vector<SpTriple>& v = (spec.dust_species >= 1) ? out.dust : omit;
    v.emplace_back("MgSiO3_dust", SpLUT::MgSiO3_dust, DustSpLUT::MgSiO3_dust);
    v.emplace_back("AC_dust", SpLUT::AC_dust, DustSpLUT::AC_dust);
  }

  {
    std::vector<SpTriple>& v = (spec.dust_species >= 2) ? out.dust : omit;
    v.emplace_back("SiM_dust", SpLUT::SiM_dust, DustSpLUT::SiM_dust);
    v.emplace_back("FeM_dust", SpLUT::FeM_dust, DustSpLUT::FeM_dust);
    v.emplace_back("Mg2SiO4_dust", SpLUT::Mg2SiO4_dust,
                   DustSpLUT::Mg2SiO4_dust);
    v.emplace_back("Fe3O4_dust", SpLUT::Fe3O4_dust, DustSpLUT::Fe3O4_dust);
    v.emplace_back("SiO2_dust", SpLUT::SiO2_dust, DustSpLUT::SiO2_dust);
    v.emplace_back("MgO_dust", SpLUT::MgO_dust, DustSpLUT::MgO_dust);
    v.emplace_back("FeS_dust", SpLUT::FeS_dust, DustSpLUT::FeS_dust);
    v.emplace_back("Al2O3_dust", SpLUT::Al2O3_dust, DustSpLUT::Al2O3_dust);
  }

  {
    std::vector<SpTriple>& v = (spec.dust_species >= 3) ? out.dust : omit;
    v.emplace_back("ref_org_dust", SpLUT::ref_org_dust,
                   DustSpLUT::ref_org_dust);
    v.emplace_back("vol_org_dust", SpLUT::vol_org_dust,
                   DustSpLUT::vol_org_dust);
    v.emplace_back("H2O_ice_dust", SpLUT::H2O_ice_dust,
                   DustSpLUT::H2O_ice_dust);
  }

  // let's use entries of omit to populate out.omit_set
  for (const SpTriple& triple : omit) {
    out.omit_set.insert(triple.name);
  }
  return out;
}

// a helper function to help with constrrucing the species maps for the tests
Expected<GRIMPL_NS::SpInitializeInfo, Error> infer_species_maps(
    ParameterSpec spec) {
  using GRIMPL_NS::GrainSpeciesInfo;
  std::optional<GrainSpeciesInfo> maybe_grainsp_info;
  const GrainSpeciesInfo* grain_info = nullptr;
  if (spec.dust_species > 0) {
    using GRIMPL_NS::GrainSpeciesInfo;
    GRIMPL_NS::Expected<GrainSpeciesInfo, GRIMPL_NS::Error> rslt =
        GrainSpeciesInfo::create(spec.dust_species);
    if (!rslt.has_value()) {
      return GRIMPL_NS::Unexpected(
          rslt.error().context_literal("unable to init GrainSpeciesInfo"));
    }
    maybe_grainsp_info = std::move(rslt).value();
    grain_info = &maybe_grainsp_info.value();
  }

  return GRIMPL_NS::infer_species_maps(spec.primordial_chemistry,
                                       spec.metal_chemistry, grain_info);
}

class SpeciesInfo : public testing::TestWithParam<ParameterSpec> {
protected:
  GRIMPL_NS::SpInitializeInfo* species_info = nullptr;

  void SetUp() override {
    // we perform setup here (rather than in a constructor) since it lets us
    // fail the test
    ParameterSpec spec = GetParam();
    Expected<GRIMPL_NS::SpInitializeInfo, Error> tmp = infer_species_maps(spec);
    if (!tmp.has_value()) {
      GTEST_FAIL() << "unable to build species maps\n\n"
                   << "  for: " << testing::PrintToString(spec) << '\n'
                   << "  error:\n"
                   << tmp.error().to_string();
    }
    species_info = new GRIMPL_NS::SpInitializeInfo(std::move(tmp).value());
  }

  ~SpeciesInfo() override {
    delete species_info;  // reminder: it's ok to delete a nullptr
  }
};

TEST_P(SpeciesInfo, LUTConsistency) {
  using GRIMPL_NS::SpKind;

  const GRIMPL_NS::FrozenKeyIdxBiMap name_map = species_info->name_map;
  const GRIMPL_NS::PartMap<SpKind> kind_map = species_info->kind_map;

  // prepare for comparisons agains our expectations
  ParameterSpec spec = GetParam();
  ExpectedSpData expected = get_expected_sp_data(spec);

  using CmpProps =
      std::tuple<const std::vector<SpTriple>&, SpKind, const char*>;
  CmpProps comparisons[3] = {
      CmpProps{expected.primordial, SpKind::PRIMORDIAL, "PrimordialSpLUT"},
      CmpProps{expected.metal, SpKind::METAL, "MetalSpLUT"},
      CmpProps{expected.dust, SpKind::DUST, "DustSpLUT"}};

  // check all expected species names
  for (const CmpProps& cmp_props : comparisons) {
    const std::vector<SpTriple>& v = std::get<0>(cmp_props);
    SpKind sp_kind = std::get<1>(cmp_props);
    const char* kindlut_name = std::get<2>(cmp_props);

    // first, check size of partition associated with sp_kind
    {
      int expected_size = static_cast<int>(v.size());
      GRIMPL_NS::IndexInterval1D idx_interval = kind_map.part_bounds(sp_kind);
      int actual_size = idx_interval.stop - idx_interval.start;
      EXPECT_EQ(actual_size, expected_size)
          << "the size of the partition associated with "
          << testing::PrintToString(sp_kind) << " species is wrong";
    }

    // second, lets iterate over all of the entries
    for (const SpTriple& entry : v) {
      std::optional<uint16_t> maybe_sp_idx = name_map.find(entry.name.c_str());

      if (!maybe_sp_idx.has_value()) {
        ADD_FAILURE() << "name_map is missing the species named \""
                      << entry.name << '"';
        continue;
      }
      int sp_idx = static_cast<int>(maybe_sp_idx.value());

      // lets check that SpLUT::species matches sp_idx
      // -> once we delete SpLUT, we can obviously remove this check
      if (sp_idx != entry.splut_idx) {
        const char* name_cstr = name_map.inverse_find(entry.splut_idx);
        std::string reason;

        if (name_cstr == nullptr) {
          reason = "(name_map doesn't have a species at the latter index)";
        } else {
          reason = ("(name_map associates the latter index with the \"" +
                    std::string(name_cstr) + "\" species)");
        }

        ADD_FAILURE() << "index associated with the \"" << entry.name
                      << "\" species "
                      << "is " << sp_idx << ". It should be " << entry.splut_idx
                      << ", the value of SpLUT::" << entry.name << ' '
                      << reason;
      }

      // let's check consistency with (Primordial|Metal|Dust)SpLUT
      std::optional<GRIMPL_NS::partmap::IdxSearch<SpKind>> maybe_rslt =
          kind_map.search_idx(sp_idx);
      if (!maybe_rslt.has_value()) {
        ADD_FAILURE() << "the index " << sp_idx << " (associated with the \""
                      << entry.name
                      << "\" species) isn't associated with any SpKind. It "
                      << "should be associated with "
                      << testing::PrintToString(sp_kind);
        continue;
      }
      GRIMPL_NS::partmap::IdxSearch<SpKind> search_rslt = maybe_rslt.value();

      if (search_rslt.pd != sp_kind) {
        ADD_FAILURE() << "the SpKind associated with the \"" << entry.name
                      << "\" species should be "
                      << testing::PrintToString(sp_kind) << ", not "
                      << testing::PrintToString(search_rslt.pd);
      } else if (search_rslt.start_offset != entry.kindlut_idx) {
        ADD_FAILURE()
            << "the species kind-specific index associated with the \""
            << entry.name << "\" species is " << search_rslt.start_offset
            << ". It should be " << entry.kindlut_idx << ", the value of "
            << kindlut_name << "::" << entry.name;
      }
    }
  }

  // confirm that name_map doesn't include all explicitly omitted names
  for (const std::string& name : expected.omit_set) {
    EXPECT_EQ(name_map.find(name.c_str()), std::nullopt)
        << "name_map should NOT include the \"" << name << "\" species";
  }
}

// this test checks inferred species info
TEST_P(SpeciesInfo, PropertyCheckFreeElectron) {
  using GRIMPL_NS::ChemSpeciesProps;

  // check the entry for free electrons
  std::optional<uint16_t> maybe_e_id = species_info->name_map.find("e");
  if (!maybe_e_id.has_value()) {
    GTEST_FAIL() << "unable to lookup id associated with 'e'";
  }
  int electron_id = maybe_e_id.value();
  const ChemSpeciesProps& prop = species_info->chem_species_props[electron_id];
  EXPECT_EQ(prop.elementary_charge, -1);
  EXPECT_FALSE(prop.composed_of_nuclides());
  EXPECT_EQ(prop.mass_factor, GRIMPL_NS::mass_factor::e);
}

TEST_P(SpeciesInfo, PropertyCheckSingleNuclide) {
  using GRIMPL_NS::ChemSpeciesProps;
  // let's spot check a few normal cases and then explicitly check each special
  // species that is composed of nuclides

  struct SingleNuclideChecks {
    const char* name;
    const char* single_nuclide;
    int charge;
  };

  const SingleNuclideChecks cases_to_check[] = {
      {"HI", "H", 0},
      {"HeIII", "He", 2},
      {"DII", "D", 1},
      {"CI", "C", 0},

      // special primordial cases
      {"HM", "H", -1},
      {"DM", "D", -1},

      // special metal cases
      {"Mg", "Mg", 0},
      {"Al", "Al", 0},
      {"S", "S", 0},
      {"Fe", "Fe", 0},
  };

  const GRIMPL_NS::NuclideModel& nuclide_model = species_info->nuclide_model;

  for (const SingleNuclideChecks& c : cases_to_check) {
    std::optional<uint16_t> maybe_species_id =
        species_info->name_map.find(c.name);
    if (!maybe_species_id) {
      continue;
    }
    int species_id = maybe_species_id.value();
    const ChemSpeciesProps& prop = species_info->chem_species_props[species_id];

    EXPECT_TRUE(prop.composed_of_nuclides())
        << "properties of \"" << c.name
        << "\" should indicate it's composed of nuclides";
    EXPECT_EQ(prop.elementary_charge, c.charge)
        << "properties of \"" << c.name << "\" reports wrong charge";

    std::optional<uint16_t> maybe_nuclide_id =
        nuclide_model.symbol_map().find(c.single_nuclide);
    if (!maybe_nuclide_id.has_value()) {
      GTEST_FAIL() << "unable to lookup nuclide id associated with the \""
                   << c.single_nuclide << "\" nuclide";
    }
    uint16_t nuclide_id = maybe_nuclide_id.value();

    EXPECT_EQ(prop.mass_factor, nuclide_model.get(nuclide_id).mass_factor)
        << "properties of \"" << c.name << "\" report the wrong mass_factor";
    if (prop.nuclide_id_coef_pairs.size() != 1) {
      ADD_FAILURE()
          << "properties of \"" << c.name << "\" indicate it's composed of "
          << prop.nuclide_id_coef_pairs.size()
          << " different types of nuclides. It's only composed of 1 type";
    } else if (prop.nuclide_id_coef_pairs[0].coef != 1 ||
               prop.nuclide_id_coef_pairs[0].id != nuclide_id) {
      std::string tmp = nuclide_model.symbol_map().inverse_find(
          prop.nuclide_id_coef_pairs[0].id);
      ADD_FAILURE() << "properties of \"" << c.name
                    << "\" indicate it's composed of "
                    << prop.nuclide_id_coef_pairs[0].coef << " \"" << tmp
                    << "\" nuclide(s). It's actually composed of 1 \""
                    << c.single_nuclide << "\" nuclide";
    }
  }
}

// we are intentionally trying to be very rigorous about covering as many cases
// as we possibly can (any errors that might arise here would be quite annoying
// to debug if we don't catch them here)
static constexpr std::array<ParameterSpec, 20> param_spec_arr({
    // {primordial_chemistry, metal_chemistry, dust_species}
    {1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 1, 3}, {2, 0, 0}, {2, 1, 0},
    {2, 1, 1}, {2, 1, 2}, {2, 1, 3}, {3, 0, 0}, {3, 1, 0}, {3, 1, 1}, {3, 1, 2},
    {3, 1, 3}, {4, 0, 0}, {4, 1, 0}, {4, 1, 1}, {4, 1, 2}, {4, 1, 3},
});
INSTANTIATE_TEST_SUITE_P(,  // First argument deliberately left empty
                         SpeciesInfo, testing::ValuesIn(param_spec_arr));
