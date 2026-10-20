//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of parsing and RxnEqnProps objects
///
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>

#include "chem_model/counter.hpp"
#include "chem_model/parse.hpp"
#include "chem_model/rxn_eqn_props.hpp"
#include "gmock/gmock.h"
#include "grtestutils/googletest/chem_model_tools.hpp"
#include "gtest/gtest.h"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/expected.hpp"
#include "support/status_reporting.hpp"

using ExpectedRxnEqnComps =
    GRIMPL_NS::Expected<GRIMPL_NS::parse::RxnEqnComps, GRIMPL_NS::parse::Err>;
using GRIMPL_NS::RxnEqnProps;
using GRIMPL_NS::parse::SViewCoefPair;
using CreationRslt = GRIMPL_NS::Expected<RxnEqnProps, GRIMPL_NS::Error>;
using GRIMPL_NS::EqnSide;
using ::grtest::HasRxnEqnProps;

TEST(ParseRxnEqnComps, EmptyString) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("").has_value());
}

TEST(ParseRxnEqnCompsErr, MissingArrow) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("HI + e HII + 2e").has_value());
}

TEST(ParseRxnEqnCompsErr, MissingLHS) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("-> HII + 2e").has_value());
}

TEST(ParseRxnEqnCompsErr, MissingRHS) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("HI + e ->").has_value());
}

TEST(ParseRxnEqnCompsErr, MissingLeftPlus) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("HI e -> HII + 2e").has_value());
}

TEST(ParseRxnEqnCompsErr, MissingRightPlus) {
  EXPECT_FALSE(GRIMPL_NS::parse::rxn_eqn("HI + e -> HII 2e").has_value());
}

// attempt to build a "rich" species name map object (i.e. enabling as many
// Grackle species as is easily possible) or immediately abort the executable
// -> if this triggers an abort, that will signal that there's a programming
//    oversight that needs to be fixed
static GRIMPL_NS::SpInitializeInfo rich_SpInitializeInfo_or_abort_() {
  grtest::ParameterSpec spec{
      .primordial_chemistry = 4, .metal_chemistry = 1, .dust_species = 3};
  GRIMPL_NS::Expected<GRIMPL_NS::SpInitializeInfo, GRIMPL_NS::Error> tmp =
      grtest::infer_species_maps(spec);
  if (!tmp.has_value()) {
    GR_INTERNAL_ERROR(
        "something changed somewhere and this function "
        "needs to be updated");
  }
  return std::move(tmp).value();
}

static GRIMPL_NS::FrozenKeyIdxBiMap rich_species_name_map_or_abort_() {
  return rich_SpInitializeInfo_or_abort_().name_map;
}

// we need to disable the following lint to avoid spurious warnings pertaining
// to the testing framework
// NOLINTBEGIN(clang-analyzer-cplusplus.NewDeleteLeaks)

// we are going to try to parse the same equation multiple times with different
// amounts of whitespace
class RxnEqnPropsK1Variants : public testing::TestWithParam<std::string> {};

TEST_P(RxnEqnPropsK1Variants, k1) {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  std::string s = GetParam();  // <- get the formula that will be parsed
  CreationRslt rslt = RxnEqnProps::from_string(s, name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(), HasRxnEqnProps(name_map, {{"HI", 1}, {"e", 1}},
                                           {{"HII", 1}, {"e", 2}}));
}

INSTANTIATE_TEST_SUITE_P(,  // First argument deliberately left empty
                         RxnEqnPropsK1Variants,
                         testing::Values("HI + e -> HII + 2e",
                                         "HI + e -> HII + 2 e",
                                         "HI  +  e  ->  HII  +  2e",
                                         "HI+e -> HII + 2e", "HI + e -> HII+2e",
                                         "HI + e-> HII + 2e",
                                         "HI + e ->HII + 2e"));

TEST(RxnEqnProps, K22) {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  CreationRslt rslt = RxnEqnProps::from_string("3HI -> H2I + HI", name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(),
              HasRxnEqnProps(name_map, {{"HI", 3}}, {{"H2I", 1}, {"HI", 1}}));
}

// confirm that things work with a photon
TEST(RxnEqnProps, K7) {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  CreationRslt rslt =
      RxnEqnProps::from_string("HI + e -> HM + photon", name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(),
              HasRxnEqnProps(name_map, {{"HI", 1}, {"e", 1}}, {{"HM", 1}},
                             /*has_photon=*/EqnSide::RIGHT));
}

TEST(RxnEqnProps, K24) {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  CreationRslt rslt =
      RxnEqnProps::from_string("HI + photon -> HII + e", name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(),
              HasRxnEqnProps(name_map, {{"HI", 1}}, {{"HII", 1}, {"e", 1}},
                             /*has_photon=*/EqnSide::LEFT));
}

TEST(RxnEqnProps, VolOrganicDustGrowth) {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  CreationRslt rslt =
      RxnEqnProps::from_string("COI + 2H2I -> vol_org_dust", name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(), HasRxnEqnProps(name_map, {{"COI", 1}, {"H2I", 2}},
                                           {{"vol_org_dust", 1}}));
}

// the following checks that the machinery works for preparsing the reaction
// components at compile-time (I have better ideas for how to actually leverage
// this machinery in the future -- this is mainly a test to exercise the
// machinery in the meantime)
TEST(RxnEqnProps, CompileTimeParse) {
  // preparse the formula components at compile-time
  static constinit ExpectedRxnEqnComps preparsed_rxn_comps_rslt =
      GRIMPL_NS::parse::rxn_eqn("COI + 2H2I -> vol_org_dust");
  ASSERT_TRUE(preparsed_rxn_comps_rslt.has_value());

  GRIMPL_NS::FrozenKeyIdxBiMap name_map = rich_species_name_map_or_abort_();
  CreationRslt rslt =
      RxnEqnProps::from_eqn_comps(preparsed_rxn_comps_rslt.value(), name_map);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string();
  }
  EXPECT_THAT(rslt.value(), HasRxnEqnProps(name_map, {{"COI", 1}, {"H2I", 2}},
                                           {{"vol_org_dust", 1}}));
}

// NOLINTEND(clang-analyzer-cplusplus.NewDeleteLeaks)

// the following test is going to go through pretty much every standard
// collisional reaction rate in order to check whether charge and nuclide
// counts are conserved
//
// this isn't really a unit-test (since it also requires the chemical formula
// parsing to be correct and the actual reaction listings to be correct), but
// it's definitely useful! Since it's not really a unit-test it *might* make
// sense to move this test

struct RxnNameEquationPair {
  std::string_view name;
  const char* equation;
};

void PrintTo(const RxnNameEquationPair& p, std::ostream* os) {
  *os << "{name=\"" << p.name << "\", equation=\"" << p.equation << "\"}";
}

#define STRINGIFY(s) STRINGIFY_INNER(s)
#define STRINGIFY_INNER(s) #s

TEST(EqnChargeConsistency, StandardCollisional) {
  using GRIMPL_NS::Error;
  using GRIMPL_NS::Expected;
  using GRIMPL_NS::NuclideAndChargeCounter;

  // construct a list of all standard collisional equations

#define ENTRY(name, formula) {STRINGIFY(name), formula},
  RxnNameEquationPair pairs[] = {
#include "collisional_rxn_rate_members.def"
  };
#undef ENTRY

  // initialize the species_info that will be used to parse the above rxns
  GRIMPL_NS::SpInitializeInfo species_info = rich_SpInitializeInfo_or_abort_();

  for (const RxnNameEquationPair& pair : pairs) {
    if (pair.equation == nullptr) {
      continue;
    }

    // parse the formula
    CreationRslt rslt =
        RxnEqnProps::from_string(pair.equation, species_info.name_map);
    if (!rslt.has_value()) {
      ADD_FAILURE() << "problem parsing the reaction named chemical eqn.\n"
                    << "rxn: \"" << testing::PrintToString(pair) << '\n'
                    << "Error message:\n"
                    << rslt.error().to_string();
      continue;
    }
    RxnEqnProps& rxn_props = rslt.value();

    // now let's count up nuclides and charge on the left sides of the equation
    NuclideAndChargeCounter counterL(species_info.chem_species_props,
                                     species_info.kind_map);
    NuclideAndChargeCounter counterR(species_info.chem_species_props,
                                     species_info.kind_map);
    bool counting_err = false;

    EqnSide eqn_sides[] = {EqnSide::LEFT, EqnSide::RIGHT};
    for (EqnSide side : eqn_sides) {
      auto& counter = (side == EqnSide::LEFT) ? counterL : counterR;
      const RxnEqnProps::InplaceIdCoefPairVec& eqn_entries =
          (side == EqnSide::LEFT) ? rxn_props.lhs_entries
                                  : rxn_props.rhs_entries;
      Expected<void, Error> count_rslt = counter.add_contribs(eqn_entries);
      if (!count_rslt.has_value()) {
        ADD_FAILURE() << "problem counting charge and nuclides for "
                      << GRIMPL_NS::to_string(side) << " side of rxn\n"
                      << "rxn: \"" << testing::PrintToString(pair) << '\n'
                      << "Error message:\n"
                      << count_rslt.error().to_string();
        counting_err = true;
      }
    }
    if (counting_err) {
      continue;
    }

    // now let's check for nuclide conservation and charge conservation
    EXPECT_EQ(counterL.charge(), counterR.charge())
        << "the net charge must match on each side of the reaction. \n"
        << "Problem rxn: " << testing::PrintToString(pair);
    if (counterL.nuclide_counts() != counterR.nuclide_counts()) {
      std::map<std::string, int> mapL, mapR;
      for (EqnSide side : eqn_sides) {
        std::map<std::string, int>& map = (side == EqnSide::LEFT) ? mapL : mapR;
        const auto& counter = (side == EqnSide::LEFT) ? counterL : counterR;
        for (const auto& kv_pair : counter.nuclide_counts()) {
          std::string name =
              species_info.nuclide_model.symbol_map().inverse_find(
                  kv_pair.first);
          map[name] = kv_pair.second;
        }
      }
      ADD_FAILURE() << "both sides of a rxn should same counts of nuclides:\n"
                    << "  problematic rxn: " << testing::PrintToString(pair)
                    << '\n'
                    << "  LHS counts: " << testing::PrintToString(mapL) << '\n'
                    << "  RHS counts: " << testing::PrintToString(mapR) << '\n';
    }
  }
}