//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Perform some basic tests of ChemSpeciesProps
/// - it tests the class itself and
/// - the parsing logic
/// (The two concepts are quite intertwined)
///
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <sys/types.h>

#include "gmock/gmock.h"
#include "grtestutils/googletest/chem_model_tools.hpp"

#include "chem_model/chemsp_props.hpp"
#include "chem_model/parse.hpp"
#include "chem_model/nuclide_model.hpp"
#include "gtest/gtest.h"

using ExpectedFormulaComps =
    GRIMPL_NS::Expected<GRIMPL_NS::parse::FormulaComps, GRIMPL_NS::parse::Err>;
using GRIMPL_NS::ChemSpeciesProps;
using CreationRslt = GRIMPL_NS::Expected<ChemSpeciesProps, GRIMPL_NS::Error>;
using ::grtest::HasChemProps;

// first we'll just preform some basic parsing tests
// -> the parsing logic constructs a GRIMPL_NS::parse::FormulaComps object
// -> a GRIMPL_NS::parse::FormulaComps object is subsequently converted into
//    a GRIMPL_NS::ChemSpeciesProps object
// -> there's a part of me that's tempted to rewrite all of these tests so that
//    they are simply testing the factory method of GRIMPL_NS::ChemSpeciesProps,
//    but I think that's a question for another day

TEST(ParseFormulaCompsErrs, EmptyString) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("").has_value());
}

TEST(ParseFormulaCompsErrs, MissingNumeral) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("He").has_value());
}

TEST(ParseFormulaCompsErrs, LooksLikeALoneNumeral) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("I").has_value());
}

TEST(ParseFormulaCompsErrs, InvalidA) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("hI").has_value())
      << "the first letter of a nuclide symbol must be a capital letter";
}

TEST(ParseFormulaCompsErrs, InvalidB) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("heI").has_value())
      << "the first letter of a nuclide symbol must be a capital letter";
}

TEST(ParseFormulaCompsErrs, InvalidC) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("H2oII").has_value())
      << "both nuclide symbols must start with a capital letter";
}

TEST(ParseFormulaComps, TrickyFromulae) {
  // check a few potentially tricky cases where there appears to be some
  // ambiguity about whether a character is part of a Roman numeral or encodes
  // a nuclide
  // - rule of thumb is: decode the largest roman numeral (that is less than
  //   100) and then parse the rest of the string
  // - we haven't really defined semantics for when the numerals exceed 100
  // - the parsing logic knows nothing about what is a valid formula
  //   (this is by design!)

  struct CmpPack {
    std::string formula;
    int roman_numeral_val;
    std::vector<GRIMPL_NS::parse::SViewCoefPair> name_subscript_pairs;
  };

  CmpPack cmp_packs[] = {
      {"HI", 1, {{"H", 1}}},
      {"HII", 2, {{"H", 1}}},

      // each of the following formulae aren't physically meaningful, but again,
      // that's by design
      {"HIII", 3, {{"H", 1}}},
      {"HIV", 4, {{"H", 1}}},
      {"HV", 5, {{"H", 1}}},
      {"HVI", 6, {{"H", 1}}},
      {"HVII", 7, {{"H", 1}}},
      {"HVIII", 8, {{"H", 1}}},
      {"HIX", 9, {{"H", 1}}},
  };

  for (const CmpPack& pack : cmp_packs) {
    ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(pack.formula);
    if (!rslt.has_value()) {
      ADD_FAILURE()
          << "error parsing formula components from \"" << pack.formula
          << "\" (parsing logic doesn't care whether a formula is meaningful)";
      continue;
    }
    EXPECT_EQ(rslt->roman_numeral_val, pack.roman_numeral_val)
        << "parsed roman numeral value from the \"" << pack.formula
        << "\" formula is wrong";
    EXPECT_THAT(rslt->symbol_subscript_pairs,
                testing::ElementsAreArray(pack.name_subscript_pairs))
        << "parsed coefficient-subscript pairs from the \"" << pack.formula
        << "\" formula is wrong";
  }
}

constexpr bool matches_expectation_() {
  using GRIMPL_NS::parse::SViewCoefPair;
  std::string_view s = "H2OII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  return (rslt.has_value() && rslt->roman_numeral_val == 2 &&
          rslt->symbol_subscript_pairs.size() == 2 &&
          rslt->symbol_subscript_pairs[0] == SViewCoefPair("H", 2) &&
          rslt->symbol_subscript_pairs[1] == SViewCoefPair("O", 1));
}

TEST(ParseFormulaComps, CompileTimeParse) {
  using MyCheck = std::bool_constant<matches_expectation_()>;
  EXPECT_TRUE(MyCheck::value);
}

// move on to testing ChemSpeciesProps
// -----------------------------------

TEST(ChemSpeciesProps, FreeElectron) {
  const ChemSpeciesProps prop = ChemSpeciesProps::free_electron();
  EXPECT_EQ(prop.elementary_charge, -1);
  EXPECT_FALSE(prop.composed_of_nuclides());
  EXPECT_EQ(prop.mass_factor, GRIMPL_NS::mass_factor::e);
}

TEST(ChemSpeciesProps, NeutralAtomicH) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HI", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 0, {{"H", 1}}));
}

TEST(ChemSpeciesProps, Proton) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 1, {{"H", 1}}));
}

TEST(ChemSpeciesProps, SinglyIonizedHelium) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HeII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 1, {{"He", 1}}));
}

TEST(ChemSpeciesProps, DoublyIonizedHelium) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HeIII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 2, {{"He", 1}}));
}

TEST(ChemSpeciesProps, NeutralDiatomicSpecies) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2I", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 0, {{"H", 2}}));
}

TEST(ChemSpeciesProps, IonizedDiatomicSpecies) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2II", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 1, {{"H", 2}}));
}

TEST(ChemSpeciesProps, MoleculeSpeciesA) {
  // the fact there was no subscript on the first nuclide was previously a bug
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("CHI", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(),
              HasChemProps(nuclide_model, 0, {{"C", 1}, {"H", 1}}));
}

TEST(ChemSpeciesProps, MoleculeSpeciesB) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2OII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(),
              HasChemProps(nuclide_model, 1, {{"H", 2}, {"O", 1}}));
}

// the following checks that the machinery works for preparsing the formula
// components at compile-time (I have better ideas for how to actually leverage
// this machinery in the future -- this is mainly a test to exercise the
// machinery in the meantime)
TEST(ChemSpeciesProps, CompileTimeParse) {
  // preparse the formula components at compile-time
  static constinit ExpectedFormulaComps preparsed_formula_comps_rslt =
      GRIMPL_NS::parse::formula_comps("H2OII");

  ASSERT_TRUE(preparsed_formula_comps_rslt.has_value());
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula(
      preparsed_formula_comps_rslt.value(), nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(),
              HasChemProps(nuclide_model, 1, {{"H", 2}, {"O", 1}}));
}

// let's check 2 modes where we expect failures
// --------------------------------------------
TEST(ChemSpeciesPropsError, ParseErr) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt =
      ChemSpeciesProps::from_formula("Not-A-Species", nuclide_model);
  EXPECT_FALSE(rslt.has_value());
}

TEST(ChemSpeciesPropsError, NotEnoughProtons) {
  GRIMPL_NS::NuclideModel nuclide_model;
  // HIII fundamentally can't exist because that would correspond to a charge of
  // +2 and a Hydrogen atom only has 1 proton
  const CreationRslt prop_rslt =
      ChemSpeciesProps::from_formula("HIII", nuclide_model);
  EXPECT_FALSE(prop_rslt.has_value());
}