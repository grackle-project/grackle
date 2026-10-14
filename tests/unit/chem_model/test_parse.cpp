//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of parsing
///
//===----------------------------------------------------------------------===//

#include <type_traits>

#include <gtest/gtest.h>

#include "chem_model/parse.hpp"

using ExpectedFormulaComps =
    GRIMPL_NS::Expected<GRIMPL_NS::parse::FormulaComps, GRIMPL_NS::parse::Err>;
using GRIMPL_NS::parse::SViewCoefPair;

TEST(FormulaCompsErrs, EmptyString) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("").has_value());
}

TEST(FormulaCompsErrs, MissingNumeral) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("He").has_value());
}

TEST(FormulaCompsErrs, NotAName) {
  EXPECT_FALSE(GRIMPL_NS::parse::formula_comps("H2oII").has_value());
}

TEST(FormulaComps, NeutralAtomicH) {
  // this is a special case -- it should never be parsed as Hydrogen Iodide
  std::string s = "HI";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 1));
}

TEST(FormulaComps, Proton) {
  // this is a special case -- it should never be parsed as Hydrogen Iodide
  std::string s = "HII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 2);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 1));
}

TEST(FormulaComps, NeutralMonatomicSpecies) {
  std::string s = "HeI";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("He", 1));
}

TEST(FormulaComps, IonizedMonatomicSpecies) {
  std::string s = "HeIII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 3);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("He", 1));
}

TEST(FormulaComps, NeutralDiatomicSpecies) {
  std::string s = "H2I";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 2));
}

TEST(FormulaComps, IonizedDiatomicSpecies) {
  std::string s = "H2III";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 3);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 1);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 2));
}

TEST(FormulaComps, MolecularSpeciesNoSubscript) {
  // this scenario previously had a bug
  std::string s = "HDI";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 2);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 1));
  EXPECT_EQ(rslt->symbol_subscript_pairs[1], SViewCoefPair("D", 1));
}

TEST(FormulaComps, MoleculeSpecies) {
  std::string s = "H2OII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 2);
  ASSERT_EQ(rslt->symbol_subscript_pairs.size(), 2);
  EXPECT_EQ(rslt->symbol_subscript_pairs[0], SViewCoefPair("H", 2));
  EXPECT_EQ(rslt->symbol_subscript_pairs[1], SViewCoefPair("O", 1));
}

constexpr bool matches_expectation_() {
  std::string_view s = "H2OII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  return (rslt.has_value() && rslt->roman_numeral_val == 2 &&
          rslt->symbol_subscript_pairs.size() == 2 &&
          rslt->symbol_subscript_pairs[0] == SViewCoefPair("H", 2) &&
          rslt->symbol_subscript_pairs[1] == SViewCoefPair("O", 1));
}

TEST(FormulaComps, CompileTimeParse) {
  using MyCheck = std::bool_constant<matches_expectation_()>;
  EXPECT_TRUE(MyCheck::value);
}