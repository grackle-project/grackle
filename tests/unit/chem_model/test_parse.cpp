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
using GRIMPL_NS::parse::CoefSViewPair;

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
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(1, "H"));
}

TEST(FormulaComps, Proton) {
  // this is a special case -- it should never be parsed as Hydrogen Iodide
  std::string s = "HII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 2);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(1, "H"));
}

TEST(FormulaComps, NeutralMonatomicSpecies) {
  std::string s = "HeI";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(1, "He"));
}

TEST(FormulaComps, IonizedMonatomicSpecies) {
  std::string s = "HeIII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 3);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(1, "He"));
}

TEST(FormulaComps, NeutralDiatomicSpecies) {
  std::string s = "H2I";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 1);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(2, "H"));
}

TEST(FormulaComps, IonizedDiatomicSpecies) {
  std::string s = "H2III";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 3);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 1);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(2, "H"));
}

TEST(FormulaComps, MoleculeSpecies) {
  std::string s = "H2OII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_EQ(rslt->roman_numeral_val, 2);
  ASSERT_EQ(rslt->subscript_symbol_pairs.size(), 2);
  EXPECT_EQ(rslt->subscript_symbol_pairs[0], CoefSViewPair(2, "H"));
  EXPECT_EQ(rslt->subscript_symbol_pairs[1], CoefSViewPair(1, "O"));
}

constexpr bool matches_expectation_() {
  std::string_view s = "H2OII";
  ExpectedFormulaComps rslt = GRIMPL_NS::parse::formula_comps(s);
  return (rslt.has_value() && rslt->roman_numeral_val == 2 &&
          rslt->subscript_symbol_pairs.size() == 2 &&
          rslt->subscript_symbol_pairs[0] == CoefSViewPair(2, "H") &&
          rslt->subscript_symbol_pairs[1] == CoefSViewPair(1, "O"));
}

TEST(FormulaComps, CompileTimeParse) {
  using MyCheck = std::bool_constant<matches_expectation_()>;
  EXPECT_TRUE(MyCheck::value);
}