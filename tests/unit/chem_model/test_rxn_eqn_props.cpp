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

#include <type_traits>

#include <gtest/gtest.h>

#include "chem_model/parse.hpp"
#include "gmock/gmock.h"
#include "grtestutils/googletest/chem_model_tools.hpp"

using ExpectedRxnEqnComps =
    GRIMPL_NS::Expected<GRIMPL_NS::parse::RxnEqnComps, GRIMPL_NS::parse::Err>;
using GRIMPL_NS::parse::SViewCoefPair;
using grtest::NameCoefPair;
using ::testing::ElementsAre;

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

// we are going to try to parse the same equation multiple times with different
// amounts of whitespace
class ParseRxnEqnCompsK1Variants : public testing::TestWithParam<std::string> {
};

TEST_P(ParseRxnEqnCompsK1Variants, k1) {
  std::string s = GetParam();
  ExpectedRxnEqnComps rslt = GRIMPL_NS::parse::rxn_eqn(s);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string(s);
  }
  EXPECT_THAT(rslt->lhs_sp_coef_pairs,
              ElementsAre(SViewCoefPair("HI", 1), SViewCoefPair("e", 1)));
  EXPECT_THAT(rslt->rhs_sp_coef_pairs,
              ElementsAre(SViewCoefPair("HII", 1), SViewCoefPair("e", 2)));
}

INSTANTIATE_TEST_SUITE_P(,  // First argument deliberately left empty
                         ParseRxnEqnCompsK1Variants,
                         testing::Values("HI + e -> HII + 2e",
                                         "HI + e -> HII + 2 e",
                                         "HI  +  e  ->  HII  +  2e",
                                         "HI+e -> HII + 2e", "HI + e -> HII+2e",
                                         "HI + e-> HII + 2e",
                                         "HI + e ->HII + 2e"));

TEST(ParseRxnEqnComps, K22) {
  std::string s = "3HI -> H2I + HI";
  ExpectedRxnEqnComps rslt = GRIMPL_NS::parse::rxn_eqn(s);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string(s);
  }
  EXPECT_THAT(rslt->lhs_sp_coef_pairs, ElementsAre(SViewCoefPair("HI", 3)));
  EXPECT_THAT(rslt->rhs_sp_coef_pairs,
              ElementsAre(SViewCoefPair("H2I", 1), SViewCoefPair("HI", 1)));
}

TEST(ParseRxnEqnComps, RefOrganicDustGrowth) {
  std::string s = "5COI + 5CH2I + 12NI -> 10ref_org_dust";
  ExpectedRxnEqnComps rslt = GRIMPL_NS::parse::rxn_eqn(s);
  if (!rslt.has_value()) {
    GTEST_FAIL() << "problem parsing chemical eqn. Error message:\n\n"
                 << rslt.error().to_string(s);
  }
  EXPECT_THAT(rslt->lhs_sp_coef_pairs,
              ElementsAre(SViewCoefPair("COI", 5), SViewCoefPair("CH2I", 5),
                          SViewCoefPair("NI", 12)));
  EXPECT_THAT(rslt->rhs_sp_coef_pairs,
              ElementsAre(SViewCoefPair("ref_org_dust", 10)));
}

constexpr bool rxn_eqn_parse_matches_expectation_() {
  std::string s = "5COI + 5CH2I + 12NI -> 10ref_org_dust";
  ExpectedRxnEqnComps rslt = GRIMPL_NS::parse::rxn_eqn(s);
  return (rslt.has_value() && rslt->lhs_sp_coef_pairs.size() == 3 &&
          rslt->lhs_sp_coef_pairs[0] == SViewCoefPair("COI", 5) &&
          rslt->lhs_sp_coef_pairs[1] == SViewCoefPair("CH2I", 5) &&
          rslt->lhs_sp_coef_pairs[2] == SViewCoefPair("NI", 12) &&
          rslt->rhs_sp_coef_pairs.size() == 1 &&
          rslt->rhs_sp_coef_pairs[0] == SViewCoefPair("ref_org_dust", 10));
}

TEST(ParseRxnEqnComps, CompileTimeParse) {
  using MyCheck = std::bool_constant<rxn_eqn_parse_matches_expectation_()>;
  EXPECT_TRUE(MyCheck::value);
}