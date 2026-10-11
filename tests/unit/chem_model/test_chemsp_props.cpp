//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Perform some basic tests of ChemSpeciesProps
///
//===----------------------------------------------------------------------===//

// we might want to shift some tests around later. We're still figuring out
// what things will look like... ChemSpeciesProps will probably lie in a weird
// middle layer that we may want to avoid burdening down with lots of tests
// (if we do shift tests, we want to make sure that the logic that calls
// ChemSpeciesProps gets tested instead)

#include <gtest/gtest.h>
#include "chem_model/chemsp_props.hpp"
#include "chem_model/nuclide_model.hpp"

using GRIMPL_NS::ChemSpeciesProps;
using CreationRslt = GRIMPL_NS::Expected<ChemSpeciesProps, GRIMPL_NS::Error>;

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
  EXPECT_TRUE(rslt->composed_of_nuclides());
  EXPECT_EQ(rslt->elementary_charge, 0);
  EXPECT_EQ(rslt->mass_factor, GRIMPL_NS::mass_factor::H);

  ASSERT_EQ(rslt->coef_nuclide_id_pairs.size(), 1);
  EXPECT_EQ(rslt->coef_nuclide_id_pairs[0].first, 1);
  EXPECT_EQ(nuclide_model.symbol_map().inverse_find(
                rslt->coef_nuclide_id_pairs[0].second),
            "H");
}

TEST(ChemSpeciesProps, Proton) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_TRUE(rslt->composed_of_nuclides());
  EXPECT_EQ(rslt->elementary_charge, +1);
  EXPECT_EQ(rslt->mass_factor, GRIMPL_NS::mass_factor::H);

  ASSERT_EQ(rslt->coef_nuclide_id_pairs.size(), 1);
  EXPECT_EQ(rslt->coef_nuclide_id_pairs[0].first, 1);
  EXPECT_EQ(nuclide_model.symbol_map().inverse_find(
                rslt->coef_nuclide_id_pairs[0].second),
            "H");
}

TEST(ChemSpeciesProps, DoublyIonizedHelium) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("HeIII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_TRUE(rslt->composed_of_nuclides());
  EXPECT_EQ(rslt->elementary_charge, +2);
  EXPECT_EQ(rslt->mass_factor, GRIMPL_NS::mass_factor::He);

  ASSERT_EQ(rslt->coef_nuclide_id_pairs.size(), 1);
  EXPECT_EQ(rslt->coef_nuclide_id_pairs[0].first, 1);
  EXPECT_EQ(nuclide_model.symbol_map().inverse_find(
                rslt->coef_nuclide_id_pairs[0].second),
            "He");
}

TEST(ChemSpeciesProps, MoleculeSpecies) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2OII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_TRUE(rslt->composed_of_nuclides());
  EXPECT_EQ(rslt->elementary_charge, +1);
  EXPECT_EQ(rslt->mass_factor,
            2 * GRIMPL_NS::mass_factor::H + GRIMPL_NS::mass_factor::O);

  ASSERT_EQ(rslt->coef_nuclide_id_pairs.size(), 2);
  EXPECT_EQ(rslt->coef_nuclide_id_pairs[0].first, 2);
  EXPECT_EQ(nuclide_model.symbol_map().inverse_find(
                rslt->coef_nuclide_id_pairs[0].second),
            "H");
  EXPECT_EQ(rslt->coef_nuclide_id_pairs[1].first, 1);
  EXPECT_EQ(nuclide_model.symbol_map().inverse_find(
                rslt->coef_nuclide_id_pairs[1].second),
            "O");
}

// let's check 2 modes where we expect the tool to fail
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