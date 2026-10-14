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
#include <gmock/gmock.h>
#include <sys/types.h>

#include "grtestutils/googletest/chem_model_tools.hpp"

#include "chem_model/chemsp_props.hpp"
#include "chem_model/nuclide_model.hpp"

using GRIMPL_NS::ChemSpeciesProps;
using CreationRslt = GRIMPL_NS::Expected<ChemSpeciesProps, GRIMPL_NS::Error>;
using grtest::HasChemProps;

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

TEST(ChemSpeciesProps, DiatomicSpecies) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2II", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(), HasChemProps(nuclide_model, 1, {{"H", 2}}));
}

TEST(ChemSpeciesProps, MoleculeSpeciesA) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("H2OII", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(),
              HasChemProps(nuclide_model, 1, {{"H", 2}, {"O", 1}}));
}

TEST(ChemSpeciesProps, MoleculeSpeciesB) {
  GRIMPL_NS::NuclideModel nuclide_model;
  CreationRslt rslt = ChemSpeciesProps::from_formula("CHI", nuclide_model);
  ASSERT_TRUE(rslt.has_value());
  EXPECT_THAT(rslt.value(),
              HasChemProps(nuclide_model, 0, {{"C", 1}, {"H", 1}}));
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