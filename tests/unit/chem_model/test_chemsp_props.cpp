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

using ExpectedFormulaComps =
    GRIMPL_NS::Expected<GRIMPL_NS::parse::FormulaComps, GRIMPL_NS::parse::Err>;
using GRIMPL_NS::ChemSpeciesProps;

TEST(ChemSpeciesProps, FreeElectron) {
  const ChemSpeciesProps prop = ChemSpeciesProps::free_electron();
  EXPECT_EQ(prop.elementary_charge, -1);
  EXPECT_FALSE(prop.composed_of_nuclides());
  EXPECT_EQ(prop.mass_factor, GRIMPL_NS::mass_factor::e);
}
