//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// test logic pertaining to the species lut
///
//===----------------------------------------------------------------------===//

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../grtestutils/RuntimeSpLUT.hpp"
#include "LUT.hpp"

std::string prep_descr(std::string descr, bool negation) {
  return ((negation) ? "doesn't have " : "has ") + descr;
}

MATCHER(HasPC1Species,
        prep_descr("the primordial_chemistry=1 species", negation)) {
  EXPECT_EQ(arg.find("e"), std::optional<int>{SpLUT::e});
  EXPECT_EQ(arg.find("HI"), std::optional<int>{SpLUT::HI});
  EXPECT_EQ(arg.find("HII"), std::optional<int>{SpLUT::HII});
  EXPECT_EQ(arg.find("HeI"), std::optional<int>{SpLUT::HeI});
  EXPECT_EQ(arg.find("HeII"), std::optional<int>{SpLUT::HeII});
  EXPECT_EQ(arg.find("HeIII"), std::optional<int>{SpLUT::HeIII});
  return true;
}

MATCHER(HasPC2Species,
        prep_descr("the primordial_chemistry=2 species", negation)) {
  EXPECT_EQ(arg.find("HM"), std::optional<int>{SpLUT::HM});
  EXPECT_EQ(arg.find("H2I"), std::optional<int>{SpLUT::H2I});
  EXPECT_EQ(arg.find("H2II"), std::optional<int>{SpLUT::H2II});
  return true;
}

MATCHER(HasPC3Species,
        prep_descr("the primordial_chemistry=3 species", negation)) {
  EXPECT_EQ(arg.find("DI"), std::optional<int>{SpLUT::DI});
  EXPECT_EQ(arg.find("DII"), std::optional<int>{SpLUT::DII});
  EXPECT_EQ(arg.find("HDI"), std::optional<int>{SpLUT::HDI});
  return true;
}

MATCHER(HasPC4Species,
        prep_descr("the primordial_chemistry=4 species", negation)) {
  EXPECT_EQ(arg.find("DM"), std::optional<int>{SpLUT::DM});
  EXPECT_EQ(arg.find("HDII"), std::optional<int>{SpLUT::HDII});
  EXPECT_EQ(arg.find("HeHII"), std::optional<int>{SpLUT::HeHII});
  return true;
}

MATCHER(HasMetalSpecies,
        prep_descr("the metal_chemistry=1 species", negation)) {
  EXPECT_EQ(arg.find("CI"), std::optional<int>{SpLUT::CI});
  EXPECT_EQ(arg.find("CII"), std::optional<int>{SpLUT::CII});
  EXPECT_EQ(arg.find("CO"), std::optional<int>{SpLUT::CO});
  EXPECT_EQ(arg.find("CO2"), std::optional<int>{SpLUT::CO2});
  EXPECT_EQ(arg.find("OI"), std::optional<int>{SpLUT::OI});
  EXPECT_EQ(arg.find("OH"), std::optional<int>{SpLUT::OH});
  EXPECT_EQ(arg.find("H2O"), std::optional<int>{SpLUT::H2O});
  EXPECT_EQ(arg.find("O2"), std::optional<int>{SpLUT::O2});
  EXPECT_EQ(arg.find("SiI"), std::optional<int>{SpLUT::SiI});
  EXPECT_EQ(arg.find("SiOI"), std::optional<int>{SpLUT::SiOI});
  EXPECT_EQ(arg.find("SiO2I"), std::optional<int>{SpLUT::SiO2I});
  EXPECT_EQ(arg.find("CH"), std::optional<int>{SpLUT::CH});
  EXPECT_EQ(arg.find("CH2"), std::optional<int>{SpLUT::CH2});
  EXPECT_EQ(arg.find("COII"), std::optional<int>{SpLUT::COII});
  EXPECT_EQ(arg.find("OII"), std::optional<int>{SpLUT::OII});
  EXPECT_EQ(arg.find("OHII"), std::optional<int>{SpLUT::OHII});
  EXPECT_EQ(arg.find("H2OII"), std::optional<int>{SpLUT::H2OII});
  EXPECT_EQ(arg.find("H3OII"), std::optional<int>{SpLUT::H3OII});
  EXPECT_EQ(arg.find("O2II"), std::optional<int>{SpLUT::O2II});
  return true;
}

MATCHER(HasAllMetalNuclideSpecies,
        prep_descr("all the metal-nuclide tracking species", negation)) {
  EXPECT_EQ(arg.find("Mg"), std::optional<int>{SpLUT::Mg});
  EXPECT_EQ(arg.find("Al"), std::optional<int>{SpLUT::Al});
  EXPECT_EQ(arg.find("S"), std::optional<int>{SpLUT::S});
  EXPECT_EQ(arg.find("Fe"), std::optional<int>{SpLUT::Fe});
  return true;
}

MATCHER(HasDS1Species, prep_descr("the dust_species=1 species", negation)) {
  EXPECT_EQ(arg.find("MgSiO3_dust"), SpLUT::MgSiO3_dust);
  EXPECT_EQ(arg.find("AC_dust"), SpLUT::AC_dust);
  return true;
}

MATCHER(HasDS2Species, prep_descr("the dust_species=2 species", negation)) {
  EXPECT_EQ(arg.find("SiM_dust"), std::optional<int>{SpLUT::SiM_dust});
  EXPECT_EQ(arg.find("FeM_dust"), std::optional<int>{SpLUT::FeM_dust});
  EXPECT_EQ(arg.find("Mg2SiO4_dust"), std::optional<int>{SpLUT::Mg2SiO4_dust});
  EXPECT_EQ(arg.find("Fe3O4_dust"), std::optional<int>{SpLUT::Fe3O4_dust});
  EXPECT_EQ(arg.find("SiO2_dust"), std::optional<int>{SpLUT::SiO2_dust});
  EXPECT_EQ(arg.find("MgO_dust"), std::optional<int>{SpLUT::MgO_dust});
  EXPECT_EQ(arg.find("FeS_dust"), std::optional<int>{SpLUT::FeS_dust});
  EXPECT_EQ(arg.find("Al2O3_dust"), std::optional<int>{SpLUT::Al2O3_dust});
  return true;
}

MATCHER(HasDS3Species, prep_descr("the dust_species=3 species", negation)) {
  EXPECT_EQ(arg.find("ref_org_dust"), std::optional<int>{SpLUT::ref_org_dust});
  EXPECT_EQ(arg.find("vol_org_dust"), std::optional<int>{SpLUT::vol_org_dust});
  EXPECT_EQ(arg.find("H2O_ice_dust"), std::optional<int>{SpLUT::H2O_ice_dust});
  return true;
}

MATCHER(IsEmptyIdxInterval, "") {
  EXPECT_EQ(arg.start, arg.stop);
  return true;
}

TEST(RuntimeSpLUT, Empty) {
  grtest::RuntimeSpLUT lut;
  EXPECT_EQ(lut.size(), 0);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL),
              IsEmptyIdxInterval());
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());
}

TEST(RuntimeSpLUT, PrimordialChemistry1) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(1, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 6);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 6);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
}

TEST(RuntimeSpLUT, PrimordialChemistry2) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(2, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 9);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 9);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
}

TEST(RuntimeSpLUT, PrimordialChemistry3) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(3, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 12);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 12);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
}

TEST(RuntimeSpLUT, PrimordialChemistry4) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 0, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 15);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 15);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
  EXPECT_THAT(lut, HasPC4Species());
}

TEST(RuntimeSpLUT, PrimordialChemistry4MetalChemistry1) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 1, nullptr);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 34);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 34);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST), IsEmptyIdxInterval());

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
  EXPECT_THAT(lut, HasPC4Species());
  EXPECT_THAT(lut, HasMetalSpecies());
}

TEST(RuntimeSpLUT, DustSpecies1) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 1, 1);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 37);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop, 35);
  // this is a little weird for the reasons we describe below
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::DUST).start, 35);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::DUST).stop, 37);

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
  EXPECT_THAT(lut, HasPC4Species());
  EXPECT_THAT(lut, HasMetalSpecies());
  EXPECT_EQ(lut.find("Mg"), std::optional<int>{SpLUT::Mg});
  // the following is weird
  // -> we need to subtract 3 from the SpLUT values since we exclude Al, S,
  //    and Fe from the runtime table (for dust_species < 2)
  // -> this is indicative of the reasons why we must move away from
  //    compile-time tables that include both chemical and dust species
  EXPECT_EQ(lut.find("MgSiO3_dust"), SpLUT::MgSiO3_dust - 3);
  EXPECT_EQ(lut.find("AC_dust"), SpLUT::AC_dust - 3);
}

TEST(RuntimeSpLUT, DustSpecies2) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 1, 2);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 48);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop,
            SpLUT::MgSiO3_dust);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST).start,
              SpLUT::MgSiO3_dust);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST).stop, 48);

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
  EXPECT_THAT(lut, HasPC4Species());
  EXPECT_THAT(lut, HasMetalSpecies());
  EXPECT_THAT(lut, HasAllMetalNuclideSpecies());
  EXPECT_THAT(lut, HasDS1Species());
  EXPECT_THAT(lut, HasDS2Species());
}

TEST(RuntimeSpLUT, DustSpecies3) {
  std::optional<grtest::RuntimeSpLUT> lut_maybe =
      grtest::RuntimeSpLUT::create(4, 1, 3);
  grtest::RuntimeSpLUT lut = lut_maybe.value();
  EXPECT_EQ(lut.size(), 51);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).start, 0);
  EXPECT_EQ(lut.kind_interval(GRIMPL_NS::SpKind::CHEMICAL).stop,
            SpLUT::MgSiO3_dust);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST).start,
              SpLUT::MgSiO3_dust);
  EXPECT_THAT(lut.kind_interval(GRIMPL_NS::SpKind::DUST).stop, 51);

  // species-checks
  EXPECT_THAT(lut, HasPC1Species());
  EXPECT_THAT(lut, HasPC2Species());
  EXPECT_THAT(lut, HasPC3Species());
  EXPECT_THAT(lut, HasPC4Species());
  EXPECT_THAT(lut, HasMetalSpecies());
  EXPECT_THAT(lut, HasAllMetalNuclideSpecies());
  EXPECT_THAT(lut, HasDS1Species());
  EXPECT_THAT(lut, HasDS2Species());
  EXPECT_THAT(lut, HasDS3Species());
}
