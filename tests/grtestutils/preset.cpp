//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// implement logic pertaining to pre-defined configuration presets
///
//===----------------------------------------------------------------------===//

#include "./preset.hpp"
#include "./utils.hpp"

#include "grackle.h"
#include "status_reporting.h"  // GR_INTERNAL_UNREACHABLE_ERROR

static std::string to_string(const grtest::ChemPreset& preset) {
  switch (preset) {
    case grtest::ChemPreset::primchem0:
      return "pc=0";
    case grtest::ChemPreset::primchem1:
      return "pc=1";
    case grtest::ChemPreset::primchem2:
      return "pc=2";
    case grtest::ChemPreset::primchem3:
      return "pc=3";
  }

  GR_INTERNAL_UNREACHABLE_ERROR();
}

grtest::InitStatus grtest::setup_chemistry_data_from_preset(
    chemistry_data* my_chem, ChemPreset preset) {
  if (local_initialize_chemistry_parameters(my_chem) != GR_SUCCESS) {
    return InitStatus::generic_fail;
  }

  if (!grtest::set_standard_datafile(*my_chem, "CloudyData_UVB=HM2012.h5")) {
    return InitStatus::datafile_notfound;
  }

  my_chem->use_grackle = 1;  // chemistry on
  my_chem->use_isrf_field = 1;
  my_chem->with_radiative_cooling = 1;  // cooling on
  my_chem->metal_cooling = 1;           // metal cooling on
  my_chem->UVbackground = 1;            // UV background on

  switch (preset) {
    case ChemPreset::primchem0: {
      my_chem->primordial_chemistry = 0;
      my_chem->dust_chemistry = 0;
      return InitStatus::success;
    }
    case ChemPreset::primchem1: {
      my_chem->primordial_chemistry = 1;
      my_chem->dust_chemistry = 1;
      return InitStatus::success;
    }
    case ChemPreset::primchem2: {
      my_chem->primordial_chemistry = 2;
      my_chem->dust_chemistry = 1;
      return InitStatus::success;
    }
    case ChemPreset::primchem3: {
      my_chem->primordial_chemistry = 3;
      my_chem->dust_chemistry = 1;
      return InitStatus::success;
    }
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

static std::string to_string(const grtest::InitialUnitPreset& preset) {
  switch (preset) {
    case grtest::InitialUnitPreset::simple_z0:
      return "simpleUnit-z=0";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

code_units grtest::setup_initial_unit(grtest::InitialUnitPreset preset) {
  // since we return in the switch statement, the compiler should always warn
  // us if we're missing an enumeration
  switch (preset) {
    case InitialUnitPreset::simple_z0: {
      double initial_redshift = 0.;
      code_units my_units;
      my_units.comoving_coordinates = 0;  // 1 if cosmological sim, 0 if not
      my_units.density_units = 1.67e-24;
      my_units.length_units = 1.0;
      my_units.time_units = 1.0e12;
      my_units.a_units = 1.0;  // units for the expansion factor
      // Set expansion factor to 1 for non-cosmological simulation.
      my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;
      set_velocity_units(&my_units);
      return my_units;
    }
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

void grtest::PrintTo(const grtest::FullConfPreset& preset, std::ostream* os) {
  *os << "Preset{" << to_string(preset.chemistry) << ','
      << to_string(preset.unit) << '}';
}
