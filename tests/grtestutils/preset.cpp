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

#include "./param.hpp"
#include "./preset.hpp"
#include "./utils.hpp"

#include "grackle.h"
#include "status_reporting.h"  // GR_INTERNAL_UNREACHABLE_ERROR

#include <ostream>

namespace grtest {
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
    case grtest::ChemPreset::primchem4_dustspecies3:
      return "pc=4-dust_species=3";
  }

  GR_INTERNAL_UNREACHABLE_ERROR();
}

std::pair<std::vector<ParamPair>, InitStatus> get_chem_preset_vals_(
    ChemPreset preset) {
  std::vector<std::pair<std::string, ParamVal>> v = make_ParamPair_vec({
      {"use_grackle", 1},  // chemistry on
      {"use_isrf_field", 1},
      {"with_radiative_cooling", 1},  // cooling on
      {"metal_cooling", 1},           // metal cooling on
      {"UVbackground", 1},            // UV background on
  });

  std::optional<std::string> maybe_datafile =
      get_standard_datafile("CloudyData_UVB=HM2012.h5");
  if (!maybe_datafile.has_value()) {
    return {{}, InitStatus::standard_datafile_notfound};
  } else {
    v.emplace_back("grackle_data_file", maybe_datafile.value());
  }

  switch (preset) {
    case ChemPreset::primchem0: {
      v.emplace_back("primordial_chemistry", 0);
      v.emplace_back("dust_chemistry", 0);
      return {v, InitStatus::success};
    }
    case ChemPreset::primchem1: {
      v.emplace_back("primordial_chemistry", 1);
      v.emplace_back("dust_chemistry", 1);
      return {v, InitStatus::success};
    }
    case ChemPreset::primchem2: {
      v.emplace_back("primordial_chemistry", 2);
      v.emplace_back("dust_chemistry", 1);
      return {v, InitStatus::success};
    }
    case ChemPreset::primchem3: {
      v.emplace_back("primordial_chemistry", 3);
      v.emplace_back("dust_chemistry", 1);
      return {v, InitStatus::success};
    }
    case ChemPreset::primchem4_dustspecies3: {
      v.emplace_back("primordial_chemistry", 4);
      v.emplace_back("dust_chemistry", 1);
      v.emplace_back("metal_chemistry", 1);
      v.emplace_back("dust_species", 3);
      v.emplace_back("use_dust_density_field", 1);
      return {v, InitStatus::success};
    }
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

static std::string to_string(const InitialUnitPreset& preset) {
  switch (preset) {
    case grtest::InitialUnitPreset::simple_z0:
      return "simpleUnit-z=0";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

code_units setup_initial_unit(InitialUnitPreset preset) {
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

void PrintTo(const grtest::FullConfPreset& preset, std::ostream* os) {
  *os << "Preset{" << to_string(preset.chemistry) << ','
      << to_string(preset.unit) << '}';
}

}  // namespace grtest
