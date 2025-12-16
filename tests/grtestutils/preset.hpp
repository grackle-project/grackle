//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// declare some standard pre-defined configuration presets
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_PRESET_HPP
#define GRTESTUTILS_PRESET_HPP

#include "grackle.h"
#include <string>
#include <ostream>

namespace grtest {

/// this only exists so that we can determine the reason that a test fails
enum class InitStatus {
  success,
  generic_fail,
  datafile_notfound,
};

/// represents different presets for initializing chemistry_data
///
/// @note
/// In the future, we probably want to add more
enum class ChemPreset {
  primchem0,
  primchem1,
  primchem2,
  primchem3,
};

/// override the settings of my_chem based on the specified preset
InitStatus setup_chemistry_data_from_preset(chemistry_data* my_chem,
                                            ChemPreset preset);

/// Preset for constructing the code_unit instance used for initializing the
/// Grackle Solver
///
/// @note
/// In the future, we probably want to add more
enum class InitialUnitPreset {
  simple_z0,  // <- no cosmology, z=0
};

/// return a code_unit instance initialized based on the specified preset
code_units setup_initial_unit(InitialUnitPreset preset);

/// Represents the preset for creating a GrackleCtxPack
struct FullConfPreset {
  ChemPreset chemistry;
  InitialUnitPreset unit;
};

// teach googletest how to print FullConfPreset
void PrintTo(const FullConfPreset& preset, std::ostream* os);

}  // namespace grtest

#endif  // GRTESTUTILS_PRESET_HPP
