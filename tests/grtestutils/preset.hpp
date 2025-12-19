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
#include <memory>
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
  primchem4_dustspecies3,
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

/// Tracks the group of Grackle objects needed for executing API functions
///
/// The primary motivation for this object's existence is making sure that
/// the allocations get cleaned up when a test fails
///
/// @note
/// Ideally, we would only make it possible to create a fully initialized
/// instance (in that case, we would delete the default constructor), but that
/// involves a bunch more work
///
/// We choose to implement this in terms of std::unique_ptr (rather than raw
/// pointers) since it implements proper move semantics for us
class GrackleCtxPack {
  /// units used for initializing chemistry_data
  code_units initial_units_;
  /// the fully initialized chemistry_data instance
  std::unique_ptr<chemistry_data> my_chemistry_;
  /// the fully initialized chemistry_data_storage_instance
  std::unique_ptr<chemistry_data_storage> my_rates_;

public:
  /// Construct an uninitialized instance
  GrackleCtxPack() : my_chemistry_(nullptr), my_rates_(nullptr) {}

  GrackleCtxPack(GrackleCtxPack&&) = default;
  GrackleCtxPack& operator=(GrackleCtxPack&&) = default;

  // forbid copy and assignment operations...
  // -> we could re-enable them if we wanted to be able add to clone the
  //    types (or if we wanted to internally use std::shared_ptr
  GrackleCtxPack(const GrackleCtxPack&) = delete;
  GrackleCtxPack& operator=(const GrackleCtxPack&) = delete;

  ~GrackleCtxPack() {
    if (!this->is_initialized()) {
      return;
    }
    local_free_chemistry_data(this->my_chemistry_.get(), this->my_rates_.get());
    // unique_ptr destructor will handle calls to delete
  }

  bool is_initialized() const { return this->my_chemistry_ != nullptr; }

  // getter functions
  const code_units& initial_units() const { return this->initial_units_; }
  chemistry_data* my_chemistry() { return this->my_chemistry_.get(); }
  chemistry_data_storage* my_rates() { return this->my_rates_.get(); }

  /// create an initialized instance from a preset
  static GrackleCtxPack create(const FullConfPreset& preset,
                               InitStatus* status) {
    std::unique_ptr<chemistry_data> my_chemistry(new chemistry_data);
    InitStatus tmp =
        setup_chemistry_data_from_preset(my_chemistry.get(), preset.chemistry);
    if (tmp != InitStatus::success) {
      if (status != nullptr) {
        *status = tmp;
      }
      return GrackleCtxPack();  // return an unitialized instance
    }

    code_units initial_unit = setup_initial_unit(preset.unit);
    std::unique_ptr<chemistry_data_storage> my_rates(
        new chemistry_data_storage);
    if (local_initialize_chemistry_data(my_chemistry.get(), my_rates.get(),
                                        &initial_unit) != GR_SUCCESS) {
      if (status != nullptr) {
        *status = InitStatus::generic_fail;
      }
      return GrackleCtxPack();  // return an unitialized instance
    }

    GrackleCtxPack out;
    out.initial_units_ = initial_unit;
    out.my_chemistry_ = std::move(my_chemistry);
    out.my_rates_ = std::move(my_rates);
    return out;
  }
};

}  // namespace grtest

#endif  // GRTESTUTILS_PRESET_HPP
