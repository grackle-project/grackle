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
#include "./param.hpp"
#include <iosfwd>
#include <vector>

namespace grtest {

/// represents different presets for initializing chemistry_data
///
/// @note
/// In the future, we probably want to add more
enum class ChemPreset {
  DEFAULT,
  primchem0,
  primchem1,
  primchem2,
  primchem3,
  primchem4_dustspecies3,
};

void PrintTo(const ChemPreset& chem_preset, std::ostream* os);

std::string to_string(const ChemPreset& preset);

/// override the settings of my_chem based on the specified preset
std::pair<std::vector<ParamPair>, Status> get_chem_preset_vals_(
    ChemPreset preset);

/// Preset for constructing the code_unit instance used for initializing the
/// Grackle Solver
///
/// @note
/// In the future, we probably want to add more
enum class InitialUnitPreset {
  simple_z0,  // <- no cosmology, z=0
};

std::string to_string(const InitialUnitPreset& preset);

/// return a code_unit instance initialized based on the specified preset
code_units setup_initial_unit(InitialUnitPreset preset);

/// Represents a configuration for creating a GrackleCtxPack with overrides
///
/// @note
/// This type is fleshed out a class, rather than just being a simpler
/// struct-like type for the sole purpose of providing nice stringification
/// operations
class ParamConf {
  // this whole class is a little clunky...

  /// the chemistry preset
  ChemPreset c_preset_;
  /// the initial unit preset
  InitialUnitPreset u_preset_;

  // A few thoughts about overrides_:
  // 1. it's probably better if it's some kind of mapping type
  // 2. long term, we may want overrides to be able to override options from
  //    u_preset_

  /// the parameter overrides
  std::vector<ParamPair> overrides_;

  // it's only possible to default construct through a factory method
  ParamConf() = default;

  bool is_overridden_(std::string_view name) const {
    for (const ParamPair& pair : overrides_) {
      if (name == pair.first) {
        return true;
      }
    }
    return false;
  }

public:
  /// Represents a ParamConfig that is simply built from presets
  static ParamConf SimplePreset(ChemPreset chemistry, InitialUnitPreset unit) {
    ParamConf out;
    out.c_preset_ = chemistry;
    out.u_preset_ = unit;
    return out;
  }

  ParamConf(const std::optional<ChemPreset>& chemistry, InitialUnitPreset unit,
            const std::vector<ParamPair>& param_overrides)
      : c_preset_(chemistry.value_or(ChemPreset::DEFAULT)),
        u_preset_{unit},
        overrides_(param_overrides) {}

  ParamConf(const std::optional<ChemPreset>& chemistry, InitialUnitPreset unit,
            std::vector<ParamPair>&& param_overrides)
      : c_preset_(chemistry.value_or(ChemPreset::DEFAULT)),
        u_preset_{unit},
        overrides_(std::move(param_overrides)) {}

  ChemPreset chem_preset() const { return c_preset_; }
  const std::vector<ParamPair>& param_overrides() const { return overrides_; }
  InitialUnitPreset unit_preset() const { return u_preset_; }

  /// teach Googletest to print the simple description
  friend void PrintTo(const ParamConf& param_conf, std::ostream* os);

  /// this is used to help describe the parameters used in a test result message
  ///
  /// @note
  /// This accepts a template parameter Stream because we want this to work
  /// with Googletest's `::testing::AssertionResult` type. While that type
  /// is probably a subclass of `std::ostream`, that's not explicitly stated in
  /// the documentation (and I don't want things to break)
  template <typename Sink>
  void print_descr(Sink* s, bool short_descr,
                   const std::string& common_indent = "") const;

  /// prints a simple description
  std::string stringify(bool short_summary,
                        const std::string& common_indent = "") const;
};

template <typename Sink>
void ParamConf::print_descr(Sink* s, bool short_descr,
                            const std::string& common_indent) const {
  std::string c_p_str = to_string(this->c_preset_);
  std::string u_p_str = to_string(this->u_preset_);

  std::size_t n_override = this->overrides_.size();

  if (short_descr) {
    if (n_override == 0) {
      *s << common_indent << "Preset{";
      *s << c_p_str << ',' << to_string(this->u_preset_) << '}';
    } else {
      *s << common_indent << "ParamConf(";
      *s << c_p_str << ",<+" << n_override << "overrides>," << u_p_str << ')';
    }
    return;
  }

  // the rest of this function is dedicated to providing a detailed summary
  std::string level1_indent = common_indent + "  ";
  std::string level2_indent = common_indent + "    ";

  *s << common_indent << "ParamConf{\n";
  *s << level1_indent << "ChemPreset: " << c_p_str << ",\n";

  // Show non-overridden parameters corresponding to chem preset
  if (this->c_preset_ != ChemPreset::DEFAULT) {
    *s << level1_indent << "Non-Overridden Preset Params: ";
    std::pair<std::vector<ParamPair>, Status> tmp =
        get_chem_preset_vals_(this->c_preset_);
    if (tmp.second.is_err()) {
      *s << "<ERROR: " << tmp.second.to_string() << ">\n";
    } else {
      // record the pair isn't overridden
      *s << "{\n";
      for (const ParamPair& pair : tmp.first) {
        if (!this->is_overridden_(pair.first)) {
          *s << level2_indent << to_string(pair, true) << ",\n";
        }
      }
      *s << level1_indent << "},\n";
    }
  }

  // handle details about the overrides
  *s << level1_indent << "Parameter Overrides: {\n";
  for (const ParamPair& pair : this->overrides_) {
    *s << level2_indent << to_string(pair, true) << ",\n";
  }
  *s << level1_indent << "},\n";

  // finally show the unit-preset
  *s << level1_indent << "Unit Preset: " << u_p_str << '\n';

  // close the brace
  *s << common_indent << "}\n";
}

}  // namespace grtest

#endif  // GRTESTUTILS_PRESET_HPP
