//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// handles chemical species properties
///
//===----------------------------------------------------------------------===//

#ifndef CHEM_MODEL_CHEMSP_PROPS_HPP
#define CHEM_MODEL_CHEMSP_PROPS_HPP

#include <cstdint>
#include "chem_model/nuclide_model.hpp"
#include "chem_model/parse.hpp"
#include "phys_constants.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/inplace_vector.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief An internal data structure used to track chemical species
///
/// @note We could potentially track other properties like the actual particle
/// mass and whether all particles all truly have one mass or if this is really
/// a mixed ensemble (due to the presence of isotopes). In practice, those
/// other properties don't actually seem all that useful
struct ChemSpeciesProps {
  /// tracks coefficient-nuclide pairs (if any) that make up the species
  InplaceVector<std::pair<uint8_t, uint16_t>, parse::MaxSymsPerFormula>
      coef_nuclide_id_pairs;

  /// The elementary charge per particle
  ///
  /// @important
  /// This field should be ignored if the @ref ParticleTableRow::charge_purity
  /// member holds a value of @ref PropertyPurity::mixed
  int elementary_charge;

  /// total mass factor
  double mass_factor;

private:
  ChemSpeciesProps() = default;  // <- only used by factory methods

public:
  /// @brief check whether the species is composed of nuclides
  bool composed_of_nuclides() const { return coef_nuclide_id_pairs.size() > 0; }

  // Define the factory methods
  // --------------------------

  /// @brief construct an object that corresponds to free-electrons
  ///
  /// A free electron is a fairly special case
  static ChemSpeciesProps free_electron() {
    ChemSpeciesProps out;
    out.elementary_charge = -1;
    out.mass_factor = mass_factor::e;
    // out.coef_nuclide_id_pairs is left empty since electrons aren't composed
    // of nuclides
    return out;
  }

  /// @brief Primary Factory Method
  ///
  /// Constructs a @ref ChemSpeciesProps by parsing a formula
  static Expected<ChemSpeciesProps, Error> from_formula(
      const std::string_view formula, const NuclideModel& model) {
    Expected<parse::FormulaComps, parse::Err> parse_rslt =
        parse::formula_comps(formula);
    if (!parse_rslt.has_value()) {
      return Unexpected(Error::parse(formula, parse_rslt.error())
                            .context_literal("problem parsing formula"));
    }
    return ChemSpeciesProps::from_formula(parse_rslt.value(), model);
  }

  /// @brief factory method
  static Expected<ChemSpeciesProps, Error> from_formula(
      const parse::FormulaComps& comps, const NuclideModel& model) {
    ChemSpeciesProps out;
    out.mass_factor = 0.0;
    int n_protons = 0;
    for (int i = 0; i < comps.subscript_symbol_pairs.size(); i++) {
      auto [coef, name] = comps.subscript_symbol_pairs[i];

      // lookup the id associated with the nuclide specified by name
      std::optional<uint16_t> maybe_id = model.symbol_map().find(name);
      if (!maybe_id.has_value()) {
        std::string name_str(name);
        return Unexpected(Error::msgf("\"%s\" is not a known chemical species",
                                      name_str.c_str()));
      }
      uint16_t nuclide_id = maybe_id.value();

      // retrieve the nuclide's properties
      const NuclideProp& nuclide_prop = model.get(nuclide_id);

      // now, let's update total number of protons and update out
      n_protons += static_cast<int>(coef) * nuclide_prop.n_proton;
      out.mass_factor += static_cast<double>(coef) * nuclide_prop.mass_factor;
      if (!out.coef_nuclide_id_pairs.try_emplace_back(coef, nuclide_id)) {
        // shouldn't be reachable since that max capacity should be identical
        return Unexpected(Error::msgf(
            "coef_nuclide_id_pairs doesn't have adequate capacity to hold "
            "data for nuclide symbol %d",
            i));
      }
    }
    out.elementary_charge = comps.roman_numeral_val - 1;

    // perform one last sanity check (perhaps we should move it to a unit test)
    if ((out.elementary_charge - n_protons) > 0) {
      return Unexpected(Error::msgf(
          "species can't have a charge of +%d (it only contains %d protons)",
          out.elementary_charge, n_protons));
    }
    return out;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_CHEMSP_PROPS_HPP