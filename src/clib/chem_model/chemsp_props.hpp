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

struct IdCoefPair {
  uint16_t id;
  uint8_t coef;
};

/// @brief An internal data structure used to track chemical species
///
/// @note We could potentially track other properties like the actual particle
/// mass and whether all particles all truly have one mass or if this is really
/// a mixed ensemble (due to the presence of isotopes). In practice, those
/// other properties don't actually seem all that useful
struct ChemSpeciesProps {
  static constexpr std::size_t MaxSymsPerFormula = parse::MaxSymsPerFormula;

  // define the data-members:
  // ------------------------
  /// tracks coefficient-nuclide pairs (if any) that make up the species
  InplaceVector<IdCoefPair, MaxSymsPerFormula> nuclide_id_coef_pairs;

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
  bool composed_of_nuclides() const { return nuclide_id_coef_pairs.size() > 0; }

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
      std::string_view formula, const NuclideModel& model);

  /// @brief factory method
  static Expected<ChemSpeciesProps, Error> from_formula(
      const parse::FormulaComps& comps, const NuclideModel& model);
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_CHEMSP_PROPS_HPP