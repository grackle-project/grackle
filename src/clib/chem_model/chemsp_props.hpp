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
#include <map>

#include "chem_model/nuclide_model.hpp"
#include "chem_model/parse.hpp"
#include "phys_constants.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/inplace_vector.hpp"
#include "support/status_reporting.hpp"

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

// forward declaration of SpeciesBuilder
class SpeciesBuilder;

/// @brief a vector like-container for accessing a @ref ChemSpeciesProps
///
/// Essentially, this lives alongside a @ref FrozenKeyIdxBiMap that is in
/// charge of mapping (bidirectionally) between species names and integer ids.
/// This container maps between the integer and chemical species properties.
///
/// This looks pretty silly (that's because it is silly).
/// - Essentially, it implements a random-access method for accessing a
///   @ref ChemSpeciesProps object
/// - In the future, the intention is just to replace this with a
///   ``std::vector<ChemSpeciesProps>``
/// - But we can't do that until we completely remove @ref SpLUT from the
///   codebase since there may be some unused species ids (e.g. if using
///   `primordial_chemistry < 4` and `metal_chemistry==1`)
class ChemSpeciesPropsVec {
  // declare SpeciesBuilder as a friend (it constructs objects of this type)
  friend SpeciesBuilder;

  using MapType_ = std::map<std::size_t, ChemSpeciesProps>;

  // actually tracks the data
  std::map<std::size_t, ChemSpeciesProps> m_;

public:
  /// @brief Overloads array subscript operation
  ///
  /// The caller is expected to know ahead of time whether an index is valid.
  /// We return a fatal error for debugging purposes (in the future, when this
  /// type is converted to a regular std::vector<ChemSpeciesProps>, this error
  /// check will no longer happen
  ///
  /// @note For less-experienced C++ developers, this is like python's
  /// `__getitem__` magic method
  const ChemSpeciesProps& operator[](std::size_t idx) const {
    std::map<std::size_t, ChemSpeciesProps>::const_iterator tmp = m_.find(idx);
    if (tmp == m_.cend()) {
      GR_INTERNAL_ERROR("programming error: %zu is an invalid index", idx);
    }
    return tmp->second;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_CHEMSP_PROPS_HPP