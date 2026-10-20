//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// declare @ref RxnEqnProps
///
//===----------------------------------------------------------------------===//

#ifndef CHEM_MODEL_CHEM_EQN_PROPS_HPP
#define CHEM_MODEL_CHEM_EQN_PROPS_HPP

#include <cstdint>
#include "chem_model/chemsp_props.hpp"
#include "chem_model/parse.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/inplace_vector.hpp"

namespace GRIMPL_NAMESPACE_DECL {

enum class EqnSide { LEFT, RIGHT };

inline const char* to_string(EqnSide v) {
  return (v == EqnSide::LEFT) ? "left" : "right";
}

/// @brief An internal data structure used to track properties about a reaction
///     equation
struct RxnEqnProps {
  using InplaceIdCoefPairVec =
      InplaceVector<IdCoefPair, parse::MaxSpeciesPerEqnSide>;

  /// species_id, coefficient pairs for the left hand side of the equation
  InplaceIdCoefPairVec lhs_entries;
  /// species_id, coefficient pairs for the right hand side of the equation
  InplaceIdCoefPairVec rhs_entries;

  /// tracks whether "photon" was in the equation (and what side it was on)
  std::optional<EqnSide> has_photon;

private:
  RxnEqnProps() = default;  // <- only used by factory methods

public:
  /// @brief Primary Factory Method
  ///
  /// Constructs a @ref RxnEqnProps by parsing a string
  ///
  /// @param s The string to be parsed
  /// @param species_name_map Maps between species names and ids
  static Expected<RxnEqnProps, Error> from_string(
      std::string_view s, const FrozenKeyIdxBiMap& species_name_map) {
    Expected<parse::RxnEqnComps, parse::Err> parse_rslt = parse::rxn_eqn(s);
    if (!parse_rslt.has_value()) {
      return Unexpected(Error::parse(s, parse_rslt.error())
                            .context_literal("problem parsing chemical eqn"));
    }
    return RxnEqnProps::from_eqn_comps(parse_rslt.value(), species_name_map);
  }

  /// Factory Method
  static Expected<RxnEqnProps, Error> from_eqn_comps(
      const parse::RxnEqnComps& comps,
      const FrozenKeyIdxBiMap& species_name_map);

  /// @brief infers the nuclide ratios for a given species in the reaction
  ///
  /// This resulting object maps nuclide ids to coefficients. The ratio between
  /// a pair of the resulting coefficients give the abundance ratios (by
  /// number) for the associated nuclides in a given species.
  ///
  /// There are 2 strong guarantees about the returned map:
  /// 1. the map holds at least 1 item
  /// 2. all the coefficients are positive
  ///
  /// As a consequence of the first guarantee, this function will report an
  /// error for a species not composed of any nuclides (e.g. a free electron)
  ///
  /// @note
  /// The actual value/normalization of the coefficients should be considered
  /// an arbitrary implementation detail that may change based on what is most
  /// convenient (or if necessary, most expedient).
  ///
  /// Implementation Note
  /// -------------------
  /// The choice to return a ``std::map`` was made purely out of convenience.
  /// Since the number of types nuclides per species is small, it may actually
  /// be faster to track the mappings in a sequential container (e.g. a vector)
  /// and use linear search
  Expected<std::map<uint16_t, int>, Error> calc_nuclide_ratio_map(
      uint16_t target_species_id, const ChemSpeciesPropsVec& prop_vec) const;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_CHEM_EQN_PROPS_HPP