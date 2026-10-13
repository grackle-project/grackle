//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// implement @ref RxnEqnProps
///
//===----------------------------------------------------------------------===//

#include <cstdint>

#include "chem_model/chemsp_props.hpp"
#include "chem_model/parse.hpp"
#include "chem_model/rxn_eqn_props.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/inplace_vector.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<RxnEqnProps, Error> RxnEqnProps::from_eqn_comps(
    const parse::RxnEqnComps& comps,
    const FrozenKeyIdxBiMap& species_name_map) {
  // the object that we are going to build
  RxnEqnProps out;

  bool handle_lhs_choices[] = {true, false};
  for (bool is_lhs : handle_lhs_choices) {
    const InplaceVector<parse::SViewCoefPair, parse::MaxSpeciesPerEqnSide>&
        name_coef_pairs =
            (is_lhs) ? comps.lhs_sp_coef_pairs : comps.rhs_sp_coef_pairs;
    InplaceVector<IdCoefPair, parse::MaxSpeciesPerEqnSide>& id_coef_pairs =
        (is_lhs) ? out.lhs_entries : out.rhs_entries;
    const char* eqn_side_name = (is_lhs) ? "left" : "right";

    std::size_t n_entries = name_coef_pairs.size();
    for (std::size_t i = 0; i < n_entries; i++) {
      std::string_view entity_name = name_coef_pairs[i].first;
      uint8_t coef = name_coef_pairs[i].second;

      // look up the species id
      std::optional<uint16_t> maybe_id = species_name_map.find(entity_name);
      if (!maybe_id.has_value()) {
        if (entity_name == "photon") {
          return Unexpected(Error::msg_literal(
              "no support for chemical equations involving \"photon\" (yet)"));
        }
        std::string tmp(entity_name);
        return Unexpected(
            Error::msgf("\"%s\" not in species_name_map", tmp.c_str()));
      }
      uint16_t species_id = maybe_id.value();

      // sanity check:
      // - in the future, we may choose to "fix things" ourself
      // - we should consider handling this during parsing (whether we want to
      //   fix things or consider it an error, this takes polynomial time and
      //   we may want to handle it @ compile-time)
      for (std::size_t j = 0; j < i; i++) {
        if (id_coef_pairs[j].id == species_id) {
          std::string tmp(entity_name);
          return Unexpected(
              Error::msgf("\"%s\" species appears multiple times on %s side of "
                          "chemical eqn",
                          tmp.c_str(), eqn_side_name));
        }
      }

      if (!id_coef_pairs.try_emplace_back(species_id, coef)) {
        // shouldn't be reachable since that max capacity should be identical
        return Unexpected(Error::msgf(
            "the buffer that holds pairs of coefficient species-id for the "
            "%s side of the chemical equation doesn't have adequate capactiy "
            "entry with index %zu",
            eqn_side_name, i));
      }
    }
  }
  return out;
}

}  // namespace GRIMPL_NAMESPACE_DECL