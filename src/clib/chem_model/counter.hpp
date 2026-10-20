//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// declare/implement @ref NuclideAndChargeCounter
///
//===----------------------------------------------------------------------===//

#ifndef CHEM_MODEL_COUNTER_HPP
#define CHEM_MODEL_COUNTER_HPP

#include <map>
#include "chem_model/species.hpp"

namespace GRIMPL_NAMESPACE_DECL {

inline bool is_prim_or_metal_(uint16_t species_id,
                              const PartMap<SpKind>& kind_map) {
  std::optional<GRIMPL_NS::partmap::IdxSearch<SpKind>> tmp =
      kind_map.search_idx(species_id);
  if (!tmp.has_value()) {
    return false;
  }
  SpKind kind = tmp.value().pd;
  return kind == SpKind::PRIMORDIAL || kind == SpKind::METAL;
}

/// @brief A counter of nuclides and charge from elements in an equation
///
/// At the moment, this is just used for debugging purposes, but in the near
/// future, we'll start using this for a purpose in @ref SpeciesBuilder
class NuclideAndChargeCounter {
  /// tracks the total charge
  int charge_;
  /// maps nuclide ids to counts
  std::map<uint16_t, int> nuc_counts_;
  /// encodes properties about chemical species
  const ChemSpeciesPropsVec& prop_vec_;
  /// describes the kind of species that a species_id corresponds to
  const PartMap<SpKind>& kind_map_;

public:
  NuclideAndChargeCounter(const ChemSpeciesPropsVec& prop_vec,
                          const PartMap<SpKind>& kind_map)
      : charge_{0}, prop_vec_{prop_vec}, kind_map_{kind_map} {}

  const int& charge() const { return charge_; }
  const std::map<uint16_t, int>& nuclide_counts() const { return nuc_counts_; }

  /// @brief adds contributions from a species to the counter
  Expected<void, Error> add_contrib(uint16_t species_id, int factor = 1) {
    if (!is_prim_or_metal_(species_id, kind_map_)) {
      return Unexpected(
          Error::msg_literal("requires primordial or metal species"));
    }
    const ChemSpeciesProps& props = prop_vec_[species_id];
    charge_ += factor * props.elementary_charge;

    for (const IdCoefPair& pair : props.nuclide_id_coef_pairs) {
      auto [nuclide_id, subscript] = pair;
      // nuc_counts_[nuclide_id] is defaulted to 0 if it didn't exist
      nuc_counts_[nuclide_id] += factor * (int)(subscript);
    }
    return Expected<void, Error>();  // <- success!
  }

  /// @brief add contributions from multiple species-coefficent pairs
  template <typename T>
  Expected<void, Error> add_contribs(T species_id_coef_pairs, int factor = 1) {
    static_assert(std::is_same_v<typename T::value_type, IdCoefPair>,
                  "invalid species_id_coef_pairs type");
    for (const IdCoefPair& pair : species_id_coef_pairs) {
      Expected<void, Error> tmp = add_contrib(pair.id, factor * (int)pair.coef);
      if (!tmp.has_value()) {
        return Unexpected(
            tmp.error().contextf("problem for species_id = %d", (int)pair.id));
      }
    }
    return Expected<void, Error>();  // <- success!
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_COUNTER_HPP