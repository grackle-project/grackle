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

#include <cstdint>
#include <map>
#include "chem_model/species.hpp"
#include "chem_model/rxn_eqn_props.hpp"
#include "support/expected.hpp"

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
/// This is defined for debugging purposes and to help implement
/// @ref calc_nuclide_ratio_map
class NuclideAndChargeCounter {
  /// tracks the total charge
  int charge_;
  /// maps nuclide ids to counts
  std::map<uint16_t, int> n_counts_;
  /// encodes properties about chemical species
  const ChemSpeciesPropsVec& prop_vec_;
  /// describes the kind of species that a species_id corresponds to
  const PartMap<SpKind>& kind_map_;

public:
  NuclideAndChargeCounter(const ChemSpeciesPropsVec& prop_vec,
                          const PartMap<SpKind>& kind_map)
      : charge_{0}, prop_vec_{prop_vec}, kind_map_{kind_map} {}

  const int& charge() const { return charge_; }
  const std::map<uint16_t, int>& nuclide_counts() const& { return n_counts_; }
  std::map<uint16_t, int>&& nuclide_counts() && { return std::move(n_counts_); }

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
      // n_counts_[nuclide_id] is defaulted to 0 if it didn't exist
      n_counts_[nuclide_id] += factor * (int)(subscript);
    }
    return Expected<void, Error>();  // <- success!
  }

  /// @brief add contributions from multiple species-coefficent pairs
  ///
  /// The function can be configured to skip a species_id. If successful, the
  /// function indicates whether a species was skipped
  template <typename T>
  Expected<bool, Error> add_contribs(
      T species_id_coef_pairs, int factor = 1,
      std::optional<uint16_t> skip_id = std::nullopt) {
    static_assert(std::is_same_v<typename T::value_type, IdCoefPair>,
                  "invalid species_id_coef_pairs type");
    // note: no species_id may exceed kind_map.n_idx()
    uint16_t unwrapped_skip = skip_id.value_or((uint16_t)kind_map_.n_idx());
    bool encountered_skip_id = false;

    for (const IdCoefPair& pair : species_id_coef_pairs) {
      if (pair.id == unwrapped_skip) {
        encountered_skip_id = true;
        continue;
      }
      Expected<void, Error> tmp = add_contrib(pair.id, factor * (int)pair.coef);
      if (!tmp.has_value()) {
        return Unexpected(
            tmp.error().contextf("problem for species_id = %d", (int)pair.id));
      }
    }
    return encountered_skip_id;
  }
};

/// @brief infers nuclide ratios for a given species from a reaction equation
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
/// Since the number of types nuclides per species is small, it may actually be
/// faster to track the mappings in a sequential container (e.g. a vector)
/// and use linear search
inline Expected<std::map<uint16_t, int>, Error> calc_nuclide_ratio_map(
    const RxnEqnProps& eqn_props, uint16_t target_species_id,
    const ChemSpeciesPropsVec& prop_vec, const PartMap<SpKind>& kind_map) {
  NuclideAndChargeCounter counter(prop_vec, kind_map);

  // add contributions from lhs & rhs (skipping over over the target species),
  Expected<bool, Error> lhs =
      counter.add_contribs(eqn_props.lhs_entries, +1, target_species_id);
  Expected<bool, Error> rhs =
      counter.add_contribs(eqn_props.rhs_entries, -1, target_species_id);
  if (!lhs.has_value()) {
    return Unexpected(lhs.error());
  } else if (!rhs.has_value()) {
    return Unexpected(rhs.error());
  }
  bool target_on_LHS = lhs.value();
  bool target_on_RHS = rhs.value();

  // extract the map, then erase entries in map where coefficient == 0
  std::map<uint16_t, int> map(std::move(counter).nuclide_counts());
  std::erase_if(map, [](const auto& kv_pair) { return kv_pair.second == 0; });

  bool any_negative_coef = false;
  for (auto& kv_pair : map) {
    if (target_on_LHS) {
      kv_pair.second *= -1;
    }
    any_negative_coef = any_negative_coef || (kv_pair.second < 0);
  }

  if (target_on_LHS && target_on_RHS) {
    return Unexpected(Error::msg_literal("target on both sides of eqn"));
  } else if (!target_on_LHS && !target_on_RHS) {
    return Unexpected(Error::msg_literal("target not in eqn"));
  } else if (map.empty()) {
    return Unexpected(Error::msg_literal("target not composed of nuclides"));
  } else if (any_negative_coef) {
    return Unexpected(Error::msg_literal("there's a negative coefficient"));
  }
  return map;
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_COUNTER_HPP