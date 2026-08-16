//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares NuclideModel
///
//===----------------------------------------------------------------------===//
#ifndef NUCLIDE_MODEL_HPP
#define NUCLIDE_MODEL_HPP

#include <vector>

#include "../ratequery.hpp"
#include "../support/config.hpp"
#include "../support/FrozenKeyIdxBiMap.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// wraps mass number (i.e. number of protons and neutrons) for most
/// commonly occuring isotope associated with a nuclide symbol.
///
/// This ONLY exists for book-keeping purposes (i.e. to clearly express where
/// information used to construct a @ref NuclideProp originates from)
struct MostCommonMassNumber {
  int val;
};

/// @brief Aggregates properties with a single nuclide symbol
///
/// See @ref NuclideModel for more details.
///
/// @note
/// A lot of these properties only exist for explanatory purposes. The most
/// important data members include @ref n_proton, @ref primordial, and
/// @ref mass_factor.
struct NuclideProp {
  /// the number of protons
  int n_proton;
  /// tracks whether the species is primordial or not
  bool primordial;
  /// tracks whether the symbol represents a single isotope
  bool single_isotope;
  /// tracks the average mass in daltons (aka unified atomic mass units)
  double mass_Da;
  /// tracks the mass_factor
  ///
  /// keep this consistent with constants in the @ref mass_factor namespace
  double mass_factor;

  // forbid (partial) improper initialization
  NuclideProp() = delete;

  /// @brief primary constructor
  ///
  /// The @p mass_factor argument is intended to help synchronize the values
  /// with the constants in the @ref mass_factor namespace.
  ///
  /// @note
  /// This has been declared constexpr in order to allow us to construct
  /// objects at compile-time
  constexpr NuclideProp(int n_proton, bool primordial, bool single_isotope,
                        double mass_Da, double mass_factor) noexcept
      : n_proton{n_proton},
        primordial{primordial},
        single_isotope(single_isotope),
        mass_Da{mass_Da},
        mass_factor{mass_factor} {};

  /// @brief constructor that infers mass_factor from most common mass number
  constexpr NuclideProp(int n_proton, bool primordial, bool single_isotope,
                        double mass_Da, MostCommonMassNumber num) noexcept
      : NuclideProp(n_proton, primordial, single_isotope, mass_Da, num.val) {}
};

/// @brief map-like type that associates properties with each "nuclide symbol"
///        (modelled in the context of Grackle's chemical networks)
///
/// In more detail, Grackle uses standard 1 to 2 character symbols to represent
/// groups of 1 or more isotopes of a given element. Grackle only originally
/// tracked 3 nuclide symbols: `H`, `He` and `D`.
/// - `H` and `D` are noteworthy cases because these symbols each map to a
///   single kind of nuclide (Hydrogen without neutrons and Deuterium).
/// - In contrast, other symbols (`He` and newer additions) refer to groupings
///   of nuclides that are all isotopes of each other. For example, `He`
///   corresponds to all nuclides containing 2 protons, and `C` corresponds to
///   all nuclides with 6 protons.
///
/// Grackle generally acts as if the ensemble of particles associated with a
/// nuclide symbol are all identical and that each property of the particle is
/// the average of that properties for each isotope in the abundance (weighted
/// by abundance).
///
/// Instances of this type can be used to map a symbol to a unique integer
/// index. It can also be used to direclty query the other properties.
///
/// @note
/// This type is "read-only" (i.e. once it's constructed, the associated
/// functions will not mutate it). While the constructor is currently hardcoded
/// to construct an instance with Grackle's standard assumptions, part of the
/// motivation for implementing the logic as a class is so that we have
/// flexibility to use different assumptions.
class NuclideModel {
  /// maps between Nuclide symbols and id
  FrozenKeyIdxBiMap symbol_map_;
  /// pointer to the array of nuclide properties
  std::vector<NuclideProp> props_;

public:
  /// @brief Primary constructor
  ///
  /// At this time, there's no reason for us to support construction of nuclide
  /// models with custom properties (but we could revisit that in the future)
  NuclideModel();

  const FrozenKeyIdxBiMap& symbol_map() const noexcept { return symbol_map_; }

  const NuclideProp* try_get(const char* symbol) const noexcept {
    std::optional<uint16_t> tmp = symbol_map_.find(symbol);
    return tmp.has_value() ? &props_[*tmp] : nullptr;
  }

  const NuclideProp* try_get(int idx) const noexcept {
    bool valid_idx = (0 <= idx) && (idx < symbol_map_.size());
    return (valid_idx) ? &props_[idx] : nullptr;
  }

  // we don't define an analogous overload for passing a symbol, since you
  // can't tell ahead of time whether a symbol will be valid
  const NuclideProp& get(int idx) const { return props_[idx]; }

  int size() const noexcept { return symbol_map_.size(); }

  /// @brief copy relevant information to make it queryable ratequery
  ///
  /// This mainly exists for exposing this information to gracklepy. In
  /// reality, we should probably make public API functions to expose relevant
  /// mass_factor information to end-users. But, I'm not so sure we actually
  /// need to expose mass_factor values for nuclide symbols (I think we
  /// actually want to directly expose mass_factor values for chemical species).
  ///
  /// We elect to copy values to @p reg_builder (as opposed to defining recipes
  /// for dynamically accessing values stored within the current object) for 2
  /// reasons:
  /// 1. there isn't an obvious reason for us to keep NuclideProp around after
  ///    we initially configure Grackle (except perhaps that we want to query
  ///    this information)
  /// 2. NuclideModel currently organizes data as an array of structs. In order
  ///    to make the data available through a recipe, we would need to
  ///    internally repack the data as a struct of arrays (this would involve
  ///    making methods for every queryable property)
  ///
  /// @note
  /// It's an error to call this more than once for a given reg_builder
  int copy_info_to_RegBuilder(ratequery::RegBuilder& reg_builder) const;
};

}  // namespace GRIMPL_NAMESPACE_DECL
#endif  // NUCLIDE_MODEL_TABLE_HPP