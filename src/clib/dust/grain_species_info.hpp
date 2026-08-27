//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares GrainSpeciesInfo as well as other associated types and functions
///
//===----------------------------------------------------------------------===//

#ifndef GRAIN_SPECIES_INFO_HPP
#define GRAIN_SPECIES_INFO_HPP

#include "../support/FrozenKeyIdxBiMap.hpp"
#include "../support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// holds information about a single gas species that is an ingredient for
/// grain growth
struct GrainGrowthIngredient {
  /// the positive integer coefficient associated with the ingredient
  int coef;

  /// species index of the ingredient
  int species_idx;

  /// the particle mass of the ingredient species (in atomic mass units)
  ///
  /// @note
  /// It's plausible that this should actually be the particle mass in hydrogen
  /// masses, which is slightly different from atomic mass units.
  /// Unfortunately, the original code didn't specify constants with enough
  /// precision to tell the difference
  double mparticle_amu;
};

/// summarizes details about a single grain species
struct GrainSpeciesInfoEntry {
  /// the species index of the grain in the #GrainSpLUT lookup table
  int species_idx;

  /// indicates whether to use the carbonaceous or silicate coefficient table
  /// to computing contributions of the grain species to the total h2dust rate
  /// (or the rate of H2 formation)
  ///
  /// @note
  /// The choice to encode this information within each GrainSpeciesInfoEntry
  /// may be somewhat inefficient since (at the time of writing) Amorphous
  /// Carbon grains are the **ONLY** grain species that use the carbonaceous
  /// coefficient table & all other species use the silicate coefficient table
  /// (to be honest, the correctness of this seems questionable, at best). As
  /// soon as more than one grain uses the carbonaceous coefficient table (or
  /// some other non sillicate table), it will be clearly be much more optimal
  /// to track the information in this manner.
  bool h2dust_uses_carbonaceous_table;

  /// The sublimation temperature in units of Kelvin
  double sublimation_temperature;

  /// specifies the density of a single grain in units of g/cm^3
  ///
  /// @note
  /// The values are consistent with the values quoted within table 2 of
  /// [Chiaki+ 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.2659C)
  /// assuming that all grains are spherical.
  double bulk_density_cgs;

  /// The number of growth ingredients
  int n_growth_ingredients;

  /// The array of growth ingredients
  const GrainGrowthIngredient* growth_ingredients;
};

/// tracks for details about every relevant grain species for a given Grackle
/// configuration. This includes information necessary for computing growth &
/// destruction rates
///
/// Relationship with OnlyGrainSpLUT
/// --------------------------------
/// In the short term, the index of each species in the
/// @ref GrainSpeciesInfo::species_info() is dictated by the
/// order of enumerators in the OnlyGrainSpLUT enumeration.
///
/// In the medium term, we plan to entirely eliminate the OnlyGrainSpLUT
/// enumeration because all of the grain species can be treated very uniformly.
/// At the time of writing, just about every place where we would use
/// OnlyGrainSpLUT corresponds to a location where we would enumerate every
/// possible grain species and perform nearly identical operations on each
/// species. In each case, it is straight-forward to replace these blocks of
/// logic with for-loops (we just need to encode species-specific variations in
/// the calculations in species_info that have the same ordering as the
/// species). To phrase it another way, in nearly all of the places where we
/// would use OnlyGrainSpLUT, we don't need to know the grain species identity.
///
/// The exception to this is when we compute the h2dust rate. In this case
/// we need to identify AC_dust since we need to do something **slightly**
/// different from the other grains, but this is easy to work around
class GrainSpeciesInfo {
  /// number of grain species considered for the current Grackle configuration
  int n_species_;

  /// holds @ref n_species entries. Each entry holds info about a separate
  /// grain species
  GrainSpeciesInfoEntry* species_info_;

  /// maps between grain species names and the associated index. The mapping is
  /// **ALWAYS** consistent with ``OnlyGrainSpLUT``.
  ///
  /// @note
  /// An argument could be made for storing this separately from the rest of
  /// the struct since the core grackle calculations don't (or at least
  /// shouldn't) use this data structure during the calculation.
  FrozenKeyIdxBiMap name_map_;

private:  // helper methods
  /// @brief helper method to deallocate a species_info array
  ///
  /// @note We could get rid of this by converting @ref GrainSpeciesInfoEntry
  /// to a proper class, or by using std::unique_ptr and encoding the cleanup
  /// of grow_ingredients in a custom deleter (this would be a little tricky
  /// since the deleter would need to know the array's length)
  static void cleanup_array_(int n_species,
                             GrainSpeciesInfoEntry* species_info) {
    if (n_species > 0) {
      for (int gsp_idx = 0; gsp_idx < n_species; gsp_idx++) {
        if ((species_info[gsp_idx].growth_ingredients) != nullptr) {
          delete[] species_info[gsp_idx].growth_ingredients;
        }
      }
      delete[] species_info;
    }
  }

public:
  /// @brief checks whether instance is valid
  explicit operator bool() const { return n_species_ > 0; }

  /// @brief number of grain species considered in current Grackle configuration
  int n_species() const { return n_species_; }

  /// @brief returns sequence of entries describing each grain species
  const GrainSpeciesInfoEntry* species_info() const { return species_info_; }

  /// @brief returns mapping between grain species names and associated indices
  const FrozenKeyIdxBiMap& name_map() const { return name_map_; }

  /// @brief Primary Constructor
  ///
  /// It is the caller's responsibility to check whether the resulting object
  /// is valid (e.g. by checking `if (obj)`).
  ///
  /// @note
  /// In the future, we could use a factory method that returns a std::optional
  /// or a C++23's std::expected. This would let us ensure that instance of
  /// this class only exists if it's valid
  explicit GrainSpeciesInfo(int dust_species_parameter);

  // the following are disabled because the default implementations won't
  // properly handle FrozenKeyIdxBiMap (since it doesn't act like a class) or
  // species_info
  GrainSpeciesInfo(const GrainSpeciesInfo&) = delete;
  GrainSpeciesInfo(GrainSpeciesInfo&&) = delete;
  GrainSpeciesInfo& operator=(const GrainSpeciesInfo&) = delete;
  GrainSpeciesInfo& operator=(GrainSpeciesInfo&&) = delete;

  ~GrainSpeciesInfo() {
    if (n_species_ > 0) {
      GrainSpeciesInfo::cleanup_array_(n_species_, species_info_);
    }
  }
};

/// return the number of grain species
///
/// @param[in] dust_species_parameter The parameter tracked by #chemistry_data
/// @return The number of grain species. A negative value indicates that the
///     @p dust_species_parameter has an invalid value
inline int get_n_grain_species(int dust_species_parameter) {
  switch (dust_species_parameter) {
    case 0:
      return 0;
    case 1:
      return 2;
    case 2:
      return 10;
    case 3:
      return 13;
    default:
      return -1;
  }
}

/// Species the maximum number of ingredients that a grain species has
///
/// @note
/// The correctness of this constant is explicitly checked in a unit test
inline constexpr int max_ingredients_per_grain_species = 3;

}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* GRAIN_SPECIES_INFO_HPP */
