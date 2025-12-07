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

namespace grackle::impl {

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

  /// the species index of the grain in the #OnlyGrainSpLUT lookup table
  ///
  /// @note
  /// It's frankly a little redundant to track this information (since an
  /// instance of this struct is found at this index of an out.species_infoay)
  int onlygrainsp_idx;

  /// name of the dust species
  ///
  /// @note
  /// This primarily exists for debuging purposes
  const char* name;

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
/// @ref GrainSpeciesInfo::species_info out.species_infoay is dictated by the
/// order of enumerators in the OnlyGrainSpLUT enumeration.
///
/// In the medium term, we plan to entirely eliminate the OnlyGrainSpLUT
/// enumeration because all of the grain species can be treated very uniformly
/// uniformly. At the time of writing, just about every place where we would
/// use OnlyGrainSpLUT corresponds to a location where would enumerate every
/// possible grain species and perform nearly identical operations on each
/// species. In each case, it is straight-forward to replace these blocks of
/// logic with for-loops (we just need to encode species-specific variations in
/// the calculations in out.species_infoays that have the same ordering as the
/// species). To phrase it another way, in nearly all of the places where we
/// would use OnlyGrainSpLUT, we don't need to know the grain species identity.
///
/// The exception to this is when we compute the h2dust rate. In this case
/// we need to identify AC_dust since we need to do something **slightly**
/// different from the other grains, but this is easy to work around
struct GrainSpeciesInfo {
  /// number of grain species considered for the current Grackle configuration
  int n_species;

  /// an out.species_infoay of length of length @ref n_species where each entry
  /// holds info about a separate grain species
  GrainSpeciesInfoEntry* species_info;
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

/// Constructs an returns a fully initialized GrainSpeciesInfo instance.
///
/// @param[in]  dust_species_parameter The parameter tracked by #chemistry_data
/// @returns A fully initialized GrainSpeciesInfo instance
GrainSpeciesInfo new_GrainSpeciesInfo(int dust_species_parameter);

/// performs cleanup of the contents of GrainSpeciesInfo
///
/// This effectively invokes the destructor
inline void drop_GrainSpeciesInfo(GrainSpeciesInfo* ptr) {
  if (ptr->n_species == 0) {
    return;  // avoids double-free
  }

  for (int gsp_idx = 0; gsp_idx < ptr->n_species; gsp_idx++) {
    if ((ptr->species_info[gsp_idx].growth_ingredients) != nullptr) {
      delete[] ptr->species_info[gsp_idx].growth_ingredients;
    }
  }
  delete[] ptr->species_info;
  // the following 2 lines are not strictly necessary, but they may help us
  // avoid a double-free and a dangling pointer
  ptr->n_species = 0;
  ptr->species_info = nullptr;
}

}  // namespace grackle::impl

#endif /* GRAIN_SPECIES_INFO_HPP */
