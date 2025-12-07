//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines details about the dust grain species
///
//===----------------------------------------------------------------------===//

#ifndef GRAIN_SPECIES_INFO_HPP
#define GRAIN_SPECIES_INFO_HPP

#include "LUT.hpp"
#include "grackle_macros.h"
#include <cstring>  // memcpy

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

// define a few macros to help implement new_GrainSpeciesInfo
// - we undef all of these macros after they're used for cleanliness
// - we may want to consider an alternative that avoids macros in the future
#define GRIMPL_INGREDIENT_LIST_SENTINEL {-1, -1, -1.0}

// does the heavy-lifting for implementing GR_MK_GRAIN_INFO_ENTRY
inline GrainSpeciesInfoEntry mk_gsp_info_entry_helper_(
    int species_idx, int onlygrainsp_idx, const char* name,
    bool h2dust_uses_carbonaceous_table, double sublimation_temperature,
    double bulk_density_cgs, const GrainGrowthIngredient* growth_ingredients) {
  GrainGrowthIngredient* out_ingredient_ptr = nullptr;
  int n_ingredients = 0;

  if (growth_ingredients != nullptr) {
    n_ingredients = 0;
    // increment the counter until we reach the sentinel
    while (growth_ingredients[n_ingredients].coef != -1) {
      n_ingredients++;
    }
    // allocate and initialize the pointer
    out_ingredient_ptr = new GrainGrowthIngredient[n_ingredients];
    std::memcpy(out_ingredient_ptr, growth_ingredients,
                (std::size_t)n_ingredients);
  }

  return GrainSpeciesInfoEntry{species_idx,
                               onlygrainsp_idx,
                               name,
                               h2dust_uses_carbonaceous_table,
                               sublimation_temperature,
                               bulk_density_cgs,
                               n_ingredients,
                               out_ingredient_ptr};
}

/// fills the specified out.species_infoay with the number of grain species.
///
/// The out.species_infoay is expected to have space for the number of entries
/// returned by get_n_grain_species(dust_species_paramter)
///
/// @param[in]  dust_species_parameter The parameter tracked by #chemistry_data
/// @returns A fully initialized GrainSpeciesInfo instance
inline GrainSpeciesInfo new_GrainSpeciesInfo(int dust_species_parameter) {
  GrainSpeciesInfo out{-1, nullptr};  // indicates an error
  out.n_species = get_n_grain_species(dust_species_parameter);

  if (out.n_species <= 0) {
    return out;
  }

  out.species_info = new GrainSpeciesInfoEntry[out.n_species];

  // At the time of writing:
  // - we **only** use h2rate_carbonaceous_coef_table for the AC_dust
  //   species (amorphous carbon grains)
  // - we use h2rate_silicate_coef_table for **ALL** other grain species
  //   (including species that are obviously NOT silicates)
  //
  // It is not obvious at all why this is the case...
  //
  // TODO: we should add documentation clearly addressing each of the
  // following questions and replace this todo-item with a reference to the
  // part of the documentation that discusses this! The documentation
  // should make sure to address:
  // - if this choice physically motivated or a common convention? (or if
  //   it was a quick & dirty choice because Gen didn't know what to do)
  // - do we have a sense for how "wrong" the choice is? (e.g. will it
  //   generally overpredict/underpredict?)
  // - why don't we use the generic my_rates->h2dusta rate as a generic
  //   fallback for non-silicate grains instead of using h2dustS? I realize
  //   that h2dusta has very different units...

  if (dust_species_parameter > 0) {
    const GrainGrowthIngredient MgSiO3_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::Mg, 24.},
        {1, SpLUT::SiOI, 44.},
        {2, SpLUT::H2O, 18.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[0] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::MgSiO3_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::MgSiO3_dust,
        /* name = */ "MgSiO3_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1222.0,
        /* bulk_density_cgs = */ 3.20185,
        /* growth_ingredients = */ MgSiO3_dust_ingred);

    const GrainGrowthIngredient AC_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::CI, 12.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[1] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::AC_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::AC_dust,
        /* name = */ "AC_dust",
        /* h2dust_uses_carbonaceous_table = */ true,
        /* sublimation_temperature = */ 1800.0,
        /* bulk_density_cgs = */ 2.27949,
        /* growth_ingredients = */ AC_dust_ingred);
  }

  if (dust_species_parameter > 1) {
    const GrainGrowthIngredient SiM_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::SiI, 28.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[2] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::SiM_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::SiM_dust,
        /* name = */ "SiM_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 2.34118,
        /* growth_ingredients = */ SiM_dust_ingred);

    const GrainGrowthIngredient FeM_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::Fe, 56.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[3] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::FeM_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::FeM_dust,
        /* name = */ "FeM_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 7.95995,
        /* growth_ingredients = */ FeM_dust_ingred);

    const GrainGrowthIngredient Mg2SiO4_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {2, SpLUT::Mg, 24.},
        {1, SpLUT::SiOI, 44.},
        {3, SpLUT::H2O, 18.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[4] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::Mg2SiO4_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::Mg2SiO4_dust,
        /* name = */ "Mg2SiO4_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1277.0,
        /* bulk_density_cgs = */ 3.22133,
        /* growth_ingredients = */ Mg2SiO4_dust_ingred);

    const GrainGrowthIngredient Fe3O4_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {3, SpLUT::Fe, 56.},
        {4, SpLUT::H2O, 18.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[5] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::Fe3O4_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::Fe3O4_dust,
        /* name = */ "Fe3O4_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 5.25096,
        /* growth_ingredients = */ Fe3O4_dust_ingred);

    const GrainGrowthIngredient SiO2_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::SiO2I, 60.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[6] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::SiO2_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::SiO2_dust,
        /* name = */ "SiO2_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 2.66235,
        /* growth_ingredients = */ SiO2_dust_ingred);

    const GrainGrowthIngredient MgO_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::Mg, 24.},
        {1, SpLUT::H2O, 18.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[7] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::MgO_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::MgO_dust,
        /* name = */ "MgO_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 3.58157,
        /* growth_ingredients = */ MgO_dust_ingred);

    const GrainGrowthIngredient FeS_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {1, SpLUT::Fe, 56.},
        {1, SpLUT::S, 32.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[8] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::FeS_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::FeS_dust,
        /* name = */ "FeS_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 680.0,
        /* bulk_density_cgs = */ 4.87265,
        /* growth_ingredients = */ FeS_dust_ingred);

    const GrainGrowthIngredient Al2O3_dust_ingred[] = {
        // {coef, species_idx, particle mass}
        {2, SpLUT::Al, 27.},
        {3, SpLUT::H2O, 18.},
        GRIMPL_INGREDIENT_LIST_SENTINEL};
    out.species_info[9] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::Al2O3_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::Al2O3_dust,
        /* name = */ "Al2O3_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 1500.0,
        /* bulk_density_cgs = */ 4.01610,
        /* growth_ingredients = */ Al2O3_dust_ingred);
  }

  if (dust_species_parameter > 2) {
    // We do not consider the growth of refractory organics, volatile
    // organics, and water ice because their sublimation temperatures
    // are low (100-600 K). They sublimate before the growth occurs.

    out.species_info[10] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::ref_org_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::ref_org_dust,
        /* name = */ "ref_org_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 575.0,
        /* bulk_density_cgs = */ 1.5,
        /* growth_ingredients = */ nullptr);
    out.species_info[11] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::vol_org_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::vol_org_dust,
        /* name = */ "vol_org_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 375.0,
        /* bulk_density_cgs = */ 1.0,
        /* growth_ingredients = */ nullptr);
    out.species_info[12] = mk_gsp_info_entry_helper_(
        /* species_idx = */ SpLUT::H2O_ice_dust,
        /* onlygrainsp_idx = */ OnlyGrainSpLUT::H2O_ice_dust,
        /* name = */ "H2O_ice_dust",
        /* h2dust_uses_carbonaceous_table = */ false,
        /* sublimation_temperature = */ 153.0,
        /* bulk_density_cgs = */ 0.92,
        /* growth_ingredients = */ nullptr);
  }

  return out;
}

#undef GRIMPL_INGREDIENT_LIST_SENTINEL

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
