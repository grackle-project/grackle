//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements functions associated with GrainSpeciesInfo (namely, the
/// new_GrainSpeciesInfo function)
///
//===----------------------------------------------------------------------===//

#include <cstring>  // memcpy

#include "LUT.hpp"
#include "grain_species_info.hpp"

// The following logic effectively does 2 (related things):
// 1. it serves as a human-readable registry of all known grain species and
//    their associated properties.
// 2. we must be able to use this information to construct an instance of the
//    GrainSpeciesInfo data structure (so that the properties of each species
//    can be queried programatically)

// Encoding the gas-phase ingredients (for grain growth) of each grain species
// poses a bit of an obstacle. Essentially, we need a way to specify a
// variable-length list of ingredients and we should be able to infer the
// length. Thus, similar to the way a c-string has a sentinel value signalling
// the end of the list (i.e. the '\0' character), we use a sentinel for the
// hardcoded ingredient lists. This is intended as a stopgap solution!
//
// Note: the version of the lists accessible through GrainSpeciesInfo NEVER
// have the sentinel!

/// Expands to an initializer list for initializing an instance of the
/// GrainGrowthIngredient struct that is used to signal the termination of
/// a variable length ingredient list (analogous to the role that '\0' plays
/// in a c-string).
#define GRIMPL_INGREDIENT_LIST_SENTINEL {-1, -1, -1.0}

namespace {  // stuff inside an anonymous namespace is local to this file

/// Constructs a GrainSpeciesInfoEntry instance
///
/// The primary work that this function must complete is determine the number
/// of ingredients that a grain species has (if any) and then construct a
/// heap-allocated copy of the ingredient list that is stored by the returned
/// GrainSpeciesInfoEntry instance. Importantly:
/// - the @p growth_ingredients argument is terminated by the sentinel
/// - the ingredient list in the returned instance is **NOT** terminated by
///   the sentinel
grackle::impl::GrainSpeciesInfoEntry mk_gsp_info_entry_helper_(
    int species_idx, int onlygrainsp_idx, const char* name,
    bool h2dust_uses_carbonaceous_table, double sublimation_temperature,
    double bulk_density_cgs,
    const grackle::impl::GrainGrowthIngredient* growth_ingredients) {
  using grackle::impl::GrainGrowthIngredient;

  // the main "work" is to determine the number of specified growth_ingredients
  // and copy them (if there are no ingredients, we use a nullptr)
  GrainGrowthIngredient* out_ingredient_ptr = nullptr;
  int n_ingredients = 0;

  if (growth_ingredients != nullptr) {
    n_ingredients = 0;
    // increment the counter until we reach the sentinel
    while (growth_ingredients[n_ingredients].coef != -1) {
      n_ingredients++;
    }
    // allocate and initialize the pointer
    out_ingredient_ptr =
        new grackle::impl::GrainGrowthIngredient[n_ingredients];
    std::memcpy(out_ingredient_ptr, growth_ingredients,
                (std::size_t)n_ingredients);
  }

  return grackle::impl::GrainSpeciesInfoEntry{species_idx,
                                              onlygrainsp_idx,
                                              name,
                                              h2dust_uses_carbonaceous_table,
                                              sublimation_temperature,
                                              bulk_density_cgs,
                                              n_ingredients,
                                              out_ingredient_ptr};
}

}  // anonymous namespace

grackle::impl::GrainSpeciesInfo grackle::impl::new_GrainSpeciesInfo(
    int dust_species_parameter) {
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
