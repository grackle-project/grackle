//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// <add a short description>
///
//===----------------------------------------------------------------------===//
#ifndef CHEMISTRY_SPECIES_HPP
#define CHEMISTRY_SPECIES_HPP

#include "grackle_chemistry_data.h"
#include "support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {
/// Encodes the choice for the species that are only used to track the mass
/// in a subset of metal nuclides
///
/// @note
/// At the time of writing, the choice is tied to the dust model
enum class MetalNuclideMassSpChoice { NONE, ONLY_Mg, ALL };

inline MetalNuclideMassSpChoice from_dust_species_param(int dust_species) {
  switch (dust_species) {
    case 1:
      return MetalNuclideMassSpChoice::ONLY_Mg;
    case 2:
      [[fallthrough]];
    case 3:
      return MetalNuclideMassSpChoice::ALL;
    default:
      return MetalNuclideMassSpChoice::NONE;
  }
}

/// categorizes the chemical species
///
/// @note
/// In the near-future, we plan to actually construct queryable tables of
/// information about chemical species. When we do so, this enum will just
/// be used to help initialize the table
enum class ChemicalSpKind {
  /// the standard case
  STANDARD,
  /// electrons are a special case
  ELECTRON,
  /// these fields (Mg, Al, S, Fe) @b ONLY existing to accounting for the
  /// amount of the relevant nuclide in the gas-phase. Unlike other chemical
  /// species:
  /// - these @b ONLY participate in dust chemistry reactions
  /// - they don't track species with a particular charge
  METAL_NUCLIDE_MASS,
};

/// describes the signature of the callback function that is used to digest
/// information about each chemical species
///
/// As is standard for C-like callback functions, a callback function accepts
/// a function-specific `void*` context argument that is used to pass extra
/// information to the callback and track information between calls to the
/// callback
typedef void (*visit_chemsp_callback)(const char* name, ChemicalSpKind kind,
                                      void* ctx);

/// calls visit_chemsp_callback for each chemical species relevant to the
/// chemical configuration specified by `my_chemistry`
///
/// @param[in] primordial_chemistry The primordial_chemistry configuration
///     parameter
/// @param[in] metal_chemistry The metal_chemistry configuration parameter
/// @param[in] choice The choice pertaining the metal nuclide mass-tracking
///     species
/// @param[in] cb The callback function
/// @param[in,out] ctx A pointer to user-defined data that is passed into the
///     callback function
///
/// @note
/// If we are willing to more fully embrace C++, it would be a lot more
/// idiomatic to eliminate the `context` argument and make `fn` a template
/// argument
inline void visit_canonical_chem_species(int primordial_chemistry,
                                         int metal_chemistry,
                                         MetalNuclideMassSpChoice choice,
                                         visit_chemsp_callback* cb, void* ctx) {
  if (primordial_chemistry >= 1) {
    (*cb)("e", ChemicalSpKind::ELECTRON, ctx);  // free electrons
    (*cb)("HI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HeI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HeII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HeIII", ChemicalSpKind::STANDARD, ctx);
  }

  if (primordial_chemistry >= 2) {
    (*cb)("HM", ChemicalSpKind::STANDARD, ctx);
    (*cb)("H2I", ChemicalSpKind::STANDARD, ctx);
    (*cb)("H2II", ChemicalSpKind::STANDARD, ctx);
  }

  if (primordial_chemistry >= 3) {
    (*cb)("DI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("DII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HDI", ChemicalSpKind::STANDARD, ctx);
  }

  if (primordial_chemistry >= 4) {
    (*cb)("DM", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HDII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("HeHII", ChemicalSpKind::STANDARD, ctx);
  }

  if (metal_chemistry == 1) {
    (*cb)("CI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("CII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("CO", ChemicalSpKind::STANDARD, ctx);
    (*cb)("CO2", ChemicalSpKind::STANDARD, ctx);
    (*cb)("OI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("OH", ChemicalSpKind::STANDARD, ctx);
    (*cb)("H2O", ChemicalSpKind::STANDARD, ctx);
    (*cb)("O2", ChemicalSpKind::STANDARD, ctx);
    (*cb)("SiI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("SiOI", ChemicalSpKind::STANDARD, ctx);
    (*cb)("SiO2I", ChemicalSpKind::STANDARD, ctx);
    (*cb)("CH", ChemicalSpKind::STANDARD, ctx);
    (*cb)("CH2", ChemicalSpKind::STANDARD, ctx);
    (*cb)("COII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("OII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("OHII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("H2OII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("H3OII", ChemicalSpKind::STANDARD, ctx);
    (*cb)("O2II", ChemicalSpKind::STANDARD, ctx);
  }

  if (choice == MetalNuclideMassSpChoice::ONLY_Mg ||
      choice == MetalNuclideMassSpChoice::ALL) {
    (*cb)("Mg", ChemicalSpKind::METAL_NUCLIDE_MASS, ctx);
  }

  if (choice == MetalNuclideMassSpChoice::ALL) {
    (*cb)("Al", ChemicalSpKind::METAL_NUCLIDE_MASS, ctx);
    (*cb)("S", ChemicalSpKind::METAL_NUCLIDE_MASS, ctx);
    (*cb)("Fe", ChemicalSpKind::METAL_NUCLIDE_MASS, ctx);
  }
}

inline void visit_canonical_chem_species(chemistry_data* my_chemistry,
                                         visit_chemsp_callback* cb, void* ctx) {
  MetalNuclideMassSpChoice choice =
      from_dust_species_param(my_chemistry->dust_species);
  visit_canonical_chem_species(my_chemistry->primordial_chemistry,
                               my_chemistry->metal_chemistry, choice, cb, ctx);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEMISTRY_SPECIES_HPP