//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines the GrainMetalInjectPathways type
///
//===----------------------------------------------------------------------===//
#ifndef GRAIN_METAL_INJECT_PATHWAYS_HPP
#define GRAIN_METAL_INJECT_PATHWAYS_HPP

// after dealing with GH Issue #446, move this file into the dust subdirectory

#include "LUT.hpp"
#include "internal_types.hpp"
#include "visitor/common.hpp"

// this is a separate namespace because we are probably going to delete the
// contents (when we deal with GH Issue #446)
namespace grackle::impl::yields {

/// Tables of values for each metal nuclide known to Grackle. Each table
/// has an entry for each injection pathway modelled in the current
/// configuration
///
/// @note
/// When we deal with GH Issue #446, we'll almost certainly delete this. We
/// avoid using the visitor pattern here since this is very much an ephemeral
/// type introduced as a stopgap solution (and because types used to implement
/// the visitor pattern, namely ``vis::idx_range_len_multiple``) should probably
/// be renamed to have a more general-sounding name
struct MetalTables {
  double* C;
  double* O;
  double* Mg;
  double* Al;
  double* Si;
  double* S;
  double* Fe;
};

/// allocates the contents of a new MetalTable
///
/// @param nelem The number of elements in each buffer
inline MetalTables new_MetalTables(int nelem) {
  MetalTables out;

  out.C = new double[nelem];
  out.O = new double[nelem];
  out.Mg = new double[nelem];
  out.Al = new double[nelem];
  out.Si = new double[nelem];
  out.S = new double[nelem];
  out.Fe = new double[nelem];

  return out;
}

/// performs cleanup of the contents of MetalTable
///
/// This effectively invokes a destructor
inline void drop_MetalTables(MetalTables* ptr) {
  delete[] ptr->C;
  delete[] ptr->O;
  delete[] ptr->Mg;
  delete[] ptr->Al;
  delete[] ptr->Si;
  delete[] ptr->S;
  delete[] ptr->Fe;
}

}  // namespace grackle::impl::yields

namespace grackle::impl {

/// The basic premise is that this tracks tables of data pertaining to the
/// list of injection pathways for dust grains and metal species
///
/// @note
/// The organization of data is a little suboptimal. But, we can reorganize
/// later.
struct GrainMetalInjectPathways {
  /// the number of modelled injection pathways
  int n_pathways;

  /// holds 1D tables (tabulated with respect to each injection pathway) that
  /// tracks the fraction of the injected mass density corresponding to
  /// occurences of given metal nuclide in both the gas phase **AND** as part
  /// of dust grains
  ///
  /// @important
  /// This counts up occurences of nulcides (i.e. a nuclide is allowed to occur
  /// multiple times within a single molecule)
  ///
  /// @note
  /// this variables holds information that is somewhat actually redundant with
  /// other tables that we are tracking
  yields::MetalTables total_metal_nuclide_yields;

  /// holds 1D tables that tracks the fraction of the injected mass density
  /// that corresponds to the **INITIAL** occurences of given metal nuclide in
  /// the gas phase.
  ///
  /// This holds a separate table for each metal nuclide and each table has an
  /// entry for each injection pathway.
  ///
  /// @important
  /// This counts up occurences of nulcides (i.e. a nuclide is allowed to occur
  /// multiple times within a single molecule)
  yields::MetalTables gas_metal_nuclide_yields;

  /// holds 1D tables (tabulated with respect to each injection pathway) that
  /// tracks the fraction of the injected mass density that corresponds to the
  /// **INITIAL** occurences of a given grain species
  ///
  /// This holds a separate table for each grain species and every table has an
  /// entry for each injection pathway.
  GrainSpeciesCollection grain_yields;
};

/// allocates the contents of a new GrainMetalInjectPathways
///
/// @param n_pathways The number of modeled injection pathways
inline GrainMetalInjectPathways new_GrainMetalInjectPathways(int n_pathways) {
  GrainMetalInjectPathways out;
  out.n_pathways = n_pathways;
  out.total_metal_nuclide_yields = yields::new_MetalTables(n_pathways);
  out.gas_metal_nuclide_yields = yields::new_MetalTables(n_pathways);
  out.grain_yields = new_GrainSpeciesCollection(n_pathways);

  return out;
}

/// performs cleanup of the contents of GrainMetalInjectPathways
///
/// This effectively invokes a destructor
inline void drop_GrainMetalInjectPathways(GrainMetalInjectPathways* ptr) {
  yields::drop_MetalTables(&ptr->total_metal_nuclide_yields);
  yields::drop_MetalTables(&ptr->gas_metal_nuclide_yields);
  drop_GrainSpeciesCollection(&ptr->grain_yields);
}

}  // namespace grackle::impl

#endif  // GRAIN_METAL_INJECT_PATHWAYS_HPP
