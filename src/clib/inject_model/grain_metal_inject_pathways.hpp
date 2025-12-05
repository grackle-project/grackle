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
  /// occurrences of given metal nuclide in both the gas phase **AND** as part
  /// of dust grains
  ///
  /// @important
  /// This counts up occurrences of nulcides (i.e. a nuclide is allowed to
  /// occur multiple times within a single molecule)
  ///
  /// @note
  /// this variables holds information that is somewhat actually redundant with
  /// other tables that we are tracking
  yields::MetalTables total_metal_nuclide_yields;

  /// holds 1D tables that tracks the fraction of the injected mass density
  /// that corresponds to the **INITIAL** occurrences of given metal nuclide in
  /// the gas phase.
  ///
  /// This holds a separate table for each metal nuclide and each table has an
  /// entry for each injection pathway.
  ///
  /// @important
  /// This counts up occurrences of nulcides (i.e. a nuclide is allowed to
  /// occur multiple times within a single molecule)
  yields::MetalTables gas_metal_nuclide_yields;

  /// holds 1D tables (tabulated with respect to each injection pathway) that
  /// tracks the fraction of the injected mass density that corresponds to the
  /// **INITIAL** occurrences of a given grain species
  ///
  /// This holds a separate table for each grain species and every table has an
  /// entry for each injection pathway.
  GrainSpeciesCollection grain_yields;

  /// holds 1D tables (tabulated with respect to each injection pathway) that
  /// tracks the 3 coefficients, derived from the 1st, 2nd, and 3rd order
  /// moments of the initial size distribution for a given grain species.
  ///
  /// This holds a separate table for each grain species and every table has 3
  /// contiguous entries for each pathway. In other words, each table has the
  /// shape ``(n_pathways, 3)`` using numpy conventions for a C-ordered array
  /// (to pass the array to Fortran, we would say a table has the shape
  /// ``(3, n_pathways)``.
  ///
  /// For a given grain species, I'm pretty confident that the entry at:
  /// - ``[j, 0]`` specifies \f$\langle r^1 \rangle_j\f$. This is the grain
  ///   species's average **INITIAL** radius, when injected by pathway ``j``.
  ///   The value has units of cm.
  /// - ``[j, 1]`` specifies \f$\langle r^2 \rangle_j\f$. The product of this
  ///   value and π gives the grain species's average **INITIAL**
  ///   cross-section, when injected by pathway ``j``. The value has units of
  ///   centimeters squared.
  /// - ``[j, 2]`` specifies \f$\langle r^3 \rangle_j\f$. The product of this
  ///   value and (4π/3) is the grain species's average **INITIAL** volume,
  ///   when injected by pathway ``j``. The value has units of centimeters
  ///   cubed.
  ///
  /// @todo
  /// What are the units of each quantity? (the dimensionality is obvious)
  ///
  /// where \f$\langle r^2 \rangle_j=\int_0^\infty r^p \Phi_j(r)\, {\rm d}r\f$
  /// is an abbreviation for the ``p``th moment of the \f$\Phi_j(r)\f$, or the
  /// initial differential grain size distribution for pathway ``j`` (this
  /// differs for each grain species). For added context:
  /// - \f$\Phi_j(r)\f$ is normalized such that the ``p=0`` moment is 1
  /// - if the number density of the given grain species injected by pathway
  ///   ``j`` is \f$n_j\f$, then \f$n_j\Phi_j(r)\f$ gives the initial number
  ///   density of the grain species with radii between ``r`` and ``r + dr``
  ///
  /// @todo
  /// We should move the detailed descriptions of the size distribution
  /// functions to the narrative docs and refer the reader to the appropriate
  /// section of documentation
  GrainSpeciesCollection size_moments;
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
  out.size_moments = new_GrainSpeciesCollection(3 * n_pathways);

  return out;
}

/// performs cleanup of the contents of GrainMetalInjectPathways
///
/// This effectively invokes a destructor
inline void drop_GrainMetalInjectPathways(GrainMetalInjectPathways* ptr) {
  yields::drop_MetalTables(&ptr->total_metal_nuclide_yields);
  yields::drop_MetalTables(&ptr->gas_metal_nuclide_yields);
  drop_GrainSpeciesCollection(&ptr->grain_yields);
  drop_GrainSpeciesCollection(&ptr->size_moments);
}

}  // namespace grackle::impl

#endif  // GRAIN_METAL_INJECT_PATHWAYS_HPP
