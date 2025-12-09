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

#include "grackle.h"  // gr_interp_grid_props
#include "../LUT.hpp"
#include "../internal_types.hpp"
#include "../interp_table_utils.hpp"
#include "../status_reporting.h"
#include "../visitor/common.hpp"

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

/// fills elements with zeros
inline void MetalTables_zero_out(MetalTables* ptr, int nelem) {
  for (int i = 0; i < nelem; i++) {
    ptr->C [i] = 0.0;
    ptr->O [i] = 0.0;
    ptr->Mg[i] = 0.0;
    ptr->Al[i] = 0.0;
    ptr->Si[i] = 0.0;
    ptr->S [i] = 0.0;
    ptr->Fe[i] = 0.0;
  }
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

  /// Tracks 1D grid constant-spaced log10(Tdust) values used for interpolation
  /// in opacity calculations.
  ///
  /// @note
  /// Using a full #gr_interp_grid_props seems like it's a little overkill for
  /// holding values for 1D interpolation, but it may be worthwhile that this
  /// datatype is **ONLY** used for interpolation. Maybe in the future, we
  /// should make a specialized version for 1D data?
  gr_interp_grid_props log10Tdust_interp_props;

  /// the number of polynomial coefficients that will be computed from the
  /// opacity_coef_table. They are used in a polynomial of degree
  /// `(n_opac_poly_coef-1)`
  int n_opac_poly_coef;

  /// tables of values for opacity calculations.
  ///
  /// This holds a 3D array for each grain species. Each array has the shape
  /// `(n_pathways, nTd, n_opac_poly_coef)`, where
  /// - `nTd` is `log10Tdust_interp_props.dimension[0]`
  /// - the axis follows numpy conventions (i.e. `n_opac_poly_coef` is the
  ///   contiguous axis).
  ///
  /// What the data means
  /// ===================
  ///
  /// To be more concrete, consider a grain species that maps to data at an
  /// index `grsp_i` and an injection pathway that maps to data at an index
  /// `injp_i`. Then, let's consider `p`, which is given by
  /// `p = opacity_coef_table.data[grsp_i] + (injp_i * nTd * n_opac_poly_coef)`.
  /// In this scenario, `p` points to a 2D array that holds `n_opac_poly_coef`
  /// coefficients for each dust temperatures in `log10Tdust_interp_props`.
  ///
  /// In practice, we construct a vector, `vec`, of values for all `nTd` dust
  /// temperatures. The value computed for dust temperature `j` is:
  ///
  /// `vec[j] = 4πζ/3 * ∑ᵢ₌₀ ( p[j * n_opac_poly_coef + i] * δrⁱ(t) )`
  ///
  /// where:
  /// - the summation runs from `i=0` through `i = n_opac_poly_coef - 1`.
  /// - `ζ` is the "bulk density" of the considered grain species (in g/cm³).
  ///   In other words, it's the mass density of a single grain.
  /// - `δr(t)` refers to the derived "size increment", in units of cm,
  ///   computed at the current time `t`. Recall that Grackle's multi-grain
  ///   species dust model involves the initial size distribution functions for
  ///   grain species when they are initially injected. This is allowed to vary
  ///   across injection pathway and grain species (in fact, moments of these
  ///   functions are stored by `GrainMetalInjectPathways::size_moments`). If
  ///   grains undergo growth, the model assumes that the distribution function
  ///   is simply translated by the size increment, `δr` (the shape doesn't
  ///   change). Thus, a `δr` of `0` means that there's been no growth.
  ///
  /// `vec[j]` holds a quantity related to computing opacity. I, MWA, **THINK**
  /// that `vec[j]`:
  /// - it directly stores opacity, in units of `cm^2/g`.
  /// - is "per unit grain mass." If my understanding of this quantity and
  ///   related logic is indeed correct, this is noteworthy because subsequent
  ///   calculations appear use values deriving from `vec[j]` to compute
  ///   opacity measured "per unit gas mass."
  ///
  /// @todo
  /// Determine whether my statements in the preceeding paragraphs are indeed
  /// correct. More generally, what kind of weighted opacity is this? (e.g. a
  /// a Rosseland mean opacity or a Planck mean opacity?)? Does it describe a
  /// restricted wavelength range?
  ///
  /// The resulting vector of values is subsequently interpolated with respect
  /// to log10Tdust_interp_props
  ///
  /// @todo
  /// Ideally, we would move some of this description to the narrative docs and
  /// then refer the reader to the appropriate section of the documentation.
  /// (This could certainly be done for the description of the size increment)
  GrainSpeciesCollection opacity_coef_table;
};

/// allocates the contents of a new GrainMetalInjectPathways
///
/// @param[in] n_pathways Number of modelled injection pathways
/// @param[in] n_log10Tdust_vals Number of log10(Tdust) values that are relevant
///     for the opacity coefficient table.
/// @param[in] n_opac_poly_coef Number of opacity coefficients that are
///     computed from the opacity coefficient table.
inline GrainMetalInjectPathways new_GrainMetalInjectPathways(
    int n_pathways, int n_log10Tdust_vals, int n_opac_poly_coef) {

  bool err = false;

  if (n_pathways <= 0) {
    GrPrintErrMsg("n_pathways must be positive\n");
    err = true;
  } else if (n_log10Tdust_vals <= 1) {
    GrPrintErrMsg("n_log10Tdust_vals must exceed 1\n");
    err = true;
  } else if (n_opac_poly_coef != 4) {
    GrPrintErrMsg(
        "the logic that uses the opacity table is hardcoded to assume that "
        "the number of opacity polynomial coefficients is exactly 4\n");
    err = true;
  }

  GrainMetalInjectPathways out;

  if (err) {
    out.n_pathways = -1;
    return out;
  }

  out.n_pathways = n_pathways;

  out.total_metal_nuclide_yields = yields::new_MetalTables(n_pathways);
  out.gas_metal_nuclide_yields = yields::new_MetalTables(n_pathways);

  out.grain_yields = new_GrainSpeciesCollection(n_pathways);
  out.size_moments = new_GrainSpeciesCollection(3 * n_pathways);

  init_empty_interp_grid_props_(&out.log10Tdust_interp_props);
  out.log10Tdust_interp_props.rank = 1LL;
  long long n_log10Tdust_vals_LL = static_cast<long long>(n_log10Tdust_vals);
  out.log10Tdust_interp_props.dimension[0] = n_log10Tdust_vals_LL;
  out.log10Tdust_interp_props.data_size = n_log10Tdust_vals_LL;
  out.log10Tdust_interp_props.parameters[0] = new double[n_log10Tdust_vals];

  out.n_opac_poly_coef = n_opac_poly_coef;

  out.opacity_coef_table = new_GrainSpeciesCollection(
      n_pathways * n_log10Tdust_vals * n_opac_poly_coef);

  return out;
}

/// Checks whether the GrainMetalInjectPathways that is returned by
/// new_GrainMetalInjectPathways is valid
inline bool GrainMetalInjectPathways_is_valid(const GrainMetalInjectPathways* ptr) {
  return ptr->n_pathways != -1;
}

/// Queries the specified GrainMetalInjectPathways instance for the number of
/// log10(Tdust) values that are relevant for the opacity coefficient table.
///
/// This function is intended to be called whenever we construct a new
/// InternalDustPropBuf instance.
///
/// @param[in] ptr The instance to query. It is ok for this to be a `nullptr`,
///     however the result is undefined if GrainMetalInjectPathways_is_valid
///     would return false.
///
/// @note
/// Frankly, this is intended to be short-term solution until we have a chance
/// to reduce the number of places where InternalDustPropBuf is constructed
inline int GrainMetalInjectPathways_get_n_log10Tdust_vals
  (const GrainMetalInjectPathways* ptr) {
  return (ptr == nullptr) ? 0 : ptr->log10Tdust_interp_props.dimension[0];
}


/// performs cleanup of the contents of GrainMetalInjectPathways
///
/// This effectively invokes a destructor
inline void drop_GrainMetalInjectPathways(GrainMetalInjectPathways* ptr) {
  if (GrainMetalInjectPathways_is_valid(ptr)){
    yields::drop_MetalTables(&ptr->total_metal_nuclide_yields);
    yields::drop_MetalTables(&ptr->gas_metal_nuclide_yields);
    drop_GrainSpeciesCollection(&ptr->grain_yields);
    drop_GrainSpeciesCollection(&ptr->size_moments);
    free_interp_grid_props_(&ptr->log10Tdust_interp_props,
                            /* use_delete = */ true);
    drop_GrainSpeciesCollection(&ptr->opacity_coef_table);
  }
}

}  // namespace grackle::impl

#endif  // GRAIN_METAL_INJECT_PATHWAYS_HPP
