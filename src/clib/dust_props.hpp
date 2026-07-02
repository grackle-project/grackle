//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declarations related to computing dust properties.
///
/// In slightly more detail, these is a series of dust-related calculations
/// that are largely independent of the gas chemistry. These calculations all
/// share a series of buffers related to dust grain properties. The buffers are
/// filled by calc_grain_size_increment_1d. They are used to compute quantities
/// that is important for the gas
/// -> by calc_all_tdust_gasgr_1d_g
/// -> by some logic lookup_cool_rates1d (this logic should be excised)
///
/// The premise is that we should try to keep these properties as separated as
/// possible from the gas-physics
///
//===----------------------------------------------------------------------===//

#ifndef DUST_PROPS_HPP
#define DUST_PROPS_HPP

// In the near future (after GH PR #405 is merged), we plan to:
// 1. rename InternalDustPropBuf to InternalGrainPropBuf so that type's role in
//    the multi-grain-species dust model is more obvious
// 2. move this file into the dust subdirectory

#include "internal_types.hpp"
#include "LUT.hpp"

namespace grackle::impl {

struct InternalDustPropBuf {
  // a decent argument could be made that this should track grain temperature

  /// geometric cross-section per unit gas mass of each grain species
  /// - this description comes from the location where the values are computed
  /// - as I understand it:
  ///   - a geometric cross-section describes the grain's physical size
  ///   - the absorption cross-section depends on the geometric cross-section,
  ///     the dust temperature and the wavelength
  ///   - Rybicki & Lightman refer to the absorption coefficient per unit gas
  ///     mass as the "mass absorption coefficient" or "opacity coefficient"
  ///     and represent it with the variable kappa
  grackle::impl::GrainSpeciesCollection grain_sigma_per_gas_mass;

  /// buffer closely related to grain_sigma_per_gas_mass
  double* sigma_per_gas_mass_tot;

  /// holds 1D opacity tables at each zone for each grain species zone
  /// - the values are tabulated with respect to dust temperature. After
  ///   constructing the tables are used to compute opacity.
  ///   - here "opacity," is the optical cross-section per unit gas mass, for a
  ///     set of wavelengths. Based on the usage of the variable, it seems to
  ///     be for the continuum
  ///   - the quantity is called the "mass absorption coefficient" or
  ///     "opacity coefficient" by Rybicki & Lightman
  /// - we have chosen a name for this struct-member that reflects its role
  ///   of holding dynamically computed interpolation tables
  /// - if you are looking at Fortran source-code, you may be a little
  ///   surprised that this isn't called `grain_alpha`, which would have been
  ///   the "obvious" choice given that this replaces a list of variables that
  ///   all shared the "al" prefix.
  ///   - I suspect the original developer picked this prefix purely for the
  ///      convenience of having shorter variable-names and they did not
  ///      consider just how confusing this would be.
  ///   - The old convention was EXTREMELY confusing for multiple reasons:
  ///     1. The lack of documentation made it totally unclear that the
  ///        opacities, identified in code as kappa-variables, were simply the
  ///        interpolated values of these tables (you had to dig through 2
  ///        layers of source code to figure this out).
  ///     2. Even if there were documentation that alpha-variables held the
  ///        tabulated versions of kappa-variables, this would still be
  ///        confusing because kappa and alpha commonly represent different
  ///        quantities. For example, Rybicki & Lightman use alpha to debote
  ///        the linear absorption coefficient, which is the product of gas
  ///        mass density and kappa.
  ///     3. To make matters worse, the original author also wrote logic
  ///        following the conventions from Rybicki and Lightman: the opacity
  ///        coefficients (tracked by the kappa-variables) are used to compute
  ///        the linear absorption coefficient (stored in a variable called
  ///        alpha).
  /// - the fact that this holds a 1d table of values for each zone for each
  ///   species seems extremely inefficient!
  /// - If we really do need to construct the full table, then (if possible)
  ///   we should restructure the algorithm so that:
  ///   - we can reuse the buffer for different species
  ///   - we can also reuse the buffer for different zones (this would probably
  ///     involve reorganizing the tables of data for each metal-source in
  ///     order to support vectorization/memory coalescence)
  grackle::impl::GrainSpeciesCollection grain_dyntab_kappa;

  /// buffer closely related to grain_dyntab_kappa
  double* dyntab_kappa_tot;
};

/// allocates the contents of a new InternalDustPropBuf
///
/// @param nelem The number of elements in each buffer
/// @param opacity_table_size The number of elements in a dynamically computed
///    opacity table
inline InternalDustPropBuf new_InternalDustPropBuf(int nelem,
                                                   int opacity_table_size) {
  InternalDustPropBuf out;
  out.grain_sigma_per_gas_mass =
      grackle::impl::new_GrainSpeciesCollection(nelem);
  out.sigma_per_gas_mass_tot = new double[nelem];

  // this is a little bit of a hack!
  int table_len = (opacity_table_size < 1) ? 1 : opacity_table_size;
  out.grain_dyntab_kappa =
      grackle::impl::new_GrainSpeciesCollection(nelem * table_len);
  out.dyntab_kappa_tot = new double[nelem * table_len];

  return out;
}

/// performs cleanup of the contents of InternalDustPropBuf
///
/// This effectively invokes the destructor
inline void drop_InternalDustPropBuf(InternalDustPropBuf* ptr) {
  grackle::impl::drop_GrainSpeciesCollection(&ptr->grain_sigma_per_gas_mass);
  delete[] ptr->sigma_per_gas_mass_tot;

  grackle::impl::drop_GrainSpeciesCollection(&ptr->grain_dyntab_kappa);
  delete[] ptr->dyntab_kappa_tot;
}

}  // namespace grackle::impl

#endif  // DUST_PROPS_HPP
