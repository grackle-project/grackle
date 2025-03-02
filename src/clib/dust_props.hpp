// See LICENSE file for license and copyright information

/// @file dust_props.hpp
/// @brief Declarations related to computing dust properties.
///
/// In slightly more detail, these is a series of dust-related calculations
/// that are largely independent of the gas chemistry. These calculations all
/// share a series of buffers related to dust grain properties. The buffers are
/// filled by calc_grain_size_increment_1d. They are used to compute quantities
/// that is important for the gas
/// -> by calc_all_tdust_gasgr_1d_g
/// -> by some logic lookup_cool_rates1d_g (this logic should be excised)
///
/// The premise is that we should try to keep these properties as separated as
/// possible from the gas-physics

#ifndef DUST_PROPS_HPP
#define DUST_PROPS_HPP

#include <cstdlib> // malloc, free
#include "internal_types.hpp"

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
  double * sigma_per_gas_mass_tot;

};

/// allocates the contents of a new InternalDustPropBuf
///
/// @param nelem The number of elements in each buffer
inline InternalDustPropBuf new_InternalDustPropBuf(int nelem) {
  InternalDustPropBuf out;
  out.grain_sigma_per_gas_mass = grackle::impl::new_GrainSpeciesCollection(
    nelem
  );
  out.sigma_per_gas_mass_tot = (double*)std::malloc(sizeof(double)*nelem);

  return out;
}

/// performs cleanup of the contents of InternalDustPropBuf
///
/// This effectively invokes the destructor
inline void drop_InternalDustPropBuf(InternalDustPropBuf* ptr) {
  grackle::impl::drop_GrainSpeciesCollection(&ptr->grain_sigma_per_gas_mass);
  std::free(ptr->sigma_per_gas_mass_tot);
}

} // namespace grackle::impl

#endif // DUST_PROPS_HPP
