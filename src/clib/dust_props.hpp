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

#include "internal_types.hpp"

namespace grackle::impl {

struct InternalDustPropBuf{};

/// allocates the contents of a new InternalDustPropBuf
///
/// @param nelem The number of elements in each buffer
inline InternalDustPropBuf new_InternalDustPropBuf(int nelem) {
  InternalDustPropBuf out;
  return out;
}

/// performs cleanup of the contents of InternalDustPropBuf
///
/// This effectively invokes the destructor
inline void drop_InternalDustPropBuf(InternalDustPropBuf* ptr) {
  // no-op
}

} // namespace grackle::impl

#endif // DUST_PROPS_HPP
