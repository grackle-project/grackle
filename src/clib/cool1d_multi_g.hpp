//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the cool1d_multi_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool1d_multi_g function from FORTRAN to C++

#ifndef COOL1D_MULTI_G_HPP
#define COOL1D_MULTI_G_HPP

#include "grackle.h"                 // gr_float
#include "fortran_func_decls.h"      // gr_mask_int
#include "internal_units.hpp"        // InternalGrUnits
#include "internal_types.hpp"        // GrainSpeciesCollection
#include "lnT_prep.hpp"              // LnTLinInterpBuf
#include "support/index_helper.hpp"  // IndexRange

namespace grackle::impl {

/// Solve radiative cooling/heating equations
///
/// @param[out] edot 1D array to hold the computed the time derivative of the
///     internal energy in the @p idx_range. Contributions are accumulated in
///     this buffer. In other words, this function does **NOT** set elements to
///     to 0 before adding contributions.
/// @param[in] alpha_continuum 1D array holding continuum linear absorption
///     coefficients for each location in the @p idx_range. This buffer may
///     already contain non-zero contributions.
/// @param[in] tgas 1D array of gas temperatures for the @p idx_range
/// @param[in] mmw 1D array of mean molecular weights for the @p idx_range
/// @param[in]  metallicity 1D array of metallicities for the @p idx_range
/// @param[in] rhoH 1D array of Hydrogen mass densities for the @p idx_range
/// @param[in] nelec_times_mH 1D array holding the number density of electrons
///     (multiplied by the Hydrogen mass) for the @p idx_range
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] itmask_metal Specifies the general metal-focused iteration-mask
///     of the @p idx_range for this calculation. Essentially, it is used to
///     skip metal-related calculations in zones that are extremely metal-poor
///     (or even metal-free).
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] my_uvb_rates Holds precomputed photorates that depend on the UV
///     background. These rates do not include the effects of self-shielding.
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[in] idx_range Specifies the current index-range
/// @param[in] logTlininterp_buf Hold values for each location in @p idx_range
///     that are used to linearly interpolate tables with respect to the natural
///     log of @p tgas1d.
/// @param[in] coolingheating_buf Pre-allocated buffers that are used by this
///     function for scratch space (to hold quantities that directly pertain to
///     cooling/heating
///
/// @par History
/// written by: Yu Zhang, Peter Anninos and Tom Abel
/// modified1: January, 1996 by Greg Bryan; adapted to KRONOS
/// modified2: October, 1996 by GB; moved to AMR
/// modified3: February, 2003 by Robert Harkness; iteration mask
/// modified4: September, 2009 by BDS to include cloudy cooling
/// modified5: March, 2025 by Christopher Bignamini & Matthew Abruzzo; C++ port
void cool1d_multi_g(
    double* edot, double* alpha_continuum, const double* tgas,
    const double* mmw, const double* metallicity, const double* rhoH,
    const double* nelec_times_mH, const gr_mask_type* itmask,
    const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
    IndexRange idx_range, grackle::impl::LnTLinInterpBuf logTlininterp_buf,
    grackle::impl::CoolHeatScratchBuf coolingheating_buf);

};  // namespace grackle::impl

#endif /* COOL1D_MULTI_G_HPP */
