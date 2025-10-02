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

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "internal_units.h"      // InternalGrUnits
#include "internal_types.hpp"    // GrainSpeciesCollection
#include "index_helper.h"        // IndexRange

namespace grackle::impl {

/// Solve radiative cooling/heating equations
///
/// This function does a lot. It probably makes sense to break off some of
/// the functionality.
///
/// @param[in] imetal Indicates whether metals are evolved
/// @param[in] iter The current iteration (the first iteration is `1`)
/// @param[out] edot 1D array to hold the computed the time derivative of the
///     internal energy in the @p idx_range
/// @param[out] tgas 1D array to hold the computed gas temperatures in the
///     @p idx_range
/// @param[out] mmw 1D array to hold the computed mean molecular weight
///     in the @p idx_range
/// @param[in] p2d 1D array to hold the computed thermal gas pressures for
///     the @p idx_range. This is computed using the user-specified nominal
///     adiabatic index value (i.e. no attempts are made to correct for presence
///     of H2)
/// @param[out] tdust 1D array to hold the computed dust temperatures at
///     each location in the @p index range. This **ONLY** used when using
///     variants of the classic 1-field dust-model or using the variant of the
///     multi-grain-species model where all grains are configured to share a
///     single temperature.
/// @param[in]  metallicity 1D array to hold the computed metallicity for the
///     @p idx_range
/// @param[out] dust2gas Holds the computed dust-to-gas ratio at each
///     location in the index range. In other words, this holds the dust mass
///     per unit gas mass (only used in certain configuration)
/// @param[out] rhoH 1D array to hold the computed Hydrogen mass density
///     for the @p idx_range
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[out] itmask_metal
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] my_uvb_rates Holds precomputed photorates that depend on the UV
///     background. These rates do not include the effects of self-shielding.
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[in] idx_range Specifies the current index-range
/// @param[in] grain_temperatures buffers to hold individual grain species
///     temperatures. This is only used in certain configurations (i.e. when we
///     aren't using the tdust argument)
/// @param[in] logTlininterp_buf Scratch space used to temporarily hold values
///     for each location in @p idx_range with valuea that are used to linearly
///     interpolate tables with respect to the natural log of @p tgas1d. (Any
///     values previously stored here will be overwritten)
/// @param[in] cool1dmulti_buf Pre-allocated buffers that are used by this
///     function for scratch space (to hold a variety of quantities)
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
/// modified5: March, 2025 by Christopher Bignamini & Matthew Abruzzo; ported to
/// C++

void cool1d_multi_g(int imetal, int iter, double* edot, double* tgas,
                    double* mmw, double* p2d, double* tdust,
                    double* metallicity, double* dust2gas, double* rhoH,
                    gr_mask_type* itmask, gr_mask_type* itmask_metal,
                    chemistry_data* my_chemistry,
                    chemistry_data_storage* my_rates,
                    grackle_field_data* my_fields,
                    photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
                    IndexRange idx_range,
                    grackle::impl::GrainSpeciesCollection grain_temperatures,
                    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
                    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
                    grackle::impl::CoolHeatScratchBuf coolingheating_buf);

};  // namespace grackle::impl

#endif /* COOL1D_MULTI_G_HPP */
