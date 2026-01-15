//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of the rate_timestep_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// rate_timestep_g function from FORTRAN to C++

#ifndef RATE_TIMESTEP_HPP
#define RATE_TIMESTEP_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "full_rxn_rate_buf.hpp"
#include "internal_types.hpp"

namespace grackle::impl {

/// Primarily computes the rate of change in the electron and HI mass densities
/// (so that they can be used to limit the max timestep)
///
/// @todo
/// At the time of writing, this also updates edot with contributions from
/// HII, HeII, HeIII recombination heating and from H2 formation heating
/// (in the future, we should relocate this logic)
/// @param[in] dedot Time derivatives of free electron mass density
/// @param[in] HIdot Time derivatives of HI mass density
/// @param[in] anydust Indicates whether we are modelling dust
/// @param[in] h2dust Rate of H2 formation on dust grains
/// @param[in] rhoH Indicates the mass density of all Hydrogen
/// @param[in] itmask Specifies the `idx_range`'s iteration-mask for this
///    calculation
/// @param[in,out] edot Specifies the time derivative of internal energy
///    density for each location in `idx_range`
/// @param[in] chunit Conversion factor
/// @param[in] dom Unit conversion to proper number density in code units
/// @param[in] my_chemistry Holds a number of configuration parameters
/// @param[in] my_fields Specifies the field data
/// @param[in] idx_range Specifies the current index-range
/// @param[in] kcr_buf Holds various pre-computed chemical reaction rates for
///    each location in `idx_range`.
/// @param[in] kshield_buf Holds various pre-computed radiative reaction rates
/// @param[in] chemheatrates_buf Holds various pre-computed chemistry-heating
/// rates at each index-range location
/// @param[in] rxn_rate_buf Holds pre-computed reaction rates for each location
///    in `idx_range`.
/// @par History
/// written by:
/// modified1: November, 2025 by Christopher Bignamini & Matthew Abruzzo; C++
/// port
void rate_timestep_g(double* dedot, double* HIdot, gr_mask_type anydust,
                     const double* h2dust, const double* rhoH,
                     const gr_mask_type* itmask, double* edot, double chunit,
                     double dom, chemistry_data* my_chemistry,
                     grackle_field_data* my_fields, IndexRange idx_range,
                     grackle::impl::CollisionalRxnRateCollection kcr_buf,
                     grackle::impl::PhotoRxnRateCollection kshield_buf,
                     grackle::impl::ChemHeatingRates chemheatrates_buf,
                     FullRxnRateBuf rxn_rate_buf);

}  // namespace grackle::impl

#endif /* RATE_TIMESTEP_H */
