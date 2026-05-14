//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_grain_size_increment_species_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_grain_size_increment_species_1d function from FORTRAN to C++

#ifndef CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP
#define CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP

#include "fortran_func_decls.h"          // gr_mask_int
#include "grackle.h"                     // gr_float
#include "support/index_helper.hpp"      // IndexRange

namespace grackle::impl {

/// Calculates properties that are derived from the grain size increment
/// for a single grain species.
///
/// @todo
/// The subroutine's name should be more descriptive
///
/// @param[in] igrgr Flag to solve grain growth reactions
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] SN0_N Number of modelled injection pathways
/// @param[in] in, jn, kn Dimensions of the computational grid
/// @param[in] idx_range Specifies the current index-range
/// @param[in] density_data Pointer to the density field data
/// @param[in] n_selected_inj_paths Number of selected injection pathways
/// @param[in] grain_species_density Pointer to the density field data for the
/// current grain species
/// @param[in] selected_inj_path_metal_densities Pointer to repacked selected
///     injection pathway metal densities for the current @p idx_range. This is
///     a 2d array where the contiguous axis holds @p in elements and the other
///     axis has space for @p SN0_N entries (in practice only the first
///     @p n_selected_inj_paths indices along this axis are used).
/// @param[in] SN_fsp Pointer to the array of values for the initial fraction of
/// the injected mass density of a given grain species
/// @param[in] SN_r0sp_data Pointer to the table of values for the initial size
/// distribution of a given grain species (1st, 2nd, and 3rd order moments)
/// @param[in] bulk_density The bulk density of the grain species (density of a
/// single grain in units of g/cm^3)
/// @param[in] sigma_per_gas_mass Pointer to the array geometric cross-section 
/// per unit gas mass of each grain species
/// @param[in,out] kappa_data Pointer to the table of values related to opacity
/// for the current grain species
/// @param[in] gr_N Bi-dimensional array of the number of tabulated values for
/// opacity-related quantities for each grain species
/// @param[in] gr_Size Total data size deriving from tabulated values for
/// opacity-related quantities
/// @param[in] opac_coef_table_data Tables of values for opacity calculations
///
/// @par History
/// modified: February, 2026 by Christopher Bignamini & Matthew Abruzzo; port to C++
void calc_grain_size_increment_species_1d(
    int igrgr, const gr_mask_type* itmask, int SN0_N, int in, int jn, int kn,
    IndexRange idx_range, const gr_float* density_data, int n_selected_inj_paths,
    const gr_float* grain_species_density, gr_float* selected_inj_path_metal_densities,
    const double* SN_fsp, double* SN_r0sp_data, double bulk_density, double* sigma_per_gas_mass,
    double* kappa_data, int* gr_N, int gr_Size, double* opac_coef_table_data);

}  // namespace grackle::impl

#endif /* CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D_HPP */
