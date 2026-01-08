//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the calc_all_tdust_gasgr_1d_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_all_tdust_gasgr_1d_g function from FORTRAN to C++

#ifndef CALC_ALL_TDUST_GASGR_1D_G_HPP
#define CALC_ALL_TDUST_GASGR_1D_G_HPP

#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"

namespace grackle::impl {

///  Calculate all dust temperature(s) and the gas to grain heat
///  transfer rate(s) (the latter is commonly called gasgr)
///
///  An argument could be made for directly using gasgr to compute
///  contributions to edot within this routine, rather than returning
///  the gasgr value(s). But that should be reconsidered in the future.
///
///  We could significantly reduce the amount of buffer space allocated
///  within the routine.
///
/// @param[in]  trad CMB ratiation temperature
/// @param[in]  tgas 1D array to hold the computed gas temperatures
/// @param[out] tdust 1D array to hold the computed dust temperatures
/// @param[in]  metallicity 1D array to hold the computed metallicity
/// @param[out] dust2gas Holds the computed dust-to-gas ratio at each
///     location in the index range. In other words, this holds the dust mass
///     per unit gas mass (only used in certain configuration)
/// @param[in]  nh 1D array of hydrogen number densities
/// @param[in]  gasgr_tdust Array of gas-grain dust temperatures
/// @param[in]  itmask_metal Mask array indicating presence of metals
/// @param[in]  coolunit Cooling unit conversion factor
/// @param[out] gasgr Array of gas-grain heat transfer rates
/// @param[in]  myisrf Array of interstellar radiation field strengths
/// @param[in]  kptot Array of total grain opacities
/// @param[in]  my_chemistry Pointer to chemistry data structure
/// @param[in]  my_rates Pointer to chemistry data storage structure
/// @param[in]  my_fields Pointer to field data structure
/// @param[in]  idx_range Specifies the current index-range
/// @param[in]  grain_temperatures buffers to hold individual grain species
///     temperatures. This is only used in certain configurations (i.e. when we
///     aren't using the tdust argument)
/// @param[in]  gas_grainsp_heatrate Gas-grain species heat rate collection
/// @param[in]  logTlininterp_buf Scratch space used to temporarily hold values
///     for each location in @p idx_range .
/// @param[in]  internal_dust_prop_buf Internal dust properties buffer
/// @param[in]  grain_kappa Grain species opacity collection
///
/// @par History
/// modified: January, 2026 by Christopher Bignamini & Matthew Abruzzo; C++ port
void calc_all_tdust_gasgr_1d_g(
    double trad, double* tgas, double* tdust, double* metallicity,
    double* dust2gas, double* nh, double* gasgr_tdust,
    gr_mask_type* itmask_metal, double coolunit, double* gasgr, double* myisrf,
    double* kptot, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    IndexRange idx_range,
    grackle::impl::GrainSpeciesCollection grain_temperatures,
    grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
    grackle::impl::InternalDustPropBuf internal_dust_prop_buf,
    grackle::impl::GrainSpeciesCollection grain_kappa);

}  // namespace grackle::impl

#endif /* CALC_ALL_TDUST_GASGR_1D_G_HPP */
