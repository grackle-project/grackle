// See LICENSE file for license and copyright information

/// @file calc_all_tdust_gasgr_1d_g-cpp.h
/// @brief Declares signature of calc_all_tdust_gasgr_1d_g

// This file was initially generated automatically during conversion of the
// calc_all_tdust_gasgr_1d_g function from FORTRAN to C++

#ifndef CALC_ALL_TDUST_GASGR_1D_G_HPP
#define CALC_ALL_TDUST_GASGR_1D_G_HPP

#include "grackle.h"             // gr_float
#include "dust_props.hpp"
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"
#include "internal_types.hpp"

namespace grackle::impl {

void calc_all_tdust_gasgr_1d_g(
  double trad, double* tgas, double* tdust, double* metallicity,
  double* dust2gas, double* nh, double* gasgr_tdust, gr_mask_type* itmask_metal,
  double coolunit, double* gasgr, double* myisrf, 
  double* kptot, 
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf,
  grackle::impl::GrainSpeciesCollection grain_kappa);

} // namespace grackle::impl

#endif /* CALC_ALL_TDUST_GASGR_1D_G_HPP */
