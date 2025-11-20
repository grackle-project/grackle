//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the cool1d_cloudy_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool1d_cloudy_g function from FORTRAN to C++

#ifndef COOL1D_CLOUDY_G_HPP
#define COOL1D_CLOUDY_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"        // IndexRange

namespace grackle::impl {

/// Solve cloudy metal cooling/heating
///
/// Solve cloudy cooling by interpolating from the data.
///
/// @param[in] rhoH 1D array to hold the computed Hydrogen mass density for the @p idx_range
/// @param[in] metallicity 1D array to hold the computed metallicity for the @p idx_range
/// @param[in] logtem Natural log of temperature values
/// @param[out] edot 1D array to hold the computed the time derivative of the
///     internal energy in the @p idx_range
/// @param[in] comp2 Constant factor to convert cooling rates to code units
/// @param[in] dom Unit conversion to proper number density in code units
/// @param[in] zr Current redshift
/// @param[in] icmbTfloor Flag to include temperature floor from cmb
/// @param[in] iClHeat Flag to include cloudy heating
/// @param[in] iZscale Flag to scale cooling by metallicity
/// @param[in] clGridRank Rank of cloudy cooling data grid
/// @param[in] clGridDim Array containing dimensions of cloudy data
/// @param[in] clPar1, clPar2, clPar3 Arrays containing cloudy
///     grid parameter values
/// @param[in] clDataSize Total size of flattened 1D cooling data array
/// @param[in] clCooling Cloudy cooling data
/// @param[in] clHeating Cloudy heating data
/// @param[in] itmask Iteration mask
/// @param[in] my_fields Grackle field data
/// @param[in] idx_range Index range specifying the portion of the grid to
///     operate on
///
/// @par History
/// written by: Britton Smith, 2009
/// modified1: November, 2025 by Christopher Bignamini & Matthew Abruzzo; C++ port
void cool1d_cloudy_g(
  double* rhoH, double* metallicity, double* logtem, double* edot,
  double comp2, double* dom, double* zr, int* icmbTfloor, int* iClHeat,
  int* iZscale, long long* clGridRank, long long* clGridDim, double* clPar1,
  double* clPar2, double* clPar3, long long* clDataSize, double* clCooling,
  double* clHeating, gr_mask_type* itmask, grackle_field_data* my_fields,
  IndexRange idx_range
);

} // namespace grackle::impl

#endif /* COOL1D_CLOUDY_G_HPP */
