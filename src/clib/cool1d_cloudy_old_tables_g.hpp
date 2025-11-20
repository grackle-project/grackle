// See LICENSE file for license and copyright information

/// @file cool1d_cloudy_old_tables_g-cpp.h
/// @brief Declares signature of cool1d_cloudy_old_tables_g

// This file was initially generated automatically during conversion of the
// cool1d_cloudy_old_tables_g function from FORTRAN to C++

#ifndef COOL1D_CLOUDY_OLD_TABLES_G_HPP
#define COOL1D_CLOUDY_OLD_TABLES_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"        // IndexRange

namespace grackle::impl {

/// Solve cloudy cooling by interpolating from the data. This version uses tables
/// formatted for the original Cloudy cooling functionality in Enzo.
///
/// Solve cloudy metal cooling for old cloudy metal cooling tables
/// 
/// @param[in] rhoH 1D array to hold the computed Hydrogen mass density for the
/// @p idx_range
/// @param[in] metallicity 1D array to hold the computed metallicity for the @p
/// idx_range
/// @param[in] logtem Natural log of temperature values
/// @param[out] edot 1D array to hold the computed the time derivative of the
///     internal energy in the @p idx_range
/// @param[in] comp2 Constant factor to convert cooling rates to code units
/// @param[in] dom Unit conversion to proper number density in code units
/// @param[in] zr Current redshift
/// @param[in] itmask Iteration mask
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] idx_range Specifies the current index-range
///
/// @par History
/// written by: Britton Smith, 2009
/// modified1: November, 2025 by Christopher Bignamini & Matthew Abruzzo; C++
/// port
void cool1d_cloudy_old_tables_g(
  double* rhoH, double* metallicity, double* logtem, double* edot,
  double comp2, double dom, double zr, gr_mask_type* itmask,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range
);

} // namespace grackle::impl

#endif /* COOL1D_CLOUDY_OLD_TABLES_G_HPP */