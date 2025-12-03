// See LICENSE file for license and copyright information

/// @file calc_temp1d_cloudy_g-cpp.h
/// @brief Declares signature of calc_temp1d_cloudy_g

// This file was initially generated automatically during conversion of the
// calc_temp1d_cloudy_g function from FORTRAN to C++

#ifndef CALC_TEMP1D_CLOUDY_G_HPP
#define CALC_TEMP1D_CLOUDY_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"
#include "internal_units.h"

namespace grackle::impl {

/// Calculate temperature and mean molecular weight for tabulated cooling.
///
/// @param[in] rhoH 1D array to hold the computed Hydrogen mass density for the
/// @p idx_range
/// @param[in] tgas 1D array to hold the temperature values for the @p
/// idx_range
/// @param[in] mmw 1D array to hold the mean molecular weight values for the
/// @p idx_range
/// @param[in] dom Unit conversion to proper number density in code units
/// @param[in] zr Current redshift
/// @param[in] imetal Flag if metal field is active (0 = no, 1 = yes)
/// @param[in] itmask Iteration mask
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] cloudy_table Cloudy cooling table data
/// @param[in] my_fields Specifies the field data.
/// @param[in] internalu Specifies unit information
/// @param[in] idx_range Specifies the current index-range
///
/// @par History
/// written by: Britton Smith, 2015
/// modified1: November, 2025 by Christopher Bignamini & Matthew Abruzzo; C++
/// port
void calc_temp1d_cloudy_g(double* rhoH, double* tgas, double* mmw, double dom,
                          double zr, int imetal, gr_mask_type* itmask,
                          chemistry_data* my_chemistry,
                          cloudy_data cloudy_table,
                          grackle_field_data* my_fields,
                          InternalGrUnits internalu, IndexRange idx_range);

}  // namespace grackle::impl

#endif /* CALC_TEMP1D_CLOUDY_G_HPP */
