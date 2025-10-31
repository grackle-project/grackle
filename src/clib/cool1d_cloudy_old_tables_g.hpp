//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of the cool1d_cloudy_old_tables_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool1d_cloudy_old_tables_g function from FORTRAN to C++

// Solve cloudy cooling by interpolating from the data. This version uses tables
// formatted for the original Cloudy cooling functionality in Enzo.

#ifndef COOL1D_CLOUDY_OLD_TABLES_G_HPP
#define COOL1D_CLOUDY_OLD_TABLES_G_HPP

#include "grackle.h"             // gr_float
#include "index_helper.h"
#include "fortran_func_decls.h"  // gr_mask_int

namespace grackle::impl {

// TODO: check my_chemistry->UVbackground argument position
void cool1d_cloudy_old_tables_g(
  double* rhoH, double* metallicity, double* logtem, double* edot,
  double comp2, double dom, double zr, gr_mask_type* itmask,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range
);

} // namespace grackle::impl
#endif /* COOL1D_CLOUDY_OLD_TABLES_G_HPP */
