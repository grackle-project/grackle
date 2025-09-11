//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of make_consistent_g
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#ifndef MAKE_CONSISTENT_HPP
#define MAKE_CONSISTENT_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

namespace grackle::impl {

  void make_consistent(
  const int* imetal, const double* dom, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields
);

} // namespace grackle::impl

#endif /* MAKE_CONSISTENT_HPP */
