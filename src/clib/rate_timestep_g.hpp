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

#ifndef RATE_TIMESTEP_H
#define RATE_TIMESTEP_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "internal_types.hpp"

namespace grackle::impl {

void rate_timestep_g(
  double* dedot, double* HIdot, gr_mask_type anydust, double* h2dust,
  double* rhoH, gr_mask_type* itmask, double* edot, double chunit, double dom,
  chemistry_data* my_chemistry, grackle_field_data* my_fields,
  photo_rate_storage my_uvb_rates, IndexRange idx_range,
  grackle::impl::CollisionalRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
);

} // namespace grackle::impl

#endif /* RATE_TIMESTEP_H */
