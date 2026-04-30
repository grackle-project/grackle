//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the cool_multi_time function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool_multi_time_g function from FORTRAN to C++

#ifndef COOL_MULTI_TIME_HPP
#define COOL_MULTI_TIME_HPP

#include "grackle.h"  // gr_float
#include "internal_units.h"
#include "support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Solve the energy cooling equations
///
/// @par history
/// * written by: Yu Zhang, Peter Anninos and Tom Abel
/// * modified1: January, 1996 by Greg Bryan; adapted to KRONOS
/// * modified2: October, 1996 by GB; moved to AMR
/// * modified3: February, 2003 by Robert Harkness; iteration mask
/// * modified4: December, 2024 by Matthew Abruzzo; ported to C++
void cool_multi_time(gr_float* cooltime_data_, int imetal,
                     InternalGrUnits internalu, chemistry_data* my_chemistry,
                     chemistry_data_storage* my_rates,
                     grackle_field_data* my_fields,
                     photo_rate_storage my_uvb_rates);
}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* COOL_MULTI_TIME_HPP */
