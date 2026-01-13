//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the cool_multi_time_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool_multi_time_g function from FORTRAN to C++

#ifndef COOL_MULTI_TIME_G_H
#define COOL_MULTI_TIME_G_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "internal_units.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
// the following function can be called from C or C++

/// Solve the energy cooling equations
///
/// @par history
/// * written by: Yu Zhang, Peter Anninos and Tom Abel
/// * modified1: January, 1996 by Greg Bryan; adapted to KRONOS
/// * modified2: October, 1996 by GB; moved to AMR
/// * modified3: February, 2003 by Robert Harkness; iteration mask
/// * modified4: December, 2024 by Matthew Abruzzo; ported to C++
///
/// @todo
/// Once the file where this routine called is adjusted to be compiled with a
/// C++ compiler, modify this function (prototype & implementation) such that:
/// - it's not enclosed by a `extern "C"` block
/// - it's defined within a `grackle::impl` namespace
void cool_multi_time_g(gr_float* cooltime_data_, int imetal,
                       InternalGrUnits internalu, chemistry_data* my_chemistry,
                       chemistry_data_storage* my_rates,
                       grackle_field_data* my_fields,
                       photo_rate_storage my_uvb_rates);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* COOL_MULTI_TIME_G_H */
