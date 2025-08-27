//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the metal chemistry rates
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_METAL_CHEMISTRY_RATES_H
#define INITIALIZE_METAL_CHEMISTRY_RATES_H

#include "grackle.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int initialize_metal_chemistry_rates(chemistry_data* my_chemistry,
                                     chemistry_data_storage* my_rates,
                                     code_units* my_units);

int local_free_metal_chemistry_rates(chemistry_data* my_chemistry,
                                     chemistry_data_storage* my_rates);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* INITIALIZE_METAL_CHEMISTRY_RATES_H */
