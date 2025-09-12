//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the dust yields
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_DUST_YIELDS_H
#define INITIALIZE_DUST_YIELDS_H

#include "grackle.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int initialize_dust_yields(chemistry_data* my_chemistry,
                           chemistry_data_storage* my_rates,
                           code_units* my_units);

int local_free_dust_yields(chemistry_data* my_chemistry,
                           chemistry_data_storage* my_rates);

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* INITIALIZE_DUST_YIELDS_H */
