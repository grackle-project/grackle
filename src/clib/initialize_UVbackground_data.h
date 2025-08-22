//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the UV background data
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_UVBACKGROUND_DATA_H
#define INITIALIZE_UVBACKGROUND_DATA_H

#include "grackle_chemistry_data.h"

/// Initializes an empty UVBtable struct with zeros and NULLs.
void initialize_empty_UVBtable_struct(UVBtable *table);

/// Initialize UV Background data
int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates);

#endif /* INITIALIZE_UVBACKGROUND_DATA_H */
