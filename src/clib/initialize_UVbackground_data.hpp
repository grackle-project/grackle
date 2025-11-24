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

#include "grackle.h"

namespace grackle::impl {

/// Initializes an empty UVBtable struct with zeros and nullptrs.
void initialize_empty_UVBtable_struct(UVBtable* table);

/// Initialize UV Background data
int initialize_UVbackground_data(chemistry_data* my_chemistry,
                                 chemistry_data_storage* my_rates);

/// deallocate memory associated with UVBtable struct
void free_UVBtable(UVBtable* table);

}  // namespace grackle::impl

#endif /* INITIALIZE_UVBACKGROUND_DATA_H */
