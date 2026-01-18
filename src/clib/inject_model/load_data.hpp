//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function for loading injection pathway data
///
//===----------------------------------------------------------------------===//

#ifndef INJECT_MODEL_LOAD_DATA_HPP
#define INJECT_MODEL_LOAD_DATA_HPP

#include "grackle.h"

namespace grackle::impl {

/// loads the model data for the various injection pathways and update
/// @p my_rates, accordingly
///
/// @returns GR_SUCCESS if successful
int load_inject_path_data(const chemistry_data* my_chemistry,
                          chemistry_data_storage* my_rates);

}  // namespace grackle::impl

#endif /* INJECT_MODEL_LOAD_DATA_HPP */
