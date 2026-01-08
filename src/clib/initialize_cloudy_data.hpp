//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function to initialize the cloudy data
///
//===----------------------------------------------------------------------===//

#ifndef INITIALIZE_CLOUDY_DATA_H
#define INITIALIZE_CLOUDY_DATA_H

#include "grackle.h"

namespace grackle::impl {

// initialize cloudy cooling data
int initialize_cloudy_data(chemistry_data* my_chemistry,
                           chemistry_data_storage* my_rates,
                           cloudy_data* my_cloudy, const char* group_name,
                           code_units* my_units, int read_data);

int free_cloudy_data(cloudy_data* my_cloudy, chemistry_data* my_chemistry,
                     int primordial);

}  // namespace grackle::impl

#endif /* INITIALIZE_CLOUDY_DATA_H */
