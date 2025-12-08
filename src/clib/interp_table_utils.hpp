//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare utility functions related to the interpolation tables
///
//===----------------------------------------------------------------------===//

#ifndef INTERP_TABLE_UTILS_HPP
#define INTERP_TABLE_UTILS_HPP

#include "grackle.h"  // gr_interp_grid, GRACKLE_CLOUDY_TABLE_MAX_DIMENSION
#include "grackle_macros.h"  // GRACKLE_FREE

namespace grackle::impl {

/// Free memory associated with a #gr_interp_grid_props instance
inline void free_interp_grid_props_(gr_interp_grid_props* props,
                                    bool use_delete) {
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    if (use_delete) {
      if (props->parameters[i] != nullptr) {
        delete[] props->parameters[i];
        props->parameters[i] = nullptr;
      }
    } else {
      GRACKLE_FREE(props->parameters[i]);
    }
  }
}

/// Free memory associated with a #gr_interp_grid
inline void free_interp_grid_(gr_interp_grid* grid) {
  free_interp_grid_props_(&(grid->props), /* use_delete = */ false);
  GRACKLE_FREE(grid->data);
}

}  // namespace grackle::impl

#endif /* INTERP_TABLE_UTILS_HPP */
