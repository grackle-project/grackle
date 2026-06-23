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

/// @brief Encodes the grid properties of an interpolation grid (other than the
///        values being interpolated.
///
/// Ideally, we will reuse this struct to help implement the cloudy_data struct
struct gr_interp_grid_props {
  /// Rank of dataset
  ///
  /// TODO: do we need this attribute? In most cases, we know the rank
  ///       of a table ahead of time
  long long rank;

  /// Dimension of dataset.
  long long dimension[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Dataset parameter values (in the common case where there there is
  /// constant spacing, we could probably track less data).
  double* parameters[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Value of the constant paramter spacing
  double parameter_spacing[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Length of 1D flattened data grid
  long long data_size;
};

/// Encodes a full interpolation grid
struct gr_interp_grid {
  /// properties of the interpolation grid
  gr_interp_grid_props props;
  /// the actual data that gets interpolated
  double* data;
};

namespace grackle::impl {

/// initialize an empty #gr_interp_grid_props
inline void init_empty_interp_grid_props_(gr_interp_grid_props* props) {
  props->rank = 0;
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    props->dimension[i] = 0;
    props->parameters[i] = nullptr;
    props->parameter_spacing[i] = 0.0;
  }
  props->data_size = 0;
}

/// Initialize an empty #gr_interp_grid
inline void initialize_empty_interp_grid_(gr_interp_grid* grid) {
  init_empty_interp_grid_props_(&(grid->props));
  grid->data = nullptr;
}

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
