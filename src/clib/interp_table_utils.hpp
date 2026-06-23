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

#include <utility>  // std::swap, std::move, std::exchange

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

public:  // interface methods
  /// @brief default constructor (makes an empty instance)
  gr_interp_grid_props() : rank{0}, data_size(0) {
    for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
      dimension[i] = 0;
      parameters[i] = nullptr;
      parameter_spacing[i] = 0.0;
    }
  }

  // delete copy constructor and assignment (these lead to dangling pointers)
  gr_interp_grid_props(const gr_interp_grid_props&) = delete;
  gr_interp_grid_props& operator=(const gr_interp_grid_props&) = delete;

  /// @brief Move constructor
  gr_interp_grid_props(gr_interp_grid_props&& other) noexcept
      : gr_interp_grid_props() {
    *this = std::move(other);  // <- use the move assigment
  }

  /// @brief Move assignment
  gr_interp_grid_props& operator=(gr_interp_grid_props&& other) noexcept {
    // swapping contents is a fairly standard idiom
    std::swap(rank, other.rank);
    std::swap(dimension, other.dimension);
    std::swap(parameters, other.parameters);
    std::swap(parameter_spacing, other.parameter_spacing);
    std::swap(data_size, other.data_size);
    return *this;
  }
};

/// Encodes a full interpolation grid
struct gr_interp_grid {
  /// properties of the interpolation grid
  gr_interp_grid_props props;
  /// the actual data that gets interpolated
  double* data;

public:  // interface methods
  /// @brief default constructor (makes an empty instance)
  gr_interp_grid() : props(), data{nullptr} {}

  // delete copy constructor and assignment (these lead to dangling pointers)
  gr_interp_grid(const gr_interp_grid&) = delete;
  gr_interp_grid& operator=(const gr_interp_grid&) = delete;

  /// @brief Move constructor
  gr_interp_grid(gr_interp_grid&& other) noexcept
      : props(std::move(other.props)),
        data{std::exchange(other.data, nullptr)} {}

  /// @brief Move assignment
  gr_interp_grid& operator=(gr_interp_grid&& other) noexcept {
    props = std::move(other.props);
    std::swap(data, other.data);
    return *this;
  }
};

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
