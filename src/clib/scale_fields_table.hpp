// See LICENSE file for license and copyright information

/// @file scale_fields_table.hpp
/// @brief Defines the scale_fields_table function

// This file was initially generated automatically during conversion of the
// scale_fields_table_g function from FORTRAN to C++

#ifndef SCALE_FIELDS_TABLE_HPP
#define SCALE_FIELDS_TABLE_HPP

#include "grackle.h"
#include "utils-cpp.hpp"

namespace grackle::impl {

/// Scales density and metal_density (if available) by factor
inline void scale_fields_table(grackle_field_data* my_fields, double factor) {
  const int grid_start[3] = {my_fields->grid_start[0], my_fields->grid_start[1],
                             my_fields->grid_start[2]};
  const int grid_end[3] = {my_fields->grid_end[0], my_fields->grid_end[1],
                           my_fields->grid_end[2]};

  // Multiply density by factor (1/a^3 or a^3)

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  for (int k = grid_start[2]; k <= grid_end[2]; k++) {
    for (int j = grid_start[1]; j <= grid_end[1]; j++) {
      for (int i = grid_start[0]; i <= grid_end[0]; i++) {
        d(i, j, k) = d(i, j, k) * factor;
      }
    }
  }

  if (my_fields->metal_density != nullptr) {
    grackle::impl::View<gr_float***> metal(
        my_fields->metal_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    for (int k = grid_start[2]; k <= grid_end[2]; k++) {
      for (int j = grid_start[1]; j <= grid_end[1]; j++) {
        for (int i = grid_start[0]; i <= grid_end[0]; i++) {
          metal(i, j, k) = metal(i, j, k) * factor;
        }
      }
    }
  }
}

}  // namespace grackle::impl

#endif /* SCALE_FIELDS_TABLE_HPP */
