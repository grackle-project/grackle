//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines some mask-related functionality
///
//===----------------------------------------------------------------------===//

#ifndef MASK_HPP
#define MASK_HPP

#include "fortran_func_decls.h"  // gr_mask_int
#include "grackle.h"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "utils-cpp.hpp"  // View

namespace GRIMPL_NAMESPACE_DECL {
namespace mask {

/// @brief adjust @p itmask based on the temperature floor for the @p idx_range
///
/// @param[inout] itmask The iteration-mask of the @p idx_range that is
///     adjusted by this function
/// @param[in] tgas 1D array of gas temperatures for the @p idx_range
/// @param[in] idx_range Specifies the current index-range
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_fields Specifies the field data.
inline void adjust_from_Tfloor(gr_mask_type* itmask, const double* tgas,
                               IndexRange idx_range,
                               const chemistry_data* my_chemistry,
                               const grackle_field_data* my_fields) {
  if (my_chemistry->use_temperature_floor == 1) {
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] <= my_chemistry->temperature_floor_scalar) {
          itmask[i] = MASK_FALSE;
        }
      }
    }
  } else if (my_chemistry->use_temperature_floor == 2) {
    grackle::impl::View<gr_float***> Tfloor(
        my_fields->temperature_floor, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] <= Tfloor(i, idx_range.j, idx_range.k)) {
          itmask[i] = MASK_FALSE;
        }
      }
    }
  }
}

/// @brief fill up the @p itmask_metal for the @p idx_range. This mask
///     specifies where metal-related calculations should be skipped
///
/// @param[inout] itmask_metal Iteration-mask to be filled by this function for
///     the @p idx_range.
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] metallicity 1D array of metallicities for the @p idx_range
/// @param[in] imetal Indicates whether grackle is configured to evolve metals
/// @param[in] idx_range Specifies the current index-range
/// @param[in] my_chemistry holds a number of configuration parameters.
inline void fill_itmask_metal(gr_mask_type* itmask_metal,
                              const gr_mask_type* itmask,
                              const double* metallicity, int imetal,
                              IndexRange idx_range,
                              const chemistry_data* my_chemistry) {
  if (imetal == 1) {
    double min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (metallicity[i] >= min_metallicity) {
        itmask_metal[i] = itmask[i];
      } else {
        itmask_metal[i] = MASK_FALSE;
      }
    }
  } else {
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      itmask_metal[i] = MASK_FALSE;
    }
  }
}

}  // namespace mask
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // MASK_HPP
