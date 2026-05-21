//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares some helper functions
///
//===----------------------------------------------------------------------===//

#ifndef TABULATED_COMMON_HPP
#define TABULATED_COMMON_HPP

#include <array>

#include "grackle.h"
#include "../support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace tabulated_detail {

inline constexpr std::size_t MAX_RANK = GRACKLE_CLOUDY_TABLE_MAX_DIMENSION;

/// Compute the step-size along each dimension.
///
/// @todo This should probably be tracked within @ref cloudy_data
inline std::array<double, MAX_RANK> param_deltas(const cloudy_data& table) {
  std::array<double, MAX_RANK> out{};  // <- all values initialized to 0
  for (long long i = 0; i < table.grid_rank; i++) {
    const double* vals = table.grid_parameters[i];
    long long n_vals = table.grid_dimension[i];
    out[i] = (vals[n_vals - 1] - vals[0]) / (double)(n_vals - 1);
  }
  return out;
}

/// retrieve the index along the redshift dimension, most closely associated
/// with @p z from a cloudy table, \p table (using bisection)
///
/// This should **ONLY** be used with new-style cloudy tables
///
/// @param z The redshift of interest
/// @param table the cloudy table
///
/// @returns The one-indexed redshift index
///
/// @todo After PR #394 has been merged, we should try to adjust the redshift
///       index so that it is now zero-indexed. I fell pretty strongly that we
///       should also transition from using `long long` values to `int64_t`
[[gnu::always_inline]] inline long long find_zindex(double z,
                                                    const cloudy_data& table) {
  if (table.grid_rank <= 2) {
    return 1LL;
  }

  // reminder (since this looks wrong at a quick glance):
  // -> in a 1D table, axis 0 maps to temperature
  // -> in a 2D table, axis 0 maps to density & axis 1 maps to temperature
  // -> in a 3D table, axis 0 maps to density, axis 1 maps to redshift, &
  //    axis 2 maps to temperature
  const double* z_vals = table.grid_parameters[1];
  const long long n_vals = table.grid_dimension[1];
  if (z <= z_vals[0]) {
    return 1LL;
  } else if (z >= z_vals[n_vals - 2]) {
    return n_vals;
  } else if (z >= z_vals[n_vals - 3]) {
    return n_vals - 2;
  } else {
    long long zindex = 1;
    long long zhighpt = n_vals - 2;
    while ((zhighpt - zindex) > 1) {
      long long zmidpt = (long long)((zhighpt + zindex) / 2);
      if (z >= z_vals[zmidpt - 1]) {
        zindex = zmidpt;
      } else {
        zhighpt = zmidpt;
      }
    }
    return zindex;
  }
}

}  // namespace tabulated_detail
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // TABULATED_COMMON_HPP