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

}  // namespace tabulated_detail
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // TABULATED_COMMON_HPP