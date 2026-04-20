//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines logic for precomputing ln(T) and filling associated buffers with
/// values pertaining to 1D interpolation
///
//===----------------------------------------------------------------------===//

#ifndef LNT_PREP_HPP
#define LNT_PREP_HPP

#include <cmath>

#include "fortran_func_decls.h"  // gr_mask_type
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "support/config.hpp"
#include "utils-cpp.hpp"  // GRIMPL_NS::clamp

namespace GRIMPL_NAMESPACE_DECL {

/// Fills buffers tracked by \p logTlininterp_buf and returns the spacing in
/// logspace
///
/// @note the way that we return the spacing in logspace feels a little "hacky"
inline double prep_lnT_lininterp_bufs(
    LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
    const chemistry_data& my_chemistry, const gr_mask_type* itmask,
    const double* temperature) {
  // Get log values of start and end of lookup tables
  const int n_bins = my_chemistry.NumberOfTemperatureBins;
  const double logtem_start = std::log(my_chemistry.TemperatureStart);
  const double logtem_end = std::log(my_chemistry.TemperatureEnd);
  const double dlogtem = (std::log(my_chemistry.TemperatureEnd) -
                          std::log(my_chemistry.TemperatureStart)) /
                         (double)(n_bins - 1);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      logTlininterp_buf.logtem[i] = std::log(temperature[i]);
      logTlininterp_buf.logtem[i] = GRIMPL_NS::clamp(
          logTlininterp_buf.logtem[i], logtem_start, logtem_end);

      // Find the corresponding index and precompute interpolation values
      logTlininterp_buf.indixe[i] = GRIMPL_NS::clamp(
          (long long)((logTlininterp_buf.logtem[i] - logtem_start) / dlogtem) +
              1LL,
          1LL, static_cast<long long>(n_bins) - 1LL);
      logTlininterp_buf.t1[i] =
          (logtem_start + (logTlininterp_buf.indixe[i] - 1) * dlogtem);
      logTlininterp_buf.t2[i] =
          (logtem_start + (logTlininterp_buf.indixe[i]) * dlogtem);
      logTlininterp_buf.tdef[i] =
          (logTlininterp_buf.logtem[i] - logTlininterp_buf.t1[i]) /
          (logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
    }
  }
  return dlogtem;
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // LNT_PREP_HPP