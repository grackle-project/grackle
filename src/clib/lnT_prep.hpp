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

/// returns the step in ln(T) of the gas temperature grid commonly shared for
/// interpolating rates
///
/// @todo
/// Consider designing a struct, embedded within @ref opaque_storage that
/// caches this value. It probably makes sense to copy the table bounds from
/// @ref chemistry_data
inline double common_1D_rate_table_lnT_step(
    const chemistry_data& my_chemistry) {
  const int n_bins = my_chemistry.NumberOfTemperatureBins;
  const double dlogtem = (std::log(my_chemistry.TemperatureEnd) -
                          std::log(my_chemistry.TemperatureStart)) /
                         (double)(n_bins - 1);
  return dlogtem;
}

namespace detail {

template <class UnaryFn>
[[gnu::always_inline]] inline double prep_lnT_lininterp_bufs_(
    LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
    const chemistry_data& my_chemistry, const gr_mask_type* itmask,
    UnaryFn get_T_fn) {
  // Get properties of the table
  // todo: stop using long long. I suspect we could just use an int32_t. We only
  //       need larger precision if the table holds more than ~2e9 elements
  //       (i.e. several gigabytes of memory)
  const long long n_bins{my_chemistry.NumberOfTemperatureBins};
  const double logtem_start = std::log(my_chemistry.TemperatureStart);
  const double logtem_end = std::log(my_chemistry.TemperatureEnd);
  const double dlogtem = common_1D_rate_table_lnT_step(my_chemistry);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      logTlininterp_buf.logtem[i] = std::log(get_T_fn(i));
      logTlininterp_buf.logtem[i] = GRIMPL_NS::clamp(
          logTlininterp_buf.logtem[i], logtem_start, logtem_end);

      // Find the corresponding index and precompute interpolation values
      logTlininterp_buf.indixe[i] = GRIMPL_NS::clamp(
          (long long)((logTlininterp_buf.logtem[i] - logtem_start) / dlogtem) +
              1LL,
          1LL, n_bins - 1LL);
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

}  // namespace detail

/// Fills buffers tracked by \p logTlininterp_buf and returns the spacing in
/// logspace
///
/// @note the way that we return the spacing in logspace feels a little "hacky"
inline double prep_lnT_lininterp_bufs(
    LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
    const chemistry_data& my_chemistry, const gr_mask_type* itmask,
    const double* temperature) {
  auto get_T = [temperature](int i) -> double { return temperature[i]; };
  return detail::prep_lnT_lininterp_bufs_(logTlininterp_buf, idx_range,
                                          my_chemistry, itmask, get_T);
}

/// Fills buffers tracked by @p logTlininterp_buf and returns the spacing in
/// logspace. In this overload, we set each value to the arithmetic average of
/// @p cur_T and @p old_T
///
/// @note the way that we return the spacing in logspace feels a little "hacky"
inline double prep_lnT_lininterp_bufs(
    LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
    const chemistry_data& my_chemistry, const gr_mask_type* itmask,
    const double* cur_T, const double* old_T) {
  auto get_T = [cur_T, old_T](int i) -> double {
    return 0.5 * (cur_T[i] + old_T[i]);
  };
  return detail::prep_lnT_lininterp_bufs_(logTlininterp_buf, idx_range,
                                          my_chemistry, itmask, get_T);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // LNT_PREP_HPP