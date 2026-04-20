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
[[gnu::always_inline]] inline void prep_lnT_lininterp_bufs_(
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
}

}  // namespace detail

/// Fills buffers tracked by @p logTlininterp_buf. In this overload, we set each
/// value to the arithmetic average of @p cur_T and @p old_T
inline void prep_lnT_lininterp_bufs(LogTLinInterpScratchBuf& logTlininterp_buf,
                                    IndexRange idx_range,
                                    const chemistry_data& my_chemistry,
                                    const gr_mask_type* itmask,
                                    const double* cur_T, const double* old_T) {
  auto get_T = [cur_T, old_T](int i) -> double {
    return 0.5 * (cur_T[i] + old_T[i]);
  };
  detail::prep_lnT_lininterp_bufs_(logTlininterp_buf, idx_range, my_chemistry,
                                   itmask, get_T);
}

/// Encapsulates logic for computing ln(T) related values
///
/// At face value instances fulfill a simple purpose:
/// - it computes the natural log of the gas temperature; this is useful
///   because log calculations are computationally expensive and most of
///   Grackle's rates are interpolated from grids where ln(T) or log10(T)
///   varies along an axis
/// - it also fills buffers within a @ref LogTLinInterpScratchBuf instance with
///   interpolation variable values that are used with the gas temperature grid
///   commonly shared for interpolating rates
///
/// This type fulfills its purpose, with two noteworthy "features"
/// - it applies "dampening" (i.e. the value passed to the natural log function
///   is the average of the provided Temperature and the temperature recorded
///   during a previous subcycle)
/// - it clamps the computed ln(T) based on the bounds of the gas temperature
///   grid commonly shared for interpolating rates
///
/// It's worth emphasizing that this class was designed to encode existing
/// behavior. As we'll discuss below, the extra "features" should probably be
/// decoupled from the act of computing ln(T) since they lead to inconsistencies
///
/// Usage Notes
/// ===========
/// At the time of writing, this class is unique among the entities that are
/// temporarily constructed within the Grackle routines that contain buffers.
/// Unlike other entities, the buffer in an instance of this type must be
/// preserved between subcycles. This extra complexity is inherent to the
/// damping feature.
///
/// An important invariant of this type is that the contained buffer must
/// outlive an instance of this type. The choice to do this (rather than say
/// embed a std::vector) was done to ease the transition to GPUs, but I could be
/// convinced to do something different.
///
/// Implementation Notes
/// ====================
/// I chose to implement this as a more traditional C++ class mostly out of
/// convenience. It could easily be converted to a more C-like style
///
/// Improving Design
/// ================
/// The design could absolutely be improved. This section discusses a number
/// questions and concerns that should be addressed in any refactoring
///
/// The way that this logic computes ln(T) and the way that temperature values
/// are use inside of Grackle is inherently inconsistent. Historically speaking:
/// - Recall: we applied damping and clamping to the ln(T) values and then used
///   those values in most places. There's nothing inherently wrong with this.
/// - Inconsistencies arose because we still continued to use the raw undamped
///   and unclamped temperature values (that would now be passed to
///   @ref prep_damped_lnT_lininterp_bufs) in a few places.
/// **To be fair,** these inconsistencies probably had a negligible impact
/// on Grackle's results.
///
/// However, the number and significance of these inconsistencies have been
/// amplified with by the introduction of the Multi-Grain-Inject-Pathway Dust
/// model. Specifically, the number of places where the temperature value is
/// undamped and unclamped has grown dramatically.
///
/// Perhaps most significantly, @ref cool1d_multi_g (mostly) used damped &
/// clamped ln(T) values while @ref lookup_cool_rates1d started using entirely
/// undamped (but still clamped) values. This was probably an oversight.
///
/// For the sake of this discussion, let's assume that we've merged in changes
/// to ensure that the exact same ln(T) values are used in both functions (at
/// least when they are invoked by @ref solve_rate_cool)
///
/// We still need to resolve a number of questions:
///
/// What is the nature of dampening?
/// --------------------------------
/// At some level its a numerical hack (that shouldn't really matter for
/// infintesimal timesteps), right?
///
/// Should we be treating it as an intrinsic part of the temperature
/// calculation? In other words should all references the gas temperature
/// in any part of the chemistry/energy integration reflect the impact of
/// dampening?
///
/// **OR, is dampening stepper-specific?** The different "steppers" all do
/// different things:
/// - Gauss-Seidel step: everything is damped
/// - decoupled energy-chemistry Newton-Raphson step: energy evolution uses
///   dampening and chemistry uses undamped
/// - coupled energy-chemistry Newton-Raphson step: every is undamped
///
/// In principle, I am not opposed to different steppers doing different things,
/// but it seems like a bad idea if we allow this in a hybrid scheme where the
/// stepping scheme varies between zones during a single subcycle (or varies
/// for a single zone across subcycles).
///
/// What is the nature of the clamped-bounds?
/// -----------------------------------------
/// This is typically, much less important since the default bounds are quite
/// extreme. But while we should try to make a decision about while we are
/// thinking about this general topic.
///
/// What about derivatives?
/// -----------------------
/// This whole discussion becomes MUCH more serious when you start to think
/// about the fact the Newton-Raphson stepper needs to compute partial
/// derivatives with respect to various quantities. Under the current
/// implementation the effect of whether dampening is modelled is handled
/// implicitly (via finite differences). Thus, it may not seem like a very
/// important choice right now.
///
/// However, we should probably make a choices compatible with a future where
/// we directly compute the parital derivatives.
/// - This approach would involve chain rule to explicitly compute temperature
///   derivatives. It's worth noting that such an approach wouldn't actually be
///   that difficult, it would probably be cheaper than finite differences, and
///   it would be definitely be more robust than finite differences (i.e. it
///   entirely sidesteps issues of picking the proper step size).
/// - It's obvious that under this approach, the use of dampening would dampen
///   the temperature derivative
///
/// @note
/// Under the current implementation, the choice to include dampening (or not)
/// is only relevant when energy and chemistry are coupled. When energy and
/// chemistry are operator split, the impact of species abundances on
/// temperature is ignored. In the future, this SHOULD absolutely change.
/// First, it should be relatively cheap if we are explicitly computing the
/// values (rather than finite differences). Moreover, it could be important.
/// Consider a cell with high H2 abundance that had a bunch of thermal energy
/// injected since the previous cycle (e.g. due to a shock). The impact on
/// temperature from rapid H2 changes (due to impacts on mmw & gamma) could be
/// significant
class LnTPreparer {
  double* old_T_;

public:
  /// Construct an instance that wraps an externally manged buffer for tracking
  /// Temperature values
  explicit LnTPreparer(double* old_T_buf) : old_T_{old_T_buf} {};

  /// record the current temperature values
  ///
  /// these values will be used for damping in subsequent calls to
  /// @ref prep_damped_lnT_lininterp_bufs
  void record_T(IndexRange idx_range, const gr_mask_type* itmask,
                const double* cur_T) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        this->old_T_[i] = cur_T[i];
      }
    }
  }

  /// Fills buffers tracked by \p logTlininterp_buf
  ///
  /// The provided temperature values are used directly (there is no damping,
  /// but clamping is still performed)
  static void prep_undamped_lnT_lininterp_bufs(
      LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
      const chemistry_data& my_chemistry, const gr_mask_type* itmask,
      const double* temperature) {
    auto get_T = [temperature](int i) -> double { return temperature[i]; };
    detail::prep_lnT_lininterp_bufs_(logTlininterp_buf, idx_range, my_chemistry,
                                     itmask, get_T);
  }

  /// Fills buffers tracked by @p logTlininterp_buf
  ///
  /// It compute the natural log of the damped temperature value (i.e. the
  /// arithmetic mean of @p cur_T and the previously recorded value)
  void prep_damped_lnT_lininterp_bufs(
      LogTLinInterpScratchBuf& logTlininterp_buf, IndexRange idx_range,
      const chemistry_data& my_chemistry, const gr_mask_type* itmask,
      const double* cur_T) const {
    const double* old_T = old_T_;
    auto get_T = [cur_T, old_T](int i) -> double {
      return 0.5 * (cur_T[i] + old_T[i]);
    };
    detail::prep_lnT_lininterp_bufs_(logTlininterp_buf, idx_range, my_chemistry,
                                     itmask, get_T);
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // LNT_PREP_HPP