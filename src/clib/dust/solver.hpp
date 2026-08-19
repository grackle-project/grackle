//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/implements @ref DustSolver
///
//===----------------------------------------------------------------------===//

#ifndef DUST_SOLVER_HPP
#define DUST_SOLVER_HPP

#include "../full_rxn_rate_buf.hpp"
#include "../internal_types.hpp"
#include "../lnT_prep.hpp"
#include "./multi_grain_species/dust_props.hpp"
#include "../support/config.hpp"
#include "../support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Provides a standard interface for computing the contributions of
///     dust chemistry to the overall chemistry solve
///
/// The basic premise is that:
/// - while performing the full chemistry calculation, the rest of grackle
///   doesn't need to know anything at all about the dust model beyond this
///   object (this is somewhat aspirational, right now make_consistent and
///   the step_rate_ functions still require some knowledge, but we're moving
///   away from that)
/// - objects of this type are immutable. After they are created, they have no
///   mutable state (any scratch buffers will need to be passed in)
class DustSolver {
  // in the future, we'll add some configuration data
  // - example: when using the Chiaki multi-grain-species-growth model, we may
  //   track injection pathway information as part of this type (after we
  //   decouple metal-chemistry from injection pathway information)

public:
  /// @brief default constructor
  DustSolver() = default;

  /// Look-up rate for H2 formation on dust & (in certain configurations) the
  /// grain growth rates for each location in the index-range.
  ///
  /// > [!note]
  /// > This function should not be invoked when we aren't using any dust model
  ///
  /// @param[in] idx_range Specifies the current index-range
  /// @param[in] tdust Precomputed dust temperatures at each location in the
  ///     index range. This **ONLY** holds meaningful values when using variants
  ///     of the classic 1-field dust-model or using variant of the
  ///     multi-grain-species model where all grains are configured to share a
  ///     single temperature.
  /// @param[in] dust2gas Holds the dust-to-gas ratio at each location in the
  ///     index range. In other words, this holds the dust mass per unit gas
  ///     mass (only used in certain configuration)
  /// @param[out] h2dust Buffer that gets filled with the rate for forming
  ///     molecular hydrogen on dust grains. **THIS IS ALWAYS FILLED**
  /// @param[in] dom a standard quantity used throughout the codebase
  /// @param[in] itmask_metal Specifies the iteration-mask for the @p idx_range
  ///     performing metal and dust calculations.
  /// @param[in] my_chemistry holds a number of configuration parameters.
  /// @param[in] my_rates Holds assorted rate data and other internal
  ///     configuration info.
  /// @param[in] my_fields Specifies the field data.
  /// @param[in] grain_temperatures individual grain species temperatures. This
  ///     is only used in certain configurations (i.e. when we aren't using the
  ///     tdust argument)
  /// @param[in] logTlininterp_buf Specifies precomputed arrays of values (for
  ///    each location in the index range) that are used to linearly interpolate
  ///    tables with respect to logT (the natural log of the gas temperature).
  /// @param[out] rxn_rate_buf output buffers to be filled with computed
  /// reaction
  ///    rates for @p idx_range
  /// @param[inout] internal_dust_prop_scratch_buf Scratch space used to hold
  ///     temporary grain species properties (only used in certain
  ///     configurations)
  void lookup_dust_rates1d(
      IndexRange idx_range, const double* tdust, const double* dust2gas,
      double dom, const gr_mask_type* itmask_metal,
      chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
      grackle_field_data* my_fields, GrainSpeciesCollection grain_temperatures,
      LnTLinInterpBuf logTlininterp_buf, FullRxnRateBuf rxn_rate_buf,
      InternalDustPropBuf internal_dust_prop_scratch_buf);
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_SOLVER_HPP