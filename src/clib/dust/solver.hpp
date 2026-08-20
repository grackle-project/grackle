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
#include "../internal_units.hpp"
#include "../lnT_prep.hpp"
#include "./multi_grain_species/dust_props.hpp"
#include "../support/config.hpp"
#include "../support/index_helper.hpp"
#include "fortran_func_decls.h"

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

  /// @brief indicates whether this solver computes any dust reaction rates
  ///
  /// TODO: consider converting to an instance method that doesn't need to be
  ///       passed any arguments
  static bool any_dust_rxn_rates(const chemistry_data& my_chemistry) noexcept {
    return my_chemistry.dust_chemistry > 0;
  }

  /// @brief Look-up rate for H2 formation on dust & (in certain configurations)
  ///        the grain growth rates for each location in the index-range.
  ///
  /// Importantly, this function must be provided with dust temperature(s) and
  /// a dust2gas buffer that were computed by a prior call to
  /// @ref calc_Tdust_and_chem_contrib.
  ///
  /// This interface function **ONLY** exists to aid in the calculation of
  /// numerical derivatives of reaction rates, under the (potentially
  /// questionable) assumption of constant dust temperatures. In other words,
  /// this function can be eliminated if/when @ref calc_Tdust_and_chem_contrib
  /// becomes capable of computing analytic derivatives.
  ///
  /// > [!note]
  /// > This function should not be invoked unless @ref any_dust_rxn_rates
  /// > returns ``true``.
  ///
  /// @param[in] idx_range Specifies the current index-range
  /// @param[in] tdust Precomputed dust temperatures at each location in the
  ///     index range. This **ONLY** holds meaningful values when using variants
  ///     of the classic 1-field dust-model or using variant of the
  ///     multi-grain-species model where all grains are configured to share a
  ///     single temperature.
  /// @param[in] dust2gas Holds the dust-to-gas ratio at each location in the
  ///     index range. In other words, this holds the dust mass per unit gas
  ///     mass (only used in certain configuration). As a basic rule of thumb,
  ///     this should **always** be passed a buffer that was previously passed
  ///     the same buffer that was used in the call to
  ///     @ref calc_Tdust_and_chem_contrib that computed dust temperature(s).
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
  ///    reaction rates for @p idx_range
  /// @param[inout] internal_dust_prop_scratch_buf Scratch space used to hold
  ///     temporary grain species properties (only used in certain
  ///     configurations)
  void lookup_dust_rxn_rates1d(
      IndexRange idx_range, const double* tdust, const double* dust2gas,
      double dom, const gr_mask_type* itmask_metal,
      chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
      grackle_field_data* my_fields, GrainSpeciesCollection grain_temperatures,
      LnTLinInterpBuf logTlininterp_buf, FullRxnRateBuf rxn_rate_buf,
      InternalDustPropBuf internal_dust_prop_scratch_buf) const;

  /// @brief compute Tdust and other dust contributions to the chemistry network
  ///
  /// This is intended to be the main entry-point for the dust model.
  ///
  /// This function writes results that fall into 2 categories of values:
  /// 1. the primary outputs that directly affect chemistry
  ///    - these include the @p edot, @p alpha_continuum_buf, and
  ///      @p rxn_rate_buf buffers
  ///    - these buffers are allowed to be passed ``nullptr`` to indicate that
  ///      we want skip calculation of a given quantity.
  /// 2. the secondary outputs:
  ///    - these include @p dust2gas and the dust temperature buffers
  ///    - strictly speaking, a chemistry solver shouldn't need to know
  ///      anything about any of these quantities. These are primarily provided
  ///      for the purpose of using @ref lookup_dust_rxn_rates1d if you need
  ///      to numerically estimating partial derivatives for a fixed dust
  ///      temperature
  ///
  /// Thoughts For The Future
  /// -----------------------
  /// We should give some thought to:
  /// - handling scratch space. Currently, we use a LOT of scratch space.
  ///   Comments in the implementation discuss strategies for reducing usage.
  ///   But, we are still going to need some scratch space. I suspect the best
  ///   strategy (for gpu compatibility) is to provide a method that expresses
  ///   how much scratch space is needed, as a function of the idx_range length
  ///   and just pass a buffer with that much space to this function
  /// - whether we need to actually consider @p dust2gass to be an output
  ///   buffer. A fairly compelling case could be made that we should just
  ///   allocate as a local variable (from scratch space) and recompute it when
  ///   we need it. Technically speaking, this introduces additional operations,
  ///   but it's not very expensive (especially compared to the rest of the
  ///   dust-related calculations). Plus it would dramatically simplify
  ///   bookkeeping.
  /// - computing analytic partial derivatives (i.e. with respect to
  ///   temperature, chemical species, and any dust grain species). By doing
  ///   this:
  ///   - we can eliminate @ref lookup_dust_rxn_rates1d (this will simplify a
  ///     bunch of bookeeping)
  ///   - in places where we previously would have used
  ///     @ref lookup_dust_rxn_rates1d, we will stop requiring the implicit
  ///     assumption that dust temperature is constant
  ///   Computing the analytic derivatives should be pretty straight-forward
  ///   (and it will probably be faster than estimating them from finite
  ///   differences). The primary challenge is deciding on an appropriate format
  ///
  /// @param[out] edot 1D array to hold the computed the time derivative of the
  ///     internal energy in the @p idx_range. Contributions are accumulated in
  ///     this buffer. In other words, this function does **NOT** set elements
  ///     to 0 before adding contributions.
  /// @param[out] alpha_continuum_buf buffer to which linear absorption
  ///     coefficients from dust are added (each element is updated in place
  ///     with the sum of its existing value and the contribution from dust). In
  ///     certain configurations this is not actually updated.
  /// @param[out] rxn_rate_buf output buffers to be filled with computed
  ///    reaction rates for @p idx_range
  /// @param[out] dust2gas Holds the computed dust-to-gas ratio at each
  ///     location in the index range. In other words, this holds the dust mass
  ///     per unit gas mass (only used in certain configuration)
  /// @param[out] tdust, grain_temperatures dust temperatures may be written
  ///     to one of these variables, based on configuration
  /// @param[in] tgas 1d array of gas temperature
  /// @param[in] rhoH 1D array of Hydrogen mass densities for the @p idx_range
  /// @param[in] nelec_times_mH 1D array holding the number density of electrons
  ///     (multiplied by the Hydrogen mass) for the @p idx_range
  /// @param[in] metallicity 1d array of metallicities
  /// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
  ///     for this calculation.
  /// @param[in] itmask_metal Specifies the metal/dust-specific iteration-mask
  ///     of the @p idx_range for this calculation.
  /// @param[in] my_chemistry holds a number of configuration parameters.
  /// @param[in] my_rates Holds assorted rate data and other internal
  ///     configuration info.
  /// @param[in] my_fields Specifies the field data.
  /// @param[in] internalu Specifies Grackle's internal unit-system
  /// @param[in] idx_range Specifies the current index-range
  /// @param[in] logTlininterp_buf hold values for each location in @p idx_range
  ///     that are used to linearly interpolate tables with respect to the
  ///     natural log of @p tgas.
  void calc_Tdust_and_chem_contrib(
      double* edot, double* alpha_continuum, FullRxnRateBuf* rxn_rate_buf,
      double* dust2gas, double* tdust,
      GrainSpeciesCollection grain_temperatures, const double* tgas,
      const double* rhoH, const double* nelec_times_mH,
      const double* metallicity, const gr_mask_type* itmask,
      const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
      chemistry_data_storage* my_rates, grackle_field_data* my_fields,
      InternalGrUnits internalu, IndexRange idx_range,
      LnTLinInterpBuf logTlininterp_buf) const;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_SOLVER_HPP