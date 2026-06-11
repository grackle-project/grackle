//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares logic pertaining to solving dust chemistry
///
//===----------------------------------------------------------------------===//

#ifndef DUST_SOLVER_HPP
#define DUST_SOLVER_HPP

#include "../field_adaptor.hpp"
#include "../full_rxn_rate_buf.hpp"
#include "../internal_types.hpp"
#include "../internal_units.hpp"
#include "../lnT_prep.hpp"
#include "./multi_grain_species/dust_props.hpp"
#include "../support/config.hpp"
#include "../support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {

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
///     index range. In other words, this holds the dust mass per unit gas mass
///     (only used in certain configuration)
/// @param[out] h2dust Buffer that gets filled with the rate for forming
///     molecular hydrogen on dust grains. **THIS IS ALWAYS FILLED**
/// @param[in] dom a standard quantity used throughout the codebase
/// @param[in] itmask_metal Specifies the iteration-mask for the @p idx_range
///     performing metal and dust calculations.
/// @param[in] dt See the warning at the end of the docstring
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] sp_densities Specifies the densities of the various species
///     that Grackle evolves (if any) in a format that allows the values to be
///     accessed with the index lookup table. Wherever possible, data should be
///     be accessed through this argument, rather than with @p my_fields
/// @param[in] grain_temperatures individual grain species temperatures. This
///     is only used in certain configurations (i.e. when we aren't using the
///     tdust argument)
/// @param[in] logTlininterp_buf Specifies precomputed arrays of values (for
///    each location in the index range) that are used to linearly interpolate
///    tables with respect to logT (the natural log of the gas temperature).
/// @param[out] rxn_rate_buf output buffers to be filled with computed reaction
///    rates for @p idx_range
/// @param[inout] internal_dust_prop_scratch_buf Scratch space used to hold
///     temporary grain species properties (only used in certain configurations)
///
/// > [!important]
/// > TODO: The role of the `dt` argument **MUST** be clarified! It is passed
/// > different values in different areas of the codebase!!!!
/// > - `solve_rate_cool_g` passes in the value of the total timestep that the
/// >   chemistry is evolved. This is the traditional meaning of `dt`
/// > - the time derivative calculation within `step_rate_newton_raphson`
/// >   passes the timestep of the current subcycle (effectively the whole
/// >   function is only being called for a single element idx_range)
/// >
/// > Internally, this arg only appears to be used to determine dust grain
/// > destruction rate.
/// > - the dust destruction rate is 0 for all temperatures below some
/// >   threshold (the threshold depends on the grain species)
/// > - above the threshold, the destruction rate is essentially the current
/// >   grain density divided by the value of the `dt` argument
/// >
/// > If you think about it:
/// > - I'd argue that setting `dt` to the whole timestep that we are evolving
/// >   the zone over is blatantly wrong. It violates the principle that you
/// >   should get consistent results whether you invoke grackle 100 separate
/// >   times or just 1 time. (The amount of dust heating would change)
/// > - setting `dt` to the current subcycle timestep makes a lot more sense
/// >   (and is the only logical choice)
/// >   - It is roughly equivalent to saying that dust is immediately destroyed
/// >     once the gas reaches a threshold temperature.
/// >   - the model is overly simplistic since dust grains can survive for
/// >     quite in ionized gas (see for example
/// >     https://ui.adsabs.harvard.edu/abs/2024ApJ...974...81R/abstract)
/// >
/// > If we stick with this instantaneous destruction model, then all
/// > dust-grain related heating and cooling should probably assume that the
/// > dust-grain density is already 0.
void lookup_dust_rates1d(IndexRange idx_range, const double* tdust,
                         const double* dust2gas, double dom,
                         const gr_mask_type* itmask_metal, double dt,
                         chemistry_data* my_chemistry,
                         chemistry_data_storage* my_rates,
                         grackle_field_data* my_fields,
                         SpeciesMultiView<const gr_float> sp_densities,
                         GrainSpeciesCollection grain_temperatures,
                         LnTLinInterpBuf logTlininterp_buf,
                         FullRxnRateBuf rxn_rate_buf,
                         InternalDustPropBuf internal_dust_prop_scratch_buf);

/// this is a helper function that handles all dust contributions pertaining
/// to cool1d_multi_g
///
/// @param[out] edot 1D array to hold the computed the time derivative of the
///     internal energy in the @p idx_range. Contributions are accumulated in
///     this buffer. In other words, this function does **NOT** set elements to
///     to 0 before adding contributions.
/// @param[out] dust2gas Holds the computed dust-to-gas ratio at each
///     location in the index range. In other words, this holds the dust mass
///     per unit gas mass (only used in certain configuration)
/// @param[out] tdust, grain_temperatures dust temperatures may be written
///     to one of these variables, based on configuration
/// @param[out] alpha_continuum buffer to which linear absorption
///     coefficients from dust are added (each element is updated in place with
///     the sum of its existing value and the contribution from dust). In
///     certain configurations this is not actually updated.
/// @param[in] tgas 1d array of gas temperature
/// @param[in] rhoH 1D array of Hydrogen mass densities for the @p idx_range
/// @param[in] nelec_times_mH 1D array holding the number density of electrons
///     (multiplied by the Hydrogen mass) for the @p idx_range
/// @param[in] metallicity 1d array of metallicities
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] itmask_metal Specifies the metal/dust-specific iteration-mask of
///     the @p idx_range for this calculation.
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] sp_densities Specifies the densities of the various species
///     that Grackle evolves (if any) in a format that allows the values to be
///     accessed with the index lookup table. Wherever possible, data should be
///     be accessed through this argument, rather than with @p my_fields
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[in] idx_range Specifies the current index-range
/// @param[in] logTlininterp_buf hold values for each location in @p idx_range
///     that are used to linearly interpolate tables with respect to the
///     natural log of @p tgas.
///
/// @note
/// In some sense, this is a step towards factoring out all of the dust logic.
/// - we need to be careful with this logic to avoid making logic harder to
///   follow.
void handle_dust_cooling_contributions(
    double* edot, double* dust2gas, double* tdust,
    GrainSpeciesCollection grain_temperatures, double* alpha_continuum,
    const double* tgas, const double* rhoH, const double* nelec_times_mH,
    const double* metallicity, const gr_mask_type* itmask,
    const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    const SpeciesMultiView<const gr_float> sp_densities,
    InternalGrUnits internalu, IndexRange idx_range,
    LnTLinInterpBuf logTlininterp_buf);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_SOLVER_HPP