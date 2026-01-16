//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Miscellaneous functionality
///
//===----------------------------------------------------------------------===//
#ifndef GRACKLE_DUST_MISC_HPP
#define GRACKLE_DUST_MISC_HPP

#include "grackle.h"
#include "fortran_func_decls.h"  // gr_mask_type
#include "fortran_func_wrappers.hpp"
#include "index_helper.h"  // IndexHelper
#include "utils-cpp.hpp"   // View
#include "internal_types.hpp"
#include "dust_props.hpp"

namespace grackle::impl {

/// Calculate various dust related properties
///
/// @todo
/// All the following items should wait until after the transcription PRs for
/// the wrapped Fortran routines are merged:
/// - Stop passing @p cool1dmulti_buf and instead pass the `mynh` and
///   `gasgr_tdust` members directly. The former is unmodified while the latter
///   is used as a scratch buffer
/// - Make all non-modified `double*` arguments `const double*`
/// - rename the @p comp2 argument so that it's known as Tcmb
///
/// @param[in] anydust Whether dust chemistry is enabled
/// @param[in] tgas 1d array of gas temperature
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] itmask_metal Specifies the metal/dust-specific iteration-mask of
///     the @p idx_range for this calculation.
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[in] idx_range Specifies the current index-range
/// @param[in] logTlininterp_buf hold values for each location in @p idx_range
///     that are used to linearly interpolate tables with respect to the
///     natural log of @p tgas.
/// @param[in] comp2 Holds the CMB temperature at the current redshift
/// @param[out] dust2gas Holds the computed dust-to-gas ratio at each
///     location in the index range. In other words, this holds the dust mass
///     per unit gas mass (only used in certain configuration)
/// @param[out] tdust, grain_temperatures dust temperatures may be written
///     to one of these variables, based on configuration
/// @param[out] gasgr, gas_grainsp_heatrate Grain/gas energy transfer rates may
///     be written to one of these variables, based on configuration
/// @param[out] kappa_tot, grain_kappa Opacity-related information may be
///     written to one of these variables, based on configuration
/// @param[in,out] cool1dmulti_buf The mynh member holds a 1D array of number
///     densities while the gasgr_tdust member is just a scratch buffer
/// @param[in,out] myisrf a scratch buffer that may be used to temporarily
///     record the interstellar radiation field
/// @param[in,out] internal_dust_prop_buf Holds scratch-space for holding
///     grain-specific information
inline void dust_related_props(
    gr_mask_type anydust, double* tgas, double* metallicity,
    const gr_mask_type* itmask, gr_mask_type* itmask_metal,
    chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
    grackle_field_data* my_fields, InternalGrUnits internalu,
    IndexRange idx_range, LogTLinInterpScratchBuf logTlininterp_buf,
    double comp2, double* dust2gas, double* tdust,
    GrainSpeciesCollection grain_temperatures, double* gasgr,
    GrainSpeciesCollection gas_grainsp_heatrate, double* kappa_tot,
    GrainSpeciesCollection grain_kappa, Cool1DMultiScratchBuf cool1dmulti_buf,
    double* myisrf, InternalDustPropBuf internal_dust_prop_buf) {
  // get relevant unit values
  double dom = internalu_calc_dom_(internalu);
  double coolunit = internalu.coolunit;

  // Compute grain size increment
  if ((my_chemistry->use_dust_density_field > 0) &&
      (my_chemistry->dust_species > 0)) {
    grackle::impl::fortran_wrapper::calc_grain_size_increment_1d(
        dom, idx_range, itmask_metal, my_chemistry, my_rates, my_fields,
        internal_dust_prop_buf);
  }

  // Calculate dust to gas ratio AND interstellar radiation field
  // -> an earlier version of this logic would store values @ indices
  //    where `itmask_metal(i) .ne. MASK_FALSE`
  // -> this was undesirable, b/c these quantities are required for
  //    photo-electric heating, which can occur when
  //    `itmask_metal(i) .eq. MASK_FALSE` (we can revisit this choice
  //    later). Moreover, in most cases, these calculations will be
  //    faster when there is no branching

  if ((anydust != MASK_FALSE) || (my_chemistry->photoelectric_heating > 0)) {
    grackle::impl::View<const gr_float***> d(
        my_fields->density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<const gr_float***> dust(
        my_fields->dust_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    if (my_chemistry->use_dust_density_field > 0) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        // REMINDER: use of `itmask` over `itmask_metal` is
        //   currently required by Photo-electric heating
        if (itmask[i] != MASK_FALSE) {
          // it may be faster to remove this branching
          dust2gas[i] = dust(i, idx_range.j, idx_range.k) /
                        d(i, idx_range.j, idx_range.k);
        }
      }
    } else {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        dust2gas[i] = my_chemistry->local_dust_to_gas_ratio * metallicity[i];
      }
    }
  }

  if ((anydust != MASK_FALSE) || (my_chemistry->photoelectric_heating > 1)) {
    if (my_chemistry->use_isrf_field > 0) {
      grackle::impl::View<const gr_float***> isrf_habing(
          my_fields->isrf_habing, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        myisrf[i] = isrf_habing(i, idx_range.j, idx_range.k);
      }
    } else {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        myisrf[i] = my_chemistry->interstellar_radiation_field;
      }
    }
  }

  // compute dust temperature and cooling due to dust
  if (anydust != MASK_FALSE) {
    grackle::impl::fortran_wrapper::calc_all_tdust_gasgr_1d_g(
        comp2, tgas, tdust, metallicity, dust2gas, cool1dmulti_buf.mynh,
        cool1dmulti_buf.gasgr_tdust, itmask_metal, coolunit, gasgr, myisrf,
        kappa_tot, my_chemistry, my_rates, my_fields, idx_range,
        grain_temperatures, gas_grainsp_heatrate, grain_kappa,
        logTlininterp_buf, internal_dust_prop_buf);
  }
}
}  // namespace grackle::impl

#endif  // GRACKLE_DUST_MISC_HPP
