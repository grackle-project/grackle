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
/// @param[in] anydust Whether dust chemistry is enabled
/// @param[in] tgas gas temperature
inline void dust_related_props(
    gr_mask_type anydust, double* tgas, double* metallicity,
    const gr_mask_type* itmask, gr_mask_type* itmask_metal,
    chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
    grackle_field_data* my_fields, InternalGrUnits internalu,
    IndexRange idx_range, LogTLinInterpScratchBuf logTlininterp_buf,
    Cool1DMultiScratchBuf cool1dmulti_buf, double comp2, double* dust2gas,
    double* myisrf, double* tdust, GrainSpeciesCollection grain_temperatures,
    double* gasgr, GrainSpeciesCollection gas_grainsp_heatrate,
    GrainSpeciesCollection grain_kappa, double* kappa_tot,
    InternalDustPropBuf internal_dust_prop_buf) {
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
    // TODO: trad -> comp2
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
