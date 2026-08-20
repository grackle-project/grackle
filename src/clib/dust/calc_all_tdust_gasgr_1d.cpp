//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_all_tdust_gasgr_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_all_tdust_gasgr_1d_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "LUT.hpp"
#include "grackle.h"

#include "calc_tdust_1d.hpp"
#include "calc_all_tdust_gasgr_1d.hpp"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "opaque_storage.hpp"
#include "support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void calc_all_tdust_gasgr_1d(
    double trad, const double* tgas, double* tdust, const double* metallicity,
    const double* dust2gas, const double* nh, double* gasgr_tdust,
    const gr_mask_type* itmask_metal, double coolunit, double* gasgr,
    const double* myisrf, double* kptot, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    IndexRange idx_range, GrainSpeciesCollection grain_temperatures,
    GrainSpeciesCollection gas_grainsp_heatrate,
    LnTLinInterpBuf logTlininterp_buf,
    InternalDustPropBuf internal_dust_prop_buf,
    GrainSpeciesCollection grain_kappa) {
  const double mh_local_var = mh_grflt;
  int i;

  const bool single_species_dust_model = my_chemistry->dust_chemistry == 1;

  int n_grain_species = 0;
  if (!single_species_dust_model) {
    GrainSpeciesInfo* gsp_info = my_rates->opaque_storage->grain_species_info;
    GRIMPL_REQUIRE(gsp_info != nullptr, "sanity check!");
    n_grain_species = gsp_info->n_species;
  }

  // Cooling/heating slice locals

  std::vector<double> gasgr_tbufs[OnlyGrainSpLUT::NUM_ENTRIES];
  std::vector<double> gisrf_bufs[OnlyGrainSpLUT::NUM_ENTRIES];

  for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
    gasgr_tbufs[gsp_idx].resize(my_fields->grid_dimension[0]);
    gisrf_bufs[gsp_idx].resize(my_fields->grid_dimension[0]);
  }
  std::vector<double> mygisrf(my_fields->grid_dimension[0]);

  grackle::impl::GrainMetalInjectPathways* inject_pathway_props =
      my_rates->opaque_storage->inject_pathway_props;

  double dlog10Tdust = 0.0;
  double* log10Tdust_vals = nullptr;

  // NOTE: gr_N is a historical name
  // -> it is pretty uninformative and should be changed!
  int gr_N[2] = {0, 0};
  if (inject_pathway_props != nullptr) {
    dlog10Tdust =
        inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0];
    log10Tdust_vals =
        inject_pathway_props->log10Tdust_interp_props.parameters[0];

    gr_N[0] = inject_pathway_props->n_opac_poly_coef;
    gr_N[1] = static_cast<int>(
        inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  };

  if (single_species_dust_model) {
    // Calculate heating from interstellar radiation field
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        mygisrf[i] = my_rates->gamma_isrf *
                     my_chemistry->local_dust_to_gas_ratio / dust2gas[i] *
                     metallicity[i];
      }
    }

    // Look up gas to grain heat transfer rates
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        gasgr[i] =
            my_rates->gas_grain[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->gas_grain[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->gas_grain[logTlininterp_buf.indixe[i] - 1]);

        gasgr_tdust[i] =
            (dust2gas[i] / metallicity[i]) * gasgr[i] * coolunit / mh_local_var;
      }
    }

    // Compute dust temperature
    calc_tdust_1d(tdust, tgas, nh, gasgr_tdust, mygisrf.data(), myisrf,
                  itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                  dlog10Tdust, log10Tdust_vals,
                  internal_dust_prop_buf.dyntab_kappa_tot, kptot, true,
                  idx_range);

  } else {
    for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
      const double* sigma_per_gas_mass =
          internal_dust_prop_buf.grain_sigma_per_gas_mass.data[gsp_idx];

      // Calculate heating from interstellar radiation field
      double* gisrf_buf = gisrf_bufs[gsp_idx].data();
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          gisrf_buf[i] = my_rates->gamma_isrf2 * sigma_per_gas_mass[i];
        }
      }

      // Look up gas to grain heat transfer rates
      double* gasgr_sp = gas_grainsp_heatrate.data[gsp_idx];
      double* gasgr_tdust_sp = gasgr_tbufs[gsp_idx].data();
      double fac = coolunit / mh_local_var;
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          double fv2k =
              my_rates->gas_grain2[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->gas_grain2[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->gas_grain2[logTlininterp_buf.indixe[i] - 1]);

          gasgr_sp[i] = fv2k * sigma_per_gas_mass[i];
          gasgr_tdust_sp[i] = gasgr_sp[i] * fac;
        }
      }

      // Compute dust temperature
      calc_tdust_1d(grain_temperatures.data[gsp_idx], tgas, nh, gasgr_tdust_sp,
                    gisrf_buf, myisrf, itmask_metal, trad,
                    my_fields->grid_dimension[0], gr_N[1], dlog10Tdust,
                    log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa.data[gsp_idx],
                    grain_kappa.data[gsp_idx], false, idx_range);
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL