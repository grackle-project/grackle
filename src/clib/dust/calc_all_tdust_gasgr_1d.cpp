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
  double fv2k, fac;
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

  // Calculate heating from interstellar radiation field
  //  -> this is ONLY used when `itmask_metal .eq. MASK_TRUE`

  if (single_species_dust_model) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        mygisrf[i] = my_rates->gamma_isrf *
                     my_chemistry->local_dust_to_gas_ratio / dust2gas[i] *
                     metallicity[i];
      }
    }
  } else {
    for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
      double* gisrf_buf = gisrf_bufs[gsp_idx].data();
      const double* sigma_per_gas_mass =
          internal_dust_prop_buf.grain_sigma_per_gas_mass.data[gsp_idx];
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          gisrf_buf[i] = my_rates->gamma_isrf2 * sigma_per_gas_mass[i];
        }
      }
    }
  }

  // --- Gas to grain heat transfer ---

  // Look up gas/grain heat transfer rates

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask_metal[i] != MASK_FALSE) {
      if (single_species_dust_model) {
        gasgr[i] =
            my_rates->gas_grain[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->gas_grain[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->gas_grain[logTlininterp_buf.indixe[i] - 1]);

        gasgr_tdust[i] =
            (dust2gas[i] / metallicity[i]) * gasgr[i] * coolunit / mh_local_var;

      } else if (my_chemistry->dust_chemistry == 2) {
        fv2k = my_rates->gas_grain2[logTlininterp_buf.indixe[i] - 1] +
               logTlininterp_buf.tdef[i] *
                   (my_rates->gas_grain2[logTlininterp_buf.indixe[i] + 1 - 1] -
                    my_rates->gas_grain2[logTlininterp_buf.indixe[i] - 1]);

        fac = coolunit / mh_local_var;

        if (my_chemistry->dust_species > 0) {
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::MgSiO3_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::AC_dust][i];
        }
        if (my_chemistry->dust_species > 1) {
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::SiM_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::FeM_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::Mg2SiO4_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::Fe3O4_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::SiO2_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::MgO_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::FeS_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::Al2O3_dust][i];
        }
        if (my_chemistry->dust_species > 2) {
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::ref_org_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::vol_org_dust][i];
          gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] =
              fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::H2O_ice_dust][i];
        }

        if (my_chemistry->dust_species > 0) {
          gasgr_tbufs[OnlyGrainSpLUT::MgSiO3_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::AC_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] * fac;
        }
        if (my_chemistry->dust_species > 1) {
          gasgr_tbufs[OnlyGrainSpLUT::SiM_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::FeM_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::Mg2SiO4_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::Fe3O4_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::SiO2_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::MgO_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::FeS_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::Al2O3_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] * fac;
        }
        if (my_chemistry->dust_species > 2) {
          gasgr_tbufs[OnlyGrainSpLUT::ref_org_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::vol_org_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] * fac;
          gasgr_tbufs[OnlyGrainSpLUT::H2O_ice_dust][i] =
              gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] * fac;
        }
      }
    }
  }

  // Compute dust temperature

  if (single_species_dust_model) {
    calc_tdust_1d(tdust, tgas, nh, gasgr_tdust, mygisrf.data(), myisrf,
                  itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                  dlog10Tdust, log10Tdust_vals,
                  internal_dust_prop_buf.dyntab_kappa_tot, kptot, true,
                  idx_range);
  } else {
    if (my_chemistry->dust_species > 0) {
      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::MgSiO3_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::MgSiO3_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::MgSiO3_dust],
                    grain_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::AC_dust], tgas, nh,
                    gasgr_tbufs[OnlyGrainSpLUT::AC_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::AC_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::AC_dust],
                    grain_kappa.data[OnlyGrainSpLUT::AC_dust], false,
                    idx_range);
    }

    if (my_chemistry->dust_species > 1) {
      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], tgas, nh,
                    gasgr_tbufs[OnlyGrainSpLUT::SiM_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::SiM_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::SiM_dust],
                    grain_kappa.data[OnlyGrainSpLUT::SiM_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], tgas, nh,
                    gasgr_tbufs[OnlyGrainSpLUT::FeM_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::FeM_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::FeM_dust],
                    grain_kappa.data[OnlyGrainSpLUT::FeM_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::Mg2SiO4_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::Mg2SiO4_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::Mg2SiO4_dust],
                    grain_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::Fe3O4_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::Fe3O4_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::Fe3O4_dust],
                    grain_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::SiO2_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::SiO2_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::SiO2_dust],
                    grain_kappa.data[OnlyGrainSpLUT::SiO2_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], tgas, nh,
                    gasgr_tbufs[OnlyGrainSpLUT::MgO_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::MgO_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::MgO_dust],
                    grain_kappa.data[OnlyGrainSpLUT::MgO_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], tgas, nh,
                    gasgr_tbufs[OnlyGrainSpLUT::FeS_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::FeS_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::FeS_dust],
                    grain_kappa.data[OnlyGrainSpLUT::FeS_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::Al2O3_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::Al2O3_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::Al2O3_dust],
                    grain_kappa.data[OnlyGrainSpLUT::Al2O3_dust], false,
                    idx_range);
    }

    if (my_chemistry->dust_species > 2) {
      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::ref_org_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::ref_org_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::ref_org_dust],
                    grain_kappa.data[OnlyGrainSpLUT::ref_org_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::vol_org_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::vol_org_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::vol_org_dust],
                    grain_kappa.data[OnlyGrainSpLUT::vol_org_dust], false,
                    idx_range);

      calc_tdust_1d(grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], tgas,
                    nh, gasgr_tbufs[OnlyGrainSpLUT::H2O_ice_dust].data(),
                    gisrf_bufs[OnlyGrainSpLUT::H2O_ice_dust].data(), myisrf,
                    itmask_metal, trad, my_fields->grid_dimension[0], gr_N[1],
                    dlog10Tdust, log10Tdust_vals,
                    internal_dust_prop_buf.grain_dyntab_kappa
                        .data[OnlyGrainSpLUT::H2O_ice_dust],
                    grain_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], false,
                    idx_range);
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL