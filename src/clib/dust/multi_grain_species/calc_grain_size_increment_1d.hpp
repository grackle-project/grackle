//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of calc_grain_size_increment_1d
///
//===----------------------------------------------------------------------===//

#ifndef CALC_GRAIN_SIZE_INCREMENT_1D_HPP
#define CALC_GRAIN_SIZE_INCREMENT_1D_HPP

// This file was initially generated automatically during conversion of the
// calc_grain_size_increment_1d function from FORTRAN to C++

#include <cstdio>
#include <limits>
#include <vector>

#include "grackle.h"
#include "dust_props.hpp"
#include "dust/grain_species_info.hpp"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "inject_model/inject_path_field_pack.hpp"
#include "LUT.hpp"
#include "utils-cpp.hpp"
#include "utils-field.hpp"

namespace grackle::impl {

/// For each grain species, compute quantities pertaining to the size
/// distribution along the specified index range.
///
/// In slightly more detail, our model assumes that the a grain species's
/// current differential size distribution is the weighted sum of the
/// initial differential size distributions from each of the injection
/// pathways that is translated by a value known as the "size increment."
/// Internally, this function computes the size increment based on the amount
/// of amount of mass that has been transferred to the grain species (compared
/// to the initial injected amount) via grain growth. For context, if there
/// hasn't been growth, the size increment is zero.
///
/// This function computes the size increment as an intermediate quantity.
/// Then, it uses that information to compute
/// - the current cross-section
/// - a 1D table of opacity-related values (for various possible dust
///   Temperatures).
///
/// @todo
/// Consider writing some narrative documentation about the detailed model of
/// the multi-grain species dust model and revising this docstring so that its
/// less dense (and points the reader to the narrative docs for more details)
///
/// @param[in] dom a standard quantity used throughout the codebase
/// @param[in] idx_range Specifies the current index-range
/// @param[in] itmask_metal Specifies the `idx_range`'s iteration-mask
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] grain_species_info holds information about each grain species
/// @param[in] inject_pathway_props holds data about the modelled injection
///     pathways for all of the grain species.
/// @param[in] my_fields specifies the field data
/// @param[in,out] internal_dust_prop_buf Holds dust-specific information that
///     gets updated by this function
inline void calc_grain_size_increment_1d(
    double dom, IndexRange idx_range, const gr_mask_type* itmask,
    const chemistry_data* my_chemistry,
    const GrainSpeciesInfo* grain_species_info,
    const GrainMetalInjectPathways* inject_pathway_props,
    grackle_field_data* my_fields,
    grackle::impl::InternalDustPropBuf internal_dust_prop_buf) {
  const int n_pathways = inject_pathway_props->n_pathways;
  const int n_log10Tdust_vals = static_cast<int>(
      inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  const int n_opac_poly_coef = inject_pathway_props->n_opac_poly_coef;

  // Step 1: Identify the set of pathways for which the max ratio b/t the
  // pathway's metal density and the total metal density exceeds a threshold.
  //
  // -> the max ratio is taken along idx_range
  // -> we repack the selected metal densities into a temporary buffer
  //
  // TODO: refactor calc_grain_size_increment_species_1d to directly accept the
  // list of selected injection pathway indices so that we can skip the process
  // of repacking.
  // -> repacking is mostly done for historical reasons (i.e. this logic was
  //    originally written in a Fortran dialect that didn't use pointers)
  //
  // TODO: for bitwise reproducibility, we should really be doing the same
  // threshold checking within calc_grain_size_increment_species_1d.
  // -> We want to ensure that a Grackle calculation **NEVER** gives
  //    different results whether it process a bunch of zones or one zone at
  //    a time.
  // -> yes, for a low enough threshold, there won't be a physically meaningful
  //    difference in the results, but we don't want Grackle to give slightly
  //    different numbers based on how its called for a few reasons:
  //    1. It can be EXTREMELY frustrating for someone trying to debug an error
  //       in a simulation code.
  //    2. It can be annoying for constructing automated Grackle tests
  //    3. Slight differences in a calculation can be a fantasitc heuristic for
  //       discovering memory errors and other bugs. But, that requires our
  //       logic to provide bitwise reproducible results.

  constexpr gr_float threshold = 0.01;
  constexpr int max_num_pathways = inj_model_input::N_Injection_Pathways;

  // to be filled with the indices of selected injection pathways
  int selected_inj_path_idx_l[max_num_pathways];

  // to be updated with the number of selected injection pathways
  int n_selected_inj_paths = 0;

  // to be filled with the metal densities for the selected injection paths
  std::vector<gr_float> repacked_inj_path_metal_densities(
      my_fields->grid_dimension[0] * n_pathways);

  // do the work
  {
    grackle::impl::View<gr_float**> SN_metal(
        repacked_inj_path_metal_densities.data(), my_fields->grid_dimension[0],
        inject_pathway_props->n_pathways);

    grackle::impl::View<const gr_float***> metal(
        const_cast<const gr_float*>(my_fields->metal_density),
        my_fields->grid_dimension[0], my_fields->grid_dimension[1],
        my_fields->grid_dimension[2]);

    InjectPathFieldPack inject_path_metal_densities =
        setup_InjectPathFieldPack(my_chemistry, my_fields);

    int start = inject_path_metal_densities.start_idx;
    int stop = inject_path_metal_densities.stop_idx;

    // make arrays
    for (int count = start; count < stop; count++) {
      // when my_chemistry->multi_metals == 0, inj_path_metal_dens wraps
      // the same pointer as `metal`

      grackle::impl::View<const gr_float***> inj_path_metal_dens(
          inject_path_metal_densities.fields[count],
          my_fields->grid_dimension[0], my_fields->grid_dimension[1],
          my_fields->grid_dimension[2]);

      // calculate the max ratio between inj_path_metal_dens and metal
      gr_float max_ratio = std::numeric_limits<gr_float>::lowest();
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        gr_float cur_ratio = (inj_path_metal_dens(i, idx_range.j, idx_range.k) /
                              metal(i, idx_range.j, idx_range.k));
        max_ratio = std::fmax(cur_ratio, max_ratio);
      }

      if (max_ratio > threshold) {
        selected_inj_path_idx_l[n_selected_inj_paths] = count;
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          SN_metal(i, n_selected_inj_paths) =
              inj_path_metal_dens(i, idx_range.j, idx_range.k);
        }
        n_selected_inj_paths++;
      }
    }
  }

  // Step 2: Compute the output quantities for each individual species

  // NOTE: gr_N and gr_Size are historical names
  // -> they are pretty uninformative and should be changed!
  int gr_N[2] = {n_opac_poly_coef, n_log10Tdust_vals};
  int gr_Size = gr_N[0] * gr_N[1];

  std::vector<double> repacked_yields(n_pathways);

  std::vector<double> repacked_size_moments_data_(n_pathways * 3);
  grackle::impl::View<double**> repacked_size_moments(
      repacked_size_moments_data_.data(), 3, n_pathways);

  std::vector<double> repacked_opac_table_data_(n_pathways * gr_Size);
  grackle::impl::View<double**> repacked_opac_table(
      repacked_opac_table_data_.data(), gr_Size, n_pathways);

  grackle::impl::SpeciesLUTFieldAdaptor field_data_adaptor{*my_fields};

  // loop over grain species
  for (int grsp_i = 0; grsp_i < grain_species_info->n_species; grsp_i++) {
    // repack the selected injection pathways for the current grain species
    grackle::impl::View<double**> orig_size_moments(
        inject_pathway_props->size_moments.data[grsp_i], 3, n_pathways);
    grackle::impl::View<double**> orig_opac_table(
        inject_pathway_props->opacity_coef_table.data[grsp_i], gr_Size,
        n_pathways);

    for (int iSN = 0; iSN < n_selected_inj_paths; iSN++) {
      int iSN0 = selected_inj_path_idx_l[iSN];
      repacked_yields[iSN] =
          inject_pathway_props->grain_yields.data[grsp_i][iSN0];
      for (int idx = 0; idx < 3; idx++) {
        repacked_size_moments(idx, iSN) = orig_size_moments(idx, iSN0);
      }
      for (int idx = 0; idx < gr_Size; idx++) {
        repacked_opac_table(idx, iSN) = orig_opac_table(idx, iSN0);
      }
    }

    // now, actually calculate the size increment and subsequent quantities
    const GrainSpeciesInfoEntry& cur_grsp_info =
        grain_species_info->species_info[grsp_i];
    double bulk_density = cur_grsp_info.bulk_density_cgs;
    const gr_float* grsp_density =
        field_data_adaptor.get_ptr_dynamic(cur_grsp_info.species_idx);

    FORTRAN_NAME(calc_grain_size_increment_species_1d)(
        &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
        &my_fields->grid_dimension[0], &my_fields->grid_dimension[1],
        &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end,
        &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
        &n_selected_inj_paths, grsp_density,
        repacked_inj_path_metal_densities.data(), repacked_yields.data(),
        repacked_size_moments.data(), &bulk_density,
        internal_dust_prop_buf.grain_sigma_per_gas_mass.data[grsp_i],
        internal_dust_prop_buf.grain_dyntab_kappa.data[grsp_i], gr_N, &gr_Size,
        &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0],
        inject_pathway_props->log10Tdust_interp_props.parameters[0],
        repacked_opac_table.data());
  }

  // step 3: calculate the total cross-section and the total opacity table
  // (i.e. that include contributions from all grain species)

  // todo: can we skip this when my_chemistry->use_multiple_dust_temperatures
  //   is not 0?

  // todo: clean up to avoid explicitly mention grain species names
  //   -> doing this will probably require an update to the gold standard
  grackle::impl::View<double**> alSiM(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFeM(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMg2SiO4(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::Mg2SiO4_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMgSiO3(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::MgSiO3_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFe3O4(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::Fe3O4_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alAC(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alSiO2D(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMgO(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFeS(
      internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alAl2O3(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::Al2O3_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alreforg(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::ref_org_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alvolorg(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::vol_org_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alH2Oice(
      internal_dust_prop_buf.grain_dyntab_kappa
          .data[OnlyGrainSpLUT::H2O_ice_dust],
      gr_N[1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> altot(internal_dust_prop_buf.dyntab_kappa_tot,
                                      gr_N[1], my_fields->grid_dimension[0]);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      if (my_chemistry->dust_species > 0) {
        internal_dust_prop_buf.sigma_per_gas_mass_tot[i] =
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::MgSiO3_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::AC_dust][i];
      }
      if (my_chemistry->dust_species > 1) {
        internal_dust_prop_buf.sigma_per_gas_mass_tot[i] =
            internal_dust_prop_buf.sigma_per_gas_mass_tot[i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::SiM_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::FeM_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::Mg2SiO4_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::Fe3O4_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::SiO2_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::MgO_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::FeS_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::Al2O3_dust][i];
      }
      if (my_chemistry->dust_species > 2) {
        internal_dust_prop_buf.sigma_per_gas_mass_tot[i] =
            internal_dust_prop_buf.sigma_per_gas_mass_tot[i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::ref_org_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::vol_org_dust][i] +
            internal_dust_prop_buf.grain_sigma_per_gas_mass
                .data[OnlyGrainSpLUT::H2O_ice_dust][i];
      }

      for (int idx = 0; idx < gr_N[1]; idx++) {
        if (my_chemistry->dust_species > 0) {
          altot(idx, i) = alMgSiO3(idx, i) + alAC(idx, i);
        }
        if (my_chemistry->dust_species > 1) {
          altot(idx, i) = altot(idx, i) + alSiM(idx, i) + alFeM(idx, i) +
                          alMg2SiO4(idx, i) + alFe3O4(idx, i) +
                          alSiO2D(idx, i) + alMgO(idx, i) + alFeS(idx, i) +
                          alAl2O3(idx, i);
        }
        if (my_chemistry->dust_species > 2) {
          altot(idx, i) = altot(idx, i) + alreforg(idx, i) + alvolorg(idx, i) +
                          alH2Oice(idx, i);
        }
      }
    }
  }
}

}  // namespace grackle::impl

#endif  // CALC_GRAIN_SIZE_INCREMENT_1D_HPP
