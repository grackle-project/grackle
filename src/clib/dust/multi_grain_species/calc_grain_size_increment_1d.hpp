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
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf
)
{
  const int n_pathways = inject_pathway_props->n_pathways;
  const int n_log10Tdust_vals = static_cast<int>(
      inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  const int n_opac_poly_coef = inject_pathway_props->n_opac_poly_coef;

  // NOTE: gr_N and gr_Size are historical names
  // -> they are pretty uninformative and should be changed!
  int gr_N[2] = {n_opac_poly_coef, n_log10Tdust_vals};
  int gr_Size = gr_N[0] * gr_N[1];


  grackle::impl::View<const gr_float***> metal(
      const_cast<const gr_float*>(my_fields->metal_density),
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

  grackle::impl::SpeciesLUTFieldAdaptor field_data_adaptor{*my_fields};

  // table
  grackle::impl::View<double**> SN0_r0SiM(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::SiM_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0FeM(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::FeM_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0Mg2SiO4(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Mg2SiO4_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0MgSiO3(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::MgSiO3_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0Fe3O4(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Fe3O4_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0AC(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::AC_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0SiO2D(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::SiO2_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0MgO(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::MgO_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0FeS(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::FeS_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0Al2O3(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Al2O3_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0reforg(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::ref_org_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0volorg(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::vol_org_dust], 3, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_r0H2Oice(inject_pathway_props->size_moments.data[OnlyGrainSpLUT::H2O_ice_dust], 3, inject_pathway_props->n_pathways);
  // opacity table
  grackle::impl::View<double**> SN0_kpSiM(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::SiM_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpFeM(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::FeM_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpMg2SiO4(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::Mg2SiO4_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpMgSiO3(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::MgSiO3_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpFe3O4(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::Fe3O4_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpAC(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::AC_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpSiO2D(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::SiO2_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpMgO(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::MgO_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpFeS(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::FeS_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpAl2O3(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::Al2O3_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpreforg(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::ref_org_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpvolorg(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::vol_org_dust], gr_Size, inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN0_kpH2Oice(inject_pathway_props->opacity_coef_table.data[OnlyGrainSpLUT::H2O_ice_dust], gr_Size, inject_pathway_props->n_pathways);
  // out
  grackle::impl::View<double**> alSiM(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFeM(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMg2SiO4(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMgSiO3(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFe3O4(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alAC(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alSiO2D(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alMgO(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alFeS(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alAl2O3(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Al2O3_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alreforg(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::ref_org_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alvolorg(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::vol_org_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> alH2Oice(internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], gr_N[2-1], my_fields->grid_dimension[0]);
  grackle::impl::View<double**> altot(internal_dust_prop_buf.dyntab_kappa_tot, gr_N[2-1], my_fields->grid_dimension[0]);
  // array
  int iSN, iSN0;
  std::vector<int> SN_i(inject_pathway_props->n_pathways);
  std::vector<gr_float> SN_metal_data_(my_fields->grid_dimension[0] * inject_pathway_props->n_pathways);
  grackle::impl::View<gr_float**> SN_metal(SN_metal_data_.data(), my_fields->grid_dimension[0], inject_pathway_props->n_pathways);

  // local
  int i, idx;

  InjectPathFieldPack inject_path_metal_densities = setup_InjectPathFieldPack(
    my_chemistry, my_fields);

  int start = inject_path_metal_densities.start_idx;
  int stop = inject_path_metal_densities.stop_idx;

  // make arrays
  int nSN = 0;
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
      gr_float cur_ratio = (
          inj_path_metal_dens(i, idx_range.j, idx_range.k) /
          metal(i, idx_range.j, idx_range.k));
      max_ratio = std::fmax(cur_ratio, max_ratio);
    }

    if (max_ratio > 0.01) {
      nSN++;

      SN_i[nSN-1] =  count + 1;
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        SN_metal(i-1,nSN-1) = inj_path_metal_dens(i-1,idx_range.jp1-1,idx_range.kp1-1);
      }
    }
  }

  // allocate some buffers
  std::vector<double> repacked_yields(n_pathways);

  std::vector<double> repacked_size_moments_data_(n_pathways * 3);
  grackle::impl::View<double**> repacked_size_moments(
      repacked_size_moments_data_.data(), 3, n_pathways);

  std::vector<double> repacked_opac_table_data_(n_pathways * gr_Size);
  grackle::impl::View<double**> repacked_opac_table(
      repacked_opac_table_data_.data(), gr_Size, n_pathways);

  // loop over grain species
  for (int grsp_i = 0; grsp_i < grain_species_info->n_species; grsp_i++) {

    // here, we repack the injection pathway for the current grain species
    grackle::impl::View<double**> orig_size_moments(
        inject_pathway_props->size_moments.data[grsp_i], 3, n_pathways);
    grackle::impl::View<double**> orig_opac_table(
        inject_pathway_props->opacity_coef_table.data[grsp_i], gr_Size,
        n_pathways);

    for (iSN = 1; iSN<=(nSN); iSN++) {
      iSN0 = SN_i[iSN-1];
      repacked_yields[iSN-1] = inject_pathway_props->grain_yields.data[grsp_i][iSN0-1];
      for (idx = 1; idx<=(3); idx++) {
        repacked_size_moments(idx-1,iSN-1) = orig_size_moments(idx-1,iSN0-1);
      }
      for (idx = 1; idx<=(gr_Size); idx++) {
        repacked_opac_table(idx-1,iSN-1) = orig_opac_table(idx-1,iSN0-1);
      }
    }

    // actually calculate the size increment and subsequent quantities
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
        &nSN, grsp_density, SN_metal.data(),
        repacked_yields.data(),
        //reduced_inject_paths.grain_yields.data[grsp_i],
        repacked_size_moments.data(),
        //reduced_inject_paths.size_moments.data[grsp_i],
        &bulk_density,
        internal_dust_prop_buf.grain_sigma_per_gas_mass.data[grsp_i],
        internal_dust_prop_buf.grain_dyntab_kappa.data[grsp_i],
        gr_N, &gr_Size,
        &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0],
        inject_pathway_props->log10Tdust_interp_props.parameters[0], 
        //reduced_inject_paths.opacity_coef_table.data[grsp_i],
        repacked_opac_table.data()
        );
  }

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
            
      if (my_chemistry->dust_species > 0)  {
        internal_dust_prop_buf.sigma_per_gas_mass_tot  [i-1] = internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust]    [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust]        [i-1];
      }
      if (my_chemistry->dust_species > 1)  {
        internal_dust_prop_buf.sigma_per_gas_mass_tot  [i-1] = internal_dust_prop_buf.sigma_per_gas_mass_tot       [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust]       [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust]       [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust]   [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust]     [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust]     [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust]       [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust]       [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust]     [i-1];
      }
      if (my_chemistry->dust_species > 2)  {
        internal_dust_prop_buf.sigma_per_gas_mass_tot  [i-1] = internal_dust_prop_buf.sigma_per_gas_mass_tot[i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust]    [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust]    [i-1]
                   + internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust]    [i-1];
      }
            
      for (idx = 1; idx<=(gr_N [ 2-1 ]); idx++) {
        if (my_chemistry->dust_species > 0)  {
          altot(idx-1,i-1) = alMgSiO3  (idx-1,i-1)
                       + alAC      (idx-1,i-1);
        }
        if (my_chemistry->dust_species > 1)  {
          altot(idx-1,i-1) = altot     (idx-1,i-1)
                       + alSiM     (idx-1,i-1)
                       + alFeM     (idx-1,i-1)
                       + alMg2SiO4 (idx-1,i-1)
                       + alFe3O4   (idx-1,i-1)
                       + alSiO2D   (idx-1,i-1)
                       + alMgO     (idx-1,i-1)
                       + alFeS     (idx-1,i-1)
                       + alAl2O3   (idx-1,i-1);
        }
        if (my_chemistry->dust_species > 2)  {
          altot(idx-1,i-1) = altot     (idx-1,i-1)
                       + alreforg  (idx-1,i-1)
                       + alvolorg  (idx-1,i-1)
                       + alH2Oice  (idx-1,i-1);
        }
      }

    }
  }
}

}  // namespace grackle::impl


#endif  // CALC_GRAIN_SIZE_INCREMENT_1D_HPP
