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

// This file was initially generated automatically during conversion of the
// calc_grain_size_increment_1d function from FORTRAN to C++

#include <cstdio>
#include <limits>
#include <vector>

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "inject_model/inject_path_field_pack.hpp"
#include "LUT.hpp"
#include "phys_constants.h"
#include "utils-cpp.hpp"

namespace grackle::impl {

inline void calc_grain_size_increment_1d(
  double dom, IndexRange idx_range,
  const gr_mask_type* itmask,  int* gr_N, int gr_Size,
  const chemistry_data* my_chemistry, grackle_field_data* my_fields,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf,
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props
)
{

  grackle::impl::View<const gr_float***> metal(
      const_cast<const gr_float*>(my_fields->metal_density),
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

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
  std::vector<double> SN_fSiM(inject_pathway_props->n_pathways);
  std::vector<double> SN_fFeM(inject_pathway_props->n_pathways);
  std::vector<double> SN_fMg2SiO4(inject_pathway_props->n_pathways);
  std::vector<double> SN_fMgSiO3(inject_pathway_props->n_pathways);
  std::vector<double> SN_fFe3O4(inject_pathway_props->n_pathways);
  std::vector<double> SN_fAC(inject_pathway_props->n_pathways);
  std::vector<double> SN_fSiO2D(inject_pathway_props->n_pathways);
  std::vector<double> SN_fMgO(inject_pathway_props->n_pathways);
  std::vector<double> SN_fFeS(inject_pathway_props->n_pathways);
  std::vector<double> SN_fAl2O3(inject_pathway_props->n_pathways);
  std::vector<double> SN_freforg(inject_pathway_props->n_pathways);
  std::vector<double> SN_fvolorg(inject_pathway_props->n_pathways);
  std::vector<double> SN_fH2Oice(inject_pathway_props->n_pathways);
  std::vector<double> SN_r0SiM_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0SiM(SN_r0SiM_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0FeM_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0FeM(SN_r0FeM_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0Mg2SiO4_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0Mg2SiO4(SN_r0Mg2SiO4_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0MgSiO3_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0MgSiO3(SN_r0MgSiO3_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0Fe3O4_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0Fe3O4(SN_r0Fe3O4_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0AC_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0AC(SN_r0AC_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0SiO2D_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0SiO2D(SN_r0SiO2D_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0MgO_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0MgO(SN_r0MgO_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0FeS_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0FeS(SN_r0FeS_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0Al2O3_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0Al2O3(SN_r0Al2O3_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0reforg_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0reforg(SN_r0reforg_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0volorg_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0volorg(SN_r0volorg_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_r0H2Oice_data_(3 * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_r0H2Oice(SN_r0H2Oice_data_.data(), 3, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpSiM_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpSiM(SN_kpSiM_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpFeM_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpFeM(SN_kpFeM_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpMg2SiO4_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpMg2SiO4(SN_kpMg2SiO4_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpMgSiO3_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpMgSiO3(SN_kpMgSiO3_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpFe3O4_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpFe3O4(SN_kpFe3O4_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpAC_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpAC(SN_kpAC_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpSiO2D_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpSiO2D(SN_kpSiO2D_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpMgO_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpMgO(SN_kpMgO_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpFeS_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpFeS(SN_kpFeS_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpAl2O3_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpAl2O3(SN_kpAl2O3_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpreforg_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpreforg(SN_kpreforg_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpvolorg_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpvolorg(SN_kpvolorg_data_.data(), gr_Size, inject_pathway_props->n_pathways);
  std::vector<double> SN_kpH2Oice_data_(gr_Size * inject_pathway_props->n_pathways);
  grackle::impl::View<double**> SN_kpH2Oice(SN_kpH2Oice_data_.data(), gr_Size, inject_pathway_props->n_pathways);

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



         
  for (iSN = 1; iSN<=(nSN); iSN++) {
    iSN0 = SN_i[iSN-1];
    if ( my_chemistry->dust_species > 0 )  {
      SN_fMgSiO3     [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::MgSiO3_dust]     [iSN0-1];
      SN_fAC         [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::AC_dust]         [iSN0-1];

      for (idx = 1; idx<=(3); idx++) {
        SN_r0MgSiO3  (idx-1,iSN-1) = SN0_r0MgSiO3  (idx-1,iSN0-1);
        SN_r0AC      (idx-1,iSN-1) = SN0_r0AC      (idx-1,iSN0-1);
      }
      for (idx = 1; idx<=(gr_Size); idx++) {
        SN_kpMgSiO3  (idx-1,iSN-1) = SN0_kpMgSiO3  (idx-1,iSN0-1);
        SN_kpAC      (idx-1,iSN-1) = SN0_kpAC      (idx-1,iSN0-1);
      }
    }
    if ( my_chemistry->dust_species > 1 )  {
      SN_fSiM        [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::SiM_dust]        [iSN0-1];
      SN_fFeM        [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::FeM_dust]        [iSN0-1];
      SN_fMg2SiO4    [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::Mg2SiO4_dust]    [iSN0-1];
      SN_fFe3O4      [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::Fe3O4_dust]      [iSN0-1];
      SN_fSiO2D      [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::SiO2_dust]      [iSN0-1];
      SN_fMgO        [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::MgO_dust]        [iSN0-1];
      SN_fFeS        [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::FeS_dust]        [iSN0-1];
      SN_fAl2O3      [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::Al2O3_dust]      [iSN0-1];

      for (idx = 1; idx<=(3); idx++) {
        SN_r0SiM     (idx-1,iSN-1) = SN0_r0SiM     (idx-1,iSN0-1);
        SN_r0FeM     (idx-1,iSN-1) = SN0_r0FeM     (idx-1,iSN0-1);
        SN_r0Mg2SiO4 (idx-1,iSN-1) = SN0_r0Mg2SiO4 (idx-1,iSN0-1);
        SN_r0Fe3O4   (idx-1,iSN-1) = SN0_r0Fe3O4   (idx-1,iSN0-1);
        SN_r0SiO2D   (idx-1,iSN-1) = SN0_r0SiO2D   (idx-1,iSN0-1);
        SN_r0MgO     (idx-1,iSN-1) = SN0_r0MgO     (idx-1,iSN0-1);
        SN_r0FeS     (idx-1,iSN-1) = SN0_r0FeS     (idx-1,iSN0-1);
        SN_r0Al2O3   (idx-1,iSN-1) = SN0_r0Al2O3   (idx-1,iSN0-1);
      }
      for (idx = 1; idx<=(gr_Size); idx++) {
        SN_kpSiM     (idx-1,iSN-1) = SN0_kpSiM     (idx-1,iSN0-1);
        SN_kpFeM     (idx-1,iSN-1) = SN0_kpFeM     (idx-1,iSN0-1);
        SN_kpMg2SiO4 (idx-1,iSN-1) = SN0_kpMg2SiO4 (idx-1,iSN0-1);
        SN_kpFe3O4   (idx-1,iSN-1) = SN0_kpFe3O4   (idx-1,iSN0-1);
        SN_kpSiO2D   (idx-1,iSN-1) = SN0_kpSiO2D   (idx-1,iSN0-1);
        SN_kpMgO     (idx-1,iSN-1) = SN0_kpMgO     (idx-1,iSN0-1);
        SN_kpFeS     (idx-1,iSN-1) = SN0_kpFeS     (idx-1,iSN0-1);
        SN_kpAl2O3   (idx-1,iSN-1) = SN0_kpAl2O3   (idx-1,iSN0-1);
      }
    }
    if ( my_chemistry->dust_species > 2 )  {
      SN_freforg     [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::ref_org_dust]     [iSN0-1];
      SN_fvolorg     [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::vol_org_dust]     [iSN0-1];
      SN_fH2Oice     [iSN-1] = inject_pathway_props->grain_yields.data[OnlyGrainSpLUT::H2O_ice_dust]     [iSN0-1];

      for (idx = 1; idx<=(3); idx++) {
        SN_r0reforg  (idx-1,iSN-1) = SN0_r0reforg  (idx-1,iSN0-1);
        SN_r0volorg  (idx-1,iSN-1) = SN0_r0volorg  (idx-1,iSN0-1);
        SN_r0H2Oice  (idx-1,iSN-1) = SN0_r0H2Oice  (idx-1,iSN0-1);
      }
      for (idx = 1; idx<=(gr_Size); idx++) {
        SN_kpreforg  (idx-1,iSN-1) = SN0_kpreforg  (idx-1,iSN0-1);
        SN_kpvolorg  (idx-1,iSN-1) = SN0_kpvolorg  (idx-1,iSN0-1);
        SN_kpH2Oice  (idx-1,iSN-1) = SN0_kpH2Oice  (idx-1,iSN0-1);
      }
    }
  }

  double bulk_densities[OnlyGrainSpLUT::NUM_ENTRIES] = {
    sMgSiO3,
    sAC,
    sSiM,
    sFeM,
    sMg2SiO4,
    sFe3O4,
    sSiO2D,
    sMgO,
    sFeS,
    sAl2O3,
    sreforg,
    svolorg,
    sH2Oice,
  };


  // ! calculate size increment

  if (my_chemistry->dust_species > 0)  {
    // !    write(*,*) 'MgSiO3'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->MgSiO3_dust_density,   SN_metal.data(), SN_fMgSiO3.data(),   SN_r0MgSiO3.data(),
             &bulk_densities[OnlyGrainSpLUT::MgSiO3_dust],   internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust],   alMgSiO3.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpMgSiO3.data()
                );

    // !    write(*,*) 'AC'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->AC_dust_density,       SN_metal.data(), SN_fAC.data(),       SN_r0AC.data(),
             &bulk_densities[OnlyGrainSpLUT::AC_dust],       internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust],       alAC.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpAC.data()
                );
  }

  if (my_chemistry->dust_species > 1)  {
    // !    write(*,*) 'SiM'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->SiM_dust_density,      SN_metal.data(), SN_fSiM.data(),      SN_r0SiM.data(),
             &bulk_densities[OnlyGrainSpLUT::SiM_dust],      internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust],      alSiM.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpSiM.data()
                );

    // !    write(*,*) 'FeM'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->FeM_dust_density,      SN_metal.data(), SN_fFeM.data(),      SN_r0FeM.data(),
             &bulk_densities[OnlyGrainSpLUT::FeM_dust],      internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust],      alFeM.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpFeM.data()
                );

    // !    write(*,*) 'Mg2SiO4'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->Mg2SiO4_dust_density,  SN_metal.data(), SN_fMg2SiO4.data(),  SN_r0Mg2SiO4.data(),
             &bulk_densities[OnlyGrainSpLUT::Mg2SiO4_dust],  internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust],  alMg2SiO4.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpMg2SiO4.data()
                );

    // !    write(*,*) 'Fe3O4'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->Fe3O4_dust_density,    SN_metal.data(), SN_fFe3O4.data(),    SN_r0Fe3O4.data(),
             &bulk_densities[OnlyGrainSpLUT::Fe3O4_dust],    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust],    alFe3O4.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpFe3O4.data()
                );

    // !    write(*,*) 'SiO2D'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->SiO2_dust_density,    SN_metal.data(), SN_fSiO2D.data(),    SN_r0SiO2D.data(),
             &bulk_densities[OnlyGrainSpLUT::SiO2_dust],    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust],    alSiO2D.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpSiO2D.data()
                );

    // !    write(*,*) 'MgO'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->MgO_dust_density,      SN_metal.data(), SN_fMgO.data(),      SN_r0MgO.data(),
             &bulk_densities[OnlyGrainSpLUT::MgO_dust],      internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust],      alMgO.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpMgO.data()
                );

    // !    write(*,*) 'FeS'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->FeS_dust_density,      SN_metal.data(), SN_fFeS.data(),      SN_r0FeS.data(),
             &bulk_densities[OnlyGrainSpLUT::FeS_dust],      internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust],      alFeS.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpFeS.data()
                );

    // !    write(*,*) 'Al2O3'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->Al2O3_dust_density,    SN_metal.data(), SN_fAl2O3.data(),    SN_r0Al2O3.data(),
             &bulk_densities[OnlyGrainSpLUT::Al2O3_dust],    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust],    alAl2O3.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpAl2O3.data()
                );
  }

  if (my_chemistry->dust_species > 2)  {
    // !    write(*,*) 'reforg'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->ref_org_dust_density,   SN_metal.data(), SN_freforg.data(),   SN_r0reforg.data(),
             &bulk_densities[OnlyGrainSpLUT::ref_org_dust],   internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust],   alreforg.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpreforg.data()
                );

    // !    write(*,*) 'volorg'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->vol_org_dust_density,   SN_metal.data(), SN_fvolorg.data(),   SN_r0volorg.data(),
             &bulk_densities[OnlyGrainSpLUT::vol_org_dust],   internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust],   alvolorg.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpvolorg.data()
                );

    // !    write(*,*) 'H2Oice'
     FORTRAN_NAME(calc_grain_size_increment_species_1d)(
              &my_chemistry->grain_growth, itmask, &inject_pathway_props->n_pathways,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
             &nSN, my_fields->H2O_ice_dust_density,   SN_metal.data(), SN_fH2Oice.data(),   SN_r0H2Oice.data(),
             &bulk_densities[OnlyGrainSpLUT::H2O_ice_dust],   internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust],   alH2Oice.data(),
             gr_N, &gr_Size, &inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0], inject_pathway_props->log10Tdust_interp_props.parameters[0], SN_kpH2Oice.data()
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

  return;
}

}  // namespace grackle::impl
