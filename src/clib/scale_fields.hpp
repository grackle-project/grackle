//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares (and in several cases also defines) field-scaling functions to
/// account for cosmology
///
//===----------------------------------------------------------------------===//

#ifndef SCALE_FIELDS_HPP
#define SCALE_FIELDS_HPP

#include "grackle.h"
#include "index_helper.h"
#include "utils-cpp.hpp"

namespace grackle::impl {

void scale_fields(int imetal, gr_float factor, chemistry_data* my_chemistry,
                  grackle_field_data* my_fields);

/// Scales density and metal_density (if available) by factor
inline void scale_fields_table(grackle_field_data* my_fields, double factor) {
  const int grid_start[3] = {my_fields->grid_start[0], my_fields->grid_start[1],
                             my_fields->grid_start[2]};
  const int grid_end[3] = {my_fields->grid_end[0], my_fields->grid_end[1],
                           my_fields->grid_end[2]};

  // Multiply density by factor (1/a^3 or a^3)

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  for (int k = grid_start[2]; k <= grid_end[2]; k++) {
    for (int j = grid_start[1]; j <= grid_end[1]; j++) {
      for (int i = grid_start[0]; i <= grid_end[0]; i++) {
        d(i, j, k) = d(i, j, k) * factor;
      }
    }
  }

  if (my_fields->metal_density != nullptr) {
    grackle::impl::View<gr_float***> metal(
        my_fields->metal_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    for (int k = grid_start[2]; k <= grid_end[2]; k++) {
      for (int j = grid_start[1]; j <= grid_end[1]; j++) {
        for (int i = grid_start[0]; i <= grid_end[0]; i++) {
          metal(i, j, k) = metal(i, j, k) * factor;
        }
      }
    }
  }
}

/// Scales fields related to computing dust temperature
inline void scale_fields_dust(chemistry_data* my_chemistry,
                              grackle_field_data* my_fields, int imetal,
                              gr_float factor) {
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_loc(
      my_fields->local_ISM_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C13(
      my_fields->ccsn13_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C20(
      my_fields->ccsn20_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C25(
      my_fields->ccsn25_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C30(
      my_fields->ccsn30_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F13(
      my_fields->fsn13_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F15(
      my_fields->fsn15_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F50(
      my_fields->fsn50_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F80(
      my_fields->fsn80_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P170(
      my_fields->pisn170_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P200(
      my_fields->pisn200_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_Y19(
      my_fields->y19_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(
      my_fields->SiM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(
      my_fields->FeM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(
      my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(
      my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(
      my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(
      my_fields->AC_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(
      my_fields->SiO2_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(
      my_fields->MgO_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(
      my_fields->FeS_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(
      my_fields->Al2O3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(
      my_fields->ref_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(
      my_fields->vol_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(
      my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  OMP_PRAGMA("omp parallel for schedule(runtime)")
  for (int t = 0; t < idx_helper.outer_ind_size; t++) {
    // construct an index-range corresponding to "i-slice"
    const IndexRange idx_range = make_idx_range_(t, &idx_helper);
    const int k = idx_range.k;
    const int j = idx_range.j;

    if (imetal == 1) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        metal(i, j, k) = metal(i, j, k) * factor;
        if (my_chemistry->multi_metals > 0) {
          metal_loc(i, j, k) = metal_loc(i, j, k) * factor;
          metal_C13(i, j, k) = metal_C13(i, j, k) * factor;
          metal_C20(i, j, k) = metal_C20(i, j, k) * factor;
          metal_C25(i, j, k) = metal_C25(i, j, k) * factor;
          metal_C30(i, j, k) = metal_C30(i, j, k) * factor;
          metal_F13(i, j, k) = metal_F13(i, j, k) * factor;
          metal_F15(i, j, k) = metal_F15(i, j, k) * factor;
          metal_F50(i, j, k) = metal_F50(i, j, k) * factor;
          metal_F80(i, j, k) = metal_F80(i, j, k) * factor;
          metal_P170(i, j, k) = metal_P170(i, j, k) * factor;
          metal_P200(i, j, k) = metal_P200(i, j, k) * factor;
          metal_Y19(i, j, k) = metal_Y19(i, j, k) * factor;
        }
      }
    }
    if (my_chemistry->use_dust_density_field == 1) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        dust(i, j, k) = dust(i, j, k) * factor;
        if ((my_chemistry->grain_growth == 1) ||
            (my_chemistry->dust_sublimation == 1)) {
          // !            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
          if (my_chemistry->dust_species > 0) {
            MgSiO3(i, j, k) = MgSiO3(i, j, k) * factor;
            AC(i, j, k) = AC(i, j, k) * factor;
          }
          if (my_chemistry->dust_species > 1) {
            SiM(i, j, k) = SiM(i, j, k) * factor;
            FeM(i, j, k) = FeM(i, j, k) * factor;
            Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * factor;
            Fe3O4(i, j, k) = Fe3O4(i, j, k) * factor;
            SiO2D(i, j, k) = SiO2D(i, j, k) * factor;
            MgO(i, j, k) = MgO(i, j, k) * factor;
            FeS(i, j, k) = FeS(i, j, k) * factor;
            Al2O3(i, j, k) = Al2O3(i, j, k) * factor;
          }
          if (my_chemistry->dust_species > 2) {
            reforg(i, j, k) = reforg(i, j, k) * factor;
            volorg(i, j, k) = volorg(i, j, k) * factor;
            H2Oice(i, j, k) = H2Oice(i, j, k) * factor;
          }
          // !            endif
        }
      }
    }
  }
}

void scale_fields_g(int imetal, gr_float factor, chemistry_data* my_chemistry,
                    grackle_field_data* my_fields);

}  // namespace grackle::impl

#endif /* SCALE_FIELDS_HPP */
