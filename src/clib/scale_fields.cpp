//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the scale_fields function
///
//===----------------------------------------------------------------------===//

#include "field_adaptor.hpp"
#include "grackle.h"
#include "scale_fields.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "utils-cpp.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void scale_inject_path_metal_densities_(grackle_field_data* my_fields,
                                        gr_float factor, int n_inj_path_ptrs,
                                        IndexRange idx_range) {
  const int* layout_left_extents = my_fields->grid_dimension;
  int offset = layoutleft_3D_index_to_1D_(layout_left_extents, 0, idx_range.j,
                                          idx_range.k);

  for (int path_idx = 0; path_idx < n_inj_path_ptrs; path_idx++) {
    gr_float* ptr = my_fields->inject_pathway_metal_density[path_idx] + offset;
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      ptr[i] = ptr[i] * factor;
    }
  }
}

void scale_fields_dust(chemistry_data* my_chemistry,
                       grackle_field_data* my_fields, gr_float factor,
                       int n_inj_path_ptrs, IndexRange idx_range,
                       SpeciesMultiView<gr_float> sp_densities,
                       const PartMap<SpKind>& sp_kind_map) {
  FortranView<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> SiM(
      my_fields->SiM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> FeM(
      my_fields->FeM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> Mg2SiO4(
      my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> MgSiO3(
      my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> Fe3O4(
      my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> AC(
      my_fields->AC_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> SiO2D(
      my_fields->SiO2_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> MgO(
      my_fields->MgO_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> FeS(
      my_fields->FeS_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> Al2O3(
      my_fields->Al2O3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> reforg(
      my_fields->ref_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> volorg(
      my_fields->vol_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> H2Oice(
      my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const bool dust_density_field_present =
      my_chemistry->dust_chemistry == 1 &&
      my_chemistry->use_dust_density_field == 1;

  const int k = idx_range.k;
  const int j = idx_range.j;

  if (my_chemistry->metal_cooling == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      metal(i, j, k) = metal(i, j, k) * factor;
    }
  }
  if (dust_density_field_present) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      dust(i, j, k) = dust(i, j, k) * factor;
    }
  } else if (my_chemistry->dust_chemistry == 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
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
    }
  }

  scale_inject_path_metal_densities_(my_fields, factor, n_inj_path_ptrs,
                                     idx_range);
}

void scale_fields(gr_float factor, chemistry_data* my_chemistry,
                  grackle_field_data* my_fields, int n_inj_path_ptrs,
                  IndexRange idx_range, SpeciesMultiView<gr_float> sp_densities,
                  const PartMap<SpKind>& sp_kind_map) {
  FortranView<gr_float***> d(my_fields->density, my_fields->grid_dimension[0],
                             my_fields->grid_dimension[1],
                             my_fields->grid_dimension[2]);
  const bool dust_density_field_present =
      my_chemistry->dust_chemistry == 1 &&
      my_chemistry->use_dust_density_field == 1;

  int j = idx_range.j;
  int k = idx_range.k;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    d(i, j, k) = d(i, j, k) * factor;
  }

  if (my_chemistry->metal_cooling == 1) {
    FortranView<gr_float***> metal(
        my_fields->metal_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      metal(i, j, k) = metal(i, j, k) * factor;
    }
  }

  if (dust_density_field_present) {
    FortranView<gr_float***> dust(
        my_fields->dust_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      dust(i, j, k) = dust(i, j, k) * factor;
    }
  }

  // scale all dynamically evolved species fields
  const SpKind sp_kinds[] = {SpKind::PRIMORDIAL, SpKind::METAL, SpKind::DUST};
  for (SpKind kind : sp_kinds) {
    // reminder: PartMap::part_bounds is designed such that when it's called
    //     for a kind that is not tracked by the kind_map, it will return an
    //     interval where start & stop are equal to each other (thus if it's
    //     used for the bounds of a for-loop, the loop will be skipped)
    IndexInterval1D bounds = sp_kind_map.part_bounds(kind);
    for (int sp_id = bounds.start; sp_id < bounds.stop; sp_id++) {
      gr_float* ptr = sp_densities.contig1d_ptr(sp_id);
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        ptr[i] = ptr[i] * factor;
      }
    }
  }

  scale_inject_path_metal_densities_(my_fields, factor, n_inj_path_ptrs,
                                     idx_range);
}

}  // namespace GRIMPL_NAMESPACE_DECL
