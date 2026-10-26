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

void scale_fields(gr_float factor, chemistry_data* my_chemistry,
                  grackle_field_data* my_fields, int n_inj_path_ptrs,
                  IndexRange idx_range, SpeciesMultiView<gr_float> sp_densities,
                  const PartMap<SpKind>& sp_kind_map) {
  FortranView<gr_float***> d(my_fields->density, my_fields->grid_dimension[0],
                             my_fields->grid_dimension[1],
                             my_fields->grid_dimension[2]);
  FortranView<gr_float***> de(
      my_fields->e_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HII(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HeI(
      my_fields->HeI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HeII(
      my_fields->HeII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HeIII(
      my_fields->HeIII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HM(
      my_fields->HM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> H2I(
      my_fields->H2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> DI(
      my_fields->DI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> DII(
      my_fields->DII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HDI(
      my_fields->HDI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> DM(
      my_fields->DM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HDII(
      my_fields->HDII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> HeHII(
      my_fields->HeHII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const bool dust_density_field_present =
      my_chemistry->dust_chemistry == 1 &&
      my_chemistry->use_dust_density_field == 1;

  int j = idx_range.j;
  int k = idx_range.k;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    d(i, j, k) = d(i, j, k) * factor;
  }

  if (my_chemistry->primordial_chemistry > 0) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      de(i, j, k) = de(i, j, k) * factor;
      HI(i, j, k) = HI(i, j, k) * factor;
      HII(i, j, k) = HII(i, j, k) * factor;
      HeI(i, j, k) = HeI(i, j, k) * factor;
      HeII(i, j, k) = HeII(i, j, k) * factor;
      HeIII(i, j, k) = HeIII(i, j, k) * factor;
    }
  }
  if (my_chemistry->primordial_chemistry > 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      HM(i, j, k) = HM(i, j, k) * factor;
      H2I(i, j, k) = H2I(i, j, k) * factor;
      H2II(i, j, k) = H2II(i, j, k) * factor;
    }
  }
  if (my_chemistry->primordial_chemistry > 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      DI(i, j, k) = DI(i, j, k) * factor;
      DII(i, j, k) = DII(i, j, k) * factor;
      HDI(i, j, k) = HDI(i, j, k) * factor;
    }
  }
  if (my_chemistry->primordial_chemistry > 3) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      DM(i, j, k) = DM(i, j, k) * factor;
      HDII(i, j, k) = HDII(i, j, k) * factor;
      HeHII(i, j, k) = HeHII(i, j, k) * factor;
    }
  }

  if (my_chemistry->metal_cooling == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      metal(i, j, k) = metal(i, j, k) * factor;
    }
  }

  if (dust_density_field_present) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      dust(i, j, k) = dust(i, j, k) * factor;
    }
  }

  // scale all dynamically evolved metal & dust species fields
  const SpKind sp_kinds[] = {SpKind::METAL, SpKind::DUST};
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
