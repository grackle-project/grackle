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
