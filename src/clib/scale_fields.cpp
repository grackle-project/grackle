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

/// @brief scale densities of dynamically evolved species
///
/// IMPORTANT: it's valid to call this function with **ANY** @p sp_id_bounds
///    object produced by a call to @ref PartMap::part_bounds (even if that
///    method was passed a partition descriptor not tracked by the underlying
///    object). For reference, this function uses @p sp_id_bounds to create
///    bounds of a for-loop and @ref PartMap::part_bounds guarantees that the
///    bounds are always valid (if the partition descriptor isn't tracked, the
///    loop exits immediately)
static void scale_dyn_sp_densities_(SpeciesMultiView<gr_float> sp_densities,
                                    gr_float factor,
                                    const IndexInterval1D& sp_id_bounds,
                                    IndexRange idx_range) {
  for (int sp_id = sp_id_bounds.start; sp_id < sp_id_bounds.stop; sp_id++) {
    gr_float* ptr = sp_densities.contig1d_ptr(sp_id);
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
  }

  // scale dynamically evolved dust species densities (if any)
  IndexInterval1D dust_sp_bounds = sp_kind_map.part_bounds(SpKind::DUST);
  // it's ok to pass bounds to scale_dyn_sp_densities_, even if ``kind``
  // isn't contained by sp_kind_map (see docstring for more details!)
  scale_dyn_sp_densities_(sp_densities, factor, dust_sp_bounds, idx_range);

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
    IndexInterval1D bounds = sp_kind_map.part_bounds(kind);
    // it's ok to pass bounds to scale_dyn_sp_densities_, even if ``kind``
    // isn't contained by sp_kind_map (see docstring for more details!)
    scale_dyn_sp_densities_(sp_densities, factor, bounds, idx_range);
  }

  scale_inject_path_metal_densities_(my_fields, factor, n_inj_path_ptrs,
                                     idx_range);
}

}  // namespace GRIMPL_NAMESPACE_DECL
