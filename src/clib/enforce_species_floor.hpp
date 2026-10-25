//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// @brief Implements enforce_species_floor
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// ceiling_species_g function from FORTRAN to C++

#include "grackle.h"                     // my_chemistry, my_fields
#include "chem_model/infer_species.hpp"  // SpKind
#include "field_adaptor.hpp"
#include "grackle_macros.h"
#include "LUT.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "support/View.hpp"
#include <cmath>  // std::fmax

#ifndef ENFORCE_SPECIES_FLOOR_HPP
#define ENFORCE_SPECIES_FLOOR_HPP

namespace GRIMPL_NAMESPACE_DECL {

/// @brief infers the species id for HeIII (returns -1, if current
///     configuration doesn't model that species)
inline int infer_species_id_HeIII_(const PartMap<SpKind>& sp_kind_map) {
  IndexInterval1D id_range = sp_kind_map.part_bounds(SpKind::PRIMORDIAL);
  int num_prim_species = id_range.stop - id_range.start;
  if (id_range.start >= 0 && num_prim_species > PrimordialSpLUT::HeIII) {
    return id_range.start + PrimordialSpLUT::HeIII;
  }
  return -1;
}

/// @brief enforces a floor for all species densities for the specified index
///     range
inline void enforce_species_floor(chemistry_data* my_chemistry,
                                  grackle_field_data* my_fields,
                                  IndexRange idx_range,
                                  SpeciesMultiView<gr_float> sp_densities,
                                  const PartMap<SpKind>& sp_kind_map) {
  FortranView<gr_float***> d(my_fields->density, my_fields->grid_dimension[0],
                             my_fields->grid_dimension[1],
                             my_fields->grid_dimension[2]);
  FortranView<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  if (my_chemistry->metal_cooling == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      int k = idx_range.k;
      int j = idx_range.j;
      metal(i, j, k) = std::fmax(metal(i, j, k), tiny_fortran_val);
      if (metal(i, j, k) > d(i, j, k)) {
        eprintf("WARNING: metal density exceeds  total density!\n");
        eprintf("i, j, k, metal, density =  %d %d %d %g %g\n", i, j, k,
                metal(i, j, k), d(i, j, k));
      }
    }
  }

  // infer the species_id for HeIII
  int species_id_HeIII = infer_species_id_HeIII_(sp_kind_map);
  const SpKind sp_kinds[] = {SpKind::PRIMORDIAL, SpKind::METAL, SpKind::DUST};
  for (SpKind kind : sp_kinds) {
    // reminder: PartMap::part_bounds is designed such that when it's called
    //     for a kind that is not tracked by the kind_map, it will return an
    //     interval where start & stop are equal to each other (thus if it's
    //     used for the bounds of a for-loop, the loop will be skipped)
    IndexInterval1D bounds = sp_kind_map.part_bounds(kind);

    for (int sp_id = bounds.start; sp_id < bounds.stop; sp_id++) {
      gr_float* ptr = sp_densities.contig1d_ptr(sp_id);
      gr_float floor_val = (sp_id == species_id_HeIII)
                               ? static_cast<gr_float>(1e-5) * tiny_fortran_val
                               : tiny_fortran_val;
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        ptr[i] = std::fmax(ptr[i], floor_val);
      }
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // ENFORCE_SPECIES_FLOOR_HPP
