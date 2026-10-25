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
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "support/View.hpp"
#include <cmath>  // std::fmax

#ifndef ENFORCE_SPECIES_FLOOR_HPP
#define ENFORCE_SPECIES_FLOOR_HPP

namespace GRIMPL_NAMESPACE_DECL {

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
  FortranView<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
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

  if (my_chemistry->primordial_chemistry > 0) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      int k = idx_range.k;
      int j = idx_range.j;
      de(i, j, k) = std::fmax(de(i, j, k), tiny_fortran_val);
      HI(i, j, k) = std::fmax(HI(i, j, k), tiny_fortran_val);
      HII(i, j, k) = std::fmax(HII(i, j, k), tiny_fortran_val);
      HeI(i, j, k) = std::fmax(HeI(i, j, k), tiny_fortran_val);
      HeII(i, j, k) = std::fmax(HeII(i, j, k), tiny_fortran_val);
      HeIII(i, j, k) =
          std::fmax(HeIII(i, j, k), (gr_float)(1e-5) * tiny_fortran_val);
    }
  }
  if (my_chemistry->primordial_chemistry > 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      int k = idx_range.k;
      int j = idx_range.j;
      HM(i, j, k) = std::fmax(HM(i, j, k), tiny_fortran_val);
      H2I(i, j, k) = std::fmax(H2I(i, j, k), tiny_fortran_val);
      H2II(i, j, k) = std::fmax(H2II(i, j, k), tiny_fortran_val);
    }
  }
  if (my_chemistry->primordial_chemistry > 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      int k = idx_range.k;
      int j = idx_range.j;
      DI(i, j, k) = std::fmax(DI(i, j, k), tiny_fortran_val);
      DII(i, j, k) = std::fmax(DII(i, j, k), tiny_fortran_val);
      HDI(i, j, k) = std::fmax(HDI(i, j, k), tiny_fortran_val);
    }
  }
  if (my_chemistry->primordial_chemistry > 3) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      int k = idx_range.k;
      int j = idx_range.j;
      DM(i, j, k) = std::fmax(DM(i, j, k), tiny_fortran_val);
      HDII(i, j, k) = std::fmax(HDII(i, j, k), tiny_fortran_val);
      HeHII(i, j, k) = std::fmax(HeHII(i, j, k), tiny_fortran_val);
    }
  }
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
        ptr[i] = std::fmax(ptr[i], tiny_fortran_val);
      }
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // ENFORCE_SPECIES_FLOOR_HPP
