//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements step_rate_gauss_seidel
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// step_rate_g function from FORTRAN to C++

#ifndef STEP_RATE_GAUSS_SEIDEL_HPP
#define STEP_RATE_GAUSS_SEIDEL_HPP

#include <cmath>
#include <cstdio>

#include "grackle.h"  // gr_float
#include "chemistry_solver_funcs.hpp"
#include "chem_model/infer_species.hpp"
#include "fortran_func_decls.h"  // gr_mask_type
#include "internal_types.hpp"
#include "LUT.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {

namespace step_gs_detail{

/// calculates electron density from other species tracked by @p spdens and
/// stores the result in the slot reserved for electron density
///
/// @note this is essentially equivalent to logic written in make_consistent
inline void update_electron_densities(IndexRange idx_range,
                                      const gr_mask_type* itmask,
                                      const gr_mask_type* itmask_metal,
                                      const chemistry_data* my_chemistry,
                                      SpeciesMultiView<gr_float> spdens) {

  // TODO: stop using SpLUT and stop accessing species unless we know that
  //       grackle is configured to actually use them
  // -> the proper way to fix this is to precalculate information about the
  //    relevant species at startup (I have a plan for that)
  gr_float* de = spdens.contig1d_ptr(SpLUT::e);
  const gr_float* HII = spdens.contig1d_ptr(SpLUT::HII);
  const gr_float* HeII = spdens.contig1d_ptr(SpLUT::HeII);
  const gr_float* HeIII = spdens.contig1d_ptr(SpLUT::HeIII);
  const gr_float* HM = spdens.contig1d_ptr(SpLUT::HM);
  const gr_float* H2II = spdens.contig1d_ptr(SpLUT::H2II);
  const gr_float* DII = spdens.contig1d_ptr(SpLUT::DII);
  const gr_float* DM = spdens.contig1d_ptr(SpLUT::DM);
  const gr_float* HDII = spdens.contig1d_ptr(SpLUT::HDII);
  const gr_float* HeHII = spdens.contig1d_ptr(SpLUT::HeHII);
  const gr_float* CII = spdens.contig1d_ptr(SpLUT::CII);
  const gr_float* COII = spdens.contig1d_ptr(SpLUT::COII);
  const gr_float* OII = spdens.contig1d_ptr(SpLUT::OII);
  const gr_float* OHII = spdens.contig1d_ptr(SpLUT::OHII);
  const gr_float* H2OII = spdens.contig1d_ptr(SpLUT::H2OII);
  const gr_float* H3OII = spdens.contig1d_ptr(SpLUT::H3OII);
  const gr_float* O2II = spdens.contig1d_ptr(SpLUT::O2II);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      de[i] = HII[i] + HeII[i]/(gr_float)(4.) +
           HeIII[i]/(gr_float)(2.);
      if (my_chemistry->primordial_chemistry > 1)
           { de[i] = de[i] - HM[i] + H2II[i]/(gr_float)(2.); }

      if (my_chemistry->primordial_chemistry > 2)
           { de[i] = de[i] + DII[i]/(gr_float)(2.); }
      if (my_chemistry->primordial_chemistry > 3)
           { de[i] = de[i] - DM[i]/(gr_float)(2.)
                + HDII[i]/(gr_float)(3.) + HeHII[i]/(gr_float)(5.); }
      if ( (my_chemistry->metal_chemistry == 1)  &&
           (itmask_metal[i] != MASK_FALSE) )
           { de[i] = de[i]
                + CII[i]/(gr_float)(12.) + COII[i]/(gr_float)(28.)
                + OII[i]/(gr_float)(16.) + OHII[i]/(gr_float)(17.)
                + H2OII[i]/(gr_float)(18.) + H3OII[i]/(gr_float)(19.)
                + O2II[i]/(gr_float)(32.); }
    }
  }
}

/// @brief overwrites species density values with values in @p species_tmpdens
///
/// For context, @p species_tmpdens holds the values from the end of the current
/// (sub-)cycle.
inline void update_from_tmpdens(
  const double* dtit, IndexRange idx_range, double* dedot_prev,
  double* HIdot_prev, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, const chemistry_data* my_chemistry,
  SpeciesMultiView<gr_float> spdens, View<const gr_float**> tmp_spdens,
  const PartMap<SpKind>& sp_kind_map
) {

  // Get pointers to the various species fields
  // ------------------------------------------
  // -> TODO: stop doing this!
  gr_float* de = spdens.contig1d_ptr(SpLUT::e);
  gr_float* HI = spdens.contig1d_ptr(SpLUT::HI);
  gr_float* HII = spdens.contig1d_ptr(SpLUT::HII);
  gr_float* HeI = spdens.contig1d_ptr(SpLUT::HeI);
  gr_float* HeII = spdens.contig1d_ptr(SpLUT::HeII);
  gr_float* HeIII = spdens.contig1d_ptr(SpLUT::HeIII);
  gr_float* HM = spdens.contig1d_ptr(SpLUT::HM);
  gr_float* H2I = spdens.contig1d_ptr(SpLUT::H2I);
  gr_float* H2II = spdens.contig1d_ptr(SpLUT::H2II);
  gr_float* DI = spdens.contig1d_ptr(SpLUT::DI);
  gr_float* DII = spdens.contig1d_ptr(SpLUT::DII);
  gr_float* HDI = spdens.contig1d_ptr(SpLUT::HDI);
  gr_float* DM = spdens.contig1d_ptr(SpLUT::DM);
  gr_float* HDII = spdens.contig1d_ptr(SpLUT::HDII);
  gr_float* HeHII = spdens.contig1d_ptr(SpLUT::HeHII);

  // Part 1: update buffers for each primordial species other than electrons
  //         that are enabled with primordial_chemistry >= 1
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      HIdot_prev[i] = std::fabs(HI[i]-tmp_spdens(SpLUT::HI, i)) /
              std::fmax((double)(dtit[i] ), tiny8);
      HI[i]    = std::fmax((gr_float)(tmp_spdens(SpLUT::HI, i) ), tiny_fortran_val);
      HII[i]   = std::fmax((gr_float)(tmp_spdens(SpLUT::HII, i) ), tiny_fortran_val);
      HeI[i]   = std::fmax((gr_float)(tmp_spdens(SpLUT::HeI, i) ), tiny_fortran_val);
      HeII[i]  = std::fmax((gr_float)(tmp_spdens(SpLUT::HeII, i) ), tiny_fortran_val);
      HeIII[i] = std::fmax((gr_float)(tmp_spdens(SpLUT::HeIII, i) ), (gr_float)(1e-5)*tiny_fortran_val);

      // this is here for historical reasons:
      // - with the way std::fmax works, I don't think it's impossible possible
      //   for HI[i] to be a NaN
      // - todo: either remove this or check if the value in the temporary
      //   buffer is a NaN
      if (std::isnan(HI[i]))  {
        OMP_PRAGMA_CRITICAL
        {
          std::printf("HUGE HI! ::  %d %d %d %g\n",
                      i, idx_range.j, idx_range.k, HI[i]);
        }
      }
    }
  }

  // Part 2: update the buffer for free electrons
  //
  // TODO: give some consideration to the case where all species have zero
  // density except for
  //   - species with neutral charges
  //   - Hminus & Dminus (H⁻ & D⁻).
  // Under the currrent charge conservation logic, the resulting electron
  // density would be negative in this scenario
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      // temporarily store the electron density from the start of the current
      // subcycle.
      dedot_prev[i] = de[i];
    }
  }

  // actually override the electron density
  update_electron_densities(idx_range, itmask, itmask_metal, my_chemistry,
                            spdens);

  // update dedot_prev to store the time derivative of the electron density
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      dedot_prev[i] = std::fabs(de[i]-dedot_prev[i])/
           std::fmax(dtit[i],tiny8);
    }
  }

  // Part 3: update buffers for all other primordial species
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {

      if (my_chemistry->primordial_chemistry > 1)  {
        HM[i]    = std::fmax((gr_float)(tmp_spdens(SpLUT::HM, i) ), tiny_fortran_val);
        H2I[i]   = std::fmax((gr_float)(tmp_spdens(SpLUT::H2I, i)), tiny_fortran_val);
        H2II[i]  = std::fmax((gr_float)(tmp_spdens(SpLUT::H2II, i) ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 2)  {
        DI[i]    = std::fmax((gr_float)(tmp_spdens(SpLUT::DI, i) ), tiny_fortran_val);
        DII[i]   = std::fmax((gr_float)(tmp_spdens(SpLUT::DII, i) ), tiny_fortran_val);
        HDI[i]   = std::fmax((gr_float)(tmp_spdens(SpLUT::HDI, i) ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 3)  {
        DM[i]    = std::fmax((gr_float)(tmp_spdens(SpLUT::DM, i) ), tiny_fortran_val);
        HDII[i]  = std::fmax((gr_float)(tmp_spdens(SpLUT::HDII, i) ), tiny_fortran_val);
        HeHII[i] = std::fmax((gr_float)(tmp_spdens(SpLUT::HeHII, i) ), tiny_fortran_val);
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
      gr_float* ptr = spdens.contig1d_ptr(sp_id);
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        // note: I'm pretty sure itmask_metal is false if itmask is false
        if ((itmask[i] != MASK_FALSE) && (itmask_metal[i] != MASK_FALSE)) {
          ptr[i] = std::fmax((gr_float)(tmp_spdens(sp_id, i)), tiny_fortran_val);
        }
      }
    }
  }
}

}  // namespace step_gs_detail

/// Uses one linearly implicit Gauss-Seidel sweep of a backward-Euler time
/// integrator to advance the rate equations by one (sub-)cycle (dtit).
inline void step_rate_gauss_seidel(
  const double* dtit, IndexRange idx_range, gr_mask_type anydust,
  const double* rhoH, double* dedot_prev,
  double* HIdot_prev, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  grackle_field_data* my_fields, SpeciesMultiView<gr_float> sp_densities,
  SpeciesCollection species_tmpdens, const FullRxnRateBuf rxn_rate_buf,
  const PartMap<SpKind>& sp_kind_map
) {
  View<gr_float**> species_tmpdens_view(species_tmpdens.data[0],
                                        MAX_EVOLVED_SPECIES_FIELDS,
                                        my_fields->grid_dimension[0]);

  // perform the Gauss-Seidel sweep to compute the species densities at the
  // end of the current timestep. The results are saved in species_tmpdens
  chemistry::species_density_updates_gauss_seidel(
      species_tmpdens_view, idx_range, dtit, anydust, rhoH, itmask,
      itmask_metal, my_chemistry, sp_densities, rxn_rate_buf, sp_kind_map);

  // update the species densities with the values in species_tmpdens
  step_gs_detail::update_from_tmpdens(
      dtit, idx_range, dedot_prev, HIdot_prev, itmask, itmask_metal,
      my_chemistry, sp_densities, species_tmpdens_view, sp_kind_map);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // STEP_RATE_GAUSS_SEIDEL_HPP
