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
#include "fortran_func_decls.h"  // gr_mask_type
#include "index_helper.h"
#include "internal_types.hpp"
#include "LUT.hpp"
#include "utils-field.hpp"

namespace grackle::impl::gauss_seidel {

/// groups some arguments together that are used with update_densities
struct DensityUpdateArgPack {
  SpeciesLUTFieldAdaptor field_adaptor;
  IndexRange idx_range;
  const gr_mask_type* itmask;
  grackle::impl::SpeciesCollection species_tmpdens;
};

/// update the specified field entry with values from species_tmpdens
///
/// this is a template function because we can have machinery to let us use
/// the SpeciesLUTFieldAdaptor machinery without incurring a runtime cost
///
/// @note
/// This is mostly intended to serve as a short-term placeholder to reduce a
/// bunch of boilerplate code while we refactor other parts of Grackle. We will
/// probably eventually remove this code
template <int lut_idx>
void update_densities(DensityUpdateArgPack& pack, gr_float floor_val) {
  grackle::impl::View<gr_float***> view(
      pack.field_adaptor.get_ptr_static<lut_idx>(),
      pack.field_adaptor.data.grid_dimension[0],
      pack.field_adaptor.data.grid_dimension[1],
      pack.field_adaptor.data.grid_dimension[2]);
  const double* src_vals = pack.species_tmpdens.data[lut_idx];

  for (int i = pack.idx_range.i_start; i < pack.idx_range.i_stop; i++) {
    if (pack.itmask[i] != MASK_FALSE) {
      view(i, pack.idx_range.j, pack.idx_range.k) =
          std::fmax((gr_float)(src_vals[i]), floor_val);
    }
  }
}

/// update my_fields->e_density using charge conservation and write the time
/// derivative to @p dedot_prev
inline void update_electron_densities(const double* dtit, IndexRange idx_range,
                                      double* dedot_prev,
                                      const gr_mask_type* itmask,
                                      const gr_mask_type* itmask_metal,
                                      const chemistry_data* my_chemistry,
                                      grackle_field_data* my_fields) {
  grackle::impl::View<gr_float***> de(
      my_fields->e_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(
      my_fields->HeII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(
      my_fields->HeIII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(
      my_fields->HM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(
      my_fields->DII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DM(
      my_fields->DM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(
      my_fields->HDII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(
      my_fields->HeHII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(
      my_fields->CII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> COII(
      my_fields->COII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OII(
      my_fields->OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OHII(
      my_fields->OHII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2OII(
      my_fields->H2OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H3OII(
      my_fields->H3OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2II(
      my_fields->O2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const int j = idx_range.j;
  const int k = idx_range.k;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      // temporarily store electron density from start of the current subcycle
      dedot_prev[i] = de(i, j, k);

      de(i, j, k) = HII(i, j, k) + HeII(i, j, k) / (gr_float)(4.) +
                    HeIII(i, j, k) / (gr_float)(2.);
      if (my_chemistry->primordial_chemistry > 1) {
        de(i, j, k) =
            de(i, j, k) - HM(i, j, k) + H2II(i, j, k) / (gr_float)(2.);
      }

      if (my_chemistry->primordial_chemistry > 2) {
        de(i, j, k) = de(i, j, k) + DII(i, j, k) / (gr_float)(2.);
      }
      if (my_chemistry->primordial_chemistry > 3) {
        de(i, j, k) = de(i, j, k) - DM(i, j, k) / (gr_float)(2.) +
                      HDII(i, j, k) / (gr_float)(3.) +
                      HeHII(i, j, k) / (gr_float)(5.);
      }
      if ((my_chemistry->metal_chemistry == 1) &&
          (itmask_metal[i] != MASK_FALSE)) {
        de(i, j, k) =
            de(i, j, k) + CII(i, j, k) / (gr_float)(12.) +
            COII(i, j, k) / (gr_float)(28.) + OII(i, j, k) / (gr_float)(16.) +
            OHII(i, j, k) / (gr_float)(17.) + H2OII(i, j, k) / (gr_float)(18.) +
            H3OII(i, j, k) / (gr_float)(19.) + O2II(i, j, k) / (gr_float)(32.);
      }

      // store the time-derivative of the electron-density in dedot
      // (don't forget that we previously stored the value from the start of
      // the current cycle within dedot_prev)
      dedot_prev[i] =
          std::fabs(de(i, j, k) - dedot_prev[i]) / std::fmax(dtit[i], tiny8);
    }
  }
}

}  // namespace grackle::impl::gauss_seidel

namespace grackle::impl {

/// This function overwrites the species density field values using the values
/// in @p species_tmpdens, which holds the values from the end of the current
/// (sub-)cycle.
inline void update_fields_from_tmpdens_gauss_seidel(
    const double* dtit, IndexRange idx_range, double* dedot_prev,
    double* HIdot_prev, const gr_mask_type* itmask,
    const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    grackle::impl::SpeciesCollection species_tmpdens) {
  grackle::impl::View<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const int j = idx_range.j;
  const int k = idx_range.k;

  // handle the primordial species
  {
    // initialize the pack variable to group a set of arguments that gets
    // repeatedly passed to a helper function
    gauss_seidel::DensityUpdateArgPack pack{SpeciesLUTFieldAdaptor{*my_fields},
                                            idx_range, itmask, species_tmpdens};

    // record the time derivative of neutral Hydrogen
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        HIdot_prev[i] =
            std::fabs(HI(i, j, k) - species_tmpdens.data[SpLUT::HI][i]) /
            std::fmax((double)(dtit[i]), tiny8);
      }
    }

    // primordial_chemistry == 1 updates
    gauss_seidel::update_densities<SpLUT::HI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::HII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::HeI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::HeII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::HeIII>(
        pack, (gr_float)(1e-5) * tiny_fortran_val);

    // Use charge conservation to determine electron fraction
    // -> in other words, we ignore species_tmpdens.data[SpLUT::e] (in
    //    practice, I think that array holds garbage values)
    // -> this uses the current value of my_fields->e_density and the updated
    //    value to actually compute the time derivative of my_fields->e_density
    //    (the result is written to dedot_prev)
    // -> IMPORTANT: by calling this function here, we are computing electron
    //    densities using only a subset of updated species densities
    //    (this is probably ok since hydrogen and helium are **SO** abundant).
    //    This is the historical behavior.
    gauss_seidel::update_electron_densities(dtit, idx_range, dedot_prev, itmask,
                                            itmask_metal, my_chemistry,
                                            my_fields);

    if (my_chemistry->primordial_chemistry > 1) {
      gauss_seidel::update_densities<SpLUT::HM>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::H2I>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::H2II>(pack, tiny_fortran_val);
    }

    if (my_chemistry->primordial_chemistry > 2) {
      gauss_seidel::update_densities<SpLUT::DI>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::DII>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::HDI>(pack, tiny_fortran_val);
    }

    if (my_chemistry->primordial_chemistry > 3) {
      gauss_seidel::update_densities<SpLUT::DM>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::HDII>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::HeHII>(pack, tiny_fortran_val);
    }
  }

  // handle metal species
  if (my_chemistry->metal_chemistry == 1) {
    // initialize the pack variable to group a set of arguments that gets
    // repeatedly passed to a helper function
    // -> we explicitly use itmask_metal instead of itmask
    gauss_seidel::DensityUpdateArgPack pack{SpeciesLUTFieldAdaptor{*my_fields},
                                            idx_range, itmask_metal,
                                            species_tmpdens};

    gauss_seidel::update_densities<SpLUT::CI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::CII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::CO>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::CO2>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::OI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::OH>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::H2O>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::O2>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::SiI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::SiOI>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::SiO2I>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::CH>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::CH2>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::COII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::OII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::OHII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::H2OII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::H3OII>(pack, tiny_fortran_val);
    gauss_seidel::update_densities<SpLUT::O2II>(pack, tiny_fortran_val);
    if ((my_chemistry->grain_growth == 1) ||
        (my_chemistry->dust_sublimation == 1)) {
      if (my_chemistry->dust_species > 0) {
        gauss_seidel::update_densities<SpLUT::Mg>(pack, tiny_fortran_val);
      }
      if (my_chemistry->dust_species > 1) {
        gauss_seidel::update_densities<SpLUT::Al>(pack, tiny_fortran_val);
        gauss_seidel::update_densities<SpLUT::S>(pack, tiny_fortran_val);
        gauss_seidel::update_densities<SpLUT::Fe>(pack, tiny_fortran_val);
      }
    }
  }

  // handle dust grain species

  if ((my_chemistry->grain_growth == 1) ||
      (my_chemistry->dust_sublimation == 1)) {
    // initialize the pack variable to group a set of arguments that gets
    // repeatedly passed to a helper function
    // -> we explicitly use itmask_metal instead of itmask
    gauss_seidel::DensityUpdateArgPack pack{SpeciesLUTFieldAdaptor{*my_fields},
                                            idx_range, itmask_metal,
                                            species_tmpdens};

    if (my_chemistry->dust_species > 0) {
      gauss_seidel::update_densities<SpLUT::MgSiO3_dust>(pack,
                                                         tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::AC_dust>(pack, tiny_fortran_val);
    }
    if (my_chemistry->dust_species > 1) {
      gauss_seidel::update_densities<SpLUT::SiM_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::FeM_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::Mg2SiO4_dust>(pack,
                                                          tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::Fe3O4_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::SiO2_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::MgO_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::FeS_dust>(pack, tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::Al2O3_dust>(pack, tiny_fortran_val);
    }
    if (my_chemistry->dust_species > 2) {
      gauss_seidel::update_densities<SpLUT::ref_org_dust>(pack,
                                                          tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::vol_org_dust>(pack,
                                                          tiny_fortran_val);
      gauss_seidel::update_densities<SpLUT::H2O_ice_dust>(pack,
                                                          tiny_fortran_val);
    }
  }

  // warn users about irregular HI densities
  // -> this is unaltered from the original version of the code
  // -> it appears that it was intended to warn users when there was a large
  //    value, but it only triggers when HI is a NaN (it won't trigger when its
  //    an inf)
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (HI(i, j, k) != HI(i, j, k)) {
      OMP_PRAGMA_CRITICAL {
        std::printf("HUGE HI! ::  %d %d %d %g\n", i, j, k, HI(i, j, k));
      }
    }
  }
}

/// Uses one linearly implicit Gauss-Seidel sweep of a backward-Euler time
/// integrator to advance the rate equations by one (sub-)cycle (dtit).
inline void step_rate_gauss_seidel(
    const double* dtit, IndexRange idx_range, gr_mask_type anydust,
    const double* h2dust, const double* rhoH, double* dedot_prev,
    double* HIdot_prev, const gr_mask_type* itmask,
    const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
    grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
    grackle::impl::GrainSpeciesCollection grain_growth_rates,
    grackle::impl::SpeciesCollection species_tmpdens,
    grackle::impl::CollisionalRxnRateCollection kcol_buf,
    grackle::impl::PhotoRxnRateCollection kshield_buf) {
  // perform the Gauss-Seidel sweep to compute the species densities at the
  // end of the current timestep. The results are saved in species_tmpdens
  grackle::impl::chemistry::species_density_updates_gauss_seidel(
      species_tmpdens, idx_range, dtit, anydust, h2dust, rhoH, itmask,
      itmask_metal, my_chemistry, my_fields, my_uvb_rates, grain_growth_rates,
      kcol_buf, kshield_buf);

  // update the entries from my_fields with the values in species_tmpdens
  update_fields_from_tmpdens_gauss_seidel(
      dtit, idx_range, dedot_prev, HIdot_prev, itmask, itmask_metal,
      my_chemistry, my_fields, species_tmpdens);
}

}  // namespace grackle::impl

#endif /* STEP_RATE_GAUSS_SEIDEL_HPP */
