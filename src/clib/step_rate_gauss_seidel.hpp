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

namespace grackle::impl {


/// This function overwrites the species density field values using the values
/// in @p species_tmpdens, which holds the values from the end of the current
/// (sub-)cycle.
inline void update_fields_from_tmpdens_gauss_seidel(
  const double* dtit, IndexRange idx_range, double* dedot_prev,
  double* HIdot_prev, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  grackle_field_data* my_fields,
  grackle::impl::SpeciesCollection species_tmpdens
) {

  // Construct views of various species fields
  // -> TODO: stop doing this!

  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(my_fields->HeI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(my_fields->HeII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(my_fields->HeIII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(my_fields->HM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DI(my_fields->DI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(my_fields->DII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(my_fields->HDI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DM(my_fields->DM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(my_fields->HDII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(my_fields->HeHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(my_fields->CI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(my_fields->CII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(my_fields->CO_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO2(my_fields->CO2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(my_fields->OI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(my_fields->OH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(my_fields->H2O_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2(my_fields->O2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiI(my_fields->SiI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiOI(my_fields->SiOI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2I(my_fields->SiO2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH(my_fields->CH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH2(my_fields->CH2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> COII(my_fields->COII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OII(my_fields->OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OHII(my_fields->OHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2OII(my_fields->H2OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H3OII(my_fields->H3OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2II(my_fields->O2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg(my_fields->Mg_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al(my_fields->Al_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> S(my_fields->S_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe(my_fields->Fe_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(my_fields->SiM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(my_fields->FeM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(my_fields->AC_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(my_fields->SiO2_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(my_fields->MgO_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(my_fields->FeS_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(my_fields->Al2O3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(my_fields->ref_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(my_fields->vol_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // --- (E) Set densities from 1D temps to 3D fields ---

  const int j = idx_range.j;
  const int k = idx_range.k;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      HIdot_prev[i] = std::fabs(HI(i,j,k)-species_tmpdens.data[SpLUT::HI][i]) /
              std::fmax((double)(dtit[i] ), tiny8);
      HI(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HI][i] ), tiny_fortran_val);
      HII(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HII][i] ), tiny_fortran_val);
      HeI(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeI][i] ), tiny_fortran_val);
      HeII(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeII][i] ), tiny_fortran_val);
      HeIII(i,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeIII][i] ), (gr_float)(1e-5)*tiny_fortran_val);

      // temporarily store the electron density from the start of the current
      // subcycle.
      dedot_prev[i] = de(i,j,k);

      // Use charge conservation to determine electron fraction
      // -> in other words, we ignore species_tmpdens.data[SpLUT::e] (in
      //    practice, I think that array holds garbage values)
      de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/(gr_float)(4.) +
           HeIII(i,j,k)/(gr_float)(2.);
      if (my_chemistry->primordial_chemistry > 1)
           { de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/(gr_float)(2.); }

      if (my_chemistry->primordial_chemistry > 2)
           { de(i,j,k) = de(i,j,k) + DII(i,j,k)/(gr_float)(2.); }
      if (my_chemistry->primordial_chemistry > 3)
           { de(i,j,k) = de(i,j,k) - DM(i,j,k)/(gr_float)(2.)
                + HDII(i,j,k)/(gr_float)(3.) + HeHII(i,j,k)/(gr_float)(5.); }
      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i] != MASK_FALSE) )
           { de(i,j,k) = de(i,j,k)
                + CII(i,j,k)/(gr_float)(12.) + COII(i,j,k)/(gr_float)(28.)
                + OII(i,j,k)/(gr_float)(16.) + OHII(i,j,k)/(gr_float)(17.)
                + H2OII(i,j,k)/(gr_float)(18.) + H3OII(i,j,k)/(gr_float)(19.)
                + O2II(i,j,k)/(gr_float)(32.); }

      // store the time-derivative of the electron-density in dedot
      // (don't forget that we previously stored the value from the start of
      // the current cycle within dedot_prev)
      dedot_prev[i] = std::fabs(de(i,j,k)-dedot_prev[i])/
           std::fmax(dtit[i],tiny8);

      if (my_chemistry->primordial_chemistry > 1)  {
        HM(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HM][i] ), tiny_fortran_val);
        H2I(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2I][i]), tiny_fortran_val);
        H2II(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2II][i] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 2)  {
        DI(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DI][i] ), tiny_fortran_val);
        DII(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DII][i] ), tiny_fortran_val);
        HDI(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDI][i] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 3)  {
        DM(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DM][i] ), tiny_fortran_val);
        HDII(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDII][i] ), tiny_fortran_val);
        HeHII(i,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeHII][i] ), tiny_fortran_val);
      }

      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i] != MASK_FALSE) )  {
        CI(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CI][i]      ), tiny_fortran_val);
        CII(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CII][i]     ), tiny_fortran_val);
        CO(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CO][i]      ), tiny_fortran_val);
        CO2(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CO2][i]     ), tiny_fortran_val);
        OI(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OI][i]      ), tiny_fortran_val);
        OH(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OH][i]      ), tiny_fortran_val);
        H2O(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2O][i]     ), tiny_fortran_val);
        O2(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2][i]      ), tiny_fortran_val);
        SiI(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiI][i]     ), tiny_fortran_val);
        SiOI(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiOI][i]    ), tiny_fortran_val);
        SiO2I(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2I][i]   ), tiny_fortran_val);
        CH(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CH][i]      ), tiny_fortran_val);
        CH2(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CH2][i]     ), tiny_fortran_val);
        COII(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::COII][i]    ), tiny_fortran_val);
        OII(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OII][i]     ), tiny_fortran_val);
        OHII(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OHII][i]    ), tiny_fortran_val);
        H2OII(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2OII][i]   ), tiny_fortran_val);
        H3OII(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H3OII][i]   ), tiny_fortran_val);
        O2II(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2II][i]    ), tiny_fortran_val);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            Mg(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg][i]      ), tiny_fortran_val);
          }
          if (my_chemistry->dust_species > 1)  {
            Al(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al][i]      ), tiny_fortran_val);
            S(i,j,k)       = std::fmax((gr_float)(species_tmpdens.data[SpLUT::S][i]       ), tiny_fortran_val);
            Fe(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe][i]      ), tiny_fortran_val);
          }
        }
      }

      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          MgSiO3(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgSiO3_dust][i]  ), tiny_fortran_val);
          AC(i,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::AC_dust][i]      ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 1)  {
          SiM(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiM_dust][i]     ), tiny_fortran_val);
          FeM(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeM_dust][i]     ), tiny_fortran_val);
          Mg2SiO4(i,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg2SiO4_dust][i] ), tiny_fortran_val);
          Fe3O4(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe3O4_dust][i]   ), tiny_fortran_val);
          SiO2D(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2_dust][i]   ), tiny_fortran_val);
          MgO(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgO_dust][i]     ), tiny_fortran_val);
          FeS(i,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeS_dust][i]     ), tiny_fortran_val);
          Al2O3(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al2O3_dust][i]   ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 2)  {
          reforg(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::ref_org_dust][i]   ), tiny_fortran_val);
          volorg(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::vol_org_dust][i]   ), tiny_fortran_val);
          H2Oice(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2O_ice_dust][i]   ), tiny_fortran_val);
        }
      }

    }
    // 

    if (HI(i,j,k) != HI(i,j,k))  {
      OMP_PRAGMA_CRITICAL
      {
        std::printf("HUGE HI! ::  %d %d %d %g\n",
                    i, j, k, HI ( i, j, k ));
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
  grackle::impl::PhotoRxnRateCollection kshield_buf
) {

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
