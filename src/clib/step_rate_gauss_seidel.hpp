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

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_type
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "LUT.hpp"
#include "utils-cpp.hpp"

namespace grackle::impl {

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
)
{

  // Construct views of various species fields
  // -----------------------------------------

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
  // -- removed line (previously just declared arg types) -- 
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

  // Radiation Fields
  grackle::impl::View<gr_float***> kphHI(my_fields->RT_HI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeI(my_fields->RT_HeI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeII(my_fields->RT_HeII_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissHDI(my_fields->RT_HDI_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphCI(my_fields->RT_CI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphOI(my_fields->RT_OI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissCO(my_fields->RT_CO_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissOH(my_fields->RT_OH_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissH2O(my_fields->RT_H2O_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // locals

  int i;
  double scoef, acoef;

  const int j = idx_range.j;
  const int k = idx_range.k;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {

        // 1) HI

        scoef  = kcol_buf.data[CollisionalRxnLUT::k2][i-1]*HII(i-1,j,k)*de(i-1,j,k);
        acoef  = kcol_buf.data[CollisionalRxnLUT::k1][i-1]*de(i-1,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k57][i-1]*HI(i-1,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k58][i-1]*HeI(i-1,j,k)/4.
               + kshield_buf.k24[i-1];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i-1,j,k); }
        species_tmpdens.data[SpLUT::HI][i-1]  = (scoef*dtit[i-1] + HI(i-1,j,k))/
             (1. + acoef*dtit[i-1]);
        if (species_tmpdens.data[SpLUT::HI][i-1] != species_tmpdens.data[SpLUT::HI][i-1])  {
          OMP_PRAGMA_CRITICAL
          {
            printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g %g %g\n",
                   i,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HI] [ i-1 ],
                   HI ( i-1, j, k ),
                   HII ( i-1, j, k ),
                   de ( i-1, j, k ),
                   kphHI ( i-1, j, k ),
                   scoef,
                   acoef,
                   dtit [ i-1 ]);
          }
          // ERROR_MESSAGE
        }

        // 2) HII
        scoef  = kcol_buf.data[CollisionalRxnLUT::k1][i-1]*species_tmpdens.data[SpLUT::HI][i-1]*de(i-1,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k57][i-1]*species_tmpdens.data[SpLUT::HI][i-1]*species_tmpdens.data[SpLUT::HI][i-1]
               + kcol_buf.data[CollisionalRxnLUT::k58][i-1]*species_tmpdens.data[SpLUT::HI][i-1]*HeI(i-1,j,k)/4.
               + kshield_buf.k24[i-1]*species_tmpdens.data[SpLUT::HI][i-1];
        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i-1,j,k)*species_tmpdens.data[SpLUT::HI][i-1]; }
        acoef  = kcol_buf.data[CollisionalRxnLUT::k2][i-1]*de (i-1,j,k);
        species_tmpdens.data[SpLUT::HII][i-1] = (scoef*dtit[i-1] + HII(i-1,j,k))/
             (1. +acoef*dtit[i-1]);
        // 
        if (species_tmpdens.data[SpLUT::HII][i-1] <= 0.)   //#####
        {
          OMP_PRAGMA_CRITICAL
          {
            printf("negative HIIp! ::  %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
                   i,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HII] [ i-1 ],
                   scoef,
                   dtit [ i-1 ],
                   HII ( i-1, j, k ),
                   acoef,
                   kcol_buf.data[CollisionalRxnLUT::k2] [ i-1 ],
                   de ( i-1, j, k ),
                   kphHI ( i-1, j, k ),
                   species_tmpdens.data[SpLUT::HI] [ i-1 ],
                   kshield_buf.k24 [ i-1 ]);
          }
        }

        // 3) Electron density

        scoef = 0.
                   + kcol_buf.data[CollisionalRxnLUT::k57][i-1]*species_tmpdens.data[SpLUT::HI][i-1]*species_tmpdens.data[SpLUT::HI][i-1]
                   + kcol_buf.data[CollisionalRxnLUT::k58][i-1]*species_tmpdens.data[SpLUT::HI][i-1]*HeI(i-1,j,k)/4.
                   + kshield_buf.k24[i-1]*HI(i-1,j,k)
                   + kshield_buf.k25[i-1]*HeII(i-1,j,k)/4.
                   + kshield_buf.k26[i-1]*HeI(i-1,j,k)/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i-1,j,k) * HI(i-1,j,k)
                  + kphHeI(i-1,j,k)  * HeI(i-1,j,k)  / 4.
                  + kphHeII(i-1,j,k) * HeII(i-1,j,k) / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i-1,j,k) * HI(i-1,j,k); }



        acoef = -(kcol_buf.data[CollisionalRxnLUT::k1][i-1]*HI(i-1,j,k)      - kcol_buf.data[CollisionalRxnLUT::k2][i-1]*HII(i-1,j,k)
                + kcol_buf.data[CollisionalRxnLUT::k3][i-1]*HeI(i-1,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k6][i-1]*HeIII(i-1,j,k)/4.
                + kcol_buf.data[CollisionalRxnLUT::k5][i-1]*HeII(i-1,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k4][i-1]*HeII(i-1,j,k)/4.);
        species_tmpdens.data[SpLUT::e][i-1]   = (scoef*dtit[i-1] + de(i-1,j,k))
                       / (1. + acoef*dtit[i-1]);

      }
    }

  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if (itmask[i-1] != MASK_FALSE)  {

      // 4) HeI

      scoef  = kcol_buf.data[CollisionalRxnLUT::k4][i-1]*HeII(i-1,j,k)*de(i-1,j,k);
      acoef  = kcol_buf.data[CollisionalRxnLUT::k3][i-1]*de(i-1,j,k)
                   + kshield_buf.k26[i-1];

      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeI(i-1,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        scoef = scoef +  4. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::k152][i-1] * HeHII(i-1,j,k) *    HI(i-1,j,k) /  5.
            + kcol_buf.data[CollisionalRxnLUT::k153][i-1] * HeHII(i-1,j,k) *    de(i-1,j,k) /  5.
            );
        acoef = acoef
            + kcol_buf.data[CollisionalRxnLUT::k148][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k149][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k150][i-1] *  H2II(i-1,j,k) /  2.;
      }
      species_tmpdens.data[SpLUT::HeI][i-1]   = ( scoef*dtit[i-1] + HeI(i-1,j,k) )
                 / ( 1. + acoef*dtit[i-1] );

      // 5) HeII

      scoef  = kcol_buf.data[CollisionalRxnLUT::k3][i-1]*species_tmpdens.data[SpLUT::HeI][i-1]*de(i-1,j,k)
             + kcol_buf.data[CollisionalRxnLUT::k6][i-1]*HeIII(i-1,j,k)*de(i-1,j,k)
             + kshield_buf.k26[i-1]*species_tmpdens.data[SpLUT::HeI][i-1];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeI(i-1,j,k)*species_tmpdens.data[SpLUT::HeI][i-1]; }

      acoef  = kcol_buf.data[CollisionalRxnLUT::k4][i-1]*de(i-1,j,k) + kcol_buf.data[CollisionalRxnLUT::k5][i-1]*de(i-1,j,k)
             + kshield_buf.k25[i-1];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeII(i-1,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        acoef = acoef
            + kcol_buf.data[CollisionalRxnLUT::k151][i-1] *    HI(i-1,j,k);
      }
      species_tmpdens.data[SpLUT::HeII][i-1]  = ( scoef*dtit[i-1] + HeII(i-1,j,k) )
                 / ( 1. + acoef*dtit[i-1] );

      // 6) HeIII

      scoef   = kcol_buf.data[CollisionalRxnLUT::k5][i-1]*species_tmpdens.data[SpLUT::HeII][i-1]*de(i-1,j,k)
              + kshield_buf.k25[i-1]*species_tmpdens.data[SpLUT::HeII][i-1];
      if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeII(i-1,j,k) * species_tmpdens.data[SpLUT::HeII][i-1]; }
      acoef   = kcol_buf.data[CollisionalRxnLUT::k6][i-1]*de(i-1,j,k);
      species_tmpdens.data[SpLUT::HeIII][i-1]  = ( scoef*dtit[i-1] + HeIII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

    }
  }

  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {

        // 1) HI
        scoef  =      kcol_buf.data[CollisionalRxnLUT::k2][i-1] * HII(i-1,j,k) * de(i-1,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k13][i-1]* HI(i-1,j,k)  * H2I(i-1,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k11][i-1]* HII(i-1,j,k) * H2I(i-1,j,k)/2.
               + 2.*kcol_buf.data[CollisionalRxnLUT::k12][i-1]* de(i-1,j,k)  * H2I(i-1,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k14][i-1]* HM(i-1,j,k)  * de(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k15][i-1]* HM(i-1,j,k)  * HI(i-1,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k16][i-1]* HM(i-1,j,k)  * HII(i-1,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k18][i-1]* H2II(i-1,j,k)* de(i-1,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k19][i-1]* H2II(i-1,j,k)* HM(i-1,j,k)/2.
               + 2.*kshield_buf.k31[i-1]   * H2I(i-1,j,k)/2.;

        acoef  =      kcol_buf.data[CollisionalRxnLUT::k1][i-1] * de(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k7][i-1] * de(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k8][i-1] * HM(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k9][i-1] * HII(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k10][i-1]* H2II(i-1,j,k)/2.
               + 2.*kcol_buf.data[CollisionalRxnLUT::k22][i-1]* std::pow(HI(i-1,j,k),2)
               +      kcol_buf.data[CollisionalRxnLUT::k57][i-1]* HI(i-1,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k58][i-1]* HeI(i-1,j,k)/4.
               + kshield_buf.k24[i-1];

        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i-1,j,k); }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if ((my_chemistry->primordial_chemistry > 2) && (my_chemistry->radiative_transfer_HDI_dissociation > 0))  {
            scoef = scoef
              + kdissHDI(i-1,j,k) * HDI(i-1,j,k)/3.0;
          }
          if ( (my_chemistry->metal_chemistry == 1)  && 
               (itmask_metal[i-1] != MASK_FALSE) )  {
            if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
              scoef = scoef
                + kdissOH (i-1,j,k) * OH(i-1,j,k) /17.0
                + kdissH2O(i-1,j,k) * H2O(i-1,j,k)/18.0;
            }
          }
        }

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i-1] != MASK_FALSE)  {
            acoef = acoef + 2. * h2dust[i-1] * rhoH[i-1];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k50][i-1] * HII(i-1,j,k) * DI(i-1,j,k)  / 2.
                + kcol_buf.data[CollisionalRxnLUT::k54][i-1] * H2I(i-1,j,k) * DI(i-1,j,k)  / 4.;
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k51][i-1] * DII(i-1,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k55][i-1] * HDI(i-1,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k131][i-1] *  HDII(i-1,j,k) *    de(i-1,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k134][i-1] *   HII(i-1,j,k) *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k135][i-1] *    HM(i-1,j,k) *    DI(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k150][i-1] *   HeI(i-1,j,k) *  H2II(i-1,j,k) /  8.
              + kcol_buf.data[CollisionalRxnLUT::k153][i-1] * HeHII(i-1,j,k) *    de(i-1,j,k) /  5.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k125][i-1] *  HDII(i-1,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k130][i-1] *   DII(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k136][i-1] *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k137][i-1] *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k151][i-1] *  HeII(i-1,j,k) /  4.
              + kcol_buf.data[CollisionalRxnLUT::k152][i-1] * HeHII(i-1,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i-1] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::kz20][i-1] *    CI(i-1,j,k) *   H2I(i-1,j,k) / 24.
              + kcol_buf.data[CollisionalRxnLUT::kz21][i-1] *    OI(i-1,j,k) *   H2I(i-1,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz22][i-1] *   HII(i-1,j,k) *    OI(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz23][i-1] *   H2I(i-1,j,k) *    CH(i-1,j,k) / 26.
              + kcol_buf.data[CollisionalRxnLUT::kz24][i-1] *   H2I(i-1,j,k) *    OH(i-1,j,k) / 34.
              + kcol_buf.data[CollisionalRxnLUT::kz26][i-1] *    OH(i-1,j,k) *    CO(i-1,j,k) / 476.
              + kcol_buf.data[CollisionalRxnLUT::kz28][i-1] *    CI(i-1,j,k) *    OH(i-1,j,k) / 204.
              + kcol_buf.data[CollisionalRxnLUT::kz32][i-1] *    OI(i-1,j,k) *    CH(i-1,j,k) / 208.
              + kcol_buf.data[CollisionalRxnLUT::kz33][i-1] *    OI(i-1,j,k) *    OH(i-1,j,k) / 272.
              + kcol_buf.data[CollisionalRxnLUT::kz34][i-1] *   HII(i-1,j,k) *    OH(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz35][i-1] *   HII(i-1,j,k) *   H2O(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz36][i-1] *   HII(i-1,j,k) *    O2(i-1,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz37][i-1] *   CII(i-1,j,k) *    OH(i-1,j,k) / 204.
              + kcol_buf.data[CollisionalRxnLUT::kz40][i-1] *   OII(i-1,j,k) *   H2I(i-1,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz41][i-1] *  OHII(i-1,j,k) *   H2I(i-1,j,k) / 34.
              + kcol_buf.data[CollisionalRxnLUT::kz42][i-1] * H2OII(i-1,j,k) *   H2I(i-1,j,k) / 36.
              + kcol_buf.data[CollisionalRxnLUT::kz46][i-1] * H2OII(i-1,j,k) *    de(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz48][i-1] * H3OII(i-1,j,k) *    de(i-1,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz49][i-1] * H3OII(i-1,j,k) *    de(i-1,j,k) / 9.5
              + kcol_buf.data[CollisionalRxnLUT::kz52][i-1] *   SiI(i-1,j,k) *    OH(i-1,j,k) / 476.
              + kcol_buf.data[CollisionalRxnLUT::kz54][i-1] *  SiOI(i-1,j,k) *    OH(i-1,j,k) / 748.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz15][i-1] *    CH(i-1,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz16][i-1] *   CH2(i-1,j,k) / 14.
              + kcol_buf.data[CollisionalRxnLUT::kz17][i-1] *    OH(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz18][i-1] *   H2O(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz19][i-1] *    O2(i-1,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz27][i-1] *    CI(i-1,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz30][i-1] *    OI(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz39][i-1] *   OII(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz43][i-1] *  COII(i-1,j,k) / 28.;
        }
        species_tmpdens.data[SpLUT::HI][i-1]  = ( scoef*dtit[i-1] + HI(i-1,j,k) ) /
                        ( 1.f + acoef*dtit[i-1] );
        if (species_tmpdens.data[SpLUT::HI][i-1] != species_tmpdens.data[SpLUT::HI][i-1])  {
          OMP_PRAGMA_CRITICAL
          {
            printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g\n",
                   i,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HI] [ i-1 ],
                   HI ( i-1, j, k ),
                   HII ( i-1, j, k ),
                   de ( i-1, j, k ),
                   H2I ( i-1, j, k ),
                   kphHI ( i-1, j, k ));
          }
        }

        // 2) HII

        scoef  =    kcol_buf.data[CollisionalRxnLUT::k1][i-1]  * HI(i-1,j,k) * de(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k10][i-1] * H2II(i-1,j,k)*HI(i-1,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k57][i-1] * HI(i-1,j,k) * HI(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k58][i-1] * HI(i-1,j,k) * HeI(i-1,j,k)/4.
               + kshield_buf.k24[i-1]*HI(i-1,j,k);

        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i-1,j,k) * HI(i-1,j,k); }

        acoef  =    kcol_buf.data[CollisionalRxnLUT::k2][i-1]  * de(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k9][i-1]  * HI(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k11][i-1] * H2I(i-1,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k16][i-1] * HM(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k17][i-1] * HM(i-1,j,k);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k51][i-1] * HI (i-1,j,k) * DII(i-1,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k52][i-1] * H2I(i-1,j,k) * DII(i-1,j,k) / 4.;
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k50][i-1] * DI (i-1,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k53][i-1] * HDI(i-1,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k125][i-1] *  HDII(i-1,j,k) *    HI(i-1,j,k) /  3.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k129][i-1] *    DI(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k134][i-1] *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k148][i-1] *   HeI(i-1,j,k) /  4.
              + kcol_buf.data[CollisionalRxnLUT::k149][i-1] *   HeI(i-1,j,k) /  4.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i-1] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::kz39][i-1] *   OII(i-1,j,k) *    HI(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz43][i-1] *  COII(i-1,j,k) *    HI(i-1,j,k) / 28.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz22][i-1] *    OI(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz34][i-1] *    OH(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz35][i-1] *   H2O(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz36][i-1] *    O2(i-1,j,k) / 32.;
        }
        species_tmpdens.data[SpLUT::HII][i-1]   = ( scoef*dtit[i-1] + HII(i-1,j,k) )
                        / ( 1. + acoef*dtit[i-1] );
        
        // 3) electrons:

        scoef =   kcol_buf.data[CollisionalRxnLUT::k8][i-1] * HM(i-1,j,k) * HI(i-1,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k15][i-1]* HM(i-1,j,k) * HI(i-1,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k17][i-1]* HM(i-1,j,k) * HII(i-1,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k57][i-1]* HI(i-1,j,k) * HI(i-1,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k58][i-1]* HI(i-1,j,k) * HeI(i-1,j,k)/4.
        // 
               + kshield_buf.k24[i-1]*species_tmpdens.data[SpLUT::HI][i-1]
               + kshield_buf.k25[i-1]*species_tmpdens.data[SpLUT::HeII][i-1]/4.
               + kshield_buf.k26[i-1]*species_tmpdens.data[SpLUT::HeI][i-1]/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i-1,j,k) * species_tmpdens.data[SpLUT::HI][i-1]
                  + kphHeI(i-1,j,k)  * species_tmpdens.data[SpLUT::HeI][i-1]  / 4.
                  + kphHeII(i-1,j,k) * species_tmpdens.data[SpLUT::HeII][i-1] / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i-1,j,k) * species_tmpdens.data[SpLUT::HI][i-1]; }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if ( (my_chemistry->metal_chemistry == 1)  && 
               (itmask_metal[i-1] != MASK_FALSE) )  {
            if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
              scoef = scoef
                + kphCI(i-1,j,k) * CI(i-1,j,k)/12.0
                + kphOI(i-1,j,k) * OI(i-1,j,k)/16.0;
            }
          }
        }

        acoef = - (kcol_buf.data[CollisionalRxnLUT::k1][i-1] *HI(i-1,j,k)    - kcol_buf.data[CollisionalRxnLUT::k2][i-1]*HII(i-1,j,k)
                +  kcol_buf.data[CollisionalRxnLUT::k3][i-1] *HeI(i-1,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k6][i-1]*HeIII(i-1,j,k)/4.
                +  kcol_buf.data[CollisionalRxnLUT::k5][i-1] *HeII(i-1,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k4][i-1]*HeII(i-1,j,k)/4.
                +  kcol_buf.data[CollisionalRxnLUT::k14][i-1]*HM(i-1,j,k)
                -  kcol_buf.data[CollisionalRxnLUT::k7][i-1] *HI(i-1,j,k)
                -  kcol_buf.data[CollisionalRxnLUT::k18][i-1]*H2II(i-1,j,k)/2.);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k56][i-1] * DI (i-1,j,k) * HM(i-1,j,k) / 2.;
          acoef = acoef
                - kcol_buf.data[CollisionalRxnLUT::k1] [i-1] * DI (i-1,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k2] [i-1] * DII(i-1,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k137][i-1] *    DM(i-1,j,k) *    HI(i-1,j,k) /  2.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k131][i-1] *  HDII(i-1,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k132][i-1] *    DI(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k153][i-1] * HeHII(i-1,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i-1] != MASK_FALSE) )  {
          scoef = scoef;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz44][i-1] *   CII(i-1,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz45][i-1] *   OII(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz46][i-1] * H2OII(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz47][i-1] * H2OII(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz48][i-1] * H3OII(i-1,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz49][i-1] * H3OII(i-1,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz50][i-1] *  O2II(i-1,j,k) / 32.;
        }
        species_tmpdens.data[SpLUT::e][i-1]  = ( scoef*dtit[i-1] + de(i-1,j,k) )
                  / ( 1. + acoef*dtit[i-1] );

        // 7) H2

        scoef = 2.*(kcol_buf.data[CollisionalRxnLUT::k8][i-1]  * HM(i-1,j,k)   * HI(i-1,j,k)
              +       kcol_buf.data[CollisionalRxnLUT::k10][i-1] * H2II(i-1,j,k) * HI(i-1,j,k)/2.
              +       kcol_buf.data[CollisionalRxnLUT::k19][i-1] * H2II(i-1,j,k) * HM(i-1,j,k)/2.
              +       kcol_buf.data[CollisionalRxnLUT::k22][i-1] * HI(i-1,j,k) * std::pow((HI(i-1,j,k)),2.));
        acoef = ( kcol_buf.data[CollisionalRxnLUT::k13][i-1]*HI(i-1,j,k) + kcol_buf.data[CollisionalRxnLUT::k11][i-1]*HII(i-1,j,k)
                + kcol_buf.data[CollisionalRxnLUT::k12][i-1]*de(i-1,j,k) )
                + kshield_buf.k29[i-1] + kshield_buf.k31[i-1];

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i-1] != MASK_FALSE)  {
            scoef = scoef + 2. * h2dust[i-1] *
                 HI(i-1,j,k) * rhoH[i-1];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef + 2. * (
                  kcol_buf.data[CollisionalRxnLUT::k53][i-1] * HDI(i-1,j,k) * HII(i-1,j,k) / 3.
                + kcol_buf.data[CollisionalRxnLUT::k55][i-1] * HDI(i-1,j,k) * HI (i-1,j,k) / 3.
                   );
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k52][i-1] * DII(i-1,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k54][i-1] * DI (i-1,j,k) / 2.;
        }
#endif
        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i-1] != MASK_FALSE) )  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::kz15][i-1] *    HI(i-1,j,k) *    CH(i-1,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz16][i-1] *    HI(i-1,j,k) *   CH2(i-1,j,k) / 14.
              + kcol_buf.data[CollisionalRxnLUT::kz17][i-1] *    HI(i-1,j,k) *    OH(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz18][i-1] *    HI(i-1,j,k) *   H2O(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz47][i-1] * H2OII(i-1,j,k) *    de(i-1,j,k) / 18.
             );
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz20][i-1] *    CI(i-1,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz21][i-1] *    OI(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz23][i-1] *    CH(i-1,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz24][i-1] *    OH(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz40][i-1] *   OII(i-1,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz41][i-1] *  OHII(i-1,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz42][i-1] * H2OII(i-1,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz51][i-1] *    CI(i-1,j,k) / 12.;
          if ((my_chemistry->grain_growth == 1)  ||  (my_chemistry->dust_sublimation == 1))  {
            if (my_chemistry->dust_species > 0)  {
              scoef = scoef + 2. *
                    grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1] * 2.;

            }
            if (my_chemistry->dust_species > 1)  {
              scoef = scoef + 2. * (
                    grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1] * 3.
                  + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i-1] * 4.
                  + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i-1]
                  + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i-1] * 3.
                );
            }
            if (my_chemistry->dust_species > 2)  {
              acoef = acoef
              + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i-1] / H2I(i-1,j,k) * 2. * 2.;
            }
          }
        }
        species_tmpdens.data[SpLUT::H2I][i-1] = ( scoef*dtit[i-1] + H2I(i-1,j,k) )
                  / ( 1. + acoef*dtit[i-1] );

        // 8) H-

        scoef = kcol_buf.data[CollisionalRxnLUT::k7][i-1] * HI(i-1,j,k) * de(i-1,j,k);
        acoef = (kcol_buf.data[CollisionalRxnLUT::k8][i-1]  + kcol_buf.data[CollisionalRxnLUT::k15][i-1])  * HI(i-1,j,k) +
                (kcol_buf.data[CollisionalRxnLUT::k16][i-1] + kcol_buf.data[CollisionalRxnLUT::k17][i-1])  * HII(i-1,j,k) +
               kcol_buf.data[CollisionalRxnLUT::k14][i-1] * de(i-1,j,k) + kcol_buf.data[CollisionalRxnLUT::k19][i-1] * H2II(i-1,j,k)/2.0f +
               my_uvb_rates.k27;
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k56][i-1] * DI (i-1,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k136][i-1] *    DM(i-1,j,k) *    HI(i-1,j,k) /  2.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k135][i-1] *    DI(i-1,j,k) /  2.;
        }
        species_tmpdens.data[SpLUT::HM][i-1] = (scoef*dtit[i-1] + HM(i-1,j,k))
             / (1.0f + acoef*dtit[i-1]);


        // 9) H2+

        species_tmpdens.data[SpLUT::H2II][i-1] = 2.*( kcol_buf.data[CollisionalRxnLUT::k9] [i-1]*species_tmpdens.data[SpLUT::HI][i-1]*species_tmpdens.data[SpLUT::HII][i-1]
                      +   kcol_buf.data[CollisionalRxnLUT::k11][i-1]*species_tmpdens.data[SpLUT::H2I][i-1]/2.*species_tmpdens.data[SpLUT::HII][i-1]
                      +   kcol_buf.data[CollisionalRxnLUT::k17][i-1]*species_tmpdens.data[SpLUT::HM][i-1]*species_tmpdens.data[SpLUT::HII][i-1]
                      + kshield_buf.k29[i-1]*species_tmpdens.data[SpLUT::H2I][i-1]
                      )
                   /  ( kcol_buf.data[CollisionalRxnLUT::k10][i-1]*species_tmpdens.data[SpLUT::HI][i-1] + kcol_buf.data[CollisionalRxnLUT::k18][i-1]*species_tmpdens.data[SpLUT::e][i-1]
                      + kcol_buf.data[CollisionalRxnLUT::k19][i-1]*species_tmpdens.data[SpLUT::HM][i-1]
                      + (kshield_buf.k28[i-1]+kshield_buf.k30[i-1])
                      );
        if (my_chemistry->primordial_chemistry > 3)  {
          species_tmpdens.data[SpLUT::H2II][i-1] = 2. * (  kcol_buf.data[CollisionalRxnLUT::k9] [i-1]*species_tmpdens.data[SpLUT::HI][i-1]*species_tmpdens.data[SpLUT::HII][i-1]
                       +   kcol_buf.data[CollisionalRxnLUT::k11][i-1]*species_tmpdens.data[SpLUT::H2I][i-1]/2.*species_tmpdens.data[SpLUT::HII][i-1]
                       +   kcol_buf.data[CollisionalRxnLUT::k17][i-1]*species_tmpdens.data[SpLUT::HM][i-1]*species_tmpdens.data[SpLUT::HII][i-1]
                       + kshield_buf.k29[i-1]*species_tmpdens.data[SpLUT::H2I][i-1]
                       + kcol_buf.data[CollisionalRxnLUT::k152][i-1]*HeHII(i-1,j,k)*species_tmpdens.data[SpLUT::HI][i-1]/5.
                       )
                    /  ( kcol_buf.data[CollisionalRxnLUT::k10][i-1]*species_tmpdens.data[SpLUT::HI][i-1] + kcol_buf.data[CollisionalRxnLUT::k18][i-1]*species_tmpdens.data[SpLUT::e][i-1]
                       + kcol_buf.data[CollisionalRxnLUT::k19][i-1]*species_tmpdens.data[SpLUT::HM][i-1]
                       + (kshield_buf.k28[i-1]+kshield_buf.k30[i-1])
                       + kcol_buf.data[CollisionalRxnLUT::k150][i-1]*species_tmpdens.data[SpLUT::HeI][i-1]/4.
                       );
        }
      }
    }
    // 
  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        
        // 1) DI
        scoef =   (       kcol_buf.data[CollisionalRxnLUT::k2][i-1] * DII(i-1,j,k) * de(i-1,j,k)
                   +      kcol_buf.data[CollisionalRxnLUT::k51][i-1]* DII(i-1,j,k) * HI(i-1,j,k)
                   + 2.*kcol_buf.data[CollisionalRxnLUT::k55][i-1]* HDI(i-1,j,k) *
                HI(i-1,j,k)/3.
                   );
        acoef  =    kcol_buf.data[CollisionalRxnLUT::k1][i-1] * de(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k50][i-1] * HII(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k54][i-1] * H2I(i-1,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k56][i-1] * HM(i-1,j,k)
               + kshield_buf.k24[i-1];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i-1,j,k); }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::k131][i-1] *  HDII(i-1,j,k) *    de(i-1,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k133][i-1] *   DII(i-1,j,k) *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k134][i-1] *   HII(i-1,j,k) *    DM(i-1,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k136][i-1] *    DM(i-1,j,k) *    HI(i-1,j,k) /  2.
              );
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k129][i-1] *   HII(i-1,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k132][i-1] *    de(i-1,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k135][i-1] *    HM(i-1,j,k);
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            scoef = scoef
              + 2. * kdissHDI(i-1,j,k) * HDI(i-1,j,k)/3.0;
          }
        }
        species_tmpdens.data[SpLUT::DI][i-1]    = ( scoef*dtit[i-1] + DI(i-1,j,k) ) /
                    ( 1. + acoef*dtit[i-1] );

        // 2) DII
        scoef =   (   kcol_buf.data[CollisionalRxnLUT::k1][i-1]  * DI(i-1,j,k) * de(i-1,j,k)
              +       kcol_buf.data[CollisionalRxnLUT::k50][i-1] * HII(i-1,j,k)* DI(i-1,j,k)
              +  2.*kcol_buf.data[CollisionalRxnLUT::k53][i-1] * HII(i-1,j,k)* HDI(i-1,j,k)/3.
              )
              + kshield_buf.k24[i-1]*DI(i-1,j,k);
        acoef = 0.;
        // ! initialize GC202002
        if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + kphHI(i-1,j,k)*DI(i-1,j,k); }
        acoef =    kcol_buf.data[CollisionalRxnLUT::k2][i-1]  * de(i-1,j,k)
              +    kcol_buf.data[CollisionalRxnLUT::k51][i-1] * HI(i-1,j,k)
              +    kcol_buf.data[CollisionalRxnLUT::k52][i-1] * H2I(i-1,j,k)/2.;
        if (my_chemistry->primordial_chemistry > 3)  {
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k130][i-1] *    HI(i-1,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k133][i-1] *    DM(i-1,j,k) /  2.;
        }
        species_tmpdens.data[SpLUT::DII][i-1]   = ( scoef*dtit[i-1] + DII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

        // 3) HDI
        scoef = 3.*(kcol_buf.data[CollisionalRxnLUT::k52][i-1] * DII(i-1,j,k)*
             H2I(i-1,j,k)/2./2.
             + kcol_buf.data[CollisionalRxnLUT::k54][i-1] * DI(i-1,j,k) * H2I(i-1,j,k)/2./2.
        // !   &           + 2._DKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
        //- ! corrected by GC202005
             +          kcol_buf.data[CollisionalRxnLUT::k56][i-1] * DI(i-1,j,k) * HM(i-1,j,k)/2.
                   );
        acoef  =    kcol_buf.data[CollisionalRxnLUT::k53][i-1] * HII(i-1,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k55][i-1] * HI(i-1,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            acoef = acoef
              + kdissHDI(i-1,j,k);
          }
        }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  3. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::k125][i-1] *  HDII(i-1,j,k) *    HI(i-1,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k137][i-1] *    DM(i-1,j,k) *    HI(i-1,j,k) /  2.
              );
        }
        species_tmpdens.data[SpLUT::HDI][i-1]   = ( scoef*dtit[i-1] + HDI(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

      }
    }
  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        
        // 1) DM
        scoef =
              kcol_buf.data[CollisionalRxnLUT::k132][i-1] *    DI(i-1,j,k) *    de(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k135][i-1] *    HM(i-1,j,k) *    DI(i-1,j,k);
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k133][i-1] *   DII(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::k134][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k136][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k137][i-1] *    HI(i-1,j,k);

        species_tmpdens.data[SpLUT::DM][i-1]    = ( scoef*dtit[i-1] + DM(i-1,j,k) ) /
                    ( 1. + acoef*dtit[i-1] );

        // 2) HDII
        scoef = 3. * (
              kcol_buf.data[CollisionalRxnLUT::k129][i-1] *    DI(i-1,j,k) *   HII(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::k130][i-1] *   DII(i-1,j,k) *    HI(i-1,j,k) /  2.
            );
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k125][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k131][i-1] *    de(i-1,j,k);

        species_tmpdens.data[SpLUT::HDII][i-1]   = ( scoef*dtit[i-1] + HDII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

        // 3) HeHII
        scoef = 5. * (
              kcol_buf.data[CollisionalRxnLUT::k148][i-1] *   HeI(i-1,j,k) *   HII(i-1,j,k) /  4.
            + kcol_buf.data[CollisionalRxnLUT::k149][i-1] *   HeI(i-1,j,k) *   HII(i-1,j,k) /  4.
            + kcol_buf.data[CollisionalRxnLUT::k150][i-1] *   HeI(i-1,j,k) *  H2II(i-1,j,k) /  8.
            + kcol_buf.data[CollisionalRxnLUT::k151][i-1] *  HeII(i-1,j,k) *    HI(i-1,j,k) /  4.
            );
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k152][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k153][i-1] *    de(i-1,j,k);

        species_tmpdens.data[SpLUT::HeHII][i-1]   = ( scoef*dtit[i-1] + HeHII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

      }
    }
  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask_metal[i-1] != MASK_FALSE)  {

        // ***** CI **********
        scoef = 0. + 12. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz15][i-1] *    HI(i-1,j,k) *    CH(i-1,j,k) / 13.
            + kcol_buf.data[CollisionalRxnLUT::kz44][i-1] *   CII(i-1,j,k) *    de(i-1,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz20][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz27][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz28][i-1] *    OH(i-1,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i-1] *    O2(i-1,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz51][i-1] *   H2I(i-1,j,k) /  2.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust]      [i-1] / CI(i-1,j,k) * 12.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            acoef = acoef
              + kphCI(i-1,j,k);
          }
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            scoef = scoef + 12. *
                kdissCO (i-1,j,k) * CO(i-1,j,k) /28.0;
          }
        }

        species_tmpdens.data[SpLUT::CI][i-1]   = ( scoef*dtit[i-1] + CI(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** CII **********
        scoef = 0. + 12. * ( 0.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz37][i-1] *    OH(i-1,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i-1] *    O2(i-1,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz44][i-1] *    de(i-1,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphCI(i-1,j,k) * CI(i-1,j,k);
          }
        }

        species_tmpdens.data[SpLUT::CII][i-1]   = ( scoef*dtit[i-1] + CII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** CO **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz28][i-1] *    CI(i-1,j,k) *    OH(i-1,j,k) / 204.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i-1] *    CI(i-1,j,k) *    O2(i-1,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i-1] *    OI(i-1,j,k) *    CH(i-1,j,k) / 208.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i-1] *   CII(i-1,j,k) *    O2(i-1,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz43][i-1] *  COII(i-1,j,k) *    HI(i-1,j,k) / 28.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz26][i-1] *    OH(i-1,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i-1] / CO(i-1,j,k) * 17. * 0.5
            + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i-1] / CO(i-1,j,k) * 17.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissCO (i-1,j,k);
          }
        }

        species_tmpdens.data[SpLUT::CO][i-1]   = ( scoef*dtit[i-1] + CO(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** CO2 **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz26][i-1] *    OH(i-1,j,k) *    CO(i-1,j,k) / 476.
           );
        acoef = 0.;

        species_tmpdens.data[SpLUT::CO2][i-1]   = ( scoef*dtit[i-1] + CO2(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** OI **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz17][i-1] *    HI(i-1,j,k) *    OH(i-1,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i-1] *    HI(i-1,j,k) *    O2(i-1,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i-1] *    OH(i-1,j,k) *    OH(i-1,j,k) / 289.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i-1] *    CI(i-1,j,k) *    O2(i-1,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz39][i-1] *   OII(i-1,j,k) *    HI(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz45][i-1] *   OII(i-1,j,k) *    de(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz47][i-1] * H2OII(i-1,j,k) *    de(i-1,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz50][i-1] *  O2II(i-1,j,k) *    de(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i-1] *   SiI(i-1,j,k) *    O2(i-1,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz21][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz22][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz30][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz31][i-1] *    OI(i-1,j,k) / 8.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i-1] *    CH(i-1,j,k) / 13.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i-1] *    OH(i-1,j,k) / 17.;
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            acoef = acoef
              + kphOI(i-1,j,k);
          }
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            scoef = scoef + 16. *
              ( kdissOH (i-1,j,k) * OH(i-1,j,k) /17.0
              + kdissCO (i-1,j,k) * CO(i-1,j,k) /28.0);
          }
        }

        species_tmpdens.data[SpLUT::OI][i-1]   = ( scoef*dtit[i-1] + OI(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** OH **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz18][i-1] *    HI(i-1,j,k) *   H2O(i-1,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i-1] *    HI(i-1,j,k) *    O2(i-1,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz21][i-1] *    OI(i-1,j,k) *   H2I(i-1,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz30][i-1] *    OI(i-1,j,k) *    HI(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz46][i-1] * H2OII(i-1,j,k) *    de(i-1,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz49][i-1] * H3OII(i-1,j,k) *    de(i-1,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz17][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz24][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i-1] *    OH(i-1,j,k) / 8.5
            + kcol_buf.data[CollisionalRxnLUT::kz26][i-1] *    CO(i-1,j,k) / 28.
            + kcol_buf.data[CollisionalRxnLUT::kz28][i-1] *    CI(i-1,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i-1] *    OI(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz34][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz37][i-1] *   CII(i-1,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i-1] *   SiI(i-1,j,k) / 28.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i-1] *  SiOI(i-1,j,k) / 44.;
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissOH (i-1,j,k);
            scoef = scoef + 17. *
                kdissH2O(i-1,j,k) * H2O(i-1,j,k)/18.0;
          }
        }

        species_tmpdens.data[SpLUT::OH][i-1]   = ( scoef*dtit[i-1] + OH(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** H2O **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz24][i-1] *   H2I(i-1,j,k) *    OH(i-1,j,k) / 34.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i-1] *    OH(i-1,j,k) *    OH(i-1,j,k) / 289.
            + kcol_buf.data[CollisionalRxnLUT::kz48][i-1] * H3OII(i-1,j,k) *    de(i-1,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz18][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz35][i-1] *   HII(i-1,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1] / H2O(i-1,j,k) * 18. * 2.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1] / H2O(i-1,j,k) * 18. * 3.
            + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i-1] / H2O(i-1,j,k) * 18. * 4.
            + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i-1] / H2O(i-1,j,k) * 18.
            + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i-1] / H2O(i-1,j,k) * 18. * 3.;
          }
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust]  [i-1] / H2O(i-1,j,k) * 18.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissH2O(i-1,j,k);
          }
        }

        species_tmpdens.data[SpLUT::H2O][i-1]   = ( scoef*dtit[i-1] + H2O(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** O2 **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz31][i-1] *    OI(i-1,j,k) *    OI(i-1,j,k) / 256.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i-1] *    OI(i-1,j,k) *    OH(i-1,j,k) / 272.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz29][i-1] *    CI(i-1,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz36][i-1] *   HII(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz38][i-1] *   CII(i-1,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i-1] *   SiI(i-1,j,k) / 28.;

        species_tmpdens.data[SpLUT::O2][i-1]   = ( scoef*dtit[i-1] + O2(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** SiI **********
        scoef = 0. + 28. * ( 0.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i-1] *    OH(i-1,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i-1] *    O2(i-1,j,k) / 32.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust]     [i-1] / SiI(i-1,j,k) * 28.;
          }
        }

        species_tmpdens.data[SpLUT::SiI][i-1]   = ( scoef*dtit[i-1] + SiI(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** SiOI **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i-1] *   SiI(i-1,j,k) *    OH(i-1,j,k) / 476.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i-1] *   SiI(i-1,j,k) *    O2(i-1,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i-1] *    OH(i-1,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1] / SiOI(i-1,j,k) * 44.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1] / SiOI(i-1,j,k) * 44.;
          }
        }

        species_tmpdens.data[SpLUT::SiOI][i-1]   = ( scoef*dtit[i-1] + SiOI(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** SiO2I **********
        scoef = 0. + 60. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i-1] *  SiOI(i-1,j,k) *    OH(i-1,j,k) / 748.
           );
        acoef = 0.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust]   [i-1] / SiO2I(i-1,j,k) * 60.;
          }
        }

        species_tmpdens.data[SpLUT::SiO2I][i-1]   = ( scoef*dtit[i-1] + SiO2I(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );

        // MINOR BUT IMPORTANT SPECIES FOR MOLECULAR FORMATION
        //- ***** CH **********
        scoef = 0. + 13. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz16][i-1] *    HI(i-1,j,k) *   CH2(i-1,j,k) / 14.
            + kcol_buf.data[CollisionalRxnLUT::kz20][i-1] *    CI(i-1,j,k) *   H2I(i-1,j,k) / 24.
            + kcol_buf.data[CollisionalRxnLUT::kz27][i-1] *    CI(i-1,j,k) *    HI(i-1,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz15][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz23][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i-1] *    OI(i-1,j,k) / 16.;

        species_tmpdens.data[SpLUT::CH][i-1]   = ( scoef*dtit[i-1] + CH(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** CH2 **********
        scoef = 0. + 14. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz23][i-1] *   H2I(i-1,j,k) *    CH(i-1,j,k) / 26.
            + kcol_buf.data[CollisionalRxnLUT::kz51][i-1] *   H2I(i-1,j,k) *    CI(i-1,j,k) / 24.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz16][i-1] *    HI(i-1,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i-1] / CH2(i-1,j,k) * 14. * 0.5;
          }
        }

        species_tmpdens.data[SpLUT::CH2][i-1]   = ( scoef*dtit[i-1] + CH2(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** COII **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz37][i-1] *   CII(i-1,j,k) *    OH(i-1,j,k) / 204.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz43][i-1] *    HI(i-1,j,k);

        species_tmpdens.data[SpLUT::COII][i-1]   = ( scoef*dtit[i-1] + COII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** OII **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz22][i-1] *   HII(i-1,j,k) *    OI(i-1,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i-1] *   CII(i-1,j,k) *    O2(i-1,j,k) / 384.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz39][i-1] *    HI(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz40][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz45][i-1] *    de(i-1,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphOI(i-1,j,k) * OI(i-1,j,k);
          }
        }

        species_tmpdens.data[SpLUT::OII][i-1]   = ( scoef*dtit[i-1] + OII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** OHII **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz34][i-1] *   HII(i-1,j,k) *    OH(i-1,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz40][i-1] *   OII(i-1,j,k) *   H2I(i-1,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz41][i-1] *   H2I(i-1,j,k) /  2.;

        species_tmpdens.data[SpLUT::OHII][i-1]   = ( scoef*dtit[i-1] + OHII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** H2OII **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz35][i-1] *   HII(i-1,j,k) *   H2O(i-1,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz41][i-1] *  OHII(i-1,j,k) *   H2I(i-1,j,k) / 34.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz42][i-1] *   H2I(i-1,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz46][i-1] *    de(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz47][i-1] *    de(i-1,j,k);

        species_tmpdens.data[SpLUT::H2OII][i-1]   = ( scoef*dtit[i-1] + H2OII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** H3OII **********
        scoef = 0. + 19. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz42][i-1] * H2OII(i-1,j,k) *   H2I(i-1,j,k) / 36.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz48][i-1] *    de(i-1,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz49][i-1] *    de(i-1,j,k);

        species_tmpdens.data[SpLUT::H3OII][i-1]   = ( scoef*dtit[i-1] + H3OII(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        // ***** O2II **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz36][i-1] *   HII(i-1,j,k) *    O2(i-1,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz50][i-1] *    de(i-1,j,k);

        species_tmpdens.data[SpLUT::O2II][i-1]   = ( scoef*dtit[i-1] + O2II(i-1,j,k) )
                   / ( 1. + acoef*dtit[i-1] );


        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            // ***** Mg **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1] / Mg(i-1,j,k) * 24.;
            if (my_chemistry->dust_species > 1)  {
              acoef = acoef
              + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1] / Mg(i-1,j,k) * 24. * 2.
              + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i-1] / Mg(i-1,j,k) * 24.;
            }

            species_tmpdens.data[SpLUT::Mg][i-1]   = ( scoef*dtit[i-1] + Mg(i-1,j,k) )
                       / ( 1. + acoef*dtit[i-1] );

          }

          if (my_chemistry->dust_species > 1)  {
            // ***** Al **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i-1] / Al(i-1,j,k) * 27. * 2.;

            species_tmpdens.data[SpLUT::Al][i-1]   = ( scoef*dtit[i-1] + Al(i-1,j,k) )
                       / ( 1. + acoef*dtit[i-1] );


            // ***** S  **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i-1] / S(i-1,j,k) * 32.;

            species_tmpdens.data[SpLUT::S][i-1]    = ( scoef*dtit[i-1] + S(i-1,j,k) )
                       / ( 1. + acoef*dtit[i-1] );


            // ***** Fe **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust]     [i-1] / Fe(i-1,j,k) * 56.
            + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i-1] / Fe(i-1,j,k) * 56. * 3.
            + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i-1] / Fe(i-1,j,k) * 56.;

            species_tmpdens.data[SpLUT::Fe][i-1]   = ( scoef*dtit[i-1] + Fe(i-1,j,k) )
                       / ( 1. + acoef*dtit[i-1] );

          }
        }

      }
    }
  }

  // --- (D4) Now do dust species ---
  if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask_metal[i-1] != MASK_FALSE)  {

        if (my_chemistry->dust_species > 0)  {
          // ***** MgSiO3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1] * 100.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::MgSiO3_dust][i-1]   = ( scoef*dtit[i-1] + MgSiO3(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** AC **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust]      [i-1] * 12.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::AC_dust][i-1]   = ( scoef*dtit[i-1] + AC(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );

        }

        if (my_chemistry->dust_species > 1)  {
          // ***** SiM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust]     [i-1] * 28.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::SiM_dust][i-1]   = ( scoef*dtit[i-1] + SiM(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** FeM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust]     [i-1] * 56.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::FeM_dust][i-1]   = ( scoef*dtit[i-1] + FeM(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** Mg2SiO4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1] * 140.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Mg2SiO4_dust][i-1]   = ( scoef*dtit[i-1] + Mg2SiO4(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** Fe3O4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i-1] * 232.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Fe3O4_dust][i-1]   = ( scoef*dtit[i-1] + Fe3O4(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** SiO2D **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust]   [i-1] * 60.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::SiO2_dust][i-1]   = ( scoef*dtit[i-1] + SiO2D(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** MgO **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i-1] * 40.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::MgO_dust][i-1]   = ( scoef*dtit[i-1] + MgO(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** FeS **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i-1] * 88.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::FeS_dust][i-1]   = ( scoef*dtit[i-1] + FeS(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** Al2O3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i-1] * 102.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Al2O3_dust][i-1]   = ( scoef*dtit[i-1] + Al2O3(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );

        }

        if (my_chemistry->dust_species > 2)  {
          // ***** reforg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i-1] * 22.68;
          acoef = 0.;

          species_tmpdens.data[SpLUT::ref_org_dust][i-1]   = ( scoef*dtit[i-1] + reforg(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** volorg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i-1] * 32.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::vol_org_dust][i-1]   = ( scoef*dtit[i-1] + volorg(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );


          // ***** H2Oice **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust]  [i-1] * 18.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::H2O_ice_dust][i-1]   = ( scoef*dtit[i-1] + H2Oice(i-1,j,k) )
                     / ( 1. + acoef*dtit[i-1] );

        }

      }
    }
  }

  // --- (E) Set densities from 1D temps to 3D fields ---

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if (itmask[i-1] != MASK_FALSE)  {
      HIdot_prev[i-1] = std::fabs(HI(i-1,j,k)-species_tmpdens.data[SpLUT::HI][i-1]) /
              std::fmax((double)(dtit[i-1] ), tiny8);
      HI(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HI][i-1] ), tiny_fortran_val);
      HII(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HII][i-1] ), tiny_fortran_val);
      HeI(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeI][i-1] ), tiny_fortran_val);
      HeII(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeII][i-1] ), tiny_fortran_val);
      HeIII(i-1,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeIII][i-1] ), (gr_float)(1e-5)*tiny_fortran_val);

      // de(i,j,k)    = dep(i)

      // Use charge conservation to determine electron fraction

      dedot_prev[i-1] = de(i-1,j,k);
      de(i-1,j,k) = HII(i-1,j,k) + HeII(i-1,j,k)/(gr_float)(4.) +
           HeIII(i-1,j,k)/(gr_float)(2.);
      if (my_chemistry->primordial_chemistry > 1)
           { de(i-1,j,k) = de(i-1,j,k) - HM(i-1,j,k) + H2II(i-1,j,k)/(gr_float)(2.); }

      if (my_chemistry->primordial_chemistry > 2)
           { de(i-1,j,k) = de(i-1,j,k) + DII(i-1,j,k)/(gr_float)(2.); }
      if (my_chemistry->primordial_chemistry > 3)
           { de(i-1,j,k) = de(i-1,j,k) - DM(i-1,j,k)/(gr_float)(2.)
                + HDII(i-1,j,k)/(gr_float)(3.) + HeHII(i-1,j,k)/(gr_float)(5.); }
      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i-1] != MASK_FALSE) )
           { de(i-1,j,k) = de(i-1,j,k)
                + CII(i-1,j,k)/(gr_float)(12.) + COII(i-1,j,k)/(gr_float)(28.)
                + OII(i-1,j,k)/(gr_float)(16.) + OHII(i-1,j,k)/(gr_float)(17.)
                + H2OII(i-1,j,k)/(gr_float)(18.) + H3OII(i-1,j,k)/(gr_float)(19.)
                + O2II(i-1,j,k)/(gr_float)(32.); }

      dedot_prev[i-1] = std::fabs(de(i-1,j,k)-dedot_prev[i-1])/
           std::fmax(dtit[i-1],tiny8);

      if (my_chemistry->primordial_chemistry > 1)  {
        HM(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HM][i-1] ), tiny_fortran_val);
        H2I(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2I][i-1]), tiny_fortran_val);
        H2II(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2II][i-1] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 2)  {
        DI(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DI][i-1] ), tiny_fortran_val);
        DII(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DII][i-1] ), tiny_fortran_val);
        HDI(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDI][i-1] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 3)  {
        DM(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DM][i-1] ), tiny_fortran_val);
        HDII(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDII][i-1] ), tiny_fortran_val);
        HeHII(i-1,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeHII][i-1] ), tiny_fortran_val);
      }

      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i-1] != MASK_FALSE) )  {
        CI(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CI][i-1]      ), tiny_fortran_val);
        CII(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CII][i-1]     ), tiny_fortran_val);
        CO(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CO][i-1]      ), tiny_fortran_val);
        CO2(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CO2][i-1]     ), tiny_fortran_val);
        OI(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OI][i-1]      ), tiny_fortran_val);
        OH(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OH][i-1]      ), tiny_fortran_val);
        H2O(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2O][i-1]     ), tiny_fortran_val);
        O2(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2][i-1]      ), tiny_fortran_val);
        SiI(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiI][i-1]     ), tiny_fortran_val);
        SiOI(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiOI][i-1]    ), tiny_fortran_val);
        SiO2I(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2I][i-1]   ), tiny_fortran_val);
        CH(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CH][i-1]      ), tiny_fortran_val);
        CH2(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CH2][i-1]     ), tiny_fortran_val);
        COII(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::COII][i-1]    ), tiny_fortran_val);
        OII(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OII][i-1]     ), tiny_fortran_val);
        OHII(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OHII][i-1]    ), tiny_fortran_val);
        H2OII(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2OII][i-1]   ), tiny_fortran_val);
        H3OII(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H3OII][i-1]   ), tiny_fortran_val);
        O2II(i-1,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2II][i-1]    ), tiny_fortran_val);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            Mg(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg][i-1]      ), tiny_fortran_val);
          }
          if (my_chemistry->dust_species > 1)  {
            Al(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al][i-1]      ), tiny_fortran_val);
            S(i-1,j,k)       = std::fmax((gr_float)(species_tmpdens.data[SpLUT::S][i-1]       ), tiny_fortran_val);
            Fe(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe][i-1]      ), tiny_fortran_val);
          }
        }
      }

      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          MgSiO3(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgSiO3_dust][i-1]  ), tiny_fortran_val);
          AC(i-1,j,k)      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::AC_dust][i-1]      ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 1)  {
          SiM(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiM_dust][i-1]     ), tiny_fortran_val);
          FeM(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeM_dust][i-1]     ), tiny_fortran_val);
          Mg2SiO4(i-1,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg2SiO4_dust][i-1] ), tiny_fortran_val);
          Fe3O4(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe3O4_dust][i-1]   ), tiny_fortran_val);
          SiO2D(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2_dust][i-1]   ), tiny_fortran_val);
          MgO(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgO_dust][i-1]     ), tiny_fortran_val);
          FeS(i-1,j,k)     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeS_dust][i-1]     ), tiny_fortran_val);
          Al2O3(i-1,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al2O3_dust][i-1]   ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 2)  {
          reforg(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::ref_org_dust][i-1]   ), tiny_fortran_val);
          volorg(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::vol_org_dust][i-1]   ), tiny_fortran_val);
          H2Oice(i-1,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2O_ice_dust][i-1]   ), tiny_fortran_val);
        }
      }

    }
    // 

    if (HI(i-1,j,k) != HI(i-1,j,k))  {
      OMP_PRAGMA_CRITICAL
      {
        printf("HUGE HI! ::  %d %d %d %g\n",
               i,
               idx_range.jp1,
               idx_range.kp1,
               HI ( i-1, j, k ));
      }
    }

  }

  return;
}

}  // namespace grackle::impl

#endif /* STEP_RATE_GAUSS_SEIDEL_HPP */
