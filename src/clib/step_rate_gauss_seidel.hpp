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

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI

        scoef  = kcol_buf.data[CollisionalRxnLUT::k2][i]*HII(i,j,k)*de(i,j,k);
        acoef  = kcol_buf.data[CollisionalRxnLUT::k1][i]*de(i,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k57][i]*HI(i,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k58][i]*HeI(i,j,k)/4.
               + kshield_buf.k24[i];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i,j,k); }
        species_tmpdens.data[SpLUT::HI][i]  = (scoef*dtit[i] + HI(i,j,k))/
             (1. + acoef*dtit[i]);
        if (species_tmpdens.data[SpLUT::HI][i] != species_tmpdens.data[SpLUT::HI][i])  {
          OMP_PRAGMA_CRITICAL
          {
            printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g %g %g\n",
                   i+1,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HI] [ i ],
                   HI ( i, j, k ),
                   HII ( i, j, k ),
                   de ( i, j, k ),
                   kphHI ( i, j, k ),
                   scoef,
                   acoef,
                   dtit [ i ]);
          }
          // ERROR_MESSAGE
        }

        // 2) HII
        scoef  = kcol_buf.data[CollisionalRxnLUT::k1][i]*species_tmpdens.data[SpLUT::HI][i]*de(i,j,k)
               + kcol_buf.data[CollisionalRxnLUT::k57][i]*species_tmpdens.data[SpLUT::HI][i]*species_tmpdens.data[SpLUT::HI][i]
               + kcol_buf.data[CollisionalRxnLUT::k58][i]*species_tmpdens.data[SpLUT::HI][i]*HeI(i,j,k)/4.
               + kshield_buf.k24[i]*species_tmpdens.data[SpLUT::HI][i];
        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i,j,k)*species_tmpdens.data[SpLUT::HI][i]; }
        acoef  = kcol_buf.data[CollisionalRxnLUT::k2][i]*de (i,j,k);
        species_tmpdens.data[SpLUT::HII][i] = (scoef*dtit[i] + HII(i,j,k))/
             (1. +acoef*dtit[i]);
        // 
        if (species_tmpdens.data[SpLUT::HII][i] <= 0.)   //#####
        {
          OMP_PRAGMA_CRITICAL
          {
            printf("negative HIIp! ::  %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
                   i+1,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HII] [ i ],
                   scoef,
                   dtit [ i ],
                   HII ( i, j, k ),
                   acoef,
                   kcol_buf.data[CollisionalRxnLUT::k2] [ i ],
                   de ( i, j, k ),
                   kphHI ( i, j, k ),
                   species_tmpdens.data[SpLUT::HI] [ i ],
                   kshield_buf.k24 [ i ]);
          }
        }

        // 3) Electron density

        scoef = 0.
                   + kcol_buf.data[CollisionalRxnLUT::k57][i]*species_tmpdens.data[SpLUT::HI][i]*species_tmpdens.data[SpLUT::HI][i]
                   + kcol_buf.data[CollisionalRxnLUT::k58][i]*species_tmpdens.data[SpLUT::HI][i]*HeI(i,j,k)/4.
                   + kshield_buf.k24[i]*HI(i,j,k)
                   + kshield_buf.k25[i]*HeII(i,j,k)/4.
                   + kshield_buf.k26[i]*HeI(i,j,k)/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k)
                  + kphHeI(i,j,k)  * HeI(i,j,k)  / 4.
                  + kphHeII(i,j,k) * HeII(i,j,k) / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k); }



        acoef = -(kcol_buf.data[CollisionalRxnLUT::k1][i]*HI(i,j,k)      - kcol_buf.data[CollisionalRxnLUT::k2][i]*HII(i,j,k)
                + kcol_buf.data[CollisionalRxnLUT::k3][i]*HeI(i,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)/4.
                + kcol_buf.data[CollisionalRxnLUT::k5][i]*HeII(i,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k4][i]*HeII(i,j,k)/4.);
        species_tmpdens.data[SpLUT::e][i]   = (scoef*dtit[i] + de(i,j,k))
                       / (1. + acoef*dtit[i]);

      }
    }

  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---

  for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {

      // 4) HeI

      scoef  = kcol_buf.data[CollisionalRxnLUT::k4][i]*HeII(i,j,k)*de(i,j,k);
      acoef  = kcol_buf.data[CollisionalRxnLUT::k3][i]*de(i,j,k)
                   + kshield_buf.k26[i];

      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeI(i,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        scoef = scoef +  4. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::k152][i] * HeHII(i,j,k) *    HI(i,j,k) /  5.
            + kcol_buf.data[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) *    de(i,j,k) /  5.
            );
        acoef = acoef
            + kcol_buf.data[CollisionalRxnLUT::k148][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k149][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k150][i] *  H2II(i,j,k) /  2.;
      }
      species_tmpdens.data[SpLUT::HeI][i]   = ( scoef*dtit[i] + HeI(i,j,k) )
                 / ( 1. + acoef*dtit[i] );

      // 5) HeII

      scoef  = kcol_buf.data[CollisionalRxnLUT::k3][i]*species_tmpdens.data[SpLUT::HeI][i]*de(i,j,k)
             + kcol_buf.data[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)*de(i,j,k)
             + kshield_buf.k26[i]*species_tmpdens.data[SpLUT::HeI][i];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeI(i,j,k)*species_tmpdens.data[SpLUT::HeI][i]; }

      acoef  = kcol_buf.data[CollisionalRxnLUT::k4][i]*de(i,j,k) + kcol_buf.data[CollisionalRxnLUT::k5][i]*de(i,j,k)
             + kshield_buf.k25[i];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeII(i,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        acoef = acoef
            + kcol_buf.data[CollisionalRxnLUT::k151][i] *    HI(i,j,k);
      }
      species_tmpdens.data[SpLUT::HeII][i]  = ( scoef*dtit[i] + HeII(i,j,k) )
                 / ( 1. + acoef*dtit[i] );

      // 6) HeIII

      scoef   = kcol_buf.data[CollisionalRxnLUT::k5][i]*species_tmpdens.data[SpLUT::HeII][i]*de(i,j,k)
              + kshield_buf.k25[i]*species_tmpdens.data[SpLUT::HeII][i];
      if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeII(i,j,k) * species_tmpdens.data[SpLUT::HeII][i]; }
      acoef   = kcol_buf.data[CollisionalRxnLUT::k6][i]*de(i,j,k);
      species_tmpdens.data[SpLUT::HeIII][i]  = ( scoef*dtit[i] + HeIII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

    }
  }

  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI
        scoef  =      kcol_buf.data[CollisionalRxnLUT::k2][i] * HII(i,j,k) * de(i,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k13][i]* HI(i,j,k)  * H2I(i,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k11][i]* HII(i,j,k) * H2I(i,j,k)/2.
               + 2.*kcol_buf.data[CollisionalRxnLUT::k12][i]* de(i,j,k)  * H2I(i,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k14][i]* HM(i,j,k)  * de(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k15][i]* HM(i,j,k)  * HI(i,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k16][i]* HM(i,j,k)  * HII(i,j,k)
               + 2.*kcol_buf.data[CollisionalRxnLUT::k18][i]* H2II(i,j,k)* de(i,j,k)/2.
               +      kcol_buf.data[CollisionalRxnLUT::k19][i]* H2II(i,j,k)* HM(i,j,k)/2.
               + 2.*kshield_buf.k31[i]   * H2I(i,j,k)/2.;

        acoef  =      kcol_buf.data[CollisionalRxnLUT::k1][i] * de(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k7][i] * de(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k8][i] * HM(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k9][i] * HII(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k10][i]* H2II(i,j,k)/2.
               + 2.*kcol_buf.data[CollisionalRxnLUT::k22][i]* std::pow(HI(i,j,k),2)
               +      kcol_buf.data[CollisionalRxnLUT::k57][i]* HI(i,j,k)
               +      kcol_buf.data[CollisionalRxnLUT::k58][i]* HeI(i,j,k)/4.
               + kshield_buf.k24[i];

        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i,j,k); }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if ((my_chemistry->primordial_chemistry > 2) && (my_chemistry->radiative_transfer_HDI_dissociation > 0))  {
            scoef = scoef
              + kdissHDI(i,j,k) * HDI(i,j,k)/3.0;
          }
          if ( (my_chemistry->metal_chemistry == 1)  && 
               (itmask_metal[i] != MASK_FALSE) )  {
            if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
              scoef = scoef
                + kdissOH (i,j,k) * OH(i,j,k) /17.0
                + kdissH2O(i,j,k) * H2O(i,j,k)/18.0;
            }
          }
        }

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i] != MASK_FALSE)  {
            acoef = acoef + 2. * h2dust[i] * rhoH[i];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k50][i] * HII(i,j,k) * DI(i,j,k)  / 2.
                + kcol_buf.data[CollisionalRxnLUT::k54][i] * H2I(i,j,k) * DI(i,j,k)  / 4.;
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k51][i] * DII(i,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k55][i] * HDI(i,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) *    de(i,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k134][i] *   HII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k135][i] *    HM(i,j,k) *    DI(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k150][i] *   HeI(i,j,k) *  H2II(i,j,k) /  8.
              + kcol_buf.data[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) *    de(i,j,k) /  5.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k130][i] *   DII(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k136][i] *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k137][i] *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k151][i] *  HeII(i,j,k) /  4.
              + kcol_buf.data[CollisionalRxnLUT::k152][i] * HeHII(i,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) *   H2I(i,j,k) / 24.
              + kcol_buf.data[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) *   H2I(i,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz22][i] *   HII(i,j,k) *    OI(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) *    CH(i,j,k) / 26.
              + kcol_buf.data[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) *    OH(i,j,k) / 34.
              + kcol_buf.data[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) *    CO(i,j,k) / 476.
              + kcol_buf.data[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) *    OH(i,j,k) / 204.
              + kcol_buf.data[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) *    CH(i,j,k) / 208.
              + kcol_buf.data[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) *    OH(i,j,k) / 272.
              + kcol_buf.data[CollisionalRxnLUT::kz34][i] *   HII(i,j,k) *    OH(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz35][i] *   HII(i,j,k) *   H2O(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz36][i] *   HII(i,j,k) *    O2(i,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) *    OH(i,j,k) / 204.
              + kcol_buf.data[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) *   H2I(i,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) *   H2I(i,j,k) / 34.
              + kcol_buf.data[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) *   H2I(i,j,k) / 36.
              + kcol_buf.data[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) *    de(i,j,k) / 9.5
              + kcol_buf.data[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) *    OH(i,j,k) / 476.
              + kcol_buf.data[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) *    OH(i,j,k) / 748.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz15][i] *    CH(i,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz16][i] *   CH2(i,j,k) / 14.
              + kcol_buf.data[CollisionalRxnLUT::kz17][i] *    OH(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz18][i] *   H2O(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz19][i] *    O2(i,j,k) / 32.
              + kcol_buf.data[CollisionalRxnLUT::kz27][i] *    CI(i,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz30][i] *    OI(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) / 28.;
        }
        species_tmpdens.data[SpLUT::HI][i]  = ( scoef*dtit[i] + HI(i,j,k) ) /
                        ( 1.f + acoef*dtit[i] );
        if (species_tmpdens.data[SpLUT::HI][i] != species_tmpdens.data[SpLUT::HI][i])  {
          OMP_PRAGMA_CRITICAL
          {
            printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g\n",
                   i+1,
                   idx_range.jp1,
                   idx_range.kp1,
                   species_tmpdens.data[SpLUT::HI] [ i ],
                   HI ( i, j, k ),
                   HII ( i, j, k ),
                   de ( i, j, k ),
                   H2I ( i, j, k ),
                   kphHI ( i, j, k ));
          }
        }

        // 2) HII

        scoef  =    kcol_buf.data[CollisionalRxnLUT::k1][i]  * HI(i,j,k) * de(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k10][i] * H2II(i,j,k)*HI(i,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k57][i] * HI(i,j,k) * HI(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k58][i] * HI(i,j,k) * HeI(i,j,k)/4.
               + kshield_buf.k24[i]*HI(i,j,k);

        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k); }

        acoef  =    kcol_buf.data[CollisionalRxnLUT::k2][i]  * de(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k9][i]  * HI(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k11][i] * H2I(i,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k16][i] * HM(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k17][i] * HM(i,j,k);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k51][i] * HI (i,j,k) * DII(i,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k52][i] * H2I(i,j,k) * DII(i,j,k) / 4.;
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k50][i] * DI (i,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k53][i] * HDI(i,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) *    HI(i,j,k) /  3.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k129][i] *    DI(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k134][i] *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k148][i] *   HeI(i,j,k) /  4.
              + kcol_buf.data[CollisionalRxnLUT::k149][i] *   HeI(i,j,k) /  4.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) *    HI(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) *    HI(i,j,k) / 28.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz22][i] *    OI(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz34][i] *    OH(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz35][i] *   H2O(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz36][i] *    O2(i,j,k) / 32.;
        }
        species_tmpdens.data[SpLUT::HII][i]   = ( scoef*dtit[i] + HII(i,j,k) )
                        / ( 1. + acoef*dtit[i] );
        
        // 3) electrons:

        scoef =   kcol_buf.data[CollisionalRxnLUT::k8][i] * HM(i,j,k) * HI(i,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k15][i]* HM(i,j,k) * HI(i,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k17][i]* HM(i,j,k) * HII(i,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k57][i]* HI(i,j,k) * HI(i,j,k)
               +  kcol_buf.data[CollisionalRxnLUT::k58][i]* HI(i,j,k) * HeI(i,j,k)/4.
        // 
               + kshield_buf.k24[i]*species_tmpdens.data[SpLUT::HI][i]
               + kshield_buf.k25[i]*species_tmpdens.data[SpLUT::HeII][i]/4.
               + kshield_buf.k26[i]*species_tmpdens.data[SpLUT::HeI][i]/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i,j,k) * species_tmpdens.data[SpLUT::HI][i]
                  + kphHeI(i,j,k)  * species_tmpdens.data[SpLUT::HeI][i]  / 4.
                  + kphHeII(i,j,k) * species_tmpdens.data[SpLUT::HeII][i] / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i,j,k) * species_tmpdens.data[SpLUT::HI][i]; }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if ( (my_chemistry->metal_chemistry == 1)  && 
               (itmask_metal[i] != MASK_FALSE) )  {
            if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
              scoef = scoef
                + kphCI(i,j,k) * CI(i,j,k)/12.0
                + kphOI(i,j,k) * OI(i,j,k)/16.0;
            }
          }
        }

        acoef = - (kcol_buf.data[CollisionalRxnLUT::k1][i] *HI(i,j,k)    - kcol_buf.data[CollisionalRxnLUT::k2][i]*HII(i,j,k)
                +  kcol_buf.data[CollisionalRxnLUT::k3][i] *HeI(i,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)/4.
                +  kcol_buf.data[CollisionalRxnLUT::k5][i] *HeII(i,j,k)/4. -
             kcol_buf.data[CollisionalRxnLUT::k4][i]*HeII(i,j,k)/4.
                +  kcol_buf.data[CollisionalRxnLUT::k14][i]*HM(i,j,k)
                -  kcol_buf.data[CollisionalRxnLUT::k7][i] *HI(i,j,k)
                -  kcol_buf.data[CollisionalRxnLUT::k18][i]*H2II(i,j,k)/2.);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf.data[CollisionalRxnLUT::k56][i] * DI (i,j,k) * HM(i,j,k) / 2.;
          acoef = acoef
                - kcol_buf.data[CollisionalRxnLUT::k1] [i] * DI (i,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k2] [i] * DII(i,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k137][i] *    DM(i,j,k) *    HI(i,j,k) /  2.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k132][i] *    DI(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz44][i] *   CII(i,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz45][i] *   OII(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) / 19.
              + kcol_buf.data[CollisionalRxnLUT::kz50][i] *  O2II(i,j,k) / 32.;
        }
        species_tmpdens.data[SpLUT::e][i]  = ( scoef*dtit[i] + de(i,j,k) )
                  / ( 1. + acoef*dtit[i] );

        // 7) H2

        scoef = 2.*(kcol_buf.data[CollisionalRxnLUT::k8][i]  * HM(i,j,k)   * HI(i,j,k)
              +       kcol_buf.data[CollisionalRxnLUT::k10][i] * H2II(i,j,k) * HI(i,j,k)/2.
              +       kcol_buf.data[CollisionalRxnLUT::k19][i] * H2II(i,j,k) * HM(i,j,k)/2.
              +       kcol_buf.data[CollisionalRxnLUT::k22][i] * HI(i,j,k) * std::pow((HI(i,j,k)),2.));
        acoef = ( kcol_buf.data[CollisionalRxnLUT::k13][i]*HI(i,j,k) + kcol_buf.data[CollisionalRxnLUT::k11][i]*HII(i,j,k)
                + kcol_buf.data[CollisionalRxnLUT::k12][i]*de(i,j,k) )
                + kshield_buf.k29[i] + kshield_buf.k31[i];

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i] != MASK_FALSE)  {
            scoef = scoef + 2. * h2dust[i] *
                 HI(i,j,k) * rhoH[i];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef + 2. * (
                  kcol_buf.data[CollisionalRxnLUT::k53][i] * HDI(i,j,k) * HII(i,j,k) / 3.
                + kcol_buf.data[CollisionalRxnLUT::k55][i] * HDI(i,j,k) * HI (i,j,k) / 3.
                   );
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k52][i] * DII(i,j,k) / 2.
                + kcol_buf.data[CollisionalRxnLUT::k54][i] * DI (i,j,k) / 2.;
        }
#endif
        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::kz15][i] *    HI(i,j,k) *    CH(i,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz16][i] *    HI(i,j,k) *   CH2(i,j,k) / 14.
              + kcol_buf.data[CollisionalRxnLUT::kz17][i] *    HI(i,j,k) *    OH(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz18][i] *    HI(i,j,k) *   H2O(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
             );
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) / 12.
              + kcol_buf.data[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz23][i] *    CH(i,j,k) / 13.
              + kcol_buf.data[CollisionalRxnLUT::kz24][i] *    OH(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) / 16.
              + kcol_buf.data[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) / 17.
              + kcol_buf.data[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) / 18.
              + kcol_buf.data[CollisionalRxnLUT::kz51][i] *    CI(i,j,k) / 12.;
          if ((my_chemistry->grain_growth == 1)  ||  (my_chemistry->dust_sublimation == 1))  {
            if (my_chemistry->dust_species > 0)  {
              scoef = scoef + 2. *
                    grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 2.;

            }
            if (my_chemistry->dust_species > 1)  {
              scoef = scoef + 2. * (
                    grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 3.
                  + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 4.
                  + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i]
                  + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i] * 3.
                );
            }
            if (my_chemistry->dust_species > 2)  {
              acoef = acoef
              + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i] / H2I(i,j,k) * 2. * 2.;
            }
          }
        }
        species_tmpdens.data[SpLUT::H2I][i] = ( scoef*dtit[i] + H2I(i,j,k) )
                  / ( 1. + acoef*dtit[i] );

        // 8) H-

        scoef = kcol_buf.data[CollisionalRxnLUT::k7][i] * HI(i,j,k) * de(i,j,k);
        acoef = (kcol_buf.data[CollisionalRxnLUT::k8][i]  + kcol_buf.data[CollisionalRxnLUT::k15][i])  * HI(i,j,k) +
                (kcol_buf.data[CollisionalRxnLUT::k16][i] + kcol_buf.data[CollisionalRxnLUT::k17][i])  * HII(i,j,k) +
               kcol_buf.data[CollisionalRxnLUT::k14][i] * de(i,j,k) + kcol_buf.data[CollisionalRxnLUT::k19][i] * H2II(i,j,k)/2.0f +
               my_uvb_rates.k27;
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          acoef = acoef
                + kcol_buf.data[CollisionalRxnLUT::k56][i] * DI (i,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf.data[CollisionalRxnLUT::k136][i] *    DM(i,j,k) *    HI(i,j,k) /  2.;
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k135][i] *    DI(i,j,k) /  2.;
        }
        species_tmpdens.data[SpLUT::HM][i] = (scoef*dtit[i] + HM(i,j,k))
             / (1.0f + acoef*dtit[i]);


        // 9) H2+

        species_tmpdens.data[SpLUT::H2II][i] = 2.*( kcol_buf.data[CollisionalRxnLUT::k9] [i]*species_tmpdens.data[SpLUT::HI][i]*species_tmpdens.data[SpLUT::HII][i]
                      +   kcol_buf.data[CollisionalRxnLUT::k11][i]*species_tmpdens.data[SpLUT::H2I][i]/2.*species_tmpdens.data[SpLUT::HII][i]
                      +   kcol_buf.data[CollisionalRxnLUT::k17][i]*species_tmpdens.data[SpLUT::HM][i]*species_tmpdens.data[SpLUT::HII][i]
                      + kshield_buf.k29[i]*species_tmpdens.data[SpLUT::H2I][i]
                      )
                   /  ( kcol_buf.data[CollisionalRxnLUT::k10][i]*species_tmpdens.data[SpLUT::HI][i] + kcol_buf.data[CollisionalRxnLUT::k18][i]*species_tmpdens.data[SpLUT::e][i]
                      + kcol_buf.data[CollisionalRxnLUT::k19][i]*species_tmpdens.data[SpLUT::HM][i]
                      + (kshield_buf.k28[i]+kshield_buf.k30[i])
                      );
        if (my_chemistry->primordial_chemistry > 3)  {
          species_tmpdens.data[SpLUT::H2II][i] = 2. * (  kcol_buf.data[CollisionalRxnLUT::k9] [i]*species_tmpdens.data[SpLUT::HI][i]*species_tmpdens.data[SpLUT::HII][i]
                       +   kcol_buf.data[CollisionalRxnLUT::k11][i]*species_tmpdens.data[SpLUT::H2I][i]/2.*species_tmpdens.data[SpLUT::HII][i]
                       +   kcol_buf.data[CollisionalRxnLUT::k17][i]*species_tmpdens.data[SpLUT::HM][i]*species_tmpdens.data[SpLUT::HII][i]
                       + kshield_buf.k29[i]*species_tmpdens.data[SpLUT::H2I][i]
                       + kcol_buf.data[CollisionalRxnLUT::k152][i]*HeHII(i,j,k)*species_tmpdens.data[SpLUT::HI][i]/5.
                       )
                    /  ( kcol_buf.data[CollisionalRxnLUT::k10][i]*species_tmpdens.data[SpLUT::HI][i] + kcol_buf.data[CollisionalRxnLUT::k18][i]*species_tmpdens.data[SpLUT::e][i]
                       + kcol_buf.data[CollisionalRxnLUT::k19][i]*species_tmpdens.data[SpLUT::HM][i]
                       + (kshield_buf.k28[i]+kshield_buf.k30[i])
                       + kcol_buf.data[CollisionalRxnLUT::k150][i]*species_tmpdens.data[SpLUT::HeI][i]/4.
                       );
        }
      }
    }
    // 
  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {
        
        // 1) DI
        scoef =   (       kcol_buf.data[CollisionalRxnLUT::k2][i] * DII(i,j,k) * de(i,j,k)
                   +      kcol_buf.data[CollisionalRxnLUT::k51][i]* DII(i,j,k) * HI(i,j,k)
                   + 2.*kcol_buf.data[CollisionalRxnLUT::k55][i]* HDI(i,j,k) *
                HI(i,j,k)/3.
                   );
        acoef  =    kcol_buf.data[CollisionalRxnLUT::k1][i] * de(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k50][i] * HII(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k54][i] * H2I(i,j,k)/2.
               +    kcol_buf.data[CollisionalRxnLUT::k56][i] * HM(i,j,k)
               + kshield_buf.k24[i];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i,j,k); }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) *    de(i,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k133][i] *   DII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k134][i] *   HII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf.data[CollisionalRxnLUT::k136][i] *    DM(i,j,k) *    HI(i,j,k) /  2.
              );
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k129][i] *   HII(i,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k132][i] *    de(i,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k135][i] *    HM(i,j,k);
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            scoef = scoef
              + 2. * kdissHDI(i,j,k) * HDI(i,j,k)/3.0;
          }
        }
        species_tmpdens.data[SpLUT::DI][i]    = ( scoef*dtit[i] + DI(i,j,k) ) /
                    ( 1. + acoef*dtit[i] );

        // 2) DII
        scoef =   (   kcol_buf.data[CollisionalRxnLUT::k1][i]  * DI(i,j,k) * de(i,j,k)
              +       kcol_buf.data[CollisionalRxnLUT::k50][i] * HII(i,j,k)* DI(i,j,k)
              +  2.*kcol_buf.data[CollisionalRxnLUT::k53][i] * HII(i,j,k)* HDI(i,j,k)/3.
              )
              + kshield_buf.k24[i]*DI(i,j,k);
        acoef = 0.;
        // ! initialize GC202002
        if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + kphHI(i,j,k)*DI(i,j,k); }
        acoef =    kcol_buf.data[CollisionalRxnLUT::k2][i]  * de(i,j,k)
              +    kcol_buf.data[CollisionalRxnLUT::k51][i] * HI(i,j,k)
              +    kcol_buf.data[CollisionalRxnLUT::k52][i] * H2I(i,j,k)/2.;
        if (my_chemistry->primordial_chemistry > 3)  {
          acoef = acoef
              + kcol_buf.data[CollisionalRxnLUT::k130][i] *    HI(i,j,k)
              + kcol_buf.data[CollisionalRxnLUT::k133][i] *    DM(i,j,k) /  2.;
        }
        species_tmpdens.data[SpLUT::DII][i]   = ( scoef*dtit[i] + DII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // 3) HDI
        scoef = 3.*(kcol_buf.data[CollisionalRxnLUT::k52][i] * DII(i,j,k)*
             H2I(i,j,k)/2./2.
             + kcol_buf.data[CollisionalRxnLUT::k54][i] * DI(i,j,k) * H2I(i,j,k)/2./2.
        // !   &           + 2._DKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
        //- ! corrected by GC202005
             +          kcol_buf.data[CollisionalRxnLUT::k56][i] * DI(i,j,k) * HM(i,j,k)/2.
                   );
        acoef  =    kcol_buf.data[CollisionalRxnLUT::k53][i] * HII(i,j,k)
               +    kcol_buf.data[CollisionalRxnLUT::k55][i] * HI(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            acoef = acoef
              + kdissHDI(i,j,k);
          }
        }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  3. * ( 0.
              + kcol_buf.data[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) *    HI(i,j,k) /  3.
              + kcol_buf.data[CollisionalRxnLUT::k137][i] *    DM(i,j,k) *    HI(i,j,k) /  2.
              );
        }
        species_tmpdens.data[SpLUT::HDI][i]   = ( scoef*dtit[i] + HDI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

      }
    }
  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {
        
        // 1) DM
        scoef =
              kcol_buf.data[CollisionalRxnLUT::k132][i] *    DI(i,j,k) *    de(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k135][i] *    HM(i,j,k) *    DI(i,j,k);
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k133][i] *   DII(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::k134][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k136][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k137][i] *    HI(i,j,k);

        species_tmpdens.data[SpLUT::DM][i]    = ( scoef*dtit[i] + DM(i,j,k) ) /
                    ( 1. + acoef*dtit[i] );

        // 2) HDII
        scoef = 3. * (
              kcol_buf.data[CollisionalRxnLUT::k129][i] *    DI(i,j,k) *   HII(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::k130][i] *   DII(i,j,k) *    HI(i,j,k) /  2.
            );
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k125][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k131][i] *    de(i,j,k);

        species_tmpdens.data[SpLUT::HDII][i]   = ( scoef*dtit[i] + HDII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // 3) HeHII
        scoef = 5. * (
              kcol_buf.data[CollisionalRxnLUT::k148][i] *   HeI(i,j,k) *   HII(i,j,k) /  4.
            + kcol_buf.data[CollisionalRxnLUT::k149][i] *   HeI(i,j,k) *   HII(i,j,k) /  4.
            + kcol_buf.data[CollisionalRxnLUT::k150][i] *   HeI(i,j,k) *  H2II(i,j,k) /  8.
            + kcol_buf.data[CollisionalRxnLUT::k151][i] *  HeII(i,j,k) *    HI(i,j,k) /  4.
            );
        acoef =
              kcol_buf.data[CollisionalRxnLUT::k152][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::k153][i] *    de(i,j,k);

        species_tmpdens.data[SpLUT::HeHII][i]   = ( scoef*dtit[i] + HeHII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

      }
    }
  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE)  {

        // ***** CI **********
        scoef = 0. + 12. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz15][i] *    HI(i,j,k) *    CH(i,j,k) / 13.
            + kcol_buf.data[CollisionalRxnLUT::kz44][i] *   CII(i,j,k) *    de(i,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz20][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz27][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz28][i] *    OH(i,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i] *    O2(i,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz51][i] *   H2I(i,j,k) /  2.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust]      [i] / CI(i,j,k) * 12.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            acoef = acoef
              + kphCI(i,j,k);
          }
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            scoef = scoef + 12. *
                kdissCO (i,j,k) * CO(i,j,k) /28.0;
          }
        }

        species_tmpdens.data[SpLUT::CI][i]   = ( scoef*dtit[i] + CI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CII **********
        scoef = 0. + 12. * ( 0.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz37][i] *    OH(i,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i] *    O2(i,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz44][i] *    de(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphCI(i,j,k) * CI(i,j,k);
          }
        }

        species_tmpdens.data[SpLUT::CII][i]   = ( scoef*dtit[i] + CII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CO **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) *    OH(i,j,k) / 204.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) *    CH(i,j,k) / 208.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) *    HI(i,j,k) / 28.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i] / CO(i,j,k) * 17. * 0.5
            + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i] / CO(i,j,k) * 17.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissCO (i,j,k);
          }
        }

        species_tmpdens.data[SpLUT::CO][i]   = ( scoef*dtit[i] + CO(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CO2 **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) *    CO(i,j,k) / 476.
           );
        acoef = 0.;

        species_tmpdens.data[SpLUT::CO2][i]   = ( scoef*dtit[i] + CO2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OI **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz17][i] *    HI(i,j,k) *    OH(i,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i] *    HI(i,j,k) *    O2(i,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) *    OH(i,j,k) / 289.
            + kcol_buf.data[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf.data[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) *    HI(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz45][i] *   OII(i,j,k) *    de(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz50][i] *  O2II(i,j,k) *    de(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) *    O2(i,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz21][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz22][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz30][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz31][i] *    OI(i,j,k) / 8.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i] *    CH(i,j,k) / 13.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i] *    OH(i,j,k) / 17.;
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            acoef = acoef
              + kphOI(i,j,k);
          }
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            scoef = scoef + 16. *
              ( kdissOH (i,j,k) * OH(i,j,k) /17.0
              + kdissCO (i,j,k) * CO(i,j,k) /28.0);
          }
        }

        species_tmpdens.data[SpLUT::OI][i]   = ( scoef*dtit[i] + OI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OH **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz18][i] *    HI(i,j,k) *   H2O(i,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i] *    HI(i,j,k) *    O2(i,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) *   H2I(i,j,k) / 32.
            + kcol_buf.data[CollisionalRxnLUT::kz30][i] *    OI(i,j,k) *    HI(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz17][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) / 8.5
            + kcol_buf.data[CollisionalRxnLUT::kz26][i] *    CO(i,j,k) / 28.
            + kcol_buf.data[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz34][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) / 28.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) / 44.;
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissOH (i,j,k);
            scoef = scoef + 17. *
                kdissH2O(i,j,k) * H2O(i,j,k)/18.0;
          }
        }

        species_tmpdens.data[SpLUT::OH][i]   = ( scoef*dtit[i] + OH(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H2O **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) *    OH(i,j,k) / 34.
            + kcol_buf.data[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) *    OH(i,j,k) / 289.
            + kcol_buf.data[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz18][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz35][i] *   HII(i,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i] / H2O(i,j,k) * 18. * 2.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / H2O(i,j,k) * 18. * 3.
            + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i] / H2O(i,j,k) * 18. * 4.
            + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i] / H2O(i,j,k) * 18.
            + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i] / H2O(i,j,k) * 18. * 3.;
          }
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust]  [i] / H2O(i,j,k) * 18.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissH2O(i,j,k);
          }
        }

        species_tmpdens.data[SpLUT::H2O][i]   = ( scoef*dtit[i] + H2O(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** O2 **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz31][i] *    OI(i,j,k) *    OI(i,j,k) / 256.
            + kcol_buf.data[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) *    OH(i,j,k) / 272.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz19][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz36][i] *   HII(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) / 12.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) / 28.;

        species_tmpdens.data[SpLUT::O2][i]   = ( scoef*dtit[i] + O2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiI **********
        scoef = 0. + 28. * ( 0.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i] *    OH(i,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i] *    O2(i,j,k) / 32.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust]     [i] / SiI(i,j,k) * 28.;
          }
        }

        species_tmpdens.data[SpLUT::SiI][i]   = ( scoef*dtit[i] + SiI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiOI **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) *    OH(i,j,k) / 476.
            + kcol_buf.data[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) *    O2(i,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i] *    OH(i,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i] / SiOI(i,j,k) * 44.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / SiOI(i,j,k) * 44.;
          }
        }

        species_tmpdens.data[SpLUT::SiOI][i]   = ( scoef*dtit[i] + SiOI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiO2I **********
        scoef = 0. + 60. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) *    OH(i,j,k) / 748.
           );
        acoef = 0.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust]   [i] / SiO2I(i,j,k) * 60.;
          }
        }

        species_tmpdens.data[SpLUT::SiO2I][i]   = ( scoef*dtit[i] + SiO2I(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // MINOR BUT IMPORTANT SPECIES FOR MOLECULAR FORMATION
        //- ***** CH **********
        scoef = 0. + 13. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz16][i] *    HI(i,j,k) *   CH2(i,j,k) / 14.
            + kcol_buf.data[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) *   H2I(i,j,k) / 24.
            + kcol_buf.data[CollisionalRxnLUT::kz27][i] *    CI(i,j,k) *    HI(i,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz15][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) / 16.;

        species_tmpdens.data[SpLUT::CH][i]   = ( scoef*dtit[i] + CH(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CH2 **********
        scoef = 0. + 14. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) *    CH(i,j,k) / 26.
            + kcol_buf.data[CollisionalRxnLUT::kz51][i] *   H2I(i,j,k) *    CI(i,j,k) / 24.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz16][i] *    HI(i,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i] / CH2(i,j,k) * 14. * 0.5;
          }
        }

        species_tmpdens.data[SpLUT::CH2][i]   = ( scoef*dtit[i] + CH2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** COII **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) *    OH(i,j,k) / 204.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz43][i] *    HI(i,j,k);

        species_tmpdens.data[SpLUT::COII][i]   = ( scoef*dtit[i] + COII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OII **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz22][i] *   HII(i,j,k) *    OI(i,j,k) / 16.
            + kcol_buf.data[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) *    O2(i,j,k) / 384.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz39][i] *    HI(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz40][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz45][i] *    de(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphOI(i,j,k) * OI(i,j,k);
          }
        }

        species_tmpdens.data[SpLUT::OII][i]   = ( scoef*dtit[i] + OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OHII **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz34][i] *   HII(i,j,k) *    OH(i,j,k) / 17.
            + kcol_buf.data[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) *   H2I(i,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz41][i] *   H2I(i,j,k) /  2.;

        species_tmpdens.data[SpLUT::OHII][i]   = ( scoef*dtit[i] + OHII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H2OII **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz35][i] *   HII(i,j,k) *   H2O(i,j,k) / 18.
            + kcol_buf.data[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) *   H2I(i,j,k) / 34.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz42][i] *   H2I(i,j,k) /  2.
            + kcol_buf.data[CollisionalRxnLUT::kz46][i] *    de(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz47][i] *    de(i,j,k);

        species_tmpdens.data[SpLUT::H2OII][i]   = ( scoef*dtit[i] + H2OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H3OII **********
        scoef = 0. + 19. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) *   H2I(i,j,k) / 36.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz48][i] *    de(i,j,k)
            + kcol_buf.data[CollisionalRxnLUT::kz49][i] *    de(i,j,k);

        species_tmpdens.data[SpLUT::H3OII][i]   = ( scoef*dtit[i] + H3OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** O2II **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf.data[CollisionalRxnLUT::kz36][i] *   HII(i,j,k) *    O2(i,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf.data[CollisionalRxnLUT::kz50][i] *    de(i,j,k);

        species_tmpdens.data[SpLUT::O2II][i]   = ( scoef*dtit[i] + O2II(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            // ***** Mg **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i] / Mg(i,j,k) * 24.;
            if (my_chemistry->dust_species > 1)  {
              acoef = acoef
              + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / Mg(i,j,k) * 24. * 2.
              + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i] / Mg(i,j,k) * 24.;
            }

            species_tmpdens.data[SpLUT::Mg][i]   = ( scoef*dtit[i] + Mg(i,j,k) )
                       / ( 1. + acoef*dtit[i] );

          }

          if (my_chemistry->dust_species > 1)  {
            // ***** Al **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i] / Al(i,j,k) * 27. * 2.;

            species_tmpdens.data[SpLUT::Al][i]   = ( scoef*dtit[i] + Al(i,j,k) )
                       / ( 1. + acoef*dtit[i] );


            // ***** S  **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i] / S(i,j,k) * 32.;

            species_tmpdens.data[SpLUT::S][i]    = ( scoef*dtit[i] + S(i,j,k) )
                       / ( 1. + acoef*dtit[i] );


            // ***** Fe **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust]     [i] / Fe(i,j,k) * 56.
            + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i] / Fe(i,j,k) * 56. * 3.
            + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i] / Fe(i,j,k) * 56.;

            species_tmpdens.data[SpLUT::Fe][i]   = ( scoef*dtit[i] + Fe(i,j,k) )
                       / ( 1. + acoef*dtit[i] );

          }
        }

      }
    }
  }

  // --- (D4) Now do dust species ---
  if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE)  {

        if (my_chemistry->dust_species > 0)  {
          // ***** MgSiO3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 100.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::MgSiO3_dust][i]   = ( scoef*dtit[i] + MgSiO3(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** AC **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust]      [i] * 12.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::AC_dust][i]   = ( scoef*dtit[i] + AC(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 1)  {
          // ***** SiM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust]     [i] * 28.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::SiM_dust][i]   = ( scoef*dtit[i] + SiM(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** FeM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust]     [i] * 56.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::FeM_dust][i]   = ( scoef*dtit[i] + FeM(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Mg2SiO4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 140.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Mg2SiO4_dust][i]   = ( scoef*dtit[i] + Mg2SiO4(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Fe3O4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 232.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Fe3O4_dust][i]   = ( scoef*dtit[i] + Fe3O4(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** SiO2D **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust]   [i] * 60.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::SiO2_dust][i]   = ( scoef*dtit[i] + SiO2D(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** MgO **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust]     [i] * 40.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::MgO_dust][i]   = ( scoef*dtit[i] + MgO(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** FeS **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust]     [i] * 88.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::FeS_dust][i]   = ( scoef*dtit[i] + FeS(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Al2O3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust]   [i] * 102.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::Al2O3_dust][i]   = ( scoef*dtit[i] + Al2O3(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 2)  {
          // ***** reforg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust]  [i] * 22.68;
          acoef = 0.;

          species_tmpdens.data[SpLUT::ref_org_dust][i]   = ( scoef*dtit[i] + reforg(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** volorg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust]  [i] * 32.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::vol_org_dust][i]   = ( scoef*dtit[i] + volorg(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** H2Oice **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust]  [i] * 18.;
          acoef = 0.;

          species_tmpdens.data[SpLUT::H2O_ice_dust][i]   = ( scoef*dtit[i] + H2Oice(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

      }
    }
  }

  // --- (E) Set densities from 1D temps to 3D fields ---

  for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      HIdot_prev[i] = std::fabs(HI(i,j,k)-species_tmpdens.data[SpLUT::HI][i]) /
              std::fmax((double)(dtit[i] ), tiny8);
      HI(i,j,k)    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HI][i] ), tiny_fortran_val);
      HII(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HII][i] ), tiny_fortran_val);
      HeI(i,j,k)   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeI][i] ), tiny_fortran_val);
      HeII(i,j,k)  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeII][i] ), tiny_fortran_val);
      HeIII(i,j,k) = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeIII][i] ), (gr_float)(1e-5)*tiny_fortran_val);

      // de(i,j,k)    = dep(i)

      // Use charge conservation to determine electron fraction

      dedot_prev[i] = de(i,j,k);
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
        printf("HUGE HI! ::  %d %d %d %g\n",
               i+1,
               idx_range.jp1,
               idx_range.kp1,
               HI ( i, j, k ));
      }
    }

  }

  return;
}

}  // namespace grackle::impl

#endif /* STEP_RATE_GAUSS_SEIDEL_HPP */
