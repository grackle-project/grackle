//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// @brief Defines chemistry reaction related functions invoked by the
///     grackle solver in order to integrate the species densities over time.
///
/// The premise is for this file to hold functions used during the calculation
/// of species density derivatives during integration and to decouple this
/// logic from the integration to the greatest extent possible. This includes:
/// - `lookup_cool_rates1d_g`, when we eventually transcribe it
/// - the logic that is duplicated between `step_rate_g` and
///   `species_density_derivatives_0d`.
///   - the ordering of calculations required by `step_rate_g` proabably can't
///     be decoupled from the derivative calculation for primoridial species
///   - it may also make sense to further divide logic by the kinds of species
///     that are affected (e.g. primordial vs grains)
///
//===----------------------------------------------------------------------===//

#ifndef CHEMISTRY_SOLVER_FUNCS_HPP
#define CHEMISTRY_SOLVER_FUNCS_HPP

#include "grackle.h"
#include "fortran_func_decls.h"  // gr_mask_type
#include "full_rxn_rate_buf.hpp"
#include "index_helper.h"
#include "internal_types.hpp"
#include "LUT.hpp"

namespace grackle::impl::chemistry {

/// Perform an implicit Gauss-Seidel sweep of a backward-Euler time integrator
/// to advance the rate equations by one (sub-)cycle (for each index in the
/// index-range that is selected by the given itmask)
///
/// @param[out] out_spdens Holds buffers that are used to store the updated
///     species densities for the @p idx_range.
/// @param[in] idx_range Specifies the current index-range
/// @param[in] dtit Specifies the timestep of the current sub-cycle for each
///     index in @p idx_range.
/// @param[in] anydust Indicates whether we are modelling dust
/// @param[in] h2dust Specifies the rate of H2 dust-formation on dust grains
///     for eacg k
/// @param[in] rhoH Indicates the mass density of all Hydrogen
/// @param[in] itmask The general iteration mask for @p idx_range.
/// @param[in] itmask_metal The iteration mask @p idx_range that specifies
///     where we should account for metal-species chemistry and grain-species
///     chemistry.
/// @param[in] my_chemistry Provides various runtime parameters (we probably
///     don't need to pass the whole thing)
/// @param[in] my_fields Specifies the current values of the field data
/// @param[in] my_uvb_rates specifies precomputed rxn rates dependent on the
///     UV background, without accounting for self-shield (we probably don't
///     need to pass the whole thing since we also pass kshield_buf)
/// @param[in] kshield_buf specifies the
///     precomputed rxn rates (depends on local physical conditions)
/// @param[in] rxn_rate_buf specifies the precomputed rxn rates (depends on
///     local physical conditions)
///
/// Refactoring Goals
/// -----------------
///
/// Right now the implementation of every rate is **very** manual, which makes
/// it very easy to forget to add a rate.
///
/// The shorter-term goal is to deduplicate as much as possible with between
/// this function & grackle::impl::chemistry::species_density_derivatives_0d.
/// Some important considerations:
/// - it will be easiest to deduplicate the metal-chemistry and grain-species
///   growth rates (we only need to introduce some light usage of templates)
/// - deduplicating other logic will involve additional templates
/// - introducing light-weight templates to deduplicate the metal-chemistry
///   and grain-species growth rate logic is well worth the cost to improve
///   the maintenance burden. Deduplicating the other logic is probably also
///   worth the cost (but it's more debatable)
///
/// The longer-term goal is to overhaul this function.
/// - In the longer term, it would make more sense to use more of a "table
///   based approach" in the regime where we have lots of chemistry.
/// - We elaborate a little more down below:
///   - we already have a table of 1d rate buffers (i.e. @p kcol_buf ).
///   - we might also track a table (with matching indices) where we track the
///     indexes associated with each reactant/product of a rate as well as
///     stoichiometric coefficients.
///   - Furthermore, we could envision tracking rate-lists for each species of
///     all creational reactions and destructive reactions.
///   - If we have all of this, then we could just iterate over each species'
///     reactions and use the tables to dynamically access reactant densities
///     by index. Essentially, we would dynamically build up the creational
///     & destructive coefficients (``acoef`` & ``scoef``) for a full row and
///     use the values to then perform the update.
///   - in practice things would be a little more complex since we have lots
///     of rates not stored in @p kcol_buf, but it's still **very** doable
///   - once this infrastructure is in place, we could replace
///     grackle::impl::chemistry::species_density_derivatives_0d
///     that iterates through all of the rates in a very different manner
///     and directly compute partial derivatives in a much more efficient
///     manner (currently, it's used with finite derivatives for computing
///     these partial derivatives)
/// - Importantly, the table-based approach is *probably* slower than the
///   manual, hand-coded approach in the limit of a small chemical network.
///   I suspect that we'll preserve the hand-coded approach for
///   primordial_chemistry == 1, probably primordial_chemistry == 2, and
///   maybe even primordial_chemistry == 3 (but, we'll obviously need to test
///   performance).
inline void species_density_updates_gauss_seidel(
  grackle::impl::SpeciesCollection out_spdens, IndexRange idx_range,
  const double* dtit, gr_mask_type anydust, const double* h2dust,
  const double* rhoH, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  const FullRxnRateBuf rxn_rate_buf
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
  const double* const* kcol_buf = FullRxnRateBuf_kcol_bufs(&rxn_rate_buf);
  const double* const* grain_growth_rates =
    FullRxnRateBuf_grain_growth_bufs(&rxn_rate_buf);

  int i;
  double scoef, acoef;

  const int j = idx_range.j;
  const int k = idx_range.k;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI

        scoef  = kcol_buf[CollisionalRxnLUT::k2][i]*HII(i,j,k)*de(i,j,k);
        acoef  = kcol_buf[CollisionalRxnLUT::k1][i]*de(i,j,k)
               + kcol_buf[CollisionalRxnLUT::k57][i]*HI(i,j,k)
               + kcol_buf[CollisionalRxnLUT::k58][i]*HeI(i,j,k)/4.
               + kshield_buf.k24[i];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i,j,k); }
        out_spdens.data[SpLUT::HI][i]  = (scoef*dtit[i] + HI(i,j,k))/
             (1. + acoef*dtit[i]);
        if (out_spdens.data[SpLUT::HI][i] != out_spdens.data[SpLUT::HI][i])  {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens.data[SpLUT::HI] [ i ],
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
        scoef  = kcol_buf[CollisionalRxnLUT::k1][i]*out_spdens.data[SpLUT::HI][i]*de(i,j,k)
               + kcol_buf[CollisionalRxnLUT::k57][i]*out_spdens.data[SpLUT::HI][i]*out_spdens.data[SpLUT::HI][i]
               + kcol_buf[CollisionalRxnLUT::k58][i]*out_spdens.data[SpLUT::HI][i]*HeI(i,j,k)/4.
               + kshield_buf.k24[i]*out_spdens.data[SpLUT::HI][i];
        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i,j,k)*out_spdens.data[SpLUT::HI][i]; }
        acoef  = kcol_buf[CollisionalRxnLUT::k2][i]*de (i,j,k);
        out_spdens.data[SpLUT::HII][i] = (scoef*dtit[i] + HII(i,j,k))/
             (1. +acoef*dtit[i]);
        // 
        if (out_spdens.data[SpLUT::HII][i] <= 0.)   //#####
        {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("negative HIIp! ::  %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens.data[SpLUT::HII] [ i ],
                   scoef,
                   dtit [ i ],
                   HII ( i, j, k ),
                   acoef,
                   kcol_buf[CollisionalRxnLUT::k2] [ i ],
                   de ( i, j, k ),
                   kphHI ( i, j, k ),
                   out_spdens.data[SpLUT::HI] [ i ],
                   kshield_buf.k24 [ i ]);
          }
        }

        // 3) Electron density

        scoef = 0.
                   + kcol_buf[CollisionalRxnLUT::k57][i]*out_spdens.data[SpLUT::HI][i]*out_spdens.data[SpLUT::HI][i]
                   + kcol_buf[CollisionalRxnLUT::k58][i]*out_spdens.data[SpLUT::HI][i]*HeI(i,j,k)/4.
                   + kshield_buf.k24[i]*HI(i,j,k)
                   + kshield_buf.k25[i]*HeII(i,j,k)/4.
                   + kshield_buf.k26[i]*HeI(i,j,k)/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k)
                  + kphHeI(i,j,k)  * HeI(i,j,k)  / 4.
                  + kphHeII(i,j,k) * HeII(i,j,k) / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k); }



        acoef = -(kcol_buf[CollisionalRxnLUT::k1][i]*HI(i,j,k)      - kcol_buf[CollisionalRxnLUT::k2][i]*HII(i,j,k)
                + kcol_buf[CollisionalRxnLUT::k3][i]*HeI(i,j,k)/4. -
             kcol_buf[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)/4.
                + kcol_buf[CollisionalRxnLUT::k5][i]*HeII(i,j,k)/4. -
             kcol_buf[CollisionalRxnLUT::k4][i]*HeII(i,j,k)/4.);
        out_spdens.data[SpLUT::e][i]   = (scoef*dtit[i] + de(i,j,k))
                       / (1. + acoef*dtit[i]);

      }
    }

  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---

  for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {

      // 4) HeI

      scoef  = kcol_buf[CollisionalRxnLUT::k4][i]*HeII(i,j,k)*de(i,j,k);
      acoef  = kcol_buf[CollisionalRxnLUT::k3][i]*de(i,j,k)
                   + kshield_buf.k26[i];

      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeI(i,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        scoef = scoef +  4. * ( 0.
            + kcol_buf[CollisionalRxnLUT::k152][i] * HeHII(i,j,k) *    HI(i,j,k) /  5.
            + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) *    de(i,j,k) /  5.
            );
        acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k148][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k149][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k150][i] *  H2II(i,j,k) /  2.;
      }
      out_spdens.data[SpLUT::HeI][i]   = ( scoef*dtit[i] + HeI(i,j,k) )
                 / ( 1. + acoef*dtit[i] );

      // 5) HeII

      scoef  = kcol_buf[CollisionalRxnLUT::k3][i]*out_spdens.data[SpLUT::HeI][i]*de(i,j,k)
             + kcol_buf[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)*de(i,j,k)
             + kshield_buf.k26[i]*out_spdens.data[SpLUT::HeI][i];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeI(i,j,k)*out_spdens.data[SpLUT::HeI][i]; }

      acoef  = kcol_buf[CollisionalRxnLUT::k4][i]*de(i,j,k) + kcol_buf[CollisionalRxnLUT::k5][i]*de(i,j,k)
             + kshield_buf.k25[i];
     
      if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { acoef = acoef + kphHeII(i,j,k); }
      if (my_chemistry->primordial_chemistry > 3)  {
        acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k151][i] *    HI(i,j,k);
      }
      out_spdens.data[SpLUT::HeII][i]  = ( scoef*dtit[i] + HeII(i,j,k) )
                 / ( 1. + acoef*dtit[i] );

      // 6) HeIII

      scoef   = kcol_buf[CollisionalRxnLUT::k5][i]*out_spdens.data[SpLUT::HeII][i]*de(i,j,k)
              + kshield_buf.k25[i]*out_spdens.data[SpLUT::HeII][i];
      if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
          { scoef = scoef + kphHeII(i,j,k) * out_spdens.data[SpLUT::HeII][i]; }
      acoef   = kcol_buf[CollisionalRxnLUT::k6][i]*de(i,j,k);
      out_spdens.data[SpLUT::HeIII][i]  = ( scoef*dtit[i] + HeIII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

    }
  }

  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI
        scoef  =      kcol_buf[CollisionalRxnLUT::k2][i] * HII(i,j,k) * de(i,j,k)
               + 2.*kcol_buf[CollisionalRxnLUT::k13][i]* HI(i,j,k)  * H2I(i,j,k)/2.
               +      kcol_buf[CollisionalRxnLUT::k11][i]* HII(i,j,k) * H2I(i,j,k)/2.
               + 2.*kcol_buf[CollisionalRxnLUT::k12][i]* de(i,j,k)  * H2I(i,j,k)/2.
               +      kcol_buf[CollisionalRxnLUT::k14][i]* HM(i,j,k)  * de(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k15][i]* HM(i,j,k)  * HI(i,j,k)
               + 2.*kcol_buf[CollisionalRxnLUT::k16][i]* HM(i,j,k)  * HII(i,j,k)
               + 2.*kcol_buf[CollisionalRxnLUT::k18][i]* H2II(i,j,k)* de(i,j,k)/2.
               +      kcol_buf[CollisionalRxnLUT::k19][i]* H2II(i,j,k)* HM(i,j,k)/2.
               + 2.*kshield_buf.k31[i]   * H2I(i,j,k)/2.;

        acoef  =      kcol_buf[CollisionalRxnLUT::k1][i] * de(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k7][i] * de(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k8][i] * HM(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k9][i] * HII(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k10][i]* H2II(i,j,k)/2.
               + 2.*kcol_buf[CollisionalRxnLUT::k22][i]* std::pow(HI(i,j,k),2)
               +      kcol_buf[CollisionalRxnLUT::k57][i]* HI(i,j,k)
               +      kcol_buf[CollisionalRxnLUT::k58][i]* HeI(i,j,k)/4.
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
                + kcol_buf[CollisionalRxnLUT::k50][i] * HII(i,j,k) * DI(i,j,k)  / 2.
                + kcol_buf[CollisionalRxnLUT::k54][i] * H2I(i,j,k) * DI(i,j,k)  / 4.;
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k51][i] * DII(i,j,k) / 2.
                + kcol_buf[CollisionalRxnLUT::k55][i] * HDI(i,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) *    de(i,j,k) /  3.
              + kcol_buf[CollisionalRxnLUT::k134][i] *   HII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k135][i] *    HM(i,j,k) *    DI(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k150][i] *   HeI(i,j,k) *  H2II(i,j,k) /  8.
              + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) *    de(i,j,k) /  5.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) /  3.
              + kcol_buf[CollisionalRxnLUT::k130][i] *   DII(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k151][i] *  HeII(i,j,k) /  4.
              + kcol_buf[CollisionalRxnLUT::k152][i] * HeHII(i,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) *   H2I(i,j,k) / 24.
              + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) *   H2I(i,j,k) / 32.
              + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII(i,j,k) *    OI(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) *    CH(i,j,k) / 26.
              + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) *    OH(i,j,k) / 34.
              + kcol_buf[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) *    CO(i,j,k) / 476.
              + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) *    OH(i,j,k) / 204.
              + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) *    CH(i,j,k) / 208.
              + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) *    OH(i,j,k) / 272.
              + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII(i,j,k) *    OH(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII(i,j,k) *   H2O(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII(i,j,k) *    O2(i,j,k) / 32.
              + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) *    OH(i,j,k) / 204.
              + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) *   H2I(i,j,k) / 32.
              + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) *   H2I(i,j,k) / 34.
              + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) *   H2I(i,j,k) / 36.
              + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
              + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) *    de(i,j,k) / 9.5
              + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) *    OH(i,j,k) / 476.
              + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) *    OH(i,j,k) / 748.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz15][i] *    CH(i,j,k) / 13.
              + kcol_buf[CollisionalRxnLUT::kz16][i] *   CH2(i,j,k) / 14.
              + kcol_buf[CollisionalRxnLUT::kz17][i] *    OH(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz18][i] *   H2O(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz19][i] *    O2(i,j,k) / 32.
              + kcol_buf[CollisionalRxnLUT::kz27][i] *    CI(i,j,k) / 12.
              + kcol_buf[CollisionalRxnLUT::kz30][i] *    OI(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) / 28.;
        }
        out_spdens.data[SpLUT::HI][i]  = ( scoef*dtit[i] + HI(i,j,k) ) /
                        ( 1.f + acoef*dtit[i] );
        if (out_spdens.data[SpLUT::HI][i] != out_spdens.data[SpLUT::HI][i])  {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens.data[SpLUT::HI] [ i ],
                   HI ( i, j, k ),
                   HII ( i, j, k ),
                   de ( i, j, k ),
                   H2I ( i, j, k ),
                   kphHI ( i, j, k ));
          }
        }

        // 2) HII

        scoef  =    kcol_buf[CollisionalRxnLUT::k1][i]  * HI(i,j,k) * de(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k10][i] * H2II(i,j,k)*HI(i,j,k)/2.
               +    kcol_buf[CollisionalRxnLUT::k57][i] * HI(i,j,k) * HI(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k58][i] * HI(i,j,k) * HeI(i,j,k)/4.
               + kshield_buf.k24[i]*HI(i,j,k);

        if (my_chemistry->use_radiative_transfer == 1)
            { scoef = scoef + kphHI(i,j,k) * HI(i,j,k); }

        acoef  =    kcol_buf[CollisionalRxnLUT::k2][i]  * de(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k9][i]  * HI(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k11][i] * H2I(i,j,k)/2.
               +    kcol_buf[CollisionalRxnLUT::k16][i] * HM(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k17][i] * HM(i,j,k);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf[CollisionalRxnLUT::k51][i] * HI (i,j,k) * DII(i,j,k) / 2.
                + kcol_buf[CollisionalRxnLUT::k52][i] * H2I(i,j,k) * DII(i,j,k) / 4.;
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k50][i] * DI (i,j,k) / 2.
                + kcol_buf[CollisionalRxnLUT::k53][i] * HDI(i,j,k) / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) *    HI(i,j,k) /  3.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k129][i] *    DI(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k134][i] *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k148][i] *   HeI(i,j,k) /  4.
              + kcol_buf[CollisionalRxnLUT::k149][i] *   HeI(i,j,k) /  4.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) *    HI(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) *    HI(i,j,k) / 28.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz22][i] *    OI(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz34][i] *    OH(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz35][i] *   H2O(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz36][i] *    O2(i,j,k) / 32.;
        }
        out_spdens.data[SpLUT::HII][i]   = ( scoef*dtit[i] + HII(i,j,k) )
                        / ( 1. + acoef*dtit[i] );
        
        // 3) electrons:

        scoef =   kcol_buf[CollisionalRxnLUT::k8][i] * HM(i,j,k) * HI(i,j,k)
               +  kcol_buf[CollisionalRxnLUT::k15][i]* HM(i,j,k) * HI(i,j,k)
               +  kcol_buf[CollisionalRxnLUT::k17][i]* HM(i,j,k) * HII(i,j,k)
               +  kcol_buf[CollisionalRxnLUT::k57][i]* HI(i,j,k) * HI(i,j,k)
               +  kcol_buf[CollisionalRxnLUT::k58][i]* HI(i,j,k) * HeI(i,j,k)/4.
        // 
               + kshield_buf.k24[i]*out_spdens.data[SpLUT::HI][i]
               + kshield_buf.k25[i]*out_spdens.data[SpLUT::HeII][i]/4.
               + kshield_buf.k26[i]*out_spdens.data[SpLUT::HeI][i]/4.;

        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0) )
            { scoef = scoef + kphHI(i,j,k) * out_spdens.data[SpLUT::HI][i]
                  + kphHeI(i,j,k)  * out_spdens.data[SpLUT::HeI][i]  / 4.
                  + kphHeII(i,j,k) * out_spdens.data[SpLUT::HeII][i] / 4.; }
        if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 1) )
            { scoef = scoef + kphHI(i,j,k) * out_spdens.data[SpLUT::HI][i]; }
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

        acoef = - (kcol_buf[CollisionalRxnLUT::k1][i] *HI(i,j,k)    - kcol_buf[CollisionalRxnLUT::k2][i]*HII(i,j,k)
                +  kcol_buf[CollisionalRxnLUT::k3][i] *HeI(i,j,k)/4. -
             kcol_buf[CollisionalRxnLUT::k6][i]*HeIII(i,j,k)/4.
                +  kcol_buf[CollisionalRxnLUT::k5][i] *HeII(i,j,k)/4. -
             kcol_buf[CollisionalRxnLUT::k4][i]*HeII(i,j,k)/4.
                +  kcol_buf[CollisionalRxnLUT::k14][i]*HM(i,j,k)
                -  kcol_buf[CollisionalRxnLUT::k7][i] *HI(i,j,k)
                -  kcol_buf[CollisionalRxnLUT::k18][i]*H2II(i,j,k)/2.);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf[CollisionalRxnLUT::k56][i] * DI (i,j,k) * HM(i,j,k) / 2.;
          acoef = acoef
                - kcol_buf[CollisionalRxnLUT::k1] [i] * DI (i,j,k) / 2.
                + kcol_buf[CollisionalRxnLUT::k2] [i] * DII(i,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM(i,j,k) *    HI(i,j,k) /  2.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) /  3.
              + kcol_buf[CollisionalRxnLUT::k132][i] *    DI(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII(i,j,k) /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          // we comment out the following line that assigns scoef to itself since
          // it has no practical impact and produces a compiler warning
          // scoef = scoef;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz44][i] *   CII(i,j,k) / 12.
              + kcol_buf[CollisionalRxnLUT::kz45][i] *   OII(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) / 19.
              + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) / 19.
              + kcol_buf[CollisionalRxnLUT::kz50][i] *  O2II(i,j,k) / 32.;
        }
        out_spdens.data[SpLUT::e][i]  = ( scoef*dtit[i] + de(i,j,k) )
                  / ( 1. + acoef*dtit[i] );

        // 7) H2

        scoef = 2.*(kcol_buf[CollisionalRxnLUT::k8][i]  * HM(i,j,k)   * HI(i,j,k)
              +       kcol_buf[CollisionalRxnLUT::k10][i] * H2II(i,j,k) * HI(i,j,k)/2.
              +       kcol_buf[CollisionalRxnLUT::k19][i] * H2II(i,j,k) * HM(i,j,k)/2.
              +       kcol_buf[CollisionalRxnLUT::k22][i] * HI(i,j,k) * std::pow((HI(i,j,k)),2.));
        acoef = ( kcol_buf[CollisionalRxnLUT::k13][i]*HI(i,j,k) + kcol_buf[CollisionalRxnLUT::k11][i]*HII(i,j,k)
                + kcol_buf[CollisionalRxnLUT::k12][i]*de(i,j,k) )
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
                  kcol_buf[CollisionalRxnLUT::k53][i] * HDI(i,j,k) * HII(i,j,k) / 3.
                + kcol_buf[CollisionalRxnLUT::k55][i] * HDI(i,j,k) * HI (i,j,k) / 3.
                   );
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k52][i] * DII(i,j,k) / 2.
                + kcol_buf[CollisionalRxnLUT::k54][i] * DI (i,j,k) / 2.;
        }
#endif
        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI(i,j,k) *    CH(i,j,k) / 13.
              + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI(i,j,k) *   CH2(i,j,k) / 14.
              + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI(i,j,k) *    OH(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI(i,j,k) *   H2O(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
             );
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) / 12.
              + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz23][i] *    CH(i,j,k) / 13.
              + kcol_buf[CollisionalRxnLUT::kz24][i] *    OH(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) / 16.
              + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) / 17.
              + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) / 18.
              + kcol_buf[CollisionalRxnLUT::kz51][i] *    CI(i,j,k) / 12.;
          if ((my_chemistry->grain_growth == 1)  ||  (my_chemistry->dust_sublimation == 1))  {
            if (my_chemistry->dust_species > 0)  {
              scoef = scoef + 2. *
                    grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 2.;

            }
            if (my_chemistry->dust_species > 1)  {
              scoef = scoef + 2. * (
                    grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 3.
                  + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 4.
                  + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i]
                  + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] * 3.
                );
            }
            if (my_chemistry->dust_species > 2)  {
              acoef = acoef
              + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] / H2I(i,j,k) * 2. * 2.;
            }
          }
        }
        out_spdens.data[SpLUT::H2I][i] = ( scoef*dtit[i] + H2I(i,j,k) )
                  / ( 1. + acoef*dtit[i] );

        // 8) H-

        scoef = kcol_buf[CollisionalRxnLUT::k7][i] * HI(i,j,k) * de(i,j,k);
        acoef = (kcol_buf[CollisionalRxnLUT::k8][i]  + kcol_buf[CollisionalRxnLUT::k15][i])  * HI(i,j,k) +
                (kcol_buf[CollisionalRxnLUT::k16][i] + kcol_buf[CollisionalRxnLUT::k17][i])  * HII(i,j,k) +
               kcol_buf[CollisionalRxnLUT::k14][i] * de(i,j,k) + kcol_buf[CollisionalRxnLUT::k19][i] * H2II(i,j,k)/2.0f +
               my_uvb_rates.k27;
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k56][i] * DI (i,j,k) / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM(i,j,k) *    HI(i,j,k) /  2.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k135][i] *    DI(i,j,k) /  2.;
        }
        out_spdens.data[SpLUT::HM][i] = (scoef*dtit[i] + HM(i,j,k))
             / (1.0f + acoef*dtit[i]);


        // 9) H2+

        out_spdens.data[SpLUT::H2II][i] = 2.*( kcol_buf[CollisionalRxnLUT::k9] [i]*out_spdens.data[SpLUT::HI][i]*out_spdens.data[SpLUT::HII][i]
                      +   kcol_buf[CollisionalRxnLUT::k11][i]*out_spdens.data[SpLUT::H2I][i]/2.*out_spdens.data[SpLUT::HII][i]
                      +   kcol_buf[CollisionalRxnLUT::k17][i]*out_spdens.data[SpLUT::HM][i]*out_spdens.data[SpLUT::HII][i]
                      + kshield_buf.k29[i]*out_spdens.data[SpLUT::H2I][i]
                      )
                   /  ( kcol_buf[CollisionalRxnLUT::k10][i]*out_spdens.data[SpLUT::HI][i] + kcol_buf[CollisionalRxnLUT::k18][i]*out_spdens.data[SpLUT::e][i]
                      + kcol_buf[CollisionalRxnLUT::k19][i]*out_spdens.data[SpLUT::HM][i]
                      + (kshield_buf.k28[i]+kshield_buf.k30[i])
                      );
        if (my_chemistry->primordial_chemistry > 3)  {
          out_spdens.data[SpLUT::H2II][i] = 2. * (  kcol_buf[CollisionalRxnLUT::k9] [i]*out_spdens.data[SpLUT::HI][i]*out_spdens.data[SpLUT::HII][i]
                       +   kcol_buf[CollisionalRxnLUT::k11][i]*out_spdens.data[SpLUT::H2I][i]/2.*out_spdens.data[SpLUT::HII][i]
                       +   kcol_buf[CollisionalRxnLUT::k17][i]*out_spdens.data[SpLUT::HM][i]*out_spdens.data[SpLUT::HII][i]
                       + kshield_buf.k29[i]*out_spdens.data[SpLUT::H2I][i]
                       + kcol_buf[CollisionalRxnLUT::k152][i]*HeHII(i,j,k)*out_spdens.data[SpLUT::HI][i]/5.
                       )
                    /  ( kcol_buf[CollisionalRxnLUT::k10][i]*out_spdens.data[SpLUT::HI][i] + kcol_buf[CollisionalRxnLUT::k18][i]*out_spdens.data[SpLUT::e][i]
                       + kcol_buf[CollisionalRxnLUT::k19][i]*out_spdens.data[SpLUT::HM][i]
                       + (kshield_buf.k28[i]+kshield_buf.k30[i])
                       + kcol_buf[CollisionalRxnLUT::k150][i]*out_spdens.data[SpLUT::HeI][i]/4.
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
        scoef =   (       kcol_buf[CollisionalRxnLUT::k2][i] * DII(i,j,k) * de(i,j,k)
                   +      kcol_buf[CollisionalRxnLUT::k51][i]* DII(i,j,k) * HI(i,j,k)
                   + 2.*kcol_buf[CollisionalRxnLUT::k55][i]* HDI(i,j,k) *
                HI(i,j,k)/3.
                   );
        acoef  =    kcol_buf[CollisionalRxnLUT::k1][i] * de(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k50][i] * HII(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k54][i] * H2I(i,j,k)/2.
               +    kcol_buf[CollisionalRxnLUT::k56][i] * HM(i,j,k)
               + kshield_buf.k24[i];
        if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + kphHI(i,j,k); }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII(i,j,k) *    de(i,j,k) /  3.
              + kcol_buf[CollisionalRxnLUT::k133][i] *   DII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k134][i] *   HII(i,j,k) *    DM(i,j,k) /  2.
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM(i,j,k) *    HI(i,j,k) /  2.
              );
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k129][i] *   HII(i,j,k)
              + kcol_buf[CollisionalRxnLUT::k132][i] *    de(i,j,k)
              + kcol_buf[CollisionalRxnLUT::k135][i] *    HM(i,j,k);
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            scoef = scoef
              + 2. * kdissHDI(i,j,k) * HDI(i,j,k)/3.0;
          }
        }
        out_spdens.data[SpLUT::DI][i]    = ( scoef*dtit[i] + DI(i,j,k) ) /
                    ( 1. + acoef*dtit[i] );

        // 2) DII
        scoef =   (   kcol_buf[CollisionalRxnLUT::k1][i]  * DI(i,j,k) * de(i,j,k)
              +       kcol_buf[CollisionalRxnLUT::k50][i] * HII(i,j,k)* DI(i,j,k)
              +  2.*kcol_buf[CollisionalRxnLUT::k53][i] * HII(i,j,k)* HDI(i,j,k)/3.
              )
              + kshield_buf.k24[i]*DI(i,j,k);
        acoef = 0.;
        // ! initialize GC202002
        if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + kphHI(i,j,k)*DI(i,j,k); }
        acoef =    kcol_buf[CollisionalRxnLUT::k2][i]  * de(i,j,k)
              +    kcol_buf[CollisionalRxnLUT::k51][i] * HI(i,j,k)
              +    kcol_buf[CollisionalRxnLUT::k52][i] * H2I(i,j,k)/2.;
        if (my_chemistry->primordial_chemistry > 3)  {
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k130][i] *    HI(i,j,k)
              + kcol_buf[CollisionalRxnLUT::k133][i] *    DM(i,j,k) /  2.;
        }
        out_spdens.data[SpLUT::DII][i]   = ( scoef*dtit[i] + DII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // 3) HDI
        scoef = 3.*(kcol_buf[CollisionalRxnLUT::k52][i] * DII(i,j,k)*
             H2I(i,j,k)/2./2.
             + kcol_buf[CollisionalRxnLUT::k54][i] * DI(i,j,k) * H2I(i,j,k)/2./2.
        // !   &           + 2._DKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
        //- ! corrected by GC202005
             +          kcol_buf[CollisionalRxnLUT::k56][i] * DI(i,j,k) * HM(i,j,k)/2.
                   );
        acoef  =    kcol_buf[CollisionalRxnLUT::k53][i] * HII(i,j,k)
               +    kcol_buf[CollisionalRxnLUT::k55][i] * HI(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
            acoef = acoef
              + kdissHDI(i,j,k);
          }
        }
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  3. * ( 0.
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII(i,j,k) *    HI(i,j,k) /  3.
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM(i,j,k) *    HI(i,j,k) /  2.
              );
        }
        out_spdens.data[SpLUT::HDI][i]   = ( scoef*dtit[i] + HDI(i,j,k) )
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
              kcol_buf[CollisionalRxnLUT::k132][i] *    DI(i,j,k) *    de(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k135][i] *    HM(i,j,k) *    DI(i,j,k);
        acoef =
              kcol_buf[CollisionalRxnLUT::k133][i] *   DII(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::k134][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k136][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k137][i] *    HI(i,j,k);

        out_spdens.data[SpLUT::DM][i]    = ( scoef*dtit[i] + DM(i,j,k) ) /
                    ( 1. + acoef*dtit[i] );

        // 2) HDII
        scoef = 3. * (
              kcol_buf[CollisionalRxnLUT::k129][i] *    DI(i,j,k) *   HII(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::k130][i] *   DII(i,j,k) *    HI(i,j,k) /  2.
            );
        acoef =
              kcol_buf[CollisionalRxnLUT::k125][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k131][i] *    de(i,j,k);

        out_spdens.data[SpLUT::HDII][i]   = ( scoef*dtit[i] + HDII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // 3) HeHII
        scoef = 5. * (
              kcol_buf[CollisionalRxnLUT::k148][i] *   HeI(i,j,k) *   HII(i,j,k) /  4.
            + kcol_buf[CollisionalRxnLUT::k149][i] *   HeI(i,j,k) *   HII(i,j,k) /  4.
            + kcol_buf[CollisionalRxnLUT::k150][i] *   HeI(i,j,k) *  H2II(i,j,k) /  8.
            + kcol_buf[CollisionalRxnLUT::k151][i] *  HeII(i,j,k) *    HI(i,j,k) /  4.
            );
        acoef =
              kcol_buf[CollisionalRxnLUT::k152][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::k153][i] *    de(i,j,k);

        out_spdens.data[SpLUT::HeHII][i]   = ( scoef*dtit[i] + HeHII(i,j,k) )
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
            + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI(i,j,k) *    CH(i,j,k) / 13.
            + kcol_buf[CollisionalRxnLUT::kz44][i] *   CII(i,j,k) *    de(i,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz20][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz27][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    OH(i,j,k) / 17.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    O2(i,j,k) / 32.
            + kcol_buf[CollisionalRxnLUT::kz51][i] *   H2I(i,j,k) /  2.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::AC_dust]      [i] / CI(i,j,k) * 12.;
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

        out_spdens.data[SpLUT::CI][i]   = ( scoef*dtit[i] + CI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CII **********
        scoef = 0. + 12. * ( 0.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz37][i] *    OH(i,j,k) / 17.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *    O2(i,j,k) / 32.
            + kcol_buf[CollisionalRxnLUT::kz44][i] *    de(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphCI(i,j,k) * CI(i,j,k);
          }
        }

        out_spdens.data[SpLUT::CII][i]   = ( scoef*dtit[i] + CII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CO **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) *    OH(i,j,k) / 204.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) *    CH(i,j,k) / 208.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII(i,j,k) *    HI(i,j,k) / 28.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] / CO(i,j,k) * 17. * 0.5
            + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] / CO(i,j,k) * 17.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissCO (i,j,k);
          }
        }

        out_spdens.data[SpLUT::CO][i]   = ( scoef*dtit[i] + CO(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CO2 **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    OH(i,j,k) *    CO(i,j,k) / 476.
           );
        acoef = 0.;

        out_spdens.data[SpLUT::CO2][i]   = ( scoef*dtit[i] + CO2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OI **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI(i,j,k) *    OH(i,j,k) / 17.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI(i,j,k) *    O2(i,j,k) / 32.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) *    OH(i,j,k) / 289.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) *    O2(i,j,k) / 384.
            + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII(i,j,k) *    HI(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz45][i] *   OII(i,j,k) *    de(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
            + kcol_buf[CollisionalRxnLUT::kz50][i] *  O2II(i,j,k) *    de(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) *    O2(i,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz21][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz30][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz31][i] *    OI(i,j,k) / 8.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    CH(i,j,k) / 13.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OH(i,j,k) / 17.;
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

        out_spdens.data[SpLUT::OI][i]   = ( scoef*dtit[i] + OI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OH **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI(i,j,k) *   H2O(i,j,k) / 18.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI(i,j,k) *    O2(i,j,k) / 32.
            + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI(i,j,k) *   H2I(i,j,k) / 32.
            + kcol_buf[CollisionalRxnLUT::kz30][i] *    OI(i,j,k) *    HI(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII(i,j,k) *    de(i,j,k) / 18.
            + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) / 8.5
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    CO(i,j,k) / 28.
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI(i,j,k) / 12.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) / 12.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) / 28.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) / 44.;
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissOH (i,j,k);
            scoef = scoef + 17. *
                kdissH2O(i,j,k) * H2O(i,j,k)/18.0;
          }
        }

        out_spdens.data[SpLUT::OH][i]   = ( scoef*dtit[i] + OH(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H2O **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I(i,j,k) *    OH(i,j,k) / 34.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OH(i,j,k) *    OH(i,j,k) / 289.
            + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII(i,j,k) *    de(i,j,k) / 19.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII(i,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / H2O(i,j,k) * 18. * 2.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / H2O(i,j,k) * 18. * 3.
            + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] / H2O(i,j,k) * 18. * 4.
            + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] / H2O(i,j,k) * 18.
            + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] / H2O(i,j,k) * 18. * 3.;
          }
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust]  [i] / H2O(i,j,k) * 18.;
          }
        }
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
            acoef = acoef
              + kdissH2O(i,j,k);
          }
        }

        out_spdens.data[SpLUT::H2O][i]   = ( scoef*dtit[i] + H2O(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** O2 **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz31][i] *    OI(i,j,k) *    OI(i,j,k) / 256.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI(i,j,k) *    OH(i,j,k) / 272.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI(i,j,k) / 12.
            + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) / 12.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) / 28.;

        out_spdens.data[SpLUT::O2][i]   = ( scoef*dtit[i] + O2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiI **********
        scoef = 0. + 28. * ( 0.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *    OH(i,j,k) / 17.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *    O2(i,j,k) / 32.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::SiM_dust]     [i] / SiI(i,j,k) * 28.;
          }
        }

        out_spdens.data[SpLUT::SiI][i]   = ( scoef*dtit[i] + SiI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiOI **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI(i,j,k) *    OH(i,j,k) / 476.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI(i,j,k) *    O2(i,j,k) / 896.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *    OH(i,j,k) / 17.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / SiOI(i,j,k) * 44.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / SiOI(i,j,k) * 44.;
          }
        }

        out_spdens.data[SpLUT::SiOI][i]   = ( scoef*dtit[i] + SiOI(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** SiO2I **********
        scoef = 0. + 60. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI(i,j,k) *    OH(i,j,k) / 748.
           );
        acoef = 0.;
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust]   [i] / SiO2I(i,j,k) * 60.;
          }
        }

        out_spdens.data[SpLUT::SiO2I][i]   = ( scoef*dtit[i] + SiO2I(i,j,k) )
                   / ( 1. + acoef*dtit[i] );

        // MINOR BUT IMPORTANT SPECIES FOR MOLECULAR FORMATION
        //- ***** CH **********
        scoef = 0. + 13. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI(i,j,k) *   CH2(i,j,k) / 14.
            + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI(i,j,k) *   H2I(i,j,k) / 24.
            + kcol_buf[CollisionalRxnLUT::kz27][i] *    CI(i,j,k) *    HI(i,j,k) / 12.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI(i,j,k) / 16.;

        out_spdens.data[SpLUT::CH][i]   = ( scoef*dtit[i] + CH(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** CH2 **********
        scoef = 0. + 14. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I(i,j,k) *    CH(i,j,k) / 26.
            + kcol_buf[CollisionalRxnLUT::kz51][i] *   H2I(i,j,k) *    CI(i,j,k) / 24.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI(i,j,k);
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] / CH2(i,j,k) * 14. * 0.5;
          }
        }

        out_spdens.data[SpLUT::CH2][i]   = ( scoef*dtit[i] + CH2(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** COII **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII(i,j,k) *    OH(i,j,k) / 204.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz43][i] *    HI(i,j,k);

        out_spdens.data[SpLUT::COII][i]   = ( scoef*dtit[i] + COII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OII **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII(i,j,k) *    OI(i,j,k) / 16.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII(i,j,k) *    O2(i,j,k) / 384.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz39][i] *    HI(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz40][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz45][i] *    de(i,j,k);
        if (my_chemistry->use_radiative_transfer == 1)  {
          if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
            scoef = scoef
              + kphOI(i,j,k) * OI(i,j,k);
          }
        }

        out_spdens.data[SpLUT::OII][i]   = ( scoef*dtit[i] + OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** OHII **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII(i,j,k) *    OH(i,j,k) / 17.
            + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII(i,j,k) *   H2I(i,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz41][i] *   H2I(i,j,k) /  2.;

        out_spdens.data[SpLUT::OHII][i]   = ( scoef*dtit[i] + OHII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H2OII **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII(i,j,k) *   H2O(i,j,k) / 18.
            + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII(i,j,k) *   H2I(i,j,k) / 34.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz42][i] *   H2I(i,j,k) /  2.
            + kcol_buf[CollisionalRxnLUT::kz46][i] *    de(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz47][i] *    de(i,j,k);

        out_spdens.data[SpLUT::H2OII][i]   = ( scoef*dtit[i] + H2OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** H3OII **********
        scoef = 0. + 19. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII(i,j,k) *   H2I(i,j,k) / 36.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz48][i] *    de(i,j,k)
            + kcol_buf[CollisionalRxnLUT::kz49][i] *    de(i,j,k);

        out_spdens.data[SpLUT::H3OII][i]   = ( scoef*dtit[i] + H3OII(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        // ***** O2II **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII(i,j,k) *    O2(i,j,k) / 32.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz50][i] *    de(i,j,k);

        out_spdens.data[SpLUT::O2II][i]   = ( scoef*dtit[i] + O2II(i,j,k) )
                   / ( 1. + acoef*dtit[i] );


        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            // ***** Mg **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / Mg(i,j,k) * 24.;
            if (my_chemistry->dust_species > 1)  {
              acoef = acoef
              + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / Mg(i,j,k) * 24. * 2.
              + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] / Mg(i,j,k) * 24.;
            }

            out_spdens.data[SpLUT::Mg][i]   = ( scoef*dtit[i] + Mg(i,j,k) )
                       / ( 1. + acoef*dtit[i] );

          }

          if (my_chemistry->dust_species > 1)  {
            // ***** Al **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] / Al(i,j,k) * 27. * 2.;

            out_spdens.data[SpLUT::Al][i]   = ( scoef*dtit[i] + Al(i,j,k) )
                       / ( 1. + acoef*dtit[i] );


            // ***** S  **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] / S(i,j,k) * 32.;

            out_spdens.data[SpLUT::S][i]    = ( scoef*dtit[i] + S(i,j,k) )
                       / ( 1. + acoef*dtit[i] );


            // ***** Fe **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::FeM_dust]     [i] / Fe(i,j,k) * 56.
            + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] / Fe(i,j,k) * 56. * 3.
            + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] / Fe(i,j,k) * 56.;

            out_spdens.data[SpLUT::Fe][i]   = ( scoef*dtit[i] + Fe(i,j,k) )
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
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 100.;
          acoef = 0.;

          out_spdens.data[SpLUT::MgSiO3_dust][i]   = ( scoef*dtit[i] + MgSiO3(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** AC **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::AC_dust]      [i] * 12.;
          acoef = 0.;

          out_spdens.data[SpLUT::AC_dust][i]   = ( scoef*dtit[i] + AC(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 1)  {
          // ***** SiM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::SiM_dust]     [i] * 28.;
          acoef = 0.;

          out_spdens.data[SpLUT::SiM_dust][i]   = ( scoef*dtit[i] + SiM(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** FeM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::FeM_dust]     [i] * 56.;
          acoef = 0.;

          out_spdens.data[SpLUT::FeM_dust][i]   = ( scoef*dtit[i] + FeM(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Mg2SiO4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 140.;
          acoef = 0.;

          out_spdens.data[SpLUT::Mg2SiO4_dust][i]   = ( scoef*dtit[i] + Mg2SiO4(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Fe3O4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 232.;
          acoef = 0.;

          out_spdens.data[SpLUT::Fe3O4_dust][i]   = ( scoef*dtit[i] + Fe3O4(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** SiO2D **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust]   [i] * 60.;
          acoef = 0.;

          out_spdens.data[SpLUT::SiO2_dust][i]   = ( scoef*dtit[i] + SiO2D(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** MgO **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] * 40.;
          acoef = 0.;

          out_spdens.data[SpLUT::MgO_dust][i]   = ( scoef*dtit[i] + MgO(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** FeS **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] * 88.;
          acoef = 0.;

          out_spdens.data[SpLUT::FeS_dust][i]   = ( scoef*dtit[i] + FeS(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** Al2O3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] * 102.;
          acoef = 0.;

          out_spdens.data[SpLUT::Al2O3_dust][i]   = ( scoef*dtit[i] + Al2O3(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 2)  {
          // ***** reforg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] * 22.68;
          acoef = 0.;

          out_spdens.data[SpLUT::ref_org_dust][i]   = ( scoef*dtit[i] + reforg(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** volorg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] * 32.;
          acoef = 0.;

          out_spdens.data[SpLUT::vol_org_dust][i]   = ( scoef*dtit[i] + volorg(i,j,k) )
                     / ( 1. + acoef*dtit[i] );


          // ***** H2Oice **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust]  [i] * 18.;
          acoef = 0.;

          out_spdens.data[SpLUT::H2O_ice_dust][i]   = ( scoef*dtit[i] + H2Oice(i,j,k) )
                     / ( 1. + acoef*dtit[i] );

        }

      }
    }
  }
}



/// Computes the time derivatives of the mass densities of various species
///
/// This function is specialized to operate on a single cell at a time. This is
/// mostly for historical reasons. In reality, it would be faster to operate on
/// many species at once.
///
/// @param[out] deriv This is where we store the time derivative of the
///     species mass densities. The caller must make sure this is filled with
///     zeros before calling this function
/// @param[in]  anydust
/// @param[in]  h2dust
/// @param[in]  rhoH Indicates the mass density of all Hydrogen
/// @param[in]  itmask_metal Indicates where we should account for metal
///     chemistry
/// @param[in]  my_chemistry Provides various runtime parameters (we probably
///     don't need to pass the whole thing)
/// @param[in]  my_fields Specifies field data
/// @param[in]  my_uvb_rates specifies precomputed rxn rates dependent on the
///     UV background
/// @param[in]  kshield_buf specifies the
///     precomputed rxn rates (depends on local physical conditions)
/// @param[in] rxn_rate_buf specifies the precomputed rxn rates (depends on
///     local physical conditions)
///
/// @note
/// Some of occurences of the const keyword is somewhat symbolic (to indicate
/// that we won't modify an arg or the contents of a struct). In reality,
/// C/C++ allows you to modify pointers in a `const struct`
///
/// @todo
/// We should consider how to deduplicate this logic and most of the logic in
/// `step_rate_g`
/// -> I am pretty sure Gen performed a copy-and-paste when he wrote this
/// -> it is totally unsustainable for us to maintain both implementations
/// -> A common obstacle for unifying this logic is:
///    - step_rate_g uses these calculations with a backward difference
///      formula (see 3.2 of the grackle paper)
///    - Here, we're computing the net time derivative.
///    We can DEFINITELY overcome this obstacle!
/// -> there's another obstacle that applies for unifying another part of
///    the logic:
///    - Parts (A), (B), and (C) in the original `step_rate_g` all apply
///      partial updates (again, see 3.2 of the grackle paper)
///    - notably, part (D) does not do this (but I'm not sure that there
///      is a reasonable motivation for that choice)
inline void species_density_derivatives_0d(
  grackle::impl::SpeciesCollection deriv, gr_mask_type anydust,
  const double* h2dust, const double* rhoH, const gr_mask_type* itmask_metal,
  const chemistry_data* my_chemistry, const grackle_field_data* my_fields,
  const photo_rate_storage my_uvb_rates,
  const grackle::impl::PhotoRxnRateCollection kshield_buf,
  const FullRxnRateBuf rxn_rate_buf
) {

  // define some local variables carried over from the fortran version:
  // - the goal is to eventually remove all of these
  // - originally we only conditionally defined some of these variables, but
  //   there honestly isn't any benefit to doing that (the memory is allocated
  //   already)
  // - originally, each of these variables was a stack allocated variable that
  //   held a copy of the corresponding entry in dsp. Now, these variables are
  //   references to casted copies taken from dsp
  // - we aren't being that consistent with the historic implementation when
  //   gr_float isn't the same as double (but, the historical implementation
  //   definitely didn't handle this case properly)

  gr_float& de      = my_fields->e_density[0];
  gr_float& HI      = my_fields->HI_density[0];
  gr_float& HII     = my_fields->HII_density[0];
  gr_float& HeI     = my_fields->HeI_density[0];
  gr_float& HeII    = my_fields->HeII_density[0];
  gr_float& HeIII   = my_fields->HeIII_density[0];
  gr_float& HM      = my_fields->HM_density[0];
  gr_float& H2I     = my_fields->H2I_density[0];
  gr_float& H2II    = my_fields->H2II_density[0];
  gr_float& DI      = my_fields->DI_density[0];
  gr_float& DII     = my_fields->DII_density[0];
  gr_float& HDI     = my_fields->HDI_density[0];
  gr_float& DM      = my_fields->DM_density[0];
  gr_float& HDII    = my_fields->HDII_density[0];
  gr_float& HeHII   = my_fields->HeHII_density[0];
  gr_float& CI      = my_fields->CI_density[0];
  gr_float& CII     = my_fields->CII_density[0];
  gr_float& CO      = my_fields->CO_density[0];
  gr_float& CO2     = my_fields->CO2_density[0];
  gr_float& OI      = my_fields->OI_density[0];
  gr_float& OH      = my_fields->OH_density[0];
  gr_float& H2O     = my_fields->H2O_density[0];
  gr_float& O2      = my_fields->O2_density[0];
  gr_float& SiI     = my_fields->SiI_density[0];
  gr_float& SiOI    = my_fields->SiOI_density[0];
  gr_float& SiO2I   = my_fields->SiO2I_density[0];
  gr_float& CH      = my_fields->CH_density[0];
  gr_float& CH2     = my_fields->CH2_density[0];
  gr_float& COII    = my_fields->COII_density[0];
  gr_float& OII     = my_fields->OII_density[0];
  gr_float& OHII    = my_fields->OHII_density[0];
  gr_float& H2OII   = my_fields->H2OII_density[0];
  gr_float& H3OII   = my_fields->H3OII_density[0];
  gr_float& O2II    = my_fields->O2II_density[0];
  gr_float& Mg      = my_fields->Mg_density[0];
  gr_float& Al      = my_fields->Al_density[0];
  gr_float& S       = my_fields->S_density[0];
  gr_float& Fe      = my_fields->Fe_density[0];
  gr_float& MgSiO3  = my_fields->MgSiO3_dust_density[0];
  gr_float& AC      = my_fields->AC_dust_density[0];
  gr_float& SiM     = my_fields->SiM_dust_density[0];
  gr_float& FeM     = my_fields->FeM_dust_density[0];
  gr_float& Mg2SiO4 = my_fields->Mg2SiO4_dust_density[0];
  gr_float& Fe3O4   = my_fields->Fe3O4_dust_density[0];
  gr_float& SiO2D   = my_fields->SiO2_dust_density[0];
  gr_float& MgO     = my_fields->MgO_dust_density[0];
  gr_float& FeS     = my_fields->FeS_dust_density[0];
  gr_float& Al2O3   = my_fields->Al2O3_dust_density[0];
  gr_float& reforg  = my_fields->ref_org_dust_density[0];
  gr_float& volorg  = my_fields->vol_org_dust_density[0];
  gr_float& H2Oice  = my_fields->H2O_ice_dust_density[0];

  const double* const* kcol_buf = FullRxnRateBuf_kcol_bufs(&rxn_rate_buf);
  const double* const* grain_growth_rates =
    FullRxnRateBuf_grain_growth_bufs(&rxn_rate_buf);

  double scoef, acoef;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {




    // 1) HI

    scoef  = kcol_buf[CollisionalRxnLUT::k2][0]   *HII       *de;
    acoef  = kcol_buf[CollisionalRxnLUT::k1][0]   *de
           + kcol_buf[CollisionalRxnLUT::k57][0]   *HI
           + kcol_buf[CollisionalRxnLUT::k58][0]   *HeI       /4.
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(my_fields->RT_HI_ionization_rate); }
    deriv.data[SpLUT::HI][0] = deriv.data[SpLUT::HI][0] + (scoef - acoef * HI);














    // 2) HII
    scoef  = kcol_buf[CollisionalRxnLUT::k1][0]   *HI    *de
           + kcol_buf[CollisionalRxnLUT::k57][0]   *HI    *HI
           + kcol_buf[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
           + kshield_buf.k24[0]   *HI;
    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)       *HI; }
    acoef  = kcol_buf[CollisionalRxnLUT::k2][0]   *de;
    deriv.data[SpLUT::HII][0] = deriv.data[SpLUT::HII][0] + (scoef - acoef * HII);
















    // 3) Electron density

    scoef = 0.
               + kcol_buf[CollisionalRxnLUT::k57][0]   *HI    *HI
               + kcol_buf[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
               + kshield_buf.k24[0]   *HI
               + kshield_buf.k25[0]   *HeII       /4.
               + kshield_buf.k26[0]   *HeI       /4.;

    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI
              + *(my_fields->RT_HeI_ionization_rate)         * HeI         / 4.
              + *(my_fields->RT_HeII_ionization_rate)        * HeII        / 4.; }
    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI; }



    acoef = -(kcol_buf[CollisionalRxnLUT::k1][0]   *HI             - kcol_buf[CollisionalRxnLUT::k2][0]   *HII
            + kcol_buf[CollisionalRxnLUT::k3][0]   *HeI       /4. -
         kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            + kcol_buf[CollisionalRxnLUT::k5][0]   *HeII       /4. -
         kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       /4.);
    deriv.data[SpLUT::e][0] = deriv.data[SpLUT::e][0] + (scoef - acoef * de);





  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---




  // 4) HeI

  scoef  = kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       *de;
  acoef  = kcol_buf[CollisionalRxnLUT::k3][0]   *de
               + kshield_buf.k26[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(my_fields->RT_HeI_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    scoef = scoef +  4. * ( 0.
        + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
        + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.
       );
    acoef = acoef
        + kcol_buf[CollisionalRxnLUT::k148][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k149][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k150][0]    *  H2II        /  2.;
  }
  deriv.data[SpLUT::HeI][0] = deriv.data[SpLUT::HeI][0] + (scoef - acoef * HeI);


  // 5) HeII

  scoef  = kcol_buf[CollisionalRxnLUT::k3][0]   *HeI    *de
         + kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       *de
         + kshield_buf.k26[0]   *HeI;

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(my_fields->RT_HeI_ionization_rate)       *HeI; }

  acoef  = kcol_buf[CollisionalRxnLUT::k4][0]   *de        + kcol_buf[CollisionalRxnLUT::k5][0]   *de
         + kshield_buf.k25[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(my_fields->RT_HeII_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    acoef = acoef
        + kcol_buf[CollisionalRxnLUT::k151][0]    *    HI;
  }
  deriv.data[SpLUT::HeII][0] = deriv.data[SpLUT::HeII][0] + (scoef - acoef * HeII);


  // 6) HeIII

  scoef   = kcol_buf[CollisionalRxnLUT::k5][0]   *HeII    *de
          + kshield_buf.k25[0]   *HeII;
  if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(my_fields->RT_HeII_ionization_rate)        * HeII; }
  acoef   = kcol_buf[CollisionalRxnLUT::k6][0]   *de;
  deriv.data[SpLUT::HeIII][0] = deriv.data[SpLUT::HeIII][0] + (scoef - acoef * HeIII);





  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms




    // 1) HI
    scoef  =      kcol_buf[CollisionalRxnLUT::k2][0]    * HII        * de
           + 2.*kcol_buf[CollisionalRxnLUT::k13][0]   * HI         * H2I       /2.
           +      kcol_buf[CollisionalRxnLUT::k11][0]   * HII        * H2I       /2.
           + 2.*kcol_buf[CollisionalRxnLUT::k12][0]   * de         * H2I       /2.
           +      kcol_buf[CollisionalRxnLUT::k14][0]   * HM         * de
           +      kcol_buf[CollisionalRxnLUT::k15][0]   * HM         * HI
           + 2.*kcol_buf[CollisionalRxnLUT::k16][0]   * HM         * HII
           + 2.*kcol_buf[CollisionalRxnLUT::k18][0]   * H2II       * de       /2.
           +      kcol_buf[CollisionalRxnLUT::k19][0]   * H2II       * HM       /2.
           + 2.*kshield_buf.k31[0]      * H2I       /2.;

    acoef  =      kcol_buf[CollisionalRxnLUT::k1][0]    * de
           +      kcol_buf[CollisionalRxnLUT::k7][0]    * de
           +      kcol_buf[CollisionalRxnLUT::k8][0]    * HM
           +      kcol_buf[CollisionalRxnLUT::k9][0]    * HII
           +      kcol_buf[CollisionalRxnLUT::k10][0]   * H2II       /2.
           + 2.*kcol_buf[CollisionalRxnLUT::k22][0]   * std::pow(HI       ,2)
           +      kcol_buf[CollisionalRxnLUT::k57][0]   * HI
           +      kcol_buf[CollisionalRxnLUT::k58][0]   * HeI       /4.
           + kshield_buf.k24[0];

    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(my_fields->RT_HI_ionization_rate); }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if ((my_chemistry->primordial_chemistry > 2) && (my_chemistry->radiative_transfer_HDI_dissociation > 0))  {
        scoef = scoef
          + *(my_fields->RT_HDI_dissociation_rate)        * HDI       /3.0;
      }
      if ((my_chemistry->metal_chemistry == 1)  && 
          (itmask_metal[0] != MASK_FALSE))  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef
            + *(my_fields->RT_OH_dissociation_rate)         * OH        /17.0
            + *(my_fields->RT_H2O_dissociation_rate)        * H2O       /18.0;
        }
      }
    }

    if (anydust != MASK_FALSE)  {
      if(itmask_metal[0] != MASK_FALSE   )  {
        acoef = acoef + 2. * h2dust[0]    * rhoH[0];
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k50][0]    * HII        * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k54][0]    * H2I        * DI         / 4.;
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k51][0]    * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k55][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM        *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
          + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        /  3.
          + kcol_buf[CollisionalRxnLUT::k130][0]    *   DII        /  2.
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k151][0]    *  HeII        /  4.
          + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *    OH        *    CO        / 476.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        *    OH        / 272.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII        *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        *    OH        / 204.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 9.5
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    CH        / 13.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *   CH2        / 14.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        / 28.;
    }
    deriv.data[SpLUT::HI][0] = deriv.data[SpLUT::HI][0] + (scoef - acoef * HI);













    // 2) HII

    scoef  =    kcol_buf[CollisionalRxnLUT::k1][0]     * HI        * de
           +    kcol_buf[CollisionalRxnLUT::k10][0]    * H2II       *HI       /2.
           +    kcol_buf[CollisionalRxnLUT::k57][0]    * HI        * HI
           +    kcol_buf[CollisionalRxnLUT::k58][0]    * HI        * HeI       /4.
           + kshield_buf.k24[0]   *HI;

    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI; }

    acoef  =    kcol_buf[CollisionalRxnLUT::k2][0]     * de
           +    kcol_buf[CollisionalRxnLUT::k9][0]     * HI
           +    kcol_buf[CollisionalRxnLUT::k11][0]    * H2I       /2.
           +    kcol_buf[CollisionalRxnLUT::k16][0]    * HM
           +    kcol_buf[CollisionalRxnLUT::k17][0]    * HM;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k51][0]    * HI         * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k52][0]    * H2I        * DII        / 4.;
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k50][0]    * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k53][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k129][0]    *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k148][0]    *   HeI        /  4.
          + kcol_buf[CollisionalRxnLUT::k149][0]    *   HeI        /  4.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *    O2        / 32.;
    }
    deriv.data[SpLUT::HII][0] = deriv.data[SpLUT::HII][0] + (scoef - acoef * HII);

    
    // 3) electrons:

    scoef =   kcol_buf[CollisionalRxnLUT::k8][0]    * HM        * HI
           +  kcol_buf[CollisionalRxnLUT::k15][0]   * HM        * HI
           +  kcol_buf[CollisionalRxnLUT::k17][0]   * HM        * HII
           +  kcol_buf[CollisionalRxnLUT::k57][0]   * HI        * HI
           +  kcol_buf[CollisionalRxnLUT::k58][0]   * HI        * HeI       /4.
    // 
           + kshield_buf.k24[0]   *HI
           + kshield_buf.k25[0]   *HeII    /4.
           + kshield_buf.k26[0]   *HeI    /4.;

    if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI
              + *(my_fields->RT_HeI_ionization_rate)         * HeI      / 4.
              + *(my_fields->RT_HeII_ionization_rate)        * HeII     / 4.; }
    if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 1) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI; }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if ((my_chemistry->metal_chemistry == 1)  && 
          (itmask_metal[0] != MASK_FALSE))  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(my_fields->RT_CI_ionization_rate)        * CI       /12.0
            + *(my_fields->RT_OI_ionization_rate)        * OI       /16.0;
        }
      }
    }

    acoef = - (kcol_buf[CollisionalRxnLUT::k1][0]    *HI           - kcol_buf[CollisionalRxnLUT::k2][0]   *HII
            +  kcol_buf[CollisionalRxnLUT::k3][0]    *HeI       /4. -
         kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            +  kcol_buf[CollisionalRxnLUT::k5][0]    *HeII       /4. -
         kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       /4.
            +  kcol_buf[CollisionalRxnLUT::k14][0]   *HM
            -  kcol_buf[CollisionalRxnLUT::k7][0]    *HI
            -  kcol_buf[CollisionalRxnLUT::k18][0]   *H2II       /2.);
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k56][0]    * DI         * HM        / 2.;
      acoef = acoef
            - kcol_buf[CollisionalRxnLUT::k1][0]     * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k2][0]     * DII        / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        /  3.
          + kcol_buf[CollisionalRxnLUT::k132][0]    *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      // we comment out the following line that assigns scoef to itself since
      // it has no practical impact and produces a compiler warning
      // scoef = scoef;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        / 19.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        / 19.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *  O2II        / 32.;
    }
    deriv.data[SpLUT::e][0] = deriv.data[SpLUT::e][0] + (scoef - acoef * de);


    // 7) H2

    scoef = 2.*(kcol_buf[CollisionalRxnLUT::k8][0]     * HM          * HI
          +       kcol_buf[CollisionalRxnLUT::k10][0]    * H2II        * HI       /2.
          +       kcol_buf[CollisionalRxnLUT::k19][0]    * H2II        * HM       /2.
          +       kcol_buf[CollisionalRxnLUT::k22][0]    * HI        * std::pow((HI       ),2.));
    acoef = ( kcol_buf[CollisionalRxnLUT::k13][0]   *HI        + kcol_buf[CollisionalRxnLUT::k11][0]   *HII
            + kcol_buf[CollisionalRxnLUT::k12][0]   *de        )
            + kshield_buf.k29[0]    + kshield_buf.k31[0];

    if (anydust != MASK_FALSE)  {
      if(itmask_metal[0] != MASK_FALSE   )  {
        scoef = scoef + 2. * h2dust[0]    *
             HI        * rhoH[0];
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef + 2. * (
              kcol_buf[CollisionalRxnLUT::k53][0]    * HDI        * HII        / 3.
            + kcol_buf[CollisionalRxnLUT::k55][0]    * HDI        * HI         / 3.
               );
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k52][0]    * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k54][0]    * DI         / 2.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
         );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *    CH        / 13.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        / 17.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *    CI        / 12.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          scoef = scoef + 2. *
                grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      * 2.;

        }
        if (my_chemistry->dust_species > 1)  {
          scoef = scoef + 2. * (
                grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 3.
              + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       * 4.
              + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]
              + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       * 3.
            );
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      / H2I        * 2. * 2.;
        }
      }
    }
    deriv.data[SpLUT::H2I][0] = deriv.data[SpLUT::H2I][0] + (scoef - acoef * H2I);


    // 8) H-

    scoef = kcol_buf[CollisionalRxnLUT::k7][0]    * HI        * de;
    acoef = (kcol_buf[CollisionalRxnLUT::k8][0]     + kcol_buf[CollisionalRxnLUT::k15][0]   )  * HI        +
            (kcol_buf[CollisionalRxnLUT::k16][0]    + kcol_buf[CollisionalRxnLUT::k17][0]   )  * HII        +
            kcol_buf[CollisionalRxnLUT::k14][0]    * de        + kcol_buf[CollisionalRxnLUT::k19][0]    * H2II       /2.0f +
            my_uvb_rates.k27;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k56][0]    * DI         / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    DI        /  2.;
    }
    deriv.data[SpLUT::HM][0] = deriv.data[SpLUT::HM][0] + (scoef - acoef * HM);



    // 9) H2+

    scoef =    2.*( kcol_buf[CollisionalRxnLUT::k9][0]    *HI    *HII
                  +   kcol_buf[CollisionalRxnLUT::k11][0]   *H2I    /2.*HII
                  +   kcol_buf[CollisionalRxnLUT::k17][0]   *HM    *HII
                  + kshield_buf.k29[0]   *H2I    /2.
                  );
    acoef =         kcol_buf[CollisionalRxnLUT::k10][0]   *HI     + kcol_buf[CollisionalRxnLUT::k18][0]   *de
                  + kcol_buf[CollisionalRxnLUT::k19][0]   *HM
                  + (kshield_buf.k28[0]   +kshield_buf.k30[0]   );
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
         );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        /  4.;
    }
    deriv.data[SpLUT::H2II][0] = deriv.data[SpLUT::H2II][0] + (scoef - acoef * H2II);









  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {


    
    // 1) DI
    scoef =   (       kcol_buf[CollisionalRxnLUT::k2][0]    * DII        * de
               +      kcol_buf[CollisionalRxnLUT::k51][0]   * DII        * HI
               + 2.*kcol_buf[CollisionalRxnLUT::k55][0]   * HDI        *
            HI       /3.
               );
    acoef  =    kcol_buf[CollisionalRxnLUT::k1][0]    * de
           +    kcol_buf[CollisionalRxnLUT::k50][0]    * HII
           +    kcol_buf[CollisionalRxnLUT::k54][0]    * H2I       /2.
           +    kcol_buf[CollisionalRxnLUT::k56][0]    * HM
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(my_fields->RT_HI_ionization_rate); }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcol_buf[CollisionalRxnLUT::k133][0]    *   DII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.
          );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k129][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::k132][0]    *    de
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM;
    }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        scoef = scoef
          + 2. * *(my_fields->RT_HDI_dissociation_rate)        * HDI       /3.0;
      }
    }
    deriv.data[SpLUT::DI][0] = deriv.data[SpLUT::DI][0] + (scoef - acoef * DI);
                                                    

    // 2) DII
    scoef =   (   kcol_buf[CollisionalRxnLUT::k1][0]     * DI        * de
          +       kcol_buf[CollisionalRxnLUT::k50][0]    * HII       * DI
          +  2.*kcol_buf[CollisionalRxnLUT::k53][0]    * HII       * HDI       /3.
          )
          + kshield_buf.k24[0]   *DI;
    acoef = 0.;
    // ! initialize GC202002
    if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + *(my_fields->RT_HI_ionization_rate)       *DI; }
    acoef =    kcol_buf[CollisionalRxnLUT::k2][0]     * de
          +    kcol_buf[CollisionalRxnLUT::k51][0]    * HI
          +    kcol_buf[CollisionalRxnLUT::k52][0]    * H2I       /2.;
    if (my_chemistry->primordial_chemistry > 3)  {
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k130][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::k133][0]    *    DM        /  2.;
    }
    deriv.data[SpLUT::DII][0] = deriv.data[SpLUT::DII][0] + (scoef - acoef * DII);


    // 3) HDI
    scoef = 3.*(kcol_buf[CollisionalRxnLUT::k52][0]    * DII       *
         H2I       /2./2.
         + kcol_buf[CollisionalRxnLUT::k54][0]    * DI        * H2I       /2./2.
    // !   &           + 2._DKIND*kcol_buf[CollisionalRxnLUT::k56][0]    * DI        * HM       /2._DKIND
    //- ! corrected by GC202005
         +          kcol_buf[CollisionalRxnLUT::k56][0]    * DI        * HM       /2.
               );
    acoef  =    kcol_buf[CollisionalRxnLUT::k53][0]    * HII
           +    kcol_buf[CollisionalRxnLUT::k55][0]    * HI;
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        acoef = acoef
          + *(my_fields->RT_HDI_dissociation_rate);
      }
    }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  3. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.
          );
    }
    deriv.data[SpLUT::HDI][0] = deriv.data[SpLUT::HDI][0] + (scoef - acoef * HDI);




  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {



    // 1) DM

    scoef =
          kcol_buf[CollisionalRxnLUT::k132][0]    *    DI        *    de
        + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM        *    DI;
    acoef =
          kcol_buf[CollisionalRxnLUT::k133][0]    *   DII        /  2.
        + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k136][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k137][0]    *    HI;

    deriv.data[SpLUT::DM][0] = deriv.data[SpLUT::DM][0] + (scoef - acoef * DM);


    // 2) HDII

    scoef = 3. * (
          kcol_buf[CollisionalRxnLUT::k129][0]    *    DI        *   HII        /  2.
        + kcol_buf[CollisionalRxnLUT::k130][0]    *   DII        *    HI        /  2.
       );
    acoef =
          kcol_buf[CollisionalRxnLUT::k125][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k131][0]    *    de;

    deriv.data[SpLUT::HDII][0] = deriv.data[SpLUT::HDII][0] + (scoef - acoef * HDII);


    // 3) HeHII

    scoef = 5. * (
          kcol_buf[CollisionalRxnLUT::k148][0]    *   HeI        *   HII        /  4.
        + kcol_buf[CollisionalRxnLUT::k149][0]    *   HeI        *   HII        /  4.
        + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
        + kcol_buf[CollisionalRxnLUT::k151][0]    *  HeII        *    HI        /  4.
       );
    acoef =
          kcol_buf[CollisionalRxnLUT::k152][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k153][0]    *    de;

    deriv.data[SpLUT::HeHII][0] = deriv.data[SpLUT::HeHII][0] + (scoef - acoef * HeHII);




  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {

    if (itmask_metal[0] != MASK_FALSE   )  {

      // ***** CI **********
      scoef = 0. + 12. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *   CII        *    de        / 12.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *   H2I        /  2.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::AC_dust][0]          / CI        * 12.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          acoef = acoef
            + *(my_fields->RT_CI_ionization_rate);
        }
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef + 12. *
              *(my_fields->RT_CO_dissociation_rate)         * CO        /28.0;
        }
      }

      deriv.data[SpLUT::CI][0] = deriv.data[SpLUT::CI][0] + (scoef - acoef * CI);



      // ***** CII **********
      scoef = 0. + 12. * ( 0.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(my_fields->RT_CI_ionization_rate)        * CI;
        }
      }

      deriv.data[SpLUT::CII][0] = deriv.data[SpLUT::CII][0] + (scoef - acoef * CII);



      // ***** CO **********
      scoef = 0. + 28. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        *    O2        / 384.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      / CO        * 17. * 0.5
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      / CO        * 17.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(my_fields->RT_CO_dissociation_rate);
        }
      }

      deriv.data[SpLUT::CO][0] = deriv.data[SpLUT::CO][0] + (scoef - acoef * CO);



      // ***** CO2 **********
      scoef = 0. + 44. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *    OH        *    CO        / 476.
         );
      acoef = 0.;

      deriv.data[SpLUT::CO2][0] = deriv.data[SpLUT::CO2][0] + (scoef - acoef * CO2);



      // ***** OI **********
      scoef = 0. + 16. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *   OII        *    de        / 16.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *  O2II        *    de        / 16.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz31][0]    *    OI        / 8.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    CH        / 13.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OH        / 17.;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          acoef = acoef
            + *(my_fields->RT_OI_ionization_rate);
        }
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef + 16. *
            ( *(my_fields->RT_OH_dissociation_rate)         * OH        /17.0
            + *(my_fields->RT_CO_dissociation_rate)         * CO        /28.0);
        }
      }

      deriv.data[SpLUT::OI][0] = deriv.data[SpLUT::OI][0] + (scoef - acoef * OI);



      // ***** OH **********
      scoef = 0. + 17. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    OI        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *    OH        / 8.5
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *    CO        / 28.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        / 28.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        / 44.;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(my_fields->RT_OH_dissociation_rate);
          scoef = scoef + 17. *
              *(my_fields->RT_H2O_dissociation_rate)        * H2O       /18.0;
        }
      }

      deriv.data[SpLUT::OH][0] = deriv.data[SpLUT::OH][0] + (scoef - acoef * OH);



      // ***** H2O **********
      scoef = 0. + 18. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      / H2O        * 18. * 2.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / H2O        * 18. * 3.
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       / H2O        * 18. * 4.
          + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         / H2O        * 18.
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       / H2O        * 18. * 3.;
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust][0]      / H2O        * 18.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(my_fields->RT_H2O_dissociation_rate);
        }
      }

      deriv.data[SpLUT::H2O][0] = deriv.data[SpLUT::H2O][0] + (scoef - acoef * H2O);



      // ***** O2 **********
      scoef = 0. + 32. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz31][0]    *    OI        *    OI        / 256.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        *    OH        / 272.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        / 28.;

      deriv.data[SpLUT::O2][0] = deriv.data[SpLUT::O2][0] + (scoef - acoef * O2);



      // ***** SiI **********
      scoef = 0. + 28. * ( 0.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *    O2        / 32.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::SiM_dust][0]         / SiI        * 28.;
        }
      }

      deriv.data[SpLUT::SiI][0] = deriv.data[SpLUT::SiI][0] + (scoef - acoef * SiI);



      // ***** SiOI **********
      scoef = 0. + 44. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      / SiOI        * 44.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / SiOI        * 44.;
        }
      }

      deriv.data[SpLUT::SiOI][0] = deriv.data[SpLUT::SiOI][0] + (scoef - acoef * SiOI);



      // ***** SiO2I **********
      scoef = 0. + 60. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.
         );
      acoef = 0.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust][0]       / SiO2I        * 60.;
        }
      }

      deriv.data[SpLUT::SiO2I][0] = deriv.data[SpLUT::SiO2I][0] + (scoef - acoef * SiO2I);



      // ***** CH **********
      scoef = 0. + 13. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    CI        *    HI        / 12.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        / 16.;

      deriv.data[SpLUT::CH][0] = deriv.data[SpLUT::CH][0] + (scoef - acoef * CH);



      // ***** CH2 **********
      scoef = 0. + 14. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *   H2I        *    CI        / 24.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      / CH2        * 14. * 0.5;
        }
      }

      deriv.data[SpLUT::CH2][0] = deriv.data[SpLUT::CH2][0] + (scoef - acoef * CH2);



      // ***** COII **********
      scoef = 0. + 28. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        *    OH        / 204.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *    HI;

      deriv.data[SpLUT::COII][0] = deriv.data[SpLUT::COII][0] + (scoef - acoef * COII);



      // ***** OII **********
      scoef = 0. + 16. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        *    O2        / 384.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(my_fields->RT_OI_ionization_rate)        * OI;
        }
      }

      deriv.data[SpLUT::OII][0] = deriv.data[SpLUT::OII][0] + (scoef - acoef * OII);



      // ***** OHII **********
      scoef = 0. + 17. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *   H2I        /  2.;

      deriv.data[SpLUT::OHII][0] = deriv.data[SpLUT::OHII][0] + (scoef - acoef * OHII);



      // ***** H2OII **********
      scoef = 0. + 18. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    *    de
          + kcol_buf[CollisionalRxnLUT::kz47][0]    *    de;

      deriv.data[SpLUT::H2OII][0] = deriv.data[SpLUT::H2OII][0] + (scoef - acoef * H2OII);



      // ***** H3OII **********
      scoef = 0. + 19. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    *    de
          + kcol_buf[CollisionalRxnLUT::kz49][0]    *    de;

      deriv.data[SpLUT::H3OII][0] = deriv.data[SpLUT::H3OII][0] + (scoef - acoef * H3OII);



      // ***** O2II **********
      scoef = 0. + 32. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII        *    O2        / 32.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *    de;

      deriv.data[SpLUT::O2II][0] = deriv.data[SpLUT::O2II][0] + (scoef - acoef * O2II);



      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          // ***** Mg **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      / Mg        * 24.;
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / Mg        * 24. * 2.
            + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         / Mg        * 24.;
          }

          deriv.data[SpLUT::Mg][0] = deriv.data[SpLUT::Mg][0] + (scoef - acoef * Mg);


        }

        if (my_chemistry->dust_species > 1)  {
          // ***** Al **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       / Al        * 27. * 2.;

          deriv.data[SpLUT::Al][0] = deriv.data[SpLUT::Al][0] + (scoef - acoef * Al);



          // ***** S  **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         / S        * 32.;

          deriv.data[SpLUT::S][0] = deriv.data[SpLUT::S][0] + (scoef - acoef * S);



          // ***** Fe **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::FeM_dust][0]         / Fe        * 56.
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       / Fe        * 56. * 3.
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         / Fe        * 56.;

          deriv.data[SpLUT::Fe][0] = deriv.data[SpLUT::Fe][0] + (scoef - acoef * Fe);


        }
      }

    }

  }

  // --- (D4) Now do dust species ---
  if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {

    if (itmask_metal[0] != MASK_FALSE   )  {

      if (my_chemistry->dust_species > 0)  {
        // ***** MgSiO3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      * 100.;
        acoef = 0.;

        deriv.data[SpLUT::MgSiO3_dust][0] = deriv.data[SpLUT::MgSiO3_dust][0] + (scoef - acoef * MgSiO3);



        // ***** AC **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::AC_dust][0]          * 12.;
        acoef = 0.;

        deriv.data[SpLUT::AC_dust][0] = deriv.data[SpLUT::AC_dust][0] + (scoef - acoef * AC);


      }

      if (my_chemistry->dust_species > 1)  {
        // ***** SiM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::SiM_dust][0]         * 28.;
        acoef = 0.;

        deriv.data[SpLUT::SiM_dust][0] = deriv.data[SpLUT::SiM_dust][0] + (scoef - acoef * SiM);



        // ***** FeM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::FeM_dust][0]         * 56.;
        acoef = 0.;

        deriv.data[SpLUT::FeM_dust][0] = deriv.data[SpLUT::FeM_dust][0] + (scoef - acoef * FeM);



        // ***** Mg2SiO4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 140.;
        acoef = 0.;

        deriv.data[SpLUT::Mg2SiO4_dust][0] = deriv.data[SpLUT::Mg2SiO4_dust][0] + (scoef - acoef * Mg2SiO4);



        // ***** Fe3O4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       * 232.;
        acoef = 0.;

        deriv.data[SpLUT::Fe3O4_dust][0] = deriv.data[SpLUT::Fe3O4_dust][0] + (scoef - acoef * Fe3O4);



        // ***** SiO2D **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust][0]       * 60.;
        acoef = 0.;

        deriv.data[SpLUT::SiO2_dust][0] = deriv.data[SpLUT::SiO2_dust][0] + (scoef - acoef * SiO2D);



        // ***** MgO **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         * 40.;
        acoef = 0.;

        deriv.data[SpLUT::MgO_dust][0] = deriv.data[SpLUT::MgO_dust][0] + (scoef - acoef * MgO);



        // ***** FeS **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         * 88.;
        acoef = 0.;

        deriv.data[SpLUT::FeS_dust][0] = deriv.data[SpLUT::FeS_dust][0] + (scoef - acoef * FeS);



        // ***** Al2O3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       * 102.;
        acoef = 0.;

        deriv.data[SpLUT::Al2O3_dust][0] = deriv.data[SpLUT::Al2O3_dust][0] + (scoef - acoef * Al2O3);


      }

      if (my_chemistry->dust_species > 2)  {
        // ***** reforg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      * 22.68;
        acoef = 0.;

        deriv.data[SpLUT::ref_org_dust][0] = deriv.data[SpLUT::ref_org_dust][0] + (scoef - acoef * reforg);



        // ***** volorg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      * 32.;
        acoef = 0.;

        deriv.data[SpLUT::vol_org_dust][0] = deriv.data[SpLUT::vol_org_dust][0] + (scoef - acoef * volorg);



        // ***** H2Oice **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust][0]      * 18.;
        acoef = 0.;

        deriv.data[SpLUT::H2O_ice_dust][0] = deriv.data[SpLUT::H2O_ice_dust][0] + (scoef - acoef * H2Oice);


      }

    }

  }

}




} // namespace grackle::impl::chemistry

#endif /* CHEMISTRY_SOLVER_FUNCS_HPP */
