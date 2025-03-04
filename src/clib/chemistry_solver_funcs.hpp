// See LICENSE file for license and copyright information

/// @file chemistry_solver_funcs.hpp
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

#ifndef CHEMISTRY_SOLVER_FUNCS_HPP
#define CHEMISTRY_SOLVER_FUNCS_HPP

#include "grackle.h"
#include "internal_types.hpp"
#include "LUT.hpp"

namespace grackle::impl::chemistry {

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
/// @param[in]  grain_growth_rates, kcr_buf, kshield_buf specifies the
///     precomputed rxn rates (depends on local physical conditions)
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
  const grackle::impl::GrainSpeciesCollection grain_growth_rates,
  const grackle::impl::ColRecRxnRateCollection kcr_buf,
  const grackle::impl::PhotoRxnRateCollection kshield_buf
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

  double scoef, acoef;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {




    // 1) HI

    scoef  = kcr_buf.data[ColRecRxnLUT::k2][0]   *HII       *de;
    acoef  = kcr_buf.data[ColRecRxnLUT::k1][0]   *de
           + kcr_buf.data[ColRecRxnLUT::k57][0]   *HI
           + kcr_buf.data[ColRecRxnLUT::k58][0]   *HeI       /4.
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(my_fields->RT_HI_ionization_rate); }
    deriv.data[SpLUT::HI][0] = deriv.data[SpLUT::HI][0] + (scoef - acoef * HI);














    // 2) HII
    scoef  = kcr_buf.data[ColRecRxnLUT::k1][0]   *HI    *de
           + kcr_buf.data[ColRecRxnLUT::k57][0]   *HI    *HI
           + kcr_buf.data[ColRecRxnLUT::k58][0]   *HI    *HeI       /4.
           + kshield_buf.k24[0]   *HI;
    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)       *HI; }
    acoef  = kcr_buf.data[ColRecRxnLUT::k2][0]   *de;
    deriv.data[SpLUT::HII][0] = deriv.data[SpLUT::HII][0] + (scoef - acoef * HII);
















    // 3) Electron density

    scoef = 0.
               + kcr_buf.data[ColRecRxnLUT::k57][0]   *HI    *HI
               + kcr_buf.data[ColRecRxnLUT::k58][0]   *HI    *HeI       /4.
               + kshield_buf.k24[0]   *HI
               + kshield_buf.k25[0]   *HeII       /4.
               + kshield_buf.k26[0]   *HeI       /4.;

    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI
              + *(my_fields->RT_HeI_ionization_rate)         * HeI         / 4.
              + *(my_fields->RT_HeII_ionization_rate)        * HeII        / 4.; }
    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI; }



    acoef = -(kcr_buf.data[ColRecRxnLUT::k1][0]   *HI             - kcr_buf.data[ColRecRxnLUT::k2][0]   *HII
            + kcr_buf.data[ColRecRxnLUT::k3][0]   *HeI       /4. -
         kcr_buf.data[ColRecRxnLUT::k6][0]   *HeIII       /4.
            + kcr_buf.data[ColRecRxnLUT::k5][0]   *HeII       /4. -
         kcr_buf.data[ColRecRxnLUT::k4][0]   *HeII       /4.);
    deriv.data[SpLUT::e][0] = deriv.data[SpLUT::e][0] + (scoef - acoef * de);





  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---




  // 4) HeI

  scoef  = kcr_buf.data[ColRecRxnLUT::k4][0]   *HeII       *de;
  acoef  = kcr_buf.data[ColRecRxnLUT::k3][0]   *de
               + kshield_buf.k26[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(my_fields->RT_HeI_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    scoef = scoef +  4. * ( 0.
        + kcr_buf.data[ColRecRxnLUT::k152][0]    * HeHII        *    HI        /  5.
        + kcr_buf.data[ColRecRxnLUT::k153][0]    * HeHII        *    de        /  5.
       );
    acoef = acoef
        + kcr_buf.data[ColRecRxnLUT::k148][0]    *   HII
        + kcr_buf.data[ColRecRxnLUT::k149][0]    *   HII
        + kcr_buf.data[ColRecRxnLUT::k150][0]    *  H2II        /  2.;
  }
  deriv.data[SpLUT::HeI][0] = deriv.data[SpLUT::HeI][0] + (scoef - acoef * HeI);


  // 5) HeII

  scoef  = kcr_buf.data[ColRecRxnLUT::k3][0]   *HeI    *de
         + kcr_buf.data[ColRecRxnLUT::k6][0]   *HeIII       *de
         + kshield_buf.k26[0]   *HeI;

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(my_fields->RT_HeI_ionization_rate)       *HeI; }

  acoef  = kcr_buf.data[ColRecRxnLUT::k4][0]   *de        + kcr_buf.data[ColRecRxnLUT::k5][0]   *de
         + kshield_buf.k25[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(my_fields->RT_HeII_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    acoef = acoef
        + kcr_buf.data[ColRecRxnLUT::k151][0]    *    HI;
  }
  deriv.data[SpLUT::HeII][0] = deriv.data[SpLUT::HeII][0] + (scoef - acoef * HeII);


  // 6) HeIII

  scoef   = kcr_buf.data[ColRecRxnLUT::k5][0]   *HeII    *de
          + kshield_buf.k25[0]   *HeII;
  if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(my_fields->RT_HeII_ionization_rate)        * HeII; }
  acoef   = kcr_buf.data[ColRecRxnLUT::k6][0]   *de;
  deriv.data[SpLUT::HeIII][0] = deriv.data[SpLUT::HeIII][0] + (scoef - acoef * HeIII);





  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms




    // 1) HI
    scoef  =      kcr_buf.data[ColRecRxnLUT::k2][0]    * HII        * de
           + 2.*kcr_buf.data[ColRecRxnLUT::k13][0]   * HI         * H2I       /2.
           +      kcr_buf.data[ColRecRxnLUT::k11][0]   * HII        * H2I       /2.
           + 2.*kcr_buf.data[ColRecRxnLUT::k12][0]   * de         * H2I       /2.
           +      kcr_buf.data[ColRecRxnLUT::k14][0]   * HM         * de
           +      kcr_buf.data[ColRecRxnLUT::k15][0]   * HM         * HI
           + 2.*kcr_buf.data[ColRecRxnLUT::k16][0]   * HM         * HII
           + 2.*kcr_buf.data[ColRecRxnLUT::k18][0]   * H2II       * de       /2.
           +      kcr_buf.data[ColRecRxnLUT::k19][0]   * H2II       * HM       /2.
           + 2.*kshield_buf.k31[0]      * H2I       /2.;

    acoef  =      kcr_buf.data[ColRecRxnLUT::k1][0]    * de
           +      kcr_buf.data[ColRecRxnLUT::k7][0]    * de
           +      kcr_buf.data[ColRecRxnLUT::k8][0]    * HM
           +      kcr_buf.data[ColRecRxnLUT::k9][0]    * HII
           +      kcr_buf.data[ColRecRxnLUT::k10][0]   * H2II       /2.
           + 2.*kcr_buf.data[ColRecRxnLUT::k22][0]   * std::pow(HI       ,2)
           +      kcr_buf.data[ColRecRxnLUT::k57][0]   * HI
           +      kcr_buf.data[ColRecRxnLUT::k58][0]   * HeI       /4.
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
            + kcr_buf.data[ColRecRxnLUT::k50][0]    * HII        * DI         / 2.
            + kcr_buf.data[ColRecRxnLUT::k54][0]    * H2I        * DI         / 4.;
      acoef = acoef
            + kcr_buf.data[ColRecRxnLUT::k51][0]    * DII        / 2.
            + kcr_buf.data[ColRecRxnLUT::k55][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcr_buf.data[ColRecRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k135][0]    *    HM        *    DI        /  2.
          + kcr_buf.data[ColRecRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
          + kcr_buf.data[ColRecRxnLUT::k153][0]    * HeHII        *    de        /  5.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k125][0]    *  HDII        /  3.
          + kcr_buf.data[ColRecRxnLUT::k130][0]    *   DII        /  2.
          + kcr_buf.data[ColRecRxnLUT::k136][0]    *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k137][0]    *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k151][0]    *  HeII        /  4.
          + kcr_buf.data[ColRecRxnLUT::k152][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcr_buf.data[ColRecRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcr_buf.data[ColRecRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcr_buf.data[ColRecRxnLUT::kz26][0]    *    OH        *    CO        / 476.
          + kcr_buf.data[ColRecRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcr_buf.data[ColRecRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcr_buf.data[ColRecRxnLUT::kz33][0]    *    OI        *    OH        / 272.
          + kcr_buf.data[ColRecRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz36][0]    *   HII        *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz37][0]    *   CII        *    OH        / 204.
          + kcr_buf.data[ColRecRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
          + kcr_buf.data[ColRecRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
          + kcr_buf.data[ColRecRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz48][0]    * H3OII        *    de        / 19.
          + kcr_buf.data[ColRecRxnLUT::kz49][0]    * H3OII        *    de        / 9.5
          + kcr_buf.data[ColRecRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcr_buf.data[ColRecRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::kz15][0]    *    CH        / 13.
          + kcr_buf.data[ColRecRxnLUT::kz16][0]    *   CH2        / 14.
          + kcr_buf.data[ColRecRxnLUT::kz17][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz18][0]    *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz19][0]    *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz27][0]    *    CI        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz30][0]    *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz39][0]    *   OII        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz43][0]    *  COII        / 28.;
    }
    deriv.data[SpLUT::HI][0] = deriv.data[SpLUT::HI][0] + (scoef - acoef * HI);













    // 2) HII

    scoef  =    kcr_buf.data[ColRecRxnLUT::k1][0]     * HI        * de
           +    kcr_buf.data[ColRecRxnLUT::k10][0]    * H2II       *HI       /2.
           +    kcr_buf.data[ColRecRxnLUT::k57][0]    * HI        * HI
           +    kcr_buf.data[ColRecRxnLUT::k58][0]    * HI        * HeI       /4.
           + kshield_buf.k24[0]   *HI;

    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(my_fields->RT_HI_ionization_rate)        * HI; }

    acoef  =    kcr_buf.data[ColRecRxnLUT::k2][0]     * de
           +    kcr_buf.data[ColRecRxnLUT::k9][0]     * HI
           +    kcr_buf.data[ColRecRxnLUT::k11][0]    * H2I       /2.
           +    kcr_buf.data[ColRecRxnLUT::k16][0]    * HM
           +    kcr_buf.data[ColRecRxnLUT::k17][0]    * HM;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcr_buf.data[ColRecRxnLUT::k51][0]    * HI         * DII        / 2.
            + kcr_buf.data[ColRecRxnLUT::k52][0]    * H2I        * DII        / 4.;
      acoef = acoef
            + kcr_buf.data[ColRecRxnLUT::k50][0]    * DI         / 2.
            + kcr_buf.data[ColRecRxnLUT::k53][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::k125][0]    *  HDII        *    HI        /  3.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k129][0]    *    DI        /  2.
          + kcr_buf.data[ColRecRxnLUT::k134][0]    *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k148][0]    *   HeI        /  4.
          + kcr_buf.data[ColRecRxnLUT::k149][0]    *   HeI        /  4.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz43][0]    *  COII        *    HI        / 28.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::kz22][0]    *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz34][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz35][0]    *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz36][0]    *    O2        / 32.;
    }
    deriv.data[SpLUT::HII][0] = deriv.data[SpLUT::HII][0] + (scoef - acoef * HII);

    
    // 3) electrons:

    scoef =   kcr_buf.data[ColRecRxnLUT::k8][0]    * HM        * HI
           +  kcr_buf.data[ColRecRxnLUT::k15][0]   * HM        * HI
           +  kcr_buf.data[ColRecRxnLUT::k17][0]   * HM        * HII
           +  kcr_buf.data[ColRecRxnLUT::k57][0]   * HI        * HI
           +  kcr_buf.data[ColRecRxnLUT::k58][0]   * HI        * HeI       /4.
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

    acoef = - (kcr_buf.data[ColRecRxnLUT::k1][0]    *HI           - kcr_buf.data[ColRecRxnLUT::k2][0]   *HII
            +  kcr_buf.data[ColRecRxnLUT::k3][0]    *HeI       /4. -
         kcr_buf.data[ColRecRxnLUT::k6][0]   *HeIII       /4.
            +  kcr_buf.data[ColRecRxnLUT::k5][0]    *HeII       /4. -
         kcr_buf.data[ColRecRxnLUT::k4][0]   *HeII       /4.
            +  kcr_buf.data[ColRecRxnLUT::k14][0]   *HM
            -  kcr_buf.data[ColRecRxnLUT::k7][0]    *HI
            -  kcr_buf.data[ColRecRxnLUT::k18][0]   *H2II       /2.);
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcr_buf.data[ColRecRxnLUT::k56][0]    * DI         * HM        / 2.;
      acoef = acoef
            - kcr_buf.data[ColRecRxnLUT::k1][0]     * DI         / 2.
            + kcr_buf.data[ColRecRxnLUT::k2][0]     * DII        / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::k137][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k131][0]    *  HDII        /  3.
          + kcr_buf.data[ColRecRxnLUT::k132][0]    *    DI        /  2.
          + kcr_buf.data[ColRecRxnLUT::k153][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::kz44][0]    *   CII        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz45][0]    *   OII        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz46][0]    * H2OII        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz47][0]    * H2OII        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz48][0]    * H3OII        / 19.
          + kcr_buf.data[ColRecRxnLUT::kz49][0]    * H3OII        / 19.
          + kcr_buf.data[ColRecRxnLUT::kz50][0]    *  O2II        / 32.;
    }
    deriv.data[SpLUT::e][0] = deriv.data[SpLUT::e][0] + (scoef - acoef * de);


    // 7) H2

    scoef = 2.*(kcr_buf.data[ColRecRxnLUT::k8][0]     * HM          * HI
          +       kcr_buf.data[ColRecRxnLUT::k10][0]    * H2II        * HI       /2.
          +       kcr_buf.data[ColRecRxnLUT::k19][0]    * H2II        * HM       /2.
          +       kcr_buf.data[ColRecRxnLUT::k22][0]    * HI        * std::pow((HI       ),2.));
    acoef = ( kcr_buf.data[ColRecRxnLUT::k13][0]   *HI        + kcr_buf.data[ColRecRxnLUT::k11][0]   *HII
            + kcr_buf.data[ColRecRxnLUT::k12][0]   *de        )
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
              kcr_buf.data[ColRecRxnLUT::k53][0]    * HDI        * HII        / 3.
            + kcr_buf.data[ColRecRxnLUT::k55][0]    * HDI        * HI         / 3.
               );
      acoef = acoef
            + kcr_buf.data[ColRecRxnLUT::k52][0]    * DII        / 2.
            + kcr_buf.data[ColRecRxnLUT::k54][0]    * DI         / 2.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcr_buf.data[ColRecRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcr_buf.data[ColRecRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz47][0]    * H2OII        *    de        / 18.
         );
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::kz20][0]    *    CI        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz21][0]    *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz23][0]    *    CH        / 13.
          + kcr_buf.data[ColRecRxnLUT::kz24][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz40][0]    *   OII        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz41][0]    *  OHII        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz42][0]    * H2OII        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz51][0]    *    CI        / 12.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          scoef = scoef + 2. *
                grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      * 2.;

        }
        if (my_chemistry->dust_species > 1)  {
          scoef = scoef + 2. * (
                grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 3.
              + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       * 4.
              + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]
              + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       * 3.
            );
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      / H2I        * 2. * 2.;
        }
      }
    }
    deriv.data[SpLUT::H2I][0] = deriv.data[SpLUT::H2I][0] + (scoef - acoef * H2I);


    // 8) H-

    scoef = kcr_buf.data[ColRecRxnLUT::k7][0]    * HI        * de;
    acoef = (kcr_buf.data[ColRecRxnLUT::k8][0]     + kcr_buf.data[ColRecRxnLUT::k15][0]   )  * HI        +
            (kcr_buf.data[ColRecRxnLUT::k16][0]    + kcr_buf.data[ColRecRxnLUT::k17][0]   )  * HII        +
            kcr_buf.data[ColRecRxnLUT::k14][0]    * de        + kcr_buf.data[ColRecRxnLUT::k19][0]    * H2II       /2.0f +
            my_uvb_rates.k27;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      acoef = acoef
            + kcr_buf.data[ColRecRxnLUT::k56][0]    * DI         / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[ColRecRxnLUT::k136][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k135][0]    *    DI        /  2.;
    }
    deriv.data[SpLUT::HM][0] = deriv.data[SpLUT::HM][0] + (scoef - acoef * HM);



    // 9) H2+

    scoef =    2.*( kcr_buf.data[ColRecRxnLUT::k9][0]    *HI    *HII
                  +   kcr_buf.data[ColRecRxnLUT::k11][0]   *H2I    /2.*HII
                  +   kcr_buf.data[ColRecRxnLUT::k17][0]   *HM    *HII
                  + kshield_buf.k29[0]   *H2I    /2.
                  );
    acoef =         kcr_buf.data[ColRecRxnLUT::k10][0]   *HI     + kcr_buf.data[ColRecRxnLUT::k18][0]   *de
                  + kcr_buf.data[ColRecRxnLUT::k19][0]   *HM
                  + (kshield_buf.k28[0]   +kshield_buf.k30[0]   );
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::k152][0]    * HeHII        *    HI        /  5.
         );
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k150][0]    *   HeI        /  4.;
    }
    deriv.data[SpLUT::H2II][0] = deriv.data[SpLUT::H2II][0] + (scoef - acoef * H2II);









  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {


    
    // 1) DI
    scoef =   (       kcr_buf.data[ColRecRxnLUT::k2][0]    * DII        * de
               +      kcr_buf.data[ColRecRxnLUT::k51][0]   * DII        * HI
               + 2.*kcr_buf.data[ColRecRxnLUT::k55][0]   * HDI        *
            HI       /3.
               );
    acoef  =    kcr_buf.data[ColRecRxnLUT::k1][0]    * de
           +    kcr_buf.data[ColRecRxnLUT::k50][0]    * HII
           +    kcr_buf.data[ColRecRxnLUT::k54][0]    * H2I       /2.
           +    kcr_buf.data[ColRecRxnLUT::k56][0]    * HM
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(my_fields->RT_HI_ionization_rate); }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcr_buf.data[ColRecRxnLUT::k133][0]    *   DII        *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcr_buf.data[ColRecRxnLUT::k136][0]    *    DM        *    HI        /  2.
          );
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k129][0]    *   HII
          + kcr_buf.data[ColRecRxnLUT::k132][0]    *    de
          + kcr_buf.data[ColRecRxnLUT::k135][0]    *    HM;
    }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        scoef = scoef
          + 2. * *(my_fields->RT_HDI_dissociation_rate)        * HDI       /3.0;
      }
    }
    deriv.data[SpLUT::DI][0] = deriv.data[SpLUT::DI][0] + (scoef - acoef * DI);
                                                    

    // 2) DII
    scoef =   (   kcr_buf.data[ColRecRxnLUT::k1][0]     * DI        * de
          +       kcr_buf.data[ColRecRxnLUT::k50][0]    * HII       * DI
          +  2.*kcr_buf.data[ColRecRxnLUT::k53][0]    * HII       * HDI       /3.
          )
          + kshield_buf.k24[0]   *DI;
    acoef = 0.;
    // ! initialize GC202002
    if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + *(my_fields->RT_HI_ionization_rate)       *DI; }
    acoef =    kcr_buf.data[ColRecRxnLUT::k2][0]     * de
          +    kcr_buf.data[ColRecRxnLUT::k51][0]    * HI
          +    kcr_buf.data[ColRecRxnLUT::k52][0]    * H2I       /2.;
    if (my_chemistry->primordial_chemistry > 3)  {
      acoef = acoef
          + kcr_buf.data[ColRecRxnLUT::k130][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::k133][0]    *    DM        /  2.;
    }
    deriv.data[SpLUT::DII][0] = deriv.data[SpLUT::DII][0] + (scoef - acoef * DII);


    // 3) HDI
    scoef = 3.*(kcr_buf.data[ColRecRxnLUT::k52][0]    * DII       *
         H2I       /2./2.
         + kcr_buf.data[ColRecRxnLUT::k54][0]    * DI        * H2I       /2./2.
    // !   &           + 2._DKIND*kcr_buf.data[ColRecRxnLUT::k56][0]    * DI        * HM       /2._DKIND
    //- ! corrected by GC202005
         +          kcr_buf.data[ColRecRxnLUT::k56][0]    * DI        * HM       /2.
               );
    acoef  =    kcr_buf.data[ColRecRxnLUT::k53][0]    * HII
           +    kcr_buf.data[ColRecRxnLUT::k55][0]    * HI;
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        acoef = acoef
          + *(my_fields->RT_HDI_dissociation_rate);
      }
    }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  3. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::k125][0]    *  HDII        *    HI        /  3.
          + kcr_buf.data[ColRecRxnLUT::k137][0]    *    DM        *    HI        /  2.
          );
    }
    deriv.data[SpLUT::HDI][0] = deriv.data[SpLUT::HDI][0] + (scoef - acoef * HDI);




  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {



    // 1) DM

    scoef =
          kcr_buf.data[ColRecRxnLUT::k132][0]    *    DI        *    de
        + kcr_buf.data[ColRecRxnLUT::k135][0]    *    HM        *    DI;
    acoef =
          kcr_buf.data[ColRecRxnLUT::k133][0]    *   DII        /  2.
        + kcr_buf.data[ColRecRxnLUT::k134][0]    *   HII
        + kcr_buf.data[ColRecRxnLUT::k136][0]    *    HI
        + kcr_buf.data[ColRecRxnLUT::k137][0]    *    HI;

    deriv.data[SpLUT::DM][0] = deriv.data[SpLUT::DM][0] + (scoef - acoef * DM);


    // 2) HDII

    scoef = 3. * (
          kcr_buf.data[ColRecRxnLUT::k129][0]    *    DI        *   HII        /  2.
        + kcr_buf.data[ColRecRxnLUT::k130][0]    *   DII        *    HI        /  2.
       );
    acoef =
          kcr_buf.data[ColRecRxnLUT::k125][0]    *    HI
        + kcr_buf.data[ColRecRxnLUT::k131][0]    *    de;

    deriv.data[SpLUT::HDII][0] = deriv.data[SpLUT::HDII][0] + (scoef - acoef * HDII);


    // 3) HeHII

    scoef = 5. * (
          kcr_buf.data[ColRecRxnLUT::k148][0]    *   HeI        *   HII        /  4.
        + kcr_buf.data[ColRecRxnLUT::k149][0]    *   HeI        *   HII        /  4.
        + kcr_buf.data[ColRecRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
        + kcr_buf.data[ColRecRxnLUT::k151][0]    *  HeII        *    HI        /  4.
       );
    acoef =
          kcr_buf.data[ColRecRxnLUT::k152][0]    *    HI
        + kcr_buf.data[ColRecRxnLUT::k153][0]    *    de;

    deriv.data[SpLUT::HeHII][0] = deriv.data[SpLUT::HeHII][0] + (scoef - acoef * HeHII);




  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {

    if (itmask_metal[0] != MASK_FALSE   )  {

      // ***** CI **********
      scoef = 0. + 12. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcr_buf.data[ColRecRxnLUT::kz44][0]    *   CII        *    de        / 12.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz20][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz27][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz28][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz29][0]    *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz51][0]    *   H2I        /  2.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][0]          / CI        * 12.;
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
          + kcr_buf.data[ColRecRxnLUT::kz37][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz38][0]    *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz44][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(my_fields->RT_CI_ionization_rate)        * CI;
        }
      }

      deriv.data[SpLUT::CII][0] = deriv.data[SpLUT::CII][0] + (scoef - acoef * CII);



      // ***** CO **********
      scoef = 0. + 28. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcr_buf.data[ColRecRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcr_buf.data[ColRecRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcr_buf.data[ColRecRxnLUT::kz38][0]    *   CII        *    O2        / 384.
          + kcr_buf.data[ColRecRxnLUT::kz43][0]    *  COII        *    HI        / 28.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz26][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      / CO        * 17. * 0.5
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      / CO        * 17.;
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
          + kcr_buf.data[ColRecRxnLUT::kz26][0]    *    OH        *    CO        / 476.
         );
      acoef = 0.;

      deriv.data[SpLUT::CO2][0] = deriv.data[SpLUT::CO2][0] + (scoef - acoef * CO2);



      // ***** OI **********
      scoef = 0. + 16. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcr_buf.data[ColRecRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcr_buf.data[ColRecRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz45][0]    *   OII        *    de        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz47][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz50][0]    *  O2II        *    de        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz21][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz22][0]    *   HII
          + kcr_buf.data[ColRecRxnLUT::kz30][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz31][0]    *    OI        / 8.
          + kcr_buf.data[ColRecRxnLUT::kz32][0]    *    CH        / 13.
          + kcr_buf.data[ColRecRxnLUT::kz33][0]    *    OH        / 17.;
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
          + kcr_buf.data[ColRecRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcr_buf.data[ColRecRxnLUT::kz30][0]    *    OI        *    HI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz49][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz17][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz24][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz25][0]    *    OH        / 8.5
          + kcr_buf.data[ColRecRxnLUT::kz26][0]    *    CO        / 28.
          + kcr_buf.data[ColRecRxnLUT::kz28][0]    *    CI        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz33][0]    *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz34][0]    *   HII
          + kcr_buf.data[ColRecRxnLUT::kz37][0]    *   CII        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz52][0]    *   SiI        / 28.
          + kcr_buf.data[ColRecRxnLUT::kz54][0]    *  SiOI        / 44.;
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
          + kcr_buf.data[ColRecRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcr_buf.data[ColRecRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcr_buf.data[ColRecRxnLUT::kz48][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz18][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz35][0]    *   HII;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / H2O        * 18. * 2.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / H2O        * 18. * 3.
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       / H2O        * 18. * 4.
          + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         / H2O        * 18.
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       / H2O        * 18. * 3.;
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][0]      / H2O        * 18.;
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
          + kcr_buf.data[ColRecRxnLUT::kz31][0]    *    OI        *    OI        / 256.
          + kcr_buf.data[ColRecRxnLUT::kz33][0]    *    OI        *    OH        / 272.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz19][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz29][0]    *    CI        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz36][0]    *   HII
          + kcr_buf.data[ColRecRxnLUT::kz38][0]    *   CII        / 12.
          + kcr_buf.data[ColRecRxnLUT::kz53][0]    *   SiI        / 28.;

      deriv.data[SpLUT::O2][0] = deriv.data[SpLUT::O2][0] + (scoef - acoef * O2);



      // ***** SiI **********
      scoef = 0. + 28. * ( 0.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz52][0]    *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz53][0]    *    O2        / 32.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][0]         / SiI        * 28.;
        }
      }

      deriv.data[SpLUT::SiI][0] = deriv.data[SpLUT::SiI][0] + (scoef - acoef * SiI);



      // ***** SiOI **********
      scoef = 0. + 44. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcr_buf.data[ColRecRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz54][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / SiOI        * 44.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / SiOI        * 44.;
        }
      }

      deriv.data[SpLUT::SiOI][0] = deriv.data[SpLUT::SiOI][0] + (scoef - acoef * SiOI);



      // ***** SiO2I **********
      scoef = 0. + 60. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.
         );
      acoef = 0.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][0]       / SiO2I        * 60.;
        }
      }

      deriv.data[SpLUT::SiO2I][0] = deriv.data[SpLUT::SiO2I][0] + (scoef - acoef * SiO2I);



      // ***** CH **********
      scoef = 0. + 13. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcr_buf.data[ColRecRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcr_buf.data[ColRecRxnLUT::kz27][0]    *    CI        *    HI        / 12.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz15][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz23][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz32][0]    *    OI        / 16.;

      deriv.data[SpLUT::CH][0] = deriv.data[SpLUT::CH][0] + (scoef - acoef * CH);



      // ***** CH2 **********
      scoef = 0. + 14. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcr_buf.data[ColRecRxnLUT::kz51][0]    *   H2I        *    CI        / 24.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz16][0]    *    HI;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      / CH2        * 14. * 0.5;
        }
      }

      deriv.data[SpLUT::CH2][0] = deriv.data[SpLUT::CH2][0] + (scoef - acoef * CH2);



      // ***** COII **********
      scoef = 0. + 28. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz37][0]    *   CII        *    OH        / 204.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz43][0]    *    HI;

      deriv.data[SpLUT::COII][0] = deriv.data[SpLUT::COII][0] + (scoef - acoef * COII);



      // ***** OII **********
      scoef = 0. + 16. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcr_buf.data[ColRecRxnLUT::kz38][0]    *   CII        *    O2        / 384.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz39][0]    *    HI
          + kcr_buf.data[ColRecRxnLUT::kz40][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz45][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(my_fields->RT_OI_ionization_rate)        * OI;
        }
      }

      deriv.data[SpLUT::OII][0] = deriv.data[SpLUT::OII][0] + (scoef - acoef * OII);



      // ***** OHII **********
      scoef = 0. + 17. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcr_buf.data[ColRecRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz41][0]    *   H2I        /  2.;

      deriv.data[SpLUT::OHII][0] = deriv.data[SpLUT::OHII][0] + (scoef - acoef * OHII);



      // ***** H2OII **********
      scoef = 0. + 18. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcr_buf.data[ColRecRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz42][0]    *   H2I        /  2.
          + kcr_buf.data[ColRecRxnLUT::kz46][0]    *    de
          + kcr_buf.data[ColRecRxnLUT::kz47][0]    *    de;

      deriv.data[SpLUT::H2OII][0] = deriv.data[SpLUT::H2OII][0] + (scoef - acoef * H2OII);



      // ***** H3OII **********
      scoef = 0. + 19. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz48][0]    *    de
          + kcr_buf.data[ColRecRxnLUT::kz49][0]    *    de;

      deriv.data[SpLUT::H3OII][0] = deriv.data[SpLUT::H3OII][0] + (scoef - acoef * H3OII);



      // ***** O2II **********
      scoef = 0. + 32. * ( 0.
          + kcr_buf.data[ColRecRxnLUT::kz36][0]    *   HII        *    O2        / 32.
         );
      acoef = 0.
          + kcr_buf.data[ColRecRxnLUT::kz50][0]    *    de;

      deriv.data[SpLUT::O2II][0] = deriv.data[SpLUT::O2II][0] + (scoef - acoef * O2II);



      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          // ***** Mg **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / Mg        * 24.;
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / Mg        * 24. * 2.
            + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         / Mg        * 24.;
          }

          deriv.data[SpLUT::Mg][0] = deriv.data[SpLUT::Mg][0] + (scoef - acoef * Mg);


        }

        if (my_chemistry->dust_species > 1)  {
          // ***** Al **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       / Al        * 27. * 2.;

          deriv.data[SpLUT::Al][0] = deriv.data[SpLUT::Al][0] + (scoef - acoef * Al);



          // ***** S  **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         / S        * 32.;

          deriv.data[SpLUT::S][0] = deriv.data[SpLUT::S][0] + (scoef - acoef * S);



          // ***** Fe **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][0]         / Fe        * 56.
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       / Fe        * 56. * 3.
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         / Fe        * 56.;

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
        + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      * 100.;
        acoef = 0.;

        deriv.data[SpLUT::MgSiO3_dust][0] = deriv.data[SpLUT::MgSiO3_dust][0] + (scoef - acoef * MgSiO3);



        // ***** AC **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][0]          * 12.;
        acoef = 0.;

        deriv.data[SpLUT::AC_dust][0] = deriv.data[SpLUT::AC_dust][0] + (scoef - acoef * AC);


      }

      if (my_chemistry->dust_species > 1)  {
        // ***** SiM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][0]         * 28.;
        acoef = 0.;

        deriv.data[SpLUT::SiM_dust][0] = deriv.data[SpLUT::SiM_dust][0] + (scoef - acoef * SiM);



        // ***** FeM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][0]         * 56.;
        acoef = 0.;

        deriv.data[SpLUT::FeM_dust][0] = deriv.data[SpLUT::FeM_dust][0] + (scoef - acoef * FeM);



        // ***** Mg2SiO4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 140.;
        acoef = 0.;

        deriv.data[SpLUT::Mg2SiO4_dust][0] = deriv.data[SpLUT::Mg2SiO4_dust][0] + (scoef - acoef * Mg2SiO4);



        // ***** Fe3O4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       * 232.;
        acoef = 0.;

        deriv.data[SpLUT::Fe3O4_dust][0] = deriv.data[SpLUT::Fe3O4_dust][0] + (scoef - acoef * Fe3O4);



        // ***** SiO2D **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][0]       * 60.;
        acoef = 0.;

        deriv.data[SpLUT::SiO2_dust][0] = deriv.data[SpLUT::SiO2_dust][0] + (scoef - acoef * SiO2D);



        // ***** MgO **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         * 40.;
        acoef = 0.;

        deriv.data[SpLUT::MgO_dust][0] = deriv.data[SpLUT::MgO_dust][0] + (scoef - acoef * MgO);



        // ***** FeS **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         * 88.;
        acoef = 0.;

        deriv.data[SpLUT::FeS_dust][0] = deriv.data[SpLUT::FeS_dust][0] + (scoef - acoef * FeS);



        // ***** Al2O3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       * 102.;
        acoef = 0.;

        deriv.data[SpLUT::Al2O3_dust][0] = deriv.data[SpLUT::Al2O3_dust][0] + (scoef - acoef * Al2O3);


      }

      if (my_chemistry->dust_species > 2)  {
        // ***** reforg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      * 22.68;
        acoef = 0.;

        deriv.data[SpLUT::ref_org_dust][0] = deriv.data[SpLUT::ref_org_dust][0] + (scoef - acoef * reforg);



        // ***** volorg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      * 32.;
        acoef = 0.;

        deriv.data[SpLUT::vol_org_dust][0] = deriv.data[SpLUT::vol_org_dust][0] + (scoef - acoef * volorg);



        // ***** H2Oice **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][0]      * 18.;
        acoef = 0.;

        deriv.data[SpLUT::H2O_ice_dust][0] = deriv.data[SpLUT::H2O_ice_dust][0] + (scoef - acoef * H2Oice);


      }

    }

  }

}


} // namespace grackle::impl::chemistry

#endif /* CHEMISTRY_SOLVER_FUNCS_HPP */
