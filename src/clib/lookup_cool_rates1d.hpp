//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares signature of lookup_cool_rates1d_g
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// lookup_cool_rates1d_g function from FORTRAN to C++

#ifndef LOOKUP_COOL_RATES1D_HPP
#define LOOKUP_COOL_RATES1D_HPP

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

namespace grackle::impl {

/// This routine uses the temperature to look up the chemical rates which are
/// tabulated in a log table as a function of temperature.
///
/// > [!important]
/// > TODO: The role of the `dt` argument **MUST** be clarified! It is passed
/// > different values in different areas of the codebase!!!!
/// > - `solve_rate_cool_g` passes in the value of the total timestep that the
/// >   chemistry is evolved. This is the traditional meaning of `dt`
/// > - the time derivative calculation within `step_rate_newton_raphson`
/// >   passes the timestep of the current subcycle (effectively the whole
/// >   function is only being called for a single element idx_range)
/// >
/// > Internally, this arg only appears to be used to determine dust grain
/// > destruction rate.
/// > - the dust destruction rate is 0 for all temperatures below some
/// >   threshold (the threshold depends on the grain species)
/// > - above the threshold, the destruction rate is essentially the current
/// >   grain density divided by the value of the `dt` argument
/// >
/// > If you think about it:
/// > - I'd argue that setting `dt` to the whole timestep that we are evolving
/// >   the zone over is blatantly wrong. It violates the principle that you
/// >   should get consistent results whether you invoke grackle 100 separate
/// >   times or just 1 time. (The amount of dust heating would change)
/// > - setting `dt` to the current subcycle timestep makes a lot more sense
/// >   (and is the only logical choice)
/// >   - It is roughly equivalent to saying that dust is immediately destroyed
/// >     once the gas reaches a threshold temperature.
/// >   - the model is overly simplistic since dust grains can survive for
/// >     quite in ionized gas (see for example
/// >     https://ui.adsabs.harvard.edu/abs/2024ApJ...974...81R/abstract)
/// >
/// > If we stick with this instantaneous destruction model, then all
/// > dust-grain related heating and cooling should probably assume that the
/// > dust-grain density is already 0.
inline void lookup_cool_rates1d_g(
  IndexRange idx_range, gr_mask_type anydust, double* tgas1d, double* mmw,
  double* tdust, double* dust2gas, double* k13dd_data_, double* h2dust,
  double dom, double dx_cgs, double c_ljeans, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, int* imetal, gr_float* rhoH, double* dt,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  InternalGrUnits internalu,
  grackle::impl::GrainSpeciesCollection grain_growth_rates,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::CollisionalRxnRateCollection kcr_buf,
  grackle::impl::CollisionalRxnRateCollection kcol_rate_tables,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf) {
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  // -------------------------------------------------------------------

  // Arguments

  // Chemistry rates as a function of temperature

  grackle::impl::View<double**> k13dda(
      my_rates->k13dd, my_chemistry->NumberOfTemperatureBins, 14);
  grackle::impl::View<double**> h2dusta(
      my_rates->h2dust, my_chemistry->NumberOfTemperatureBins,
      my_chemistry->NumberOfDustTemperatureBins);

  // Density fields

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(
      my_fields->HeI_density, my_fields->grid_dimension[0],
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
  grackle::impl::View<gr_float***> H2I(
      my_fields->H2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DI(
      my_fields->DI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(
      my_fields->DII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(
      my_fields->HDI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(
      my_fields->CI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(
      my_fields->H2O_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiI(
      my_fields->SiI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiOI(
      my_fields->SiOI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2I(
      my_fields->SiO2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg(
      my_fields->Mg_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al(
      my_fields->Al_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> S(
      my_fields->S_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe(
      my_fields->Fe_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(
      my_fields->SiM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(
      my_fields->FeM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(
      my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(
      my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(
      my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(
      my_fields->AC_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(
      my_fields->SiO2_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(
      my_fields->MgO_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(
      my_fields->FeS_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(
      my_fields->Al2O3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(
      my_fields->ref_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(
      my_fields->vol_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(
      my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Radiation fields

  grackle::impl::View<gr_float***> kdissH2I(
      my_fields->RT_H2_dissociation_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // H2 self-shielding length-scale field

  grackle::impl::View<gr_float***> xH2shield(
      my_fields->H2_self_shielding_length, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Custom H2 shielding factor

  grackle::impl::View<gr_float***> f_shield_custom(
      my_fields->H2_custom_shielding_factor, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Returned rate values

  grackle::impl::View<double**> k13dd(k13dd_data_, my_fields->grid_dimension[0],
                                      14);

  // 1D temporaries (locally allocated)

  std::vector<long long> d_indixe(my_fields->grid_dimension[0]);
  std::vector<double> d_t1(my_fields->grid_dimension[0]);
  std::vector<double> d_t2(my_fields->grid_dimension[0]);
  std::vector<double> d_logtem(my_fields->grid_dimension[0]);
  std::vector<double> d_tdef(my_fields->grid_dimension[0]);
  std::vector<double> dusti1(my_fields->grid_dimension[0]);
  std::vector<double> dusti2(my_fields->grid_dimension[0]);
  std::vector<double> divrhoa(6);
  std::vector<double> f_shield_H(my_fields->grid_dimension[0]);
  std::vector<double> f_shield_He(my_fields->grid_dimension[0]);

  // Parameters

  const double everg = ev2erg_grflt;
  const double e24 = 13.6;
  const double e26 = 24.6;

  // locals

  int i, n1;
  double factor, x, logtem0, logtem9, dlogtem, nh, d_logtem0, d_logtem9,
      d_dlogtem, divrho, N_H2, f_shield, b_doppler, l_H2shield;
  double k13_CID, k13_DT;
  double k13ind;
  std::vector<double> logT(my_fields->grid_dimension[0]);
  std::vector<double> logrho(my_fields->grid_dimension[0]);
  // opacity table
  double log_kh2, log_kgg;

  // stuff related to grain-growth
  //
  // internal_dust_prop_buf holds buffers of intermediate quantities used
  // within dust-routines
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf =
      grackle::impl::new_InternalDustPropBuf(my_fields->grid_dimension[0],
                                             my_rates->gr_N[1]);
  // these are some legacy variables that referene allocations now tracked
  // within internal_dust_prop_buf.
  //
  // TODO: directly access members of internal_dust_prop_buf (and get rid of
  // these demporary variables)
  double* sgSiM = internal_dust_prop_buf.grain_sigma_per_gas_mass
                      .data[OnlyGrainSpLUT::SiM_dust];
  double* sgFeM = internal_dust_prop_buf.grain_sigma_per_gas_mass
                      .data[OnlyGrainSpLUT::FeM_dust];
  double* sgMg2SiO4 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                          .data[OnlyGrainSpLUT::Mg2SiO4_dust];
  double* sgMgSiO3 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::MgSiO3_dust];
  double* sgFe3O4 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::Fe3O4_dust];
  double* sgAC = internal_dust_prop_buf.grain_sigma_per_gas_mass
                     .data[OnlyGrainSpLUT::AC_dust];
  double* sgSiO2D = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::SiO2_dust];
  double* sgMgO = internal_dust_prop_buf.grain_sigma_per_gas_mass
                      .data[OnlyGrainSpLUT::MgO_dust];
  double* sgFeS = internal_dust_prop_buf.grain_sigma_per_gas_mass
                      .data[OnlyGrainSpLUT::FeS_dust];
  double* sgAl2O3 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::Al2O3_dust];
  double* sgreforg = internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::ref_org_dust];
  double* sgvolorg = internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::vol_org_dust];
  double* sgH2Oice = internal_dust_prop_buf.grain_sigma_per_gas_mass
                         .data[OnlyGrainSpLUT::H2O_ice_dust];
  double* sgtot = internal_dust_prop_buf.sigma_per_gas_mass_tot;

  double h2SiM, h2FeM, h2Mg2SiO4, h2MgSiO3, h2Fe3O4, h2AC, h2SiO2D, h2MgO,
      h2FeS, h2Al2O3, h2reforg, h2volorg, h2H2Oice;
  // tabulate h2 formation rate
  long long d_Size;
  std::vector<long long> d_N(2);
  double d_dTd, d_dTg;
  std::vector<double> d_Td(my_chemistry->NumberOfDustTemperatureBins);
  std::vector<double> d_Tg(my_chemistry->NumberOfTemperatureBins);
  int idratec, iratec;
  double kd;
  const double mh_local_var = mh_grflt;

  // locals for H2 self-shielding as WG+19

  double tgas_touse, ngas_touse, aWG2019;

  double nSSh, nratio;

  // this is a temporary variable to help with transcription
  // -> we do this to circumvent the transcription issue of casting a
  //    variable and passing the result by reference.
  // -> we make this a 1 element array to let us avoid some warnings
  //    from some simple analysis-routines that check that the array
  //    rank of a variable (or the lack thereof) is consistent with the
  //    declaration within the called subroutine
  std::vector<long long> nratec_single_elem_arr(1);

  // Cast nratec
  nratec_single_elem_arr[1 - 1] =
      (long long)(my_chemistry->NumberOfTemperatureBins);

  // Set log values of start and end of lookup tables

  logtem0 = std::log(my_chemistry->TemperatureStart);
  logtem9 = std::log(my_chemistry->TemperatureEnd);
  dlogtem = (std::log(my_chemistry->TemperatureEnd) -
             std::log(my_chemistry->TemperatureStart)) /
            (double)(my_chemistry->NumberOfTemperatureBins - 1);

  for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
    if (itmask[i - 1] != MASK_FALSE) {
      // Compute temp-centered temperature (and log)

      // logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
      logTlininterp_buf.logtem[i - 1] = std::log(tgas1d[i - 1]);
      logTlininterp_buf.logtem[i - 1] =
          std::fmax(logTlininterp_buf.logtem[i - 1], logtem0);
      logTlininterp_buf.logtem[i - 1] =
          std::fmin(logTlininterp_buf.logtem[i - 1], logtem9);

      // Find index into tble and precompute interpolation values

      logTlininterp_buf.indixe[i - 1] = std::fmin(
          my_chemistry->NumberOfTemperatureBins - 1,
          std::fmax(1, (long long)((logTlininterp_buf.logtem[i - 1] - logtem0) /
                                   dlogtem) +
                           1));
      logTlininterp_buf.t1[i - 1] =
          (logtem0 + (logTlininterp_buf.indixe[i - 1] - 1) * dlogtem);
      logTlininterp_buf.t2[i - 1] =
          (logtem0 + (logTlininterp_buf.indixe[i - 1]) * dlogtem);
      logTlininterp_buf.tdef[i - 1] =
          (logTlininterp_buf.logtem[i - 1] - logTlininterp_buf.t1[i - 1]) /
          (logTlininterp_buf.t2[i - 1] - logTlininterp_buf.t1[i - 1]);

      // Do linear table lookup (in log temperature)

      kcr_buf.data[CollisionalRxnLUT::k1][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k1]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k1]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k1]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k2][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k2]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k2]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k2]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k3][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k3]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k3]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k3]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k4][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k4]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k4]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k4]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k5][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k5]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k5]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k5]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k6][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k6]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k6]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k6]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k57][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k57]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k57]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k57]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
      kcr_buf.data[CollisionalRxnLUT::k58][i - 1] =
          kcol_rate_tables.data[CollisionalRxnLUT::k58]
                               [logTlininterp_buf.indixe[i - 1] - 1] +
          (kcol_rate_tables.data[CollisionalRxnLUT::k58]
                                [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
           kcol_rate_tables.data[CollisionalRxnLUT::k58]
                                [logTlininterp_buf.indixe[i - 1] - 1]) *
              logTlininterp_buf.tdef[i - 1];
    }
  }

  // Look-up for 9-species model

  if (my_chemistry->primordial_chemistry > 1) {
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        kcr_buf.data[CollisionalRxnLUT::k7][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k7]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k7]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k7]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k8][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k8]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k8]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k8]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k9][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k9]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k9]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k9]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k10][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k10]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k10]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k10]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k11][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k11]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k11]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k11]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k12][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k12]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k12]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k12]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k13][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k13]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k13]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k13]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k14][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k14]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k14]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k14]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k15][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k15]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k15]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k15]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k16][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k16]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k16]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k16]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k17][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k17]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k17]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k17]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k18][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k18]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k18]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k18]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k19][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k19]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k19]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k19]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k22][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k22]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k22]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k22]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];

        // H2 formation heating terms.

        chemheatrates_buf.n_cr_n[i - 1] =
            my_rates->n_cr_n[logTlininterp_buf.indixe[i - 1] - 1] +
            (my_rates->n_cr_n[logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             my_rates->n_cr_n[logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        chemheatrates_buf.n_cr_d1[i - 1] =
            my_rates->n_cr_d1[logTlininterp_buf.indixe[i - 1] - 1] +
            (my_rates->n_cr_d1[logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             my_rates->n_cr_d1[logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        chemheatrates_buf.n_cr_d2[i - 1] =
            my_rates->n_cr_d2[logTlininterp_buf.indixe[i - 1] - 1] +
            (my_rates->n_cr_d2[logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             my_rates->n_cr_d2[logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
      }
    }

    for (n1 = 1; n1 <= (14); n1++) {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          k13dd(i - 1, n1 - 1) =
              k13dda(logTlininterp_buf.indixe[i - 1] - 1, n1 - 1) +
              (k13dda(logTlininterp_buf.indixe[i - 1] + 1 - 1, n1 - 1) -
               k13dda(logTlininterp_buf.indixe[i - 1] - 1, n1 - 1)) *
                  logTlininterp_buf.tdef[i - 1];
        }
      }
    }
  }

  // Look-up for 12-species model

  if (my_chemistry->primordial_chemistry > 2) {
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        kcr_buf.data[CollisionalRxnLUT::k50][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k50]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k50]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k50]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k51][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k51]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k51]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k51]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k52][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k52]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k52]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k52]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k53][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k53]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k53]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k53]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k54][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k54]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k54]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k54]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k55][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k55]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k55]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k55]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k56][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k56]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k56]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k56]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
      }
    }
  }

  // Look-up for 15-species model

  if (my_chemistry->primordial_chemistry > 3) {
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        kcr_buf.data[CollisionalRxnLUT::k125][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k125]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k125]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k125]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k129][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k129]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k129]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k129]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k130][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k130]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k130]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k130]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k131][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k131]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k131]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k131]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k132][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k132]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k132]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k132]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k133][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k133]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k133]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k133]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k134][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k134]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k134]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k134]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k135][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k135]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k135]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k135]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k136][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k136]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k136]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k136]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k137][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k137]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k137]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k137]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k148][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k148]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k148]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k148]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k149][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k149]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k149]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k149]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k150][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k150]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k150]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k150]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k151][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k151]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k151]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k151]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k152][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k152]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k152]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k152]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::k153][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::k153]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::k153]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::k153]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
      }
    }
  }

  // Look-up for metal species model

  if (my_chemistry->metal_chemistry == 1) {
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        kcr_buf.data[CollisionalRxnLUT::kz15][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz15]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz15]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz15]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz16][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz16]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz16]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz16]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz17][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz17]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz17]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz17]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz18][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz18]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz18]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz18]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz19][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz19]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz19]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz19]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz20][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz20]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz20]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz20]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz21][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz21]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz21]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz21]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz22][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz22]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz22]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz22]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz23][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz23]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz23]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz23]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz24][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz24]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz24]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz24]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz25][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz25]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz25]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz25]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz26][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz26]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz26]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz26]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz27][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz27]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz27]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz27]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz28][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz28]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz28]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz28]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz29][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz29]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz29]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz29]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz30][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz30]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz30]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz30]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz31][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz31]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz31]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz31]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz32][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz32]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz32]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz32]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz33][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz33]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz33]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz33]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz34][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz34]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz34]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz34]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz35][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz35]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz35]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz35]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz36][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz36]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz36]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz36]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz37][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz37]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz37]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz37]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz38][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz38]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz38]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz38]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz39][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz39]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz39]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz39]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz40][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz40]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz40]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz40]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz41][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz41]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz41]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz41]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz42][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz42]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz42]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz42]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz43][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz43]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz43]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz43]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz44][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz44]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz44]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz44]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz45][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz45]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz45]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz45]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz46][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz46]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz46]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz46]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz47][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz47]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz47]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz47]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz48][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz48]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz48]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz48]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz49][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz49]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz49]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz49]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz50][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz50]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz50]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz50]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz51][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz51]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz51]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz51]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz52][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz52]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz52]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz52]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz53][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz53]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz53]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz53]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
        kcr_buf.data[CollisionalRxnLUT::kz54][i - 1] =
            kcol_rate_tables.data[CollisionalRxnLUT::kz54]
                                 [logTlininterp_buf.indixe[i - 1] - 1] +
            (kcol_rate_tables.data[CollisionalRxnLUT::kz54]
                                  [logTlininterp_buf.indixe[i - 1] + 1 - 1] -
             kcol_rate_tables.data[CollisionalRxnLUT::kz54]
                                  [logTlininterp_buf.indixe[i - 1] - 1]) *
                logTlininterp_buf.tdef[i - 1];
      }
    }
  }

  // Compute grain size increment

  if ((anydust != MASK_FALSE) && (my_chemistry->dust_species > 0)) {
    f_wrap::calc_grain_size_increment_1d(dom, idx_range, itmask_metal,
                                         my_chemistry, my_rates, my_fields,
                                         internal_dust_prop_buf);
  }

  // Look-up for H2 formation on dust

  if (anydust != MASK_FALSE) {
    d_logtem0 = std::log(my_chemistry->DustTemperatureStart);
    d_logtem9 = std::log(my_chemistry->DustTemperatureEnd);
    d_dlogtem = (std::log(my_chemistry->DustTemperatureEnd) -
                 std::log(my_chemistry->DustTemperatureStart)) /
                (double)(my_chemistry->NumberOfDustTemperatureBins - 1);

    if (my_chemistry->dust_species == 0) {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask_metal[i - 1] != MASK_FALSE) {
          // Assume dust melting at T > 1500 K

          if (tdust[i - 1] > my_chemistry->DustTemperatureEnd) {
            h2dust[i - 1] = tiny8;
          } else {
            // Get log dust temperature

            d_logtem[i - 1] = std::log(tdust[i - 1]);
            d_logtem[i - 1] = std::fmax(d_logtem[i - 1], d_logtem0);
            d_logtem[i - 1] = std::fmin(d_logtem[i - 1], d_logtem9);

            // Find index into table and precompute interpolation values

            d_indixe[i - 1] = std::fmin(
                my_chemistry->NumberOfDustTemperatureBins - 1,
                std::fmax(
                    1, (long long)((d_logtem[i - 1] - d_logtem0) / d_dlogtem) +
                           1));
            d_t1[i - 1] = (d_logtem0 + (d_indixe[i - 1] - 1) * d_dlogtem);
            d_t2[i - 1] = (d_logtem0 + (d_indixe[i - 1]) * d_dlogtem);
            d_tdef[i - 1] =
                (d_logtem[i - 1] - d_t1[i - 1]) / (d_t2[i - 1] - d_t1[i - 1]);

            // Get rate from 2D interpolation

            dusti1[i - 1] = h2dusta(logTlininterp_buf.indixe[i - 1] - 1,
                                    d_indixe[i - 1] - 1) +
                            (h2dusta(logTlininterp_buf.indixe[i - 1] + 1 - 1,
                                     d_indixe[i - 1] - 1) -
                             h2dusta(logTlininterp_buf.indixe[i - 1] - 1,
                                     d_indixe[i - 1] - 1)) *
                                logTlininterp_buf.tdef[i - 1];
            dusti2[i - 1] = h2dusta(logTlininterp_buf.indixe[i - 1] - 1,
                                    d_indixe[i - 1] + 1 - 1) +
                            (h2dusta(logTlininterp_buf.indixe[i - 1] + 1 - 1,
                                     d_indixe[i - 1] + 1 - 1) -
                             h2dusta(logTlininterp_buf.indixe[i - 1] - 1,
                                     d_indixe[i - 1] + 1 - 1)) *
                                logTlininterp_buf.tdef[i - 1];
            h2dust[i - 1] =
                dusti1[i - 1] + (dusti2[i - 1] - dusti1[i - 1]) * d_tdef[i - 1];

            // Multiply by dust to gas ratio

            h2dust[i - 1] = h2dust[i - 1] * dust2gas[i - 1];
          }
        }
      }

    } else {
      // create table for interpolation
      d_N[1 - 1] = (long long)(my_chemistry->NumberOfDustTemperatureBins);
      d_N[2 - 1] = (long long)(my_chemistry->NumberOfTemperatureBins);
      d_Size = d_N[1 - 1] * d_N[2 - 1];
      d_dTd = d_dlogtem;
      d_dTg = dlogtem;
      for (idratec = 1; idratec <= (my_chemistry->NumberOfDustTemperatureBins);
           idratec++) {
        d_Td[idratec - 1] = d_logtem0 + (double)(idratec - 1) * d_dlogtem;
      }
      for (iratec = 1; iratec <= (my_chemistry->NumberOfTemperatureBins);
           iratec++) {
        d_Tg[iratec - 1] = logtem0 + (double)(iratec - 1) * dlogtem;
      }

      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask_metal[i - 1] != MASK_FALSE) {
          if (my_chemistry->use_multiple_dust_temperatures == 0) {
            if (my_chemistry->dust_species > 0) {
              d_logtem[i - 1] = std::log(tdust[i - 1]);
              h2MgSiO3 = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              h2AC = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustC);

              h2dust[i - 1] = h2MgSiO3 * sgMgSiO3[i - 1] + h2AC * sgAC[i - 1];
            }
            if (my_chemistry->dust_species > 1) {
              h2dust[i - 1] =
                  h2dust[i - 1] + h2MgSiO3 * sgSiM[i - 1] +
                  h2MgSiO3 * sgFeM[i - 1] + h2MgSiO3 * sgMg2SiO4[i - 1] +
                  h2MgSiO3 * sgFe3O4[i - 1] + h2MgSiO3 * sgSiO2D[i - 1] +
                  h2MgSiO3 * sgMgO[i - 1] + h2MgSiO3 * sgFeS[i - 1] +
                  h2MgSiO3 * sgAl2O3[i - 1];
            }
            if (my_chemistry->dust_species > 2) {
              h2dust[i - 1] = h2dust[i - 1] + h2MgSiO3 * sgreforg[i - 1] +
                              h2MgSiO3 * sgvolorg[i - 1] +
                              h2MgSiO3 * sgH2Oice[i - 1];
            }

          } else {
            if (my_chemistry->dust_species > 0) {
              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i - 1]);
              h2MgSiO3 = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i - 1]);
              h2AC = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustC);
            }

            if (my_chemistry->dust_species > 1) {
              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i - 1]);
              h2SiM = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i - 1]);
              h2FeM = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1]);
              h2Mg2SiO4 = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1]);
              h2Fe3O4 = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i - 1]);
              h2SiO2D = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i - 1]);
              h2MgO = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i - 1]);
              h2FeS = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i - 1]);
              h2Al2O3 = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);
            }

            if (my_chemistry->dust_species > 2) {
              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][i - 1]);
              h2reforg = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i - 1]);
              h2volorg = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);

              d_logtem[i - 1] = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i - 1]);
              h2H2Oice = f_wrap::interpolate_2d_g(
                  d_logtem[i - 1], logTlininterp_buf.logtem[i - 1], d_N.data(),
                  d_Td.data(), d_dTd, d_Tg.data(), d_dTg, d_Size,
                  my_rates->h2dustS);
            }

            if (my_chemistry->dust_species > 0) {
              h2dust[i - 1] = h2MgSiO3 * sgMgSiO3[i - 1] + h2AC * sgAC[i - 1];
            }
            if (my_chemistry->dust_species > 1) {
              h2dust[i - 1] =
                  h2dust[i - 1] + h2SiM * sgSiM[i - 1] + h2FeM * sgFeM[i - 1] +
                  h2Mg2SiO4 * sgMg2SiO4[i - 1] + h2Fe3O4 * sgFe3O4[i - 1] +
                  h2SiO2D * sgSiO2D[i - 1] + h2MgO * sgMgO[i - 1] +
                  h2FeS * sgFeS[i - 1] + h2Al2O3 * sgAl2O3[i - 1];
            }
            if (my_chemistry->dust_species > 2) {
              h2dust[i - 1] = h2dust[i - 1] + h2reforg * sgreforg[i - 1] +
                              h2volorg * sgvolorg[i - 1] +
                              h2H2Oice * sgH2Oice[i - 1];
            }
          }
        }
      }

      // Compute grain growth rate

      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask_metal[i - 1] != MASK_FALSE) {
          if (my_chemistry->grain_growth == 1) {
            if (my_chemistry->dust_species > 0) {
              kd = f_wrap::interpolate_1d_g(
                  logTlininterp_buf.logtem[i - 1],
                  nratec_single_elem_arr.data(), d_Tg.data(), d_dTg,
                  nratec_single_elem_arr[0], my_rates->grain_growth_rate);

              grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i - 1] =
                  kd * sgMgSiO3[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  grackle::impl::fmin(
                      Mg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(24., 1.5),
                      SiOI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(44., 1.5),
                      H2O(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(18., 1.5) / 2.);
              // !             if ( idsub .eq. 1 )
              // !   &            kdMgSiO3  (i) = kdMgSiO3  (i) * ( 1.d0 -
              // !   &           sqrt(tMgSiO3 (i) / tgas1d(i))
              // !   &         * exp(-min(37.2400d4/tgas1d(i) -
              // 104.872d0, 5.d1)) / ( !   &           (     Mg
              // (i,j,k)*dom/24._DKIND * kboltz*tgas1d(i)) !   &         * (
              // SiOI(i,j,k)*dom/44._DKIND * kboltz*tgas1d(i)) !   &         *
              // (2.d0*H2O (i,j,k)*dom/18._DKIND * kboltz*tgas1d(i)) !   & /
              // (1.d0*H2I (i,j,k)*dom/ 2._DKIND * kboltz*tgas1d(i)) !   & ) )
              // !             Formulation from Nozawa et al. (2003, 2012)

              grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i - 1] =
                  kd * sgAC[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  CI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                  std::pow(12., 1.5);
            }

            if (my_chemistry->dust_species > 1) {
              grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i - 1] =
                  kd * sgSiM[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  SiI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                  std::pow(28., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i - 1] =
                  kd * sgFeM[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  Fe(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                  std::pow(56., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1] =
                  kd * sgMg2SiO4[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  grackle::impl::fmin(
                      Mg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(24., 1.5) / 2.,
                      SiOI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(44., 1.5),
                      H2O(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                          std::pow(18., 1.5) / 3.);

              grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1] =
                  kd * sgFe3O4[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  std::fmin(Fe(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(56., 1.5) / 3.,
                            H2O(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(18., 1.5) / 4.);

              grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i - 1] =
                  kd * sgSiO2D[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  SiO2I(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                  std::pow(60., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i - 1] =
                  kd * sgMgO[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  std::fmin(Mg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(24., 1.5),
                            H2O(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(18., 1.5));

              grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i - 1] =
                  kd * sgFeS[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  std::fmin(S(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(32., 1.5),
                            Fe(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(56., 1.5));

              grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i - 1] =
                  kd * sgAl2O3[i - 1] *
                  d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                  std::fmin(Al(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(27., 1.5) / 2.,
                            H2O(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                                std::pow(18., 1.5) / 3.);
            }

            // We do not consider the growth of refractory organics, volatile
            // organics, and water ice because their sublimation temperatures
            // are low (100-600 K). They sublimate before the growth occurs.
            if (my_chemistry->dust_species > 2) {
              grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i - 1] =
                  0.e0;
            }
          }

          if (my_chemistry->dust_sublimation == 1) {
            if (my_chemistry->dust_species > 0) {
              grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i - 1] = 0.e0;
            }
            if (my_chemistry->dust_species > 1) {
              grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i - 1] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i - 1] = 0.e0;
            }
            if (my_chemistry->dust_species > 2) {
              grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i - 1] =
                  0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i - 1] =
                  0.e0;
            }

            if (my_chemistry->use_multiple_dust_temperatures == 0) {
              if (my_chemistry->dust_species > 0) {
                if (tdust[i - 1] > 1222.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i - 1] =
                      (tiny8 -
                       MgSiO3(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1800.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i - 1] =
                      (tiny8 -
                       AC(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }

              if (my_chemistry->dust_species > 1) {
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i - 1] =
                      (tiny8 -
                       SiM(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i - 1] =
                      (tiny8 -
                       FeM(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1277.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1] =
                      (tiny8 -
                       Mg2SiO4(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1] =
                      (tiny8 -
                       Fe3O4(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i - 1] =
                      (tiny8 -
                       SiO2D(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i - 1] =
                      (tiny8 -
                       MgO(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 680.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i - 1] =
                      (tiny8 -
                       FeS(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i - 1] =
                      (tiny8 -
                       Al2O3(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }

              if (my_chemistry->dust_species > 2) {
                if (tdust[i - 1] > 575.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i - 1] =
                      (tiny8 -
                       reforg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 375.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i - 1] =
                      (tiny8 -
                       volorg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (tdust[i - 1] > 153.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i - 1] =
                      (tiny8 -
                       H2Oice(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }

            } else {
              if (my_chemistry->dust_species > 0) {
                if (grain_temperatures
                        .data[OnlyGrainSpLUT::MgSiO3_dust][i - 1] > 1222.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i - 1] =
                      (tiny8 -
                       MgSiO3(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i - 1] >
                    1800.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i - 1] =
                      (tiny8 -
                       AC(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }

              if (my_chemistry->dust_species > 1) {
                if (grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i - 1] =
                      (tiny8 -
                       SiM(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i - 1] =
                      (tiny8 -
                       FeM(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures
                        .data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1] > 1277.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i - 1] =
                      (tiny8 -
                       Mg2SiO4(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i - 1] =
                      (tiny8 -
                       Fe3O4(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i - 1] =
                      (tiny8 -
                       SiO2D(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i - 1] =
                      (tiny8 -
                       MgO(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i - 1] >
                    680.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i - 1] =
                      (tiny8 -
                       FeS(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i - 1] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i - 1] =
                      (tiny8 -
                       Al2O3(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }

              if (my_chemistry->dust_species > 2) {
                if (grain_temperatures
                        .data[OnlyGrainSpLUT::ref_org_dust][i - 1] > 575.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i - 1] =
                      (tiny8 -
                       reforg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures
                        .data[OnlyGrainSpLUT::vol_org_dust][i - 1] > 375.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i - 1] =
                      (tiny8 -
                       volorg(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
                if (grain_temperatures
                        .data[OnlyGrainSpLUT::H2O_ice_dust][i - 1] > 153.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i - 1] =
                      (tiny8 -
                       H2Oice(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) /
                      (*dt);
                }
              }
            }
          }
        }
      }
    }
  }

  // Include approximate self-shielding factors if requested

  for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
    if (itmask[i - 1] != MASK_FALSE) {
      kshield_buf.k24[i - 1] = my_uvb_rates.k24;
      kshield_buf.k25[i - 1] = my_uvb_rates.k25;
      kshield_buf.k26[i - 1] = my_uvb_rates.k26;
      kshield_buf.k28[i - 1] = my_uvb_rates.k28;
      kshield_buf.k29[i - 1] = my_uvb_rates.k29;
      kshield_buf.k30[i - 1] = my_uvb_rates.k30;
    }
  }

  // H2 self-shielding (Sobolev-like, spherically averaged, Wolcott-Green+ 2011)

  if (my_chemistry->primordial_chemistry > 1) {
    // If radiative transfer for LW photons have been already solved in
    // your hydro code, add kdissH2I later
    if (my_chemistry->use_radiative_transfer == 0 ||
        my_chemistry->radiative_transfer_use_H2_shielding == 1) {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          kshield_buf.k31[i - 1] = my_uvb_rates.k31;
        }
      }
    } else {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          kshield_buf.k31[i - 1] =
              my_uvb_rates.k31 +
              kdissH2I(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
        }
      }
    }

    if (my_chemistry->H2_self_shielding > 0) {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          // Calculate a Sobolev-like length assuming a 3D grid.
          if (my_chemistry->H2_self_shielding == 1) {
            divrhoa[1 - 1] =
                d(i + 1 - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrhoa[2 - 1] =
                d(i - 1 - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrhoa[3 - 1] =
                d(i - 1, idx_range.jp1 + 1 - 1, idx_range.kp1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrhoa[4 - 1] =
                d(i - 1, idx_range.jp1 - 1 - 1, idx_range.kp1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrhoa[5 - 1] =
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 + 1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrhoa[6 - 1] =
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1 - 1) -
                d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
            divrho = tiny_fortran_val;
            // Exclude directions with (drho/ds > 0)
            for (n1 = 1; n1 <= (6); n1++) {
              if (divrhoa[n1 - 1] < 0.) {
                divrho = divrho + divrhoa[n1 - 1];
              }
            }
            // (rho / divrho) is the Sobolev-like length in cell widths
            l_H2shield = std::fmin(
                dx_cgs * d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                    (2. * std::fabs(divrho)),
                internalu.xbase1);

            // User-supplied length-scale field.
          } else if (my_chemistry->H2_self_shielding == 2) {
            l_H2shield =
                xH2shield(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                internalu.uxyz;

            // Jeans Length
          } else if (my_chemistry->H2_self_shielding == 3) {
            l_H2shield = c_ljeans *
                         std::sqrt(tgas1d[i - 1] / (d(i - 1, idx_range.jp1 - 1,
                                                      idx_range.kp1 - 1) *
                                                    mmw[i - 1]));

          } else {
            l_H2shield = (gr_float)(0.);
          }

          N_H2 = dom * H2I(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
                 l_H2shield;

          // update: self-shielding following Wolcott-Green & Haiman (2019)
          // range of validity: T=100-8000 K, n<=1e7 cm^-3

          tgas_touse = std::fmax(tgas1d[i - 1], 1e2);
          tgas_touse = std::fmin(tgas_touse, 8e3);
          ngas_touse = d(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) * dom /
                       mmw[i - 1];
          ngas_touse = std::fmin(ngas_touse, 1e7);

          //_// PORT:                aWG2019 = (0.8711_DKIND *
          //_// PORT:      &              log10(tgas_touse) - 1.928_DKIND) *
          //_// PORT:      &              exp(-0.2856_DKIND * log10(ngas_touse))
          //+
          //_// PORT:      &              (-0.9639_DKIND * log10(tgas_touse)
          //+ 3.892_DKIND)

          x = 2.0e-15 * N_H2;
          b_doppler = 1e-5 * std::sqrt(2. * kboltz_grflt * tgas1d[i - 1] /
                                       (2. * mh_grflt));
          f_shield =
              0.965 / std::pow((1. + x / b_doppler), aWG2019) +
              0.035 * std::exp(-8.5e-4 * std::sqrt(1. + x)) / std::sqrt(1. + x);

          // avoid f>1
          f_shield = std::fmin(f_shield, 1.);

          kshield_buf.k31[i - 1] = f_shield * kshield_buf.k31[i - 1];
        }
      }
    }
    // If radiative transfer for LW photons have been already solved in
    // your hydro code, add kdissH2I here
    if (my_chemistry->use_radiative_transfer == 1 &&
        my_chemistry->radiative_transfer_use_H2_shielding == 1) {
      // write(*,*) 'kdissH2I included'
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          kshield_buf.k31[i - 1] =
              kshield_buf.k31[i - 1] +
              kdissH2I(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);
        }
      }
    }

    // Custom H2 shielding
    if (my_chemistry->H2_custom_shielding > 0) {
      for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
        if (itmask[i - 1] != MASK_FALSE) {
          kshield_buf.k31[i - 1] =
              f_shield_custom(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) *
              kshield_buf.k31[i - 1];
        }
      }
    }
  }

  if (my_chemistry->self_shielding_method > 0) {
    // Compute shielding factors
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        // Compute shielding factor for H
        nSSh = 6.73e-3 * std::pow((my_uvb_rates.crsHI / 2.49e-18), (-2. / 3.)) *
               std::pow((tgas1d[i - 1] / 1.0e4), (0.17)) *
               std::pow((my_uvb_rates.k24 / internalu.tbase1 / 1.0e-12),
                        (2.0 / 3.0));

        // Compute the total Hydrogen number density
        nratio = (HI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                  HII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1));
        if (my_chemistry->primordial_chemistry > 1) {
          nratio = nratio + HM(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                   H2I(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                   H2II(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1);

          if (my_chemistry->primordial_chemistry > 2) {
            nratio =
                nratio +
                0.5 * (DI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                       DII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) +
                2.0 * HDI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) / 3.0;
          }
        }

        nratio = nratio * dom / nSSh;

        f_shield_H[i - 1] =
            (0.98 * std::pow((1.0 + std::pow(nratio, (1.64))), (-2.28)) +
             0.02 * std::pow((1.0 + nratio), (-0.84)));

        // Compute shielding factor for He

        nSSh = 6.73e-3 *
               std::pow((my_uvb_rates.crsHeI / 2.49e-18), (-2. / 3.)) *
               std::pow((tgas1d[i - 1] / 1.0e4), (0.17)) *
               std::pow((my_uvb_rates.k26 / internalu.tbase1 / 1.0e-12),
                        (2.0 / 3.0));

        nratio = 0.25 *
                 (HeI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                  HeII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                  HeIII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)) *
                 dom / nSSh;

        f_shield_He[i - 1] =
            (0.98 * std::pow((1.0 + std::pow(nratio, (1.64))), (-2.28)) +
             0.02 * std::pow((1.0 + nratio), (-0.84)));
      }
    }
  }

  if (my_chemistry->self_shielding_method == 1) {
    // approximate self shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // to shield HI, while leaving HeI and HeII optically thin

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling. (rate k29)
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i - 1] = 0.;
        } else {
          kshield_buf.k24[i - 1] = kshield_buf.k24[i - 1] * f_shield_H[i - 1];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i - 1] = 0.;
        } else {
          kshield_buf.k29[i - 1] = kshield_buf.k29[i - 1] * f_shield_H[i - 1];
        }

        kshield_buf.k25[i - 1] = my_uvb_rates.k25;
        kshield_buf.k26[i - 1] = my_uvb_rates.k26;
      }
    }

  } else if (my_chemistry->self_shielding_method == 2) {
    // Better self-shielding in HI using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // approximate self shielding in HeI and HeII

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling as HI. (rate k29)

    //   Attenuate radiation rates for H2+ dissociation (30 eV)
    //   using same scaling as HeII. (rate k28 and k30)

    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i - 1] = 0.;
        } else {
          kshield_buf.k24[i - 1] = kshield_buf.k24[i - 1] * f_shield_H[i - 1];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i - 1] = 0.;
        } else {
          kshield_buf.k29[i - 1] = kshield_buf.k29[i - 1] * f_shield_H[i - 1];
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i - 1] = 0.;
        } else {
          kshield_buf.k26[i - 1] = kshield_buf.k26[i - 1] * f_shield_He[i - 1];
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i - 1] = 0.0;
        } else {
          kshield_buf.k28[i - 1] = kshield_buf.k28[i - 1] * f_shield_He[i - 1];
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i - 1] = 0.0;
        } else {
          kshield_buf.k30[i - 1] = kshield_buf.k30[i - 1] * f_shield_He[i - 1];
        }

        kshield_buf.k25[i - 1] = my_uvb_rates.k25;
      }
    }

  } else if (my_chemistry->self_shielding_method == 3) {
    // shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // in HI and HeI, but ignoring HeII heating entirely
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i - 1] = 0.;
        } else {
          kshield_buf.k24[i - 1] = kshield_buf.k24[i - 1] * f_shield_H[i - 1];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i - 1] = 0.;
        } else {
          kshield_buf.k29[i - 1] = kshield_buf.k29[i - 1] * f_shield_H[i - 1];
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i - 1] = 0.;
        } else {
          kshield_buf.k26[i - 1] = kshield_buf.k26[i - 1] * f_shield_He[i - 1];
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i - 1] = 0.0;
        } else {
          kshield_buf.k28[i - 1] = kshield_buf.k28[i - 1] * f_shield_He[i - 1];
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i - 1] = 0.0;
        } else {
          kshield_buf.k30[i - 1] = kshield_buf.k30[i - 1] * f_shield_He[i - 1];
        }

        kshield_buf.k25[i - 1] = 0.0;
      }
    }
  }

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
  // If using a high-energy radiation field, then account for
  //   effects of secondary electrons (Shull * Steenberg 1985)
  //   (see calc_rate.src)

  for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
    if (itmask[i - 1] != MASK_FALSE) {
      x = std::fmax(HII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) /
                        (HI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) +
                         HII(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1)),
                    1.0e-4);
      factor = 0.3908 * std::pow((1. - std::pow(x, 0.4092)), 1.7592);
      kshield_buf.k24[i - 1] =
          kshield_buf.k24[i - 1] +
          factor * (my_uvb_rates.piHI + 0.08 * my_uvb_rates.piHeI) /
              (e24 * everg) * internalu.coolunit * internalu.tbase1;
      factor = 0.0554 * std::pow((1. - std::pow(x, 0.4614)), 1.6660);
      kshield_buf.k26[i - 1] =
          kshield_buf.k26[i - 1] +
          factor * (my_uvb_rates.piHI / 0.08 + my_uvb_rates.piHeI) /
              (e26 * everg) * internalu.coolunit * internalu.tbase1;
    }
  }
#endif

  //   If using H2, and using the density-dependent collisional
  //     H2 dissociation rate, then replace the the density-independant
  //        k13 rate with the new one.
  // May/00: there appears to be a problem with the density-dependent
  //     collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
  if (my_chemistry->primordial_chemistry > 1 &&
      my_chemistry->three_body_rate == 0) {
    for (i = idx_range.i_start + 1; i <= (idx_range.i_end + 1); i++) {
      if (itmask[i - 1] != MASK_FALSE) {
        nh = std::fmin(HI(i - 1, idx_range.jp1 - 1, idx_range.kp1 - 1) * dom,
                       1.0e9);
        kcr_buf.data[CollisionalRxnLUT::k13][i - 1] = tiny8;
        if (tgas1d[i - 1] >= 500. && tgas1d[i - 1] < 1.0e6) {
          // Direct collisional dissociation
          k13_CID =
              k13dd(i - 1, 1 - 1) -
              k13dd(i - 1, 2 - 1) / (1. + std::pow((nh / k13dd(i - 1, 5 - 1)),
                                                   k13dd(i - 1, 7 - 1))) +
              k13dd(i - 1, 3 - 1) -
              k13dd(i - 1, 4 - 1) / (1. + std::pow((nh / k13dd(i - 1, 6 - 1)),
                                                   k13dd(i - 1, 7 - 1)));
          k13_CID = std::fmax(std::pow(10., k13_CID), tiny8);
          // Dissociative tunnelling
          k13_DT =
              k13dd(i - 1, 8 - 1) -
              k13dd(i - 1, 9 - 1) / (1. + std::pow((nh / k13dd(i - 1, 12 - 1)),
                                                   k13dd(i - 1, 14 - 1))) +
              k13dd(i - 1, 10 - 1) -
              k13dd(i - 1, 11 - 1) / (1. + std::pow((nh / k13dd(i - 1, 13 - 1)),
                                                    k13dd(i - 1, 14 - 1)));
          k13_DT = std::fmax(std::pow(10., k13_DT), tiny8);
          //
          kcr_buf.data[CollisionalRxnLUT::k13][i - 1] = k13_DT + k13_CID;
        }
      }
    }
  }
  // #define USE_PALLA_SALPETER_STAHLER1983
  // #ifdef USE_PALLA_SALPETER_STAHLER1983
  //            if (ispecies .gt. 1 .and. ithreebody .eq. 1) then
  //               do i = is+1, ie+1
  //                  if (itmask(i)) then
  //                  nh = (HI(i,j,k) + H2I(i,j,k)/2._DKIND)*dom
  //                  k13ind = 1._DKIND / (1._DKIND + nh / k13dd(i,3))
  //                  k13(i) = 10._DKIND**(
  //     &                     (1._DKIND-k13ind) * k13dd(i,2)
  //     &                             + k13ind  * k13dd(i,1) )
  //                  endif
  //               enddo
  //            endif
  // #endif
#endif  //  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE

  drop_InternalDustPropBuf(&internal_dust_prop_buf);

  return;
}

}  // namespace grackle::impl

#endif  // LOOKUP_COOL_RATES1D_HPP
