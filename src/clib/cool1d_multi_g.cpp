//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the cool1d_multi_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool1d_multi_g function from FORTRAN to C++

#include <cstdio>
#include <vector>
#include <iostream>

#include "cool1d_cloudy_g.hpp"
#include "cool1d_multi_g.hpp"
#include "grackle.h"
#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "dust_props.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

void grackle::impl::cool1d_multi_g(
    int imetal, int iter, double* edot, double* tgas, double* mmw, double* p2d,
    double* tdust, double* metallicity, double* dust2gas, double* rhoH,
    gr_mask_type* itmask, gr_mask_type* itmask_metal,
    chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
    grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
    InternalGrUnits internalu, IndexRange idx_range,
    grackle::impl::GrainSpeciesCollection grain_temperatures,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
    grackle::impl::CoolHeatScratchBuf coolingheating_buf) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(
      my_fields->internal_energy, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> de(
      my_fields->e_density, my_fields->grid_dimension[0],
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
  grackle::impl::View<gr_float***> HDI(
      my_fields->HDI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Vheat(
      my_fields->volumetric_heating_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mheat(
      my_fields->specific_heating_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Tfloor(
      my_fields->temperature_floor, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> photogamma(
      my_fields->RT_heating_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> isrf_habing(
      my_fields->isrf_habing, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(
      my_fields->CI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(
      my_fields->CII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(
      my_fields->CO_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(
      my_fields->OI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(
      my_fields->OH_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(
      my_fields->H2O_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Declare some constants:
  const double mh_local_var = mh_grflt;
  const double mu_metal = 16.;  // approx. mean molecular weight of metals

  // Locals
  int i, iZscale, mycmbTfloor;
  double dom, qq, vibl, logtem0, logtem9, dlogtem, zr, hdlte1, hdlow1, gamma2,
      x, fudge, gphdl1, dom_inv, tau, ciefudge, coolunit, tbase1, nH2, nother,
      nSSh, nratio, nssh_he, nratio_he, fSShHI, fSShHeI, pe_eps, pe_X, grbeta,
      ih2cox, min_metallicity;
  double comp1, comp2;

  // Performing heap allocations for all of the subsequent buffers within this
  // function is a major impediment to (i) CPU performance and (ii) adding GPU
  // support to Grackle. Future work should work on addressing this
  // - in the immediate short-term, we need to focus on aggregating these
  //   variables into logically organized structs. In a lot of cases, it
  //   may make sense to move the buffers into the existing
  //   Cool1DMultiScratchBuf or CoolHeatScratchBuf structs.
  // - in the longer term the goal is to refactor this logic to remove as many
  //   of these buffers as possible (without crippling cache performance on
  //   CPUs)

  std::vector<double> gaHI(my_fields->grid_dimension[0]);
  std::vector<double> gaH2(my_fields->grid_dimension[0]);
  std::vector<double> gaHe(my_fields->grid_dimension[0]);
  std::vector<double> gaHp(my_fields->grid_dimension[0]);
  std::vector<double> gael(my_fields->grid_dimension[0]);
  std::vector<double> h2lte(my_fields->grid_dimension[0]);
  std::vector<double> galdl(my_fields->grid_dimension[0]);
  // gas/grain heat transfer rate
  std::vector<double> gasgr(my_fields->grid_dimension[0]);
  // holds values of the interstellar radiation field
  std::vector<double> myisrf(my_fields->grid_dimension[0]);
  std::vector<double> cieY06(my_fields->grid_dimension[0]);

  std::vector<double> logT(my_fields->grid_dimension[0]);
  std::vector<double> logTcmb(my_fields->grid_dimension[0]);
  std::vector<double> logrho(my_fields->grid_dimension[0]);
  std::vector<double> logH(my_fields->grid_dimension[0]);
  std::vector<double> logH2I(my_fields->grid_dimension[0]);
  std::vector<double> logHDI(my_fields->grid_dimension[0]);
  std::vector<double> logH2(my_fields->grid_dimension[0]);
  std::vector<double> logCI(my_fields->grid_dimension[0]);
  std::vector<double> logCII(my_fields->grid_dimension[0]);
  std::vector<double> logOI(my_fields->grid_dimension[0]);
  std::vector<double> logCO(my_fields->grid_dimension[0]);
  std::vector<double> logOH(my_fields->grid_dimension[0]);
  std::vector<double> logH2O(my_fields->grid_dimension[0]);
  double lognhat;
  std::vector<double> logdvdr(my_fields->grid_dimension[0]);
  double log_Linv, log_Ginv, L, G;
  std::vector<double> Lpri(my_fields->grid_dimension[0]);
  std::vector<double> LH2(my_fields->grid_dimension[0]);
  std::vector<double> LCIE(my_fields->grid_dimension[0]);
  std::vector<double> LHD(my_fields->grid_dimension[0]);
  std::vector<double> LCI(my_fields->grid_dimension[0]);
  std::vector<double> LCII(my_fields->grid_dimension[0]);
  std::vector<double> LOI(my_fields->grid_dimension[0]);
  std::vector<double> LCO(my_fields->grid_dimension[0]);
  std::vector<double> LOH(my_fields->grid_dimension[0]);
  std::vector<double> LH2O(my_fields->grid_dimension[0]);
  std::vector<double> Ldst(my_fields->grid_dimension[0]);
  std::vector<double> alpha(my_fields->grid_dimension[0]);
  std::vector<double> alphad(my_fields->grid_dimension[0]);
  std::vector<double> lshield_con(my_fields->grid_dimension[0]);
  std::vector<double> tau_con(my_fields->grid_dimension[0]);
  double log_a;

  // buffers of intermediate quantities used within dust-routines (for
  // calculating quantites related to heating/cooling)
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf =
      grackle::impl::new_InternalDustPropBuf(my_fields->grid_dimension[0],
                                             my_rates->gr_N[1]);
  // opacity coefficients for each dust grain (the product of opacity
  // coefficient & gas mass density is the linear absortpion coefficient)
  grackle::impl::GrainSpeciesCollection grain_kappa =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);
  // closely related to grain_kappa
  std::vector<double> kappa_tot(my_fields->grid_dimension[0]);
  // holds the gas/grain-species heat transfer rates
  grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

  // Iteration mask

  gr_mask_type anydust, interp;
  std::vector<gr_mask_type> itmask_tab(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Set flag for dust-related options

  if ((my_chemistry->h2_on_dust > 0) || (my_chemistry->dust_chemistry > 0) ||
      (my_chemistry->dust_recombination_cooling > 0)) {
    anydust = MASK_TRUE;
  } else {
    anydust = MASK_FALSE;
  }

  // Set flag for needing interpolation variables

  if ((my_chemistry->primordial_chemistry > 0) ||
      (my_chemistry->dust_chemistry > 0)) {
    interp = MASK_TRUE;
  } else {
    interp = MASK_FALSE;
  }
  // Set log values of start and end of lookup tables

  logtem0 = std::log(my_chemistry->TemperatureStart);
  logtem9 = std::log(my_chemistry->TemperatureEnd);
  dlogtem = (std::log(my_chemistry->TemperatureEnd) -
             std::log(my_chemistry->TemperatureStart)) /
            (double)(my_chemistry->NumberOfTemperatureBins - 1);

  // Set units

  dom = internalu_calc_dom_(internalu);
  dom_inv = 1. / dom;
  tbase1 = internalu.tbase1;
  coolunit = internalu.coolunit;

  zr = 1. / (internalu.a_value * internalu.a_units) - 1.;
  fudge = 1.;

  // Set compton cooling coefficients (and temperature)

  comp1 = my_rates->comp * std::pow((1. + zr), 4);
  comp2 = 2.73 * (1. + zr);

  // multiplicative factor for including/excluding H2 cooling
  ih2cox = (double)(my_chemistry->ih2co);

  // ignore metal chemistry/cooling below this metallicity
  min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;

  // Initialize edot

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      edot[i] = 0.;
    }
  }

  // Compute Pressure

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      p2d[i] = (my_chemistry->Gamma - 1.) * d(i, idx_range.j, idx_range.k) *
               e(i, idx_range.j, idx_range.k);
    }
  }

  // Compute Temperature

  // If no chemistry, use a tabulated mean molecular weight
  // and iterate to convergence.

  if (my_chemistry->primordial_chemistry == 0) {
    // fh is H mass fraction in metal-free gas.

    if (imetal == 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          rhoH[i] = my_chemistry->HydrogenFractionByMass *
                    (d(i, idx_range.j, idx_range.k) -
                     metal(i, idx_range.j, idx_range.k));
        }
      }
    } else {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          rhoH[i] = my_chemistry->HydrogenFractionByMass *
                    d(i, idx_range.j, idx_range.k);
        }
      }
    }

    grackle::impl::fortran_wrapper::calc_temp1d_cloudy_g(
        rhoH, idx_range, tgas, mmw, dom, zr, imetal,
        my_rates->cloudy_primordial, itmask, my_chemistry, my_fields,
        internalu);

  } else {
    // Compute mean molecular weight (and temperature) directly

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        mmw[i] = (HeI(i, idx_range.j, idx_range.k) +
                  HeII(i, idx_range.j, idx_range.k) +
                  HeIII(i, idx_range.j, idx_range.k)) /
                     4. +
                 HI(i, idx_range.j, idx_range.k) +
                 HII(i, idx_range.j, idx_range.k) +
                 de(i, idx_range.j, idx_range.k);
        rhoH[i] =
            HI(i, idx_range.j, idx_range.k) + HII(i, idx_range.j, idx_range.k);
        cool1dmulti_buf.myde[i] = de(i, idx_range.j, idx_range.k);
      }
    }

    // (include molecular hydrogen, but ignore deuterium)

    if (my_chemistry->primordial_chemistry > 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          mmw[i] = mmw[i] + HM(i, idx_range.j, idx_range.k) +
                   (H2I(i, idx_range.j, idx_range.k) +
                    H2II(i, idx_range.j, idx_range.k)) /
                       2.;
          rhoH[i] = rhoH[i] + H2I(i, idx_range.j, idx_range.k) +
                    H2II(i, idx_range.j, idx_range.k);
        }
      }
    }

    // Include metal species

    if (imetal == 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          mmw[i] = mmw[i] + metal(i, idx_range.j, idx_range.k) / mu_metal;
        }
      }
    }

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        tgas[i] = std::fmax(p2d[i] * internalu.utem / mmw[i],
                            my_chemistry->TemperatureStart);
        mmw[i] = d(i, idx_range.j, idx_range.k) / mmw[i];
      }
    }

    // Correct temperature for gamma from H2

    if (my_chemistry->primordial_chemistry > 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          nH2 = 0.5 * (H2I(i, idx_range.j, idx_range.k) +
                       H2II(i, idx_range.j, idx_range.k));
          nother = (HeI(i, idx_range.j, idx_range.k) +
                    HeII(i, idx_range.j, idx_range.k) +
                    HeIII(i, idx_range.j, idx_range.k)) /
                       4. +
                   HI(i, idx_range.j, idx_range.k) +
                   HII(i, idx_range.j, idx_range.k) +
                   de(i, idx_range.j, idx_range.k);

          int iter_tgas = 0;
          double tgas_err = huge8;
          while ((iter_tgas < 100) && (tgas_err > 1.e-3)) {
            // tgas0 is used when CALCULATE_TGAS_SELF_CONSISTENTLY is defined
            [[maybe_unused]] double tgas0 = tgas[i];
            if (nH2 / nother > 1.0e-3) {
              x = 6100. / tgas[i];  // not quite self-consistent
              if (x > 10.) {
                gamma2 = 0.5 * 5.;
              } else {
                gamma2 = 0.5 * (5. + 2. * std::pow(x, 2) * std::exp(x) /
                                         std::pow((std::exp(x) - 1), 2));
              }
            } else {
              gamma2 = 2.5;
            }
            gamma2 =
                1. + (nH2 + nother) /
                         (nH2 * gamma2 + nother / (my_chemistry->Gamma - 1.));
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
            tgas[i] =
                std::fmax((gamma2 - 1.) * mmw[i] *
                              e(i, idx_range.j, idx_range.k) * internalu.utem,
                          my_chemistry->TemperatureStart);
            tgas_err = grackle::impl::dabs(tgas0 - tgas[i]) / tgas0;
            iter_tgas = iter_tgas + 1;
#else
            tgas[i] = tgas[i] * (gamma2 - 1.) / (my_chemistry->Gamma - 1.);
            iter_tgas = 101;
#endif
          }
        }
      }
    }
  }

  // Skip if below the temperature floor

  if (my_chemistry->use_temperature_floor == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] <= my_chemistry->temperature_floor_scalar) {
          edot[i] = tiny_fortran_val;
          itmask[i] = MASK_FALSE;
        }
      }
    }
  } else if (my_chemistry->use_temperature_floor == 2) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] <= Tfloor(i, idx_range.j, idx_range.k)) {
          edot[i] = tiny_fortran_val;
          itmask[i] = MASK_FALSE;
        }
      }
    }
  }

  // Calculate metallicity and H number density

  if (imetal == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        metallicity[i] = metal(i, idx_range.j, idx_range.k) /
                         d(i, idx_range.j, idx_range.k) /
                         my_chemistry->SolarMetalFractionByMass;
      }
    }
  } else {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        metallicity[i] = tiny_fortran_val;
      }
    }
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      cool1dmulti_buf.mynh[i] = rhoH[i] * dom;
    }
  }

  // If this is the first time through, just set tgasold to tgas

  if (iter == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        cool1dmulti_buf.tgasold[i] = tgas[i];
      }
    }
  }

  // Compute log densities

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      logT[i] = std::log10(tgas[i]);
      if (my_chemistry->cmb_temperature_floor == 1)
        logTcmb[i] = std::log10(comp2);
      logrho[i] = std::log10(d(i, idx_range.j, idx_range.k) * dom * mh);
      if (my_chemistry->primordial_chemistry > 0) {
        logH[i] = std::log10(HI(i, idx_range.j, idx_range.k) * dom);
        logH2[i] = std::log10(HI(i, idx_range.j, idx_range.k) * dom);
      }
      if (my_chemistry->primordial_chemistry > 1) {
        logH2[i] = std::log10((HI(i, idx_range.j, idx_range.k) +
                               H2I(i, idx_range.j, idx_range.k) / 2.0) *
                              dom);
        logH2I[i] = std::log10(H2I(i, idx_range.j, idx_range.k) * dom / 2.0);
      }
      if (my_chemistry->primordial_chemistry > 2) {
        logHDI[i] = std::log10(HDI(i, idx_range.j, idx_range.k) * dom / 3.0);
      }
      if (my_chemistry->metal_cooling == 1) {
        if (my_chemistry->metal_chemistry == 1) {
          logCI[i] = std::log10(CI(i, idx_range.j, idx_range.k) * dom / 12.0);
          logCII[i] = std::log10(CII(i, idx_range.j, idx_range.k) * dom / 12.0);
          logOI[i] = std::log10(OI(i, idx_range.j, idx_range.k) * dom / 16.0);
          logCO[i] = std::log10(CO(i, idx_range.j, idx_range.k) * dom / 28.0);
          logOH[i] = std::log10(OH(i, idx_range.j, idx_range.k) * dom / 17.0);
          logH2O[i] = std::log10(H2O(i, idx_range.j, idx_range.k) * dom / 18.0);
        }
      }

      // From Chiaki & Wise (2019), approximate dv/dr as 1/(3 * t_ff)
      logdvdr[i] = -8.79947961814e0 + 0.5e0 * logrho[i];  // km/s / cm
      lshield_con[i] = std::sqrt(
          (my_chemistry->Gamma * pi_fortran_val * kboltz_grflt * tgas[i]) /
          (GravConst_grflt * mmw[i] * mh_local_var *
           d(i, idx_range.j, idx_range.k) * dom * mh_local_var));
    }
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Compute log temperature and truncate if above/below table max/min

      logTlininterp_buf.logtem[i] =
          std::log(0.5 * (tgas[i] + cool1dmulti_buf.tgasold[i]));
      logTlininterp_buf.logtem[i] =
          std::fmax(logTlininterp_buf.logtem[i], logtem0);
      logTlininterp_buf.logtem[i] =
          std::fmin(logTlininterp_buf.logtem[i], logtem9);
    }
  }

  // Compute interpolation indices

  if (interp != MASK_FALSE) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Compute index into the table and precompute parts of linear interp

        logTlininterp_buf.indixe[i] = std::fmin(
            my_chemistry->NumberOfTemperatureBins - 1,
            std::fmax(1, (long long)((logTlininterp_buf.logtem[i] - logtem0) /
                                     dlogtem) +
                             1));
        logTlininterp_buf.t1[i] =
            (logtem0 + (logTlininterp_buf.indixe[i] - 1) * dlogtem);
        logTlininterp_buf.t2[i] =
            (logtem0 + (logTlininterp_buf.indixe[i]) * dlogtem);
        logTlininterp_buf.tdef[i] =
            (logTlininterp_buf.logtem[i] - logTlininterp_buf.t1[i]) /
            (logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
      }
    }
  }

  // --- 6 species cooling ---

  if (my_chemistry->primordial_chemistry > 0) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Lookup cooling values and do a linear temperature in log(T)

        coolingheating_buf.ceHI[i] =
            my_rates->ceHI[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ceHI[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ceHI[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ceHeI[i] =
            my_rates->ceHeI[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ceHeI[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ceHeI[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ceHeII[i] =
            my_rates->ceHeII[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ceHeII[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ceHeII[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ciHI[i] =
            my_rates->ciHI[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ciHI[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ciHI[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ciHeI[i] =
            my_rates->ciHeI[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ciHeI[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ciHeI[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ciHeIS[i] =
            my_rates->ciHeIS[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ciHeIS[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ciHeIS[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.ciHeII[i] =
            my_rates->ciHeII[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->ciHeII[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->ciHeII[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.reHII[i] =
            my_rates->reHII[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->reHII[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->reHII[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.reHeII1[i] =
            my_rates->reHeII1[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->reHeII1[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->reHeII1[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.reHeII2[i] =
            my_rates->reHeII2[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->reHeII2[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->reHeII2[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.reHeIII[i] =
            my_rates->reHeIII[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->reHeIII[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->reHeIII[logTlininterp_buf.indixe[i] - 1]);
        coolingheating_buf.brem[i] =
            my_rates->brem[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->brem[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->brem[logTlininterp_buf.indixe[i] - 1]);
      }
    }

    // Compute the cooling function

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        edot[i] = (

            // Collisional excitations

            -coolingheating_buf.ceHI[i] * HI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k)  // ce of HI
            - coolingheating_buf.ceHeI[i] * HeII(i, idx_range.j, idx_range.k) *
                  std::pow(de(i, idx_range.j, idx_range.k), 2) * dom /
                  4.  // ce of HeI
            - coolingheating_buf.ceHeII[i] * HeII(i, idx_range.j, idx_range.k) *
                  de(i, idx_range.j, idx_range.k) / 4.  // ce of HeII

            // Collisional ionizations

            - coolingheating_buf.ciHI[i] * HI(i, idx_range.j, idx_range.k) *
                  de(i, idx_range.j, idx_range.k)  // ci of HI
            - coolingheating_buf.ciHeI[i] * HeI(i, idx_range.j, idx_range.k) *
                  de(i, idx_range.j, idx_range.k) / 4.  // ci of HeI
            - coolingheating_buf.ciHeII[i] * HeII(i, idx_range.j, idx_range.k) *
                  de(i, idx_range.j, idx_range.k) / 4.  // ci of HeII
            - coolingheating_buf.ciHeIS[i] * HeII(i, idx_range.j, idx_range.k) *
                  std::pow(de(i, idx_range.j, idx_range.k), 2) * dom /
                  4.  // ci of HeIS

            // Recombinations

            - coolingheating_buf.reHII[i] * HII(i, idx_range.j, idx_range.k) *
                  de(i, idx_range.j, idx_range.k)  // re of HII
            -
            coolingheating_buf.reHeII1[i] * HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4.  // re of HeII
            -
            coolingheating_buf.reHeII2[i] * HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4.  // re of HeII
            -
            coolingheating_buf.reHeIII[i] * HeIII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4.  // re of HeIII

            // Bremsstrahlung

            - coolingheating_buf.brem[i] *
                  (HII(i, idx_range.j, idx_range.k) +
                   HeII(i, idx_range.j, idx_range.k) / 4. +
                   HeIII(i, idx_range.j, idx_range.k)) *
                  de(i, idx_range.j, idx_range.k)

        );
        Lpri[i] = edot[i];

        if (edot[i] != edot[i]) {
          OMP_PRAGMA_CRITICAL {
            eprintf("NaN in edot[1]:  %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
                    i, idx_range.j, idx_range.k, edot[i],
                    HI(i, idx_range.j, idx_range.k),
                    HII(i, idx_range.j, idx_range.k),
                    HeI(i, idx_range.j, idx_range.k),
                    HeII(i, idx_range.j, idx_range.k),
                    HeIII(i, idx_range.j, idx_range.k),
                    de(i, idx_range.j, idx_range.k),
                    d(i, idx_range.j, idx_range.k), tgas[i], p2d[i]);
          }
        }
      }
    }
  }

  // --- H2 cooling ---

  if (my_chemistry->primordial_chemistry > 1) {
    // Chiaki & Wise (2019) H2 cooling rate
    if (my_chemistry->h2_cooling_rate == 3) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          lognhat = logH2I[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH[i], my_rates->LH2.props.dimension,
              my_rates->LH2.props.parameters[0],
              my_rates->LH2.props.parameter_spacing[0],
              my_rates->LH2.props.parameters[1],
              my_rates->LH2.props.parameter_spacing[1],
              my_rates->LH2.props.parameters[2],
              my_rates->LH2.props.parameter_spacing[2],
              my_rates->LH2.props.data_size, my_rates->LH2.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH[i], my_rates->LH2.props.dimension,
                my_rates->LH2.props.parameters[0],
                my_rates->LH2.props.parameter_spacing[0],
                my_rates->LH2.props.parameters[1],
                my_rates->LH2.props.parameter_spacing[1],
                my_rates->LH2.props.parameters[2],
                my_rates->LH2.props.parameter_spacing[2],
                my_rates->LH2.props.data_size, my_rates->LH2.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LH2[i] =
              ih2cox * (G - L) / dom * H2I(i, idx_range.j, idx_range.k) / 2.e0;
          if (LH2[i] != LH2[i]) {
            LH2[i] = 0.e0;
          }
          edot[i] = edot[i] + LH2[i];
        }
      }

      // Glover & Abel (2008) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 2) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          gaHI[i] = my_rates->GAHI[logTlininterp_buf.indixe[i] - 1] +
                    logTlininterp_buf.tdef[i] *
                        (my_rates->GAHI[logTlininterp_buf.indixe[i] + 1 - 1] -
                         my_rates->GAHI[logTlininterp_buf.indixe[i] - 1]);
          gaH2[i] = my_rates->GAH2[logTlininterp_buf.indixe[i] - 1] +
                    logTlininterp_buf.tdef[i] *
                        (my_rates->GAH2[logTlininterp_buf.indixe[i] + 1 - 1] -
                         my_rates->GAH2[logTlininterp_buf.indixe[i] - 1]);
          gaHe[i] = my_rates->GAHe[logTlininterp_buf.indixe[i] - 1] +
                    logTlininterp_buf.tdef[i] *
                        (my_rates->GAHe[logTlininterp_buf.indixe[i] + 1 - 1] -
                         my_rates->GAHe[logTlininterp_buf.indixe[i] - 1]);
          gaHp[i] = my_rates->GAHp[logTlininterp_buf.indixe[i] - 1] +
                    logTlininterp_buf.tdef[i] *
                        (my_rates->GAHp[logTlininterp_buf.indixe[i] + 1 - 1] -
                         my_rates->GAHp[logTlininterp_buf.indixe[i] - 1]);
          gael[i] = my_rates->GAel[logTlininterp_buf.indixe[i] - 1] +
                    logTlininterp_buf.tdef[i] *
                        (my_rates->GAel[logTlininterp_buf.indixe[i] + 1 - 1] -
                         my_rates->GAel[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.gphdl[i] =
              my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i] +
                                                  1 - 1] -
                   my_rates
                       ->GP99HighDensityLimit[logTlininterp_buf.indixe[i] - 1]);
          h2lte[i] = my_rates->H2LTE[logTlininterp_buf.indixe[i] - 1] +
                     logTlininterp_buf.tdef[i] *
                         (my_rates->H2LTE[logTlininterp_buf.indixe[i] + 1 - 1] -
                          my_rates->H2LTE[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.cieco[i] =
              my_rates->cieco[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->cieco[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->cieco[logTlininterp_buf.indixe[i] - 1]);
        }
      }

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5 * H2I(i, idx_range.j, idx_range.k);
          nother = (HeI(i, idx_range.j, idx_range.k) +
                    HeII(i, idx_range.j, idx_range.k) +
                    HeIII(i, idx_range.j, idx_range.k)) /
                       4. +
                   HI(i, idx_range.j, idx_range.k) +
                   HII(i, idx_range.j, idx_range.k) +
                   de(i, idx_range.j, idx_range.k);
          double fH2 = nH2 / (nH2 + nother);
          fudge = std::sqrt(
              (40. *
               std::pow(
                   10.,
                   (4.8 * std::sqrt(std::max(std::log10(tgas[i]), 2.) - 2.))) /
               std::pow(fH2, 2)) /
              ((nH2 + nother) * dom));
          fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */
          // Note that this optical depth approximation comes from
          // RA04.
          if (my_chemistry->h2_optical_depth_approximation == 1) {
            fudge = std::pow(
                (0.76 * d(i, idx_range.j, idx_range.k) * dom / 8.e9), (-0.45));
            fudge = std::fmin(fudge, 1.);
          } else {
            fudge = 1.;
          }
          galdl[i] = gaHI[i] * HI(i, idx_range.j, idx_range.k) +
                     gaH2[i] * H2I(i, idx_range.j, idx_range.k) / 2. +
                     gaHe[i] * HeI(i, idx_range.j, idx_range.k) / 4. +
                     gaHp[i] * HII(i, idx_range.j, idx_range.k) +
                     gael[i] * de(i, idx_range.j, idx_range.k);
          // gphdl1 = gphdl(i)/dom
          gphdl1 = h2lte[i] / dom;
          edot[i] = edot[i] - ih2cox * fudge *
                                  H2I(i, idx_range.j, idx_range.k) * h2lte[i] /
                                  (1. + gphdl1 / galdl[i]) / (2. * dom);
        }
      }

      // Galli & Palla (1998) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          coolingheating_buf.gpldl[i] =
              my_rates->GP99LowDensityLimit[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->GP99LowDensityLimit[logTlininterp_buf.indixe[i] +
                                                 1 - 1] -
                   my_rates
                       ->GP99LowDensityLimit[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.gphdl[i] =
              my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i] +
                                                  1 - 1] -
                   my_rates
                       ->GP99HighDensityLimit[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.cieco[i] =
              my_rates->cieco[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->cieco[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->cieco[logTlininterp_buf.indixe[i] - 1]);
        }
      }

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
#define NO_OPTICAL_DEPTH_FUDGE
#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5 * H2I(i, idx_range.j, idx_range.k);
          nother = (HeI(i, idx_range.j, idx_range.k) +
                    HeII(i, idx_range.j, idx_range.k) +
                    HeIII(i, idx_range.j, idx_range.k)) /
                       4. +
                   HI(i, idx_range.j, idx_range.k) +
                   HII(i, idx_range.j, idx_range.k) +
                   de(i, idx_range.j, idx_range.k);
          double fH2 = nH2 / (nH2 + nother);
          // TODO: code duplication with line 726
          fudge = std::sqrt(
              (40. *
               std::pow(
                   10.,
                   (4.8 * std::sqrt(std::max(std::log10(tgas[i]), 2.) - 2.))) /
               std::pow(fH2, 2)) /
              ((nH2 + nother) * dom)) fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */
          // Note that this optical depth approximation comes from
          // RA04.
          if (my_chemistry->h2_optical_depth_approximation == 1) {
            fudge = std::pow(
                (0.76 * d(i, idx_range.j, idx_range.k) * dom / 8.e9), (-0.45));
            fudge = std::fmin(fudge, 1.);
          } else {
            fudge = 1.;
          }
          gphdl1 = coolingheating_buf.gphdl[i] /
                   (HI(i, idx_range.j, idx_range.k) * dom);
          edot[i] = edot[i] - ih2cox * fudge *
                                  H2I(i, idx_range.j, idx_range.k) *
                                  coolingheating_buf.gphdl[i] /
                                  (1. + gphdl1 / coolingheating_buf.gpldl[i]) /
                                  (2. * dom);
        }
      }

      // Lepp & Shull (1983) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 0) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          coolingheating_buf.hyd01k[i] =
              my_rates->hyd01k[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->hyd01k[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->hyd01k[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.h2k01[i] =
              my_rates->h2k01[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->h2k01[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->h2k01[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.vibh[i] =
              my_rates->vibh[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->vibh[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->vibh[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.roth[i] =
              my_rates->roth[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->roth[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->roth[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.rotl[i] =
              my_rates->rotl[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->rotl[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->rotl[logTlininterp_buf.indixe[i] - 1]);
          coolingheating_buf.cieco[i] =
              my_rates->cieco[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->cieco[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->cieco[logTlininterp_buf.indixe[i] - 1]);
        }
      }

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          qq = 1.2 * std::pow((HI(i, idx_range.j, idx_range.k) * dom), 0.77) +
               std::pow((H2I(i, idx_range.j, idx_range.k) * dom / 2.), 0.77);
          vibl =
              (HI(i, idx_range.j, idx_range.k) * coolingheating_buf.hyd01k[i] +
               H2I(i, idx_range.j, idx_range.k) / 2. *
                   coolingheating_buf.h2k01[i]) *
              dom * 8.18e-13;

#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5 * H2I(i, idx_range.j, idx_range.k);
          nother = (HeI(i, idx_range.j, idx_range.k) +
                    HeII(i, idx_range.j, idx_range.k) +
                    HeIII(i, idx_range.j, idx_range.k)) /
                       4. +
                   HI(i, idx_range.j, idx_range.k) +
                   HII(i, idx_range.j, idx_range.k) +
                   de(i, idx_range.j, idx_range.k);
          double fH2 = nH2 / (nH2 + nother);
          fudge = std::sqrt(
              (40. *
               std::pow(
                   10.,
                   (4.8 * std::sqrt(std::max(std::log10(tgas[i]), 2.) - 2.))) /
               std::pow(fH2, 2)) /
              ((nH2 + nother) * dom)) fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */

          edot[i] =
              edot[i] -
              ih2cox * fudge * H2I(i, idx_range.j, idx_range.k) *
                  (coolingheating_buf.vibh[i] /
                       (1. + coolingheating_buf.vibh[i] /
                                 std::fmax(vibl, tiny_fortran_val)) +
                   coolingheating_buf.roth[i] /
                       (1. + coolingheating_buf.roth[i] /
                                 std::fmax(qq * coolingheating_buf.rotl[i],
                                           tiny_fortran_val))) /
                  2. / dom;
        }
      }
    }

    // CIE
    // cooling from H2-H2 and He-H2 collisional induced emission comes
    //- with its own radiative transfer correction as discussed in
    //- Ripamonti & Abel 2003
    if (my_chemistry->cie_cooling == 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // Only calculate if H2I(i) is a substantial fraction
          if (d(i, idx_range.j, idx_range.k) * dom > 1e10) {
            ciefudge = 1.;
            tau = std::pow(((d(i, idx_range.j, idx_range.k) / 2e16) * dom),
                           2.8);  // 2e16 is in units of cm^-3
            tau = std::fmax(tau, 1.e-5);
            ciefudge = std::fmin((1. - std::exp(-tau)) / tau, 1.);
            // Matt's attempt at a second exponentialier cutoff
            tau = std::pow(((d(i, idx_range.j, idx_range.k) / 2.e18) * dom),
                           8.);  // 2e18 is in units of cm^-3
            tau = std::fmax(tau, 1.e-5);
            ciefudge = ciefudge * std::fmin((1.f - std::exp(-tau)) / tau, 1.);
            // ciefudge, which is applied to the continuum, is applied to edot
            edot[i] =
                ciefudge * (edot[i] - H2I(i, idx_range.j, idx_range.k) *
                                          (d(i, idx_range.j, idx_range.k) *
                                           coolingheating_buf.cieco[i]));
          }
        }
      }
      // CIE H2 cooling using Yoshida et al. (2006)
    } else if (my_chemistry->cie_cooling == 2) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          cieY06[i] =
              my_rates->cieY06[logTlininterp_buf.indixe[i] - 1] +
              logTlininterp_buf.tdef[i] *
                  (my_rates->cieY06[logTlininterp_buf.indixe[i] + 1 - 1] -
                   my_rates->cieY06[logTlininterp_buf.indixe[i] - 1]);
          LCIE[i] = -cieY06[i] *
                    std::pow((H2I(i, idx_range.j, idx_range.k) / 2.e0), 2);
          edot[i] = edot[i] + LCIE[i];
        }
      }
    }
  }

  // --- Cooling from HD ---

  if (my_chemistry->primordial_chemistry > 2) {
    // Chiaki & Wise (2019) HD cooling rate
    if (my_chemistry->hd_cooling_rate == 1) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          lognhat = logHDI[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH[i], my_rates->LHD.props.dimension,
              my_rates->LHD.props.parameters[0],
              my_rates->LHD.props.parameter_spacing[0],
              my_rates->LHD.props.parameters[1],
              my_rates->LHD.props.parameter_spacing[1],
              my_rates->LHD.props.parameters[2],
              my_rates->LHD.props.parameter_spacing[2],
              my_rates->LHD.props.data_size, my_rates->LHD.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH[i], my_rates->LHD.props.dimension,
                my_rates->LHD.props.parameters[0],
                my_rates->LHD.props.parameter_spacing[0],
                my_rates->LHD.props.parameters[1],
                my_rates->LHD.props.parameter_spacing[1],
                my_rates->LHD.props.parameters[2],
                my_rates->LHD.props.parameter_spacing[2],
                my_rates->LHD.props.data_size, my_rates->LHD.data);

            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LHD[i] = (G - L) / dom * HDI(i, idx_range.j, idx_range.k) / 3.e0;
          if (LHD[i] != LHD[i]) {
            LHD[i] = 0.e0;
          }
          edot[i] = edot[i] + LHD[i];
        }
      }

      // Coppola et al (2011) and Wrathmall, Gusdorf, & Flower (2007) HD cooling
      // rate
    } else if (my_chemistry->hd_cooling_rate == 0) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // CMB cooling floor
          if (tgas[i] > comp2) {
            coolingheating_buf.hdlte[i] =
                my_rates->HDlte[logTlininterp_buf.indixe[i] - 1] +
                logTlininterp_buf.tdef[i] *
                    (my_rates->HDlte[logTlininterp_buf.indixe[i] + 1 - 1] -
                     my_rates->HDlte[logTlininterp_buf.indixe[i] - 1]);
            coolingheating_buf.hdlow[i] =
                my_rates->HDlow[logTlininterp_buf.indixe[i] - 1] +
                logTlininterp_buf.tdef[i] *
                    (my_rates->HDlow[logTlininterp_buf.indixe[i] + 1 - 1] -
                     my_rates->HDlow[logTlininterp_buf.indixe[i] - 1]);
          } else {
            coolingheating_buf.hdlte[i] = tiny_fortran_val;
            coolingheating_buf.hdlow[i] = tiny_fortran_val;
          }
        }
      }

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // old (incorrect) way:
          //              hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2._DKIND)
          //              hdlow1 = max(hdlow(i), tiny)
          //              edot(i) = edot(i) - HDI(i,j,k)*
          //    .                     (hdlte1/(1._DKIND +
          //    hdlte1/hdlow1)/(2._DKIND*dom))
          // new (correct) way: (april 4, 2007)
          hdlte1 = coolingheating_buf.hdlte[i] /
                   (HI(i, idx_range.j, idx_range.k) * dom);
          hdlow1 = std::fmax(coolingheating_buf.hdlow[i], tiny_fortran_val);
          edot[i] = edot[i] -
                    HDI(i, idx_range.j, idx_range.k) *
                        (coolingheating_buf.hdlte[i] / (1. + hdlte1 / hdlow1)) /
                        (3. * dom);
        }
      }
    }
  }

  // Iteration mask for metal-rich cells
  if (imetal == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (metallicity[i] >= min_metallicity) {
        itmask_metal[i] = itmask[i];
      } else {
        itmask_metal[i] = MASK_FALSE;
      }
    }
  } else {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      itmask_metal[i] = MASK_FALSE;
    }
  }

  // Compute grain size increment
  if ((my_chemistry->use_dust_density_field > 0) &&
      (my_chemistry->dust_species > 0)) {
    grackle::impl::fortran_wrapper::calc_grain_size_increment_1d(
        dom, idx_range, itmask_metal, my_chemistry, my_rates, my_fields,
        internal_dust_prop_buf);
  }

  // Calculate dust to gas ratio AND interstellar radiation field
  // -> an earlier version of this logic would store values @ indices
  //    where `itmask_metal(i) .ne. MASK_FALSE`
  // -> this was undesirable, b/c these quantities are required for
  //    photo-electric heating, which can occur when
  //    `itmask_metal(i) .eq. MASK_FALSE` (we can revisit this choice
  //    later). Moreover, in most cases, these calculations will be
  //    faster when there is no branching

  if ((anydust != MASK_FALSE) || (my_chemistry->photoelectric_heating > 0)) {
    if (my_chemistry->use_dust_density_field > 0) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        // REMINDER: use of `itmask` over `itmask_metal` is
        //   currently required by Photo-electric heating
        if (itmask[i] != MASK_FALSE) {
          // it may be faster to remove this branching
          dust2gas[i] = dust(i, idx_range.j, idx_range.k) /
                        d(i, idx_range.j, idx_range.k);
        }
      }
    } else {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        dust2gas[i] = my_chemistry->local_dust_to_gas_ratio * metallicity[i];
      }
    }
  }

  if ((anydust != MASK_FALSE) || (my_chemistry->photoelectric_heating > 1)) {
    if (my_chemistry->use_isrf_field > 0) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        myisrf[i] = isrf_habing(i, idx_range.j, idx_range.k);
      }
    } else {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        myisrf[i] = my_chemistry->interstellar_radiation_field;
      }
    }
  }

  // compute dust temperature and cooling due to dust
  if (anydust != MASK_FALSE) {
    // TODO: trad -> comp2
    grackle::impl::fortran_wrapper::calc_all_tdust_gasgr_1d_g(
        comp2, tgas, tdust, metallicity, dust2gas, cool1dmulti_buf.mynh,
        cool1dmulti_buf.gasgr_tdust, itmask_metal, coolunit, gasgr.data(),
        myisrf.data(), kappa_tot.data(), my_chemistry, my_rates, my_fields,
        idx_range, grain_temperatures, gas_grainsp_heatrate, grain_kappa,
        logTlininterp_buf, internal_dust_prop_buf);
  }

  // Calculate dust cooling rate
  if (anydust != MASK_FALSE) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        if (my_chemistry->dust_species == 0) {
          Ldst[i] = -gasgr[i] * (tgas[i] - tdust[i]) * dust2gas[i] * rhoH[i] *
                    rhoH[i];

        } else {
          if (my_chemistry->use_multiple_dust_temperatures == 0) {
            Ldst[i] = -gasgr[i] * (tgas[i] - tdust[i]) *
                      d(i, idx_range.j, idx_range.k) * rhoH[i];
          } else {
            if (my_chemistry->dust_species > 0) {
              Ldst[i] =
                  -(gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] *
                        (tgas[i] - grain_temperatures
                                       .data[OnlyGrainSpLUT::MgSiO3_dust][i]) +
                    gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] *
                        (tgas[i] -
                         grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i])) *
                  d(i, idx_range.j, idx_range.k) * rhoH[i];
            }

            if (my_chemistry->dust_species > 1) {
              Ldst[i] =
                  Ldst[i] -
                  (gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] *
                       (tgas[i] -
                        grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] *
                       (tgas[i] -
                        grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::Mg2SiO4_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::Fe3O4_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] *
                       (tgas[i] -
                        grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] *
                       (tgas[i] -
                        grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] *
                       (tgas[i] -
                        grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::Al2O3_dust][i])) *
                      d(i, idx_range.j, idx_range.k) * rhoH[i];
            }

            if (my_chemistry->dust_species > 2) {
              Ldst[i] =
                  Ldst[i] -
                  (gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::ref_org_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::vol_org_dust][i]) +
                   gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] *
                       (tgas[i] - grain_temperatures
                                      .data[OnlyGrainSpLUT::H2O_ice_dust][i])) *
                      d(i, idx_range.j, idx_range.k) * rhoH[i];
            }
          }
        }

        edot[i] = edot[i] + Ldst[i];
      }
    }
  }

  // Compute continuum opacity

  if (my_chemistry->use_primordial_continuum_opacity == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        // ! primordial continuum opacity !!
        log_a = grackle::impl::fortran_wrapper::interpolate_2d_g(
            logrho[i], logT[i], my_rates->alphap.props.dimension,
            my_rates->alphap.props.parameters[0],
            my_rates->alphap.props.parameter_spacing[0],
            my_rates->alphap.props.parameters[1],
            my_rates->alphap.props.parameter_spacing[1],
            my_rates->alphap.props.data_size, my_rates->alphap.data);

        alpha[i] = std::pow(1.e1, log_a);
      }
    }

  } else {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        alpha[i] = 0.f;
      }
    }
  }

  // Add contributions from dust opacity to alpha, the linear absorption
  // coefficient
  //
  //  The original Fortran version of this function had the following 2
  //  comments:
  //    ! if (idspecies .eq. 0), dust opacity is overestimated at Td > 50 K
  //    ! We better not include dust opacity.
  // It's a little unclear how relevant these comments actually are.
  if ((anydust != MASK_FALSE) && (my_chemistry->dust_species > 0)) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        if (my_chemistry->use_multiple_dust_temperatures == 0) {
          // In the future, we should consider renaming `alphad`. The
          // current name is a little confusing since:
          // - the related `alpha` variable holds linear absorption
          //   coefficients (which is commonly denoted by the Greek
          //   letter alpha)
          // - in contrast, `alphad` only ever holds the sum of
          //   opacity coefficients (commonly denoted by the Greek
          //   letter kappa)

          alphad[i] = kappa_tot[i];

        } else {
          if (my_chemistry->dust_species > 0) {
            alphad[i] = grain_kappa.data[OnlyGrainSpLUT::MgSiO3_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::AC_dust][i];
          }
          if (my_chemistry->dust_species > 1) {
            alphad[i] = alphad[i] +
                        grain_kappa.data[OnlyGrainSpLUT::SiM_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::FeM_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::Fe3O4_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::SiO2_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::MgO_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::FeS_dust][i] +
                        grain_kappa.data[OnlyGrainSpLUT::Al2O3_dust][i];
          }
          if (my_chemistry->dust_species > 2) {
            alphad[i] =
                alphad[i] +
                gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] +
                gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] +
                gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i];
          }
        }

        alpha[i] = alpha[i] + alphad[i] * d(i, idx_range.j, idx_range.k) * dom *
                                  mh_local_var;
      }
    }
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      tau_con[i] = alpha[i] * lshield_con[i];
    }
  }

  // --- Compute (external) radiative heating terms ---
  // Photoionization heating

  if (my_chemistry->primordial_chemistry > 0) {
    if (my_chemistry->self_shielding_method == 0) {  // no shielding
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          edot[i] =
              edot[i] +
              (double)(my_chemistry->ipiht) *
                  (my_uvb_rates.piHI *
                       HI(i, idx_range.j, idx_range.k)  // pi of HI
                   + my_uvb_rates.piHeI * HeI(i, idx_range.j, idx_range.k) *
                         0.25  // pi of HeI
                   + my_uvb_rates.piHeII * HeII(i, idx_range.j, idx_range.k) *
                         0.25  // pi of HeII
                   ) /
                  dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 1) {
      // approximate self shielding using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // to shield HI, while leaving HeI and HeII optically thin

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          if (my_uvb_rates.k24 < tiny8) {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
                   std::pow((my_uvb_rates.crsHI / 2.49e-18), (-2. / 3.)) *
                   std::pow((tgas[i] / 1.0e4), (0.17)) *
                   std::pow((my_uvb_rates.k24 / tbase1 / 1.0e-12), (2. / 3.));
            nratio = (HI(i, idx_range.j, idx_range.k) +
                      HII(i, idx_range.j, idx_range.k)) *
                     dom / nSSh;
            fSShHI = 0.98 * std::pow((1. + std::pow(nratio, (1.64))), (-2.28)) +
                     0.02 * std::pow((1. + nratio), (-0.84));
          }

          edot[i] =
              edot[i] + (double)(my_chemistry->ipiht) *
                            (my_uvb_rates.piHI *
                                 HI(i, idx_range.j, idx_range.k) * fSShHI +
                             my_uvb_rates.piHeI *
                                 HeI(i, idx_range.j, idx_range.k) * 0.25 +
                             my_uvb_rates.piHeII *
                                 HeII(i, idx_range.j, idx_range.k) * 0.25) /
                            dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 2) {
      // Better self-shielding in HI using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // approximate self shielding in HeI and HeII

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // HI self shielding ratio
          if (my_uvb_rates.k24 < tiny8) {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
                   std::pow((my_uvb_rates.crsHI / 2.49e-18), (-2. / 3.)) *
                   std::pow((tgas[i] / 1.0e4), (0.17)) *
                   std::pow((my_uvb_rates.k24 / tbase1 / 1.0e-12), (2. / 3.));
            nratio = (HI(i, idx_range.j, idx_range.k) +
                      HII(i, idx_range.j, idx_range.k)) *
                     dom / nSSh;
            fSShHI = 0.98 * std::pow((1. + std::pow(nratio, (1.64))), (-2.28)) +
                     0.02 * std::pow((1. + nratio), (-0.84));
          }

          // HeI self shielding ratio
          if (my_uvb_rates.k26 < tiny8) {
            fSShHeI = 1.;
          } else {
            nssh_he =
                6.73e-3 *
                std::pow((my_uvb_rates.crsHeI / 2.49e-18), (-2. / 3.)) *
                std::pow((tgas[i] / 1.0e4), (0.17)) *
                std::pow((my_uvb_rates.k26 / tbase1 / 1.0e-12), (2. / 3.));
            nratio_he = 0.25 *
                        (HeI(i, idx_range.j, idx_range.k) +
                         HeII(i, idx_range.j, idx_range.k) +
                         HeIII(i, idx_range.j, idx_range.k)) *
                        dom / nssh_he;
            fSShHeI =
                0.98 * std::pow((1. + std::pow(nratio_he, (1.64))), (-2.28)) +
                0.02 * std::pow((1. + nratio_he), (-0.84));
          }

          edot[i] = edot[i] +
                    (double)(my_chemistry->ipiht) *
                        (my_uvb_rates.piHI * HI(i, idx_range.j, idx_range.k) *
                             fSShHI +
                         my_uvb_rates.piHeI * HeI(i, idx_range.j, idx_range.k) *
                             0.25 * fSShHeI +
                         my_uvb_rates.piHeII *
                             HeII(i, idx_range.j, idx_range.k) * 0.25) /
                        dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 3) {
      // shielding using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // in HI and HeI, but ignoring HeII heating entirely

      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // HI self shielding ratio
          if (my_uvb_rates.k24 < tiny8) {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
                   std::pow((my_uvb_rates.crsHI / 2.49e-18), (-2. / 3.)) *
                   std::pow((tgas[i] / 1.0e4), (0.17)) *
                   std::pow((my_uvb_rates.k24 / tbase1 / 1.0e-12), (2. / 3.));
            nratio = (HI(i, idx_range.j, idx_range.k) +
                      HII(i, idx_range.j, idx_range.k)) *
                     dom / nSSh;
            fSShHI = 0.98 * std::pow((1. + std::pow(nratio, (1.64))), (-2.28)) +
                     0.02 * std::pow((1. + nratio), (-0.84));
          }

          // HeI self shielding ratio
          if (my_uvb_rates.k26 < tiny8) {
            fSShHeI = 1.;
          } else {
            nssh_he =
                6.73e-3 *
                std::pow((my_uvb_rates.crsHeI / 2.49e-18), (-2. / 3.)) *
                std::pow((tgas[i] / 1.0e4), (0.17)) *
                std::pow((my_uvb_rates.k26 / tbase1 / 1.0e-12), (2. / 3.));
            nratio_he = 0.25 *
                        (HeI(i, idx_range.j, idx_range.k) +
                         HeII(i, idx_range.j, idx_range.k) +
                         HeIII(i, idx_range.j, idx_range.k)) *
                        dom / nssh_he;
            fSShHeI =
                0.98 * std::pow((1. + std::pow(nratio_he, (1.64))), (-2.28)) +
                0.02 * std::pow((1. + nratio_he), (-0.84));
          }

          edot[i] =
              edot[i] + (double)(my_chemistry->ipiht) *
                            (my_uvb_rates.piHI *
                                 HI(i, idx_range.j, idx_range.k) * fSShHI +
                             my_uvb_rates.piHeI *
                                 HeI(i, idx_range.j, idx_range.k) * fSShHeI) /
                            dom;

          // Ignoring HeII heating (HeII heating rate -> 0)
        }
      }
    }
  }

  // --- Cloudy primordial cooling and heating ---

  if (my_chemistry->primordial_chemistry == 0) {
    iZscale = 0;
    mycmbTfloor = 0;
    grackle::impl::cool1d_cloudy_g(
        rhoH, metallicity, logTlininterp_buf.logtem, edot,
        comp2, dom, zr, mycmbTfloor, my_chemistry->UVbackground, iZscale,
        my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension, my_rates->cloudy_primordial.grid_parameters[0],
        my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.data_size,
        *my_rates->cloudy_primordial.cooling_data, *my_rates->cloudy_primordial.heating_data, itmask, my_fields,
        idx_range);


    // Calculate electron density from mean molecular weight

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        cool1dmulti_buf.myde[i] =
            1 -
            mmw[i] * (3.0 * my_chemistry->HydrogenFractionByMass + 1.0) / 4.0;
        if (imetal == 1) {
          cool1dmulti_buf.myde[i] =
              cool1dmulti_buf.myde[i] -
              mmw[i] * metal(i, idx_range.j, idx_range.k) /
                  (d(i, idx_range.j, idx_range.k) * mu_metal);
        }
        cool1dmulti_buf.myde[i] =
            d(i, idx_range.j, idx_range.k) * cool1dmulti_buf.myde[i] / mmw[i];
        cool1dmulti_buf.myde[i] = std::fmax(cool1dmulti_buf.myde[i], 0.);
      }
    }
  }

  // Photo-electric heating by UV-irradiated dust

  if (my_chemistry->photoelectric_heating == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] > 2.e4) {
          cool1dmulti_buf.gammaha_eff[i] = 0.;
        } else {
          cool1dmulti_buf.gammaha_eff[i] = my_rates->gammah;
        }
      }
    }

    // Use eqn. 1 of Wolfire et al. (1995)
  } else if (my_chemistry->photoelectric_heating == 2) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (tgas[i] > 2.e4) {
          cool1dmulti_buf.gammaha_eff[i] = 0.;
        } else {
          // Assume constant epsilon = 0.05.
          cool1dmulti_buf.gammaha_eff[i] = my_rates->gammah * 0.05 * myisrf[i];
        }
      }
    }

    // Full calculation of epsilon (eqn. 2 of Wolfire 1995)
  } else if (my_chemistry->photoelectric_heating == 3) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        pe_X =
            myisrf[i] * dom_inv * std::sqrt(tgas[i]) / cool1dmulti_buf.myde[i];
        pe_eps = (4.9e-2 / (1. + std::pow((pe_X / 1925.), 0.73))) +
                 ((3.7e-2 * std::pow((tgas[i] / 1.e4), 0.7)) /
                  (1. + (pe_X / 5000.)));
        cool1dmulti_buf.gammaha_eff[i] = my_rates->gammah * pe_eps * myisrf[i];
      }
    }
  }

  if (my_chemistry->photoelectric_heating > 0) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        edot[i] = edot[i] + cool1dmulti_buf.gammaha_eff[i] * rhoH[i] * dom_inv *
                                dust2gas[i] /
                                my_chemistry->local_dust_to_gas_ratio;
      }
    }
  }

  // Electron recombination onto dust grains (eqn. 9 of Wolfire 1995)

  if ((my_chemistry->dust_chemistry > 0) ||
      (my_chemistry->dust_recombination_cooling > 0)) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        cool1dmulti_buf.regr[i] =
            my_rates->regr[logTlininterp_buf.indixe[i] - 1] +
            logTlininterp_buf.tdef[i] *
                (my_rates->regr[logTlininterp_buf.indixe[i] + 1 - 1] -
                 my_rates->regr[logTlininterp_buf.indixe[i] - 1]);
      }
    }

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        grbeta = 0.74 / std::pow(tgas[i], 0.068);
        edot[i] = edot[i] -
                  cool1dmulti_buf.regr[i] *
                      std::pow((myisrf[i] * dom_inv / cool1dmulti_buf.myde[i]),
                               grbeta) *
                      cool1dmulti_buf.myde[i] * rhoH[i] * dust2gas[i] /
                      my_chemistry->local_dust_to_gas_ratio;
      }
    }
  }

  // Compton cooling or heating and X-ray compton heating

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      edot[i] = edot[i]

                // Compton cooling or heating

                - comp1 * (tgas[i] - comp2) * cool1dmulti_buf.myde[i] * dom_inv

                // X-ray compton heating

                - my_uvb_rates.comp_xray * (tgas[i] - my_uvb_rates.temp_xray) *
                      cool1dmulti_buf.myde[i] * dom_inv;
    }
  }

  // Photoheating from radiative transfer

  if (my_chemistry->use_radiative_transfer == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        edot[i] = edot[i] + (double)(my_chemistry->ipiht) *
                                photogamma(i, idx_range.j, idx_range.k) /
                                coolunit * HI(i, idx_range.j, idx_range.k) /
                                dom;

        if (edot[i] != edot[i]) {
          OMP_PRAGMA_CRITICAL {
            eprintf(
                "NaN in edot[2]:  %d %d %d %g %g %g %g %g %g %g %g %g %g %g "
                "%g\n",
                i, idx_range.j, idx_range.k, edot[i],
                photogamma(i, idx_range.j, idx_range.k),
                HI(i, idx_range.j, idx_range.k),
                de(i, idx_range.j, idx_range.k), d(i, idx_range.j, idx_range.k),
                e(i, idx_range.j, idx_range.k), p2d[i], tgas[i], dom,
                internalu.urho, internalu.a_value, mh_local_var);
          }
        }
      }
    }
  }

  // --- Cloudy metal cooling and heating ---

  if (my_chemistry->metal_cooling == 1) {
    // Determine if the temperature is above the threshold to do tabulated
    // cooling.
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      itmask_tab[i] = itmask_metal[i];
      if (itmask_tab[i] != MASK_FALSE) {
        if ((my_chemistry->tabulated_cooling_minimum_temperature > 0.0e0) &&
            (tgas[i] < my_chemistry->tabulated_cooling_minimum_temperature)) {
          itmask_tab[i] = MASK_FALSE;
        }
      }
    }

    if (my_rates->cloudy_data_new == 1) {
      iZscale = 1;
      grackle::impl::cool1d_cloudy_g(
          rhoH, metallicity, logTlininterp_buf.logtem, edot,
          comp2, dom, zr, my_chemistry->cmb_temperature_floor, my_chemistry->UVbackground, iZscale,
          my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension, my_rates->cloudy_metal.grid_parameters[0],
          my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.data_size,
          *my_rates->cloudy_metal.cooling_data, *my_rates->cloudy_metal.heating_data, itmask_tab.data(), my_fields,
          idx_range);

    } else {
      FORTRAN_NAME(cool1d_cloudy_old_tables_g)(
          d.data(), de.data(), rhoH, metallicity, &my_fields->grid_dimension[0],
          &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
          &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
          logTlininterp_buf.logtem, edot, &comp2,
          &my_chemistry->primordial_chemistry, &dom, &zr,
          &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground,
          &my_chemistry->cloudy_electron_fraction_factor,
          &my_rates->cloudy_metal.grid_rank,
          my_rates->cloudy_metal.grid_dimension,
          my_rates->cloudy_metal.grid_parameters[0],
          my_rates->cloudy_metal.grid_parameters[1],
          my_rates->cloudy_metal.grid_parameters[2],
          my_rates->cloudy_metal.grid_parameters[3],
          my_rates->cloudy_metal.grid_parameters[4],
          &my_rates->cloudy_metal.data_size,
          my_rates->cloudy_metal.cooling_data,
          my_rates->cloudy_metal.heating_data, itmask_tab.data());
    }

    if (my_chemistry->metal_chemistry == 1) {
      // --- C/O fine-structure, metal molecular rotational cooling for low
      // temperatures ---

      // C/O fine-structure cooling
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          // CI
          lognhat = logCI[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH[i], my_rates->LCI.props.dimension,
              my_rates->LCI.props.parameters[0],
              my_rates->LCI.props.parameter_spacing[0],
              my_rates->LCI.props.parameters[1],
              my_rates->LCI.props.parameter_spacing[1],
              my_rates->LCI.props.parameters[2],
              my_rates->LCI.props.parameter_spacing[2],
              my_rates->LCI.props.data_size, my_rates->LCI.data);

          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH[i], my_rates->LCI.props.dimension,
                my_rates->LCI.props.parameters[0],
                my_rates->LCI.props.parameter_spacing[0],
                my_rates->LCI.props.parameters[1],
                my_rates->LCI.props.parameter_spacing[1],
                my_rates->LCI.props.parameters[2],
                my_rates->LCI.props.parameter_spacing[2],
                my_rates->LCI.props.data_size, my_rates->LCI.data);

            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LCI[i] = (G - L) / dom * CI(i, idx_range.j, idx_range.k) / 12.e0;
          if (LCI[i] != LCI[i]) {
            LCI[i] = 0.e0;
          }
          edot[i] = edot[i] + LCI[i];

          // CII
          lognhat = logCII[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH[i], my_rates->LCII.props.dimension,
              my_rates->LCII.props.parameters[0],
              my_rates->LCII.props.parameter_spacing[0],
              my_rates->LCII.props.parameters[1],
              my_rates->LCII.props.parameter_spacing[1],
              my_rates->LCII.props.parameters[2],
              my_rates->LCII.props.parameter_spacing[2],
              my_rates->LCII.props.data_size, my_rates->LCII.data);

          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH[i], my_rates->LCII.props.dimension,
                my_rates->LCII.props.parameters[0],
                my_rates->LCII.props.parameter_spacing[0],
                my_rates->LCII.props.parameters[1],
                my_rates->LCII.props.parameter_spacing[1],
                my_rates->LCII.props.parameters[2],
                my_rates->LCII.props.parameter_spacing[2],
                my_rates->LCII.props.data_size, my_rates->LCII.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LCII[i] = (G - L) / dom * CII(i, idx_range.j, idx_range.k) / 12.e0;
          if (LCII[i] != LCII[i]) {
            LCII[i] = 0.e0;
          }
          edot[i] = edot[i] + LCII[i];

          // OI
          lognhat = logOI[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH[i], my_rates->LOI.props.dimension,
              my_rates->LOI.props.parameters[0],
              my_rates->LOI.props.parameter_spacing[0],
              my_rates->LOI.props.parameters[1],
              my_rates->LOI.props.parameter_spacing[1],
              my_rates->LOI.props.parameters[2],
              my_rates->LOI.props.parameter_spacing[2],
              my_rates->LOI.props.data_size, my_rates->LOI.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH[i], my_rates->LOI.props.dimension,
                my_rates->LOI.props.parameters[0],
                my_rates->LOI.props.parameter_spacing[0],
                my_rates->LOI.props.parameters[1],
                my_rates->LOI.props.parameter_spacing[1],
                my_rates->LOI.props.parameters[2],
                my_rates->LOI.props.parameter_spacing[2],
                my_rates->LOI.props.data_size, my_rates->LOI.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LOI[i] = (G - L) / dom * OI(i, idx_range.j, idx_range.k) / 16.e0;
          if (LOI[i] != LOI[i]) {
            LOI[i] = 0.e0;
          }
          edot[i] = edot[i] + LOI[i];

          // Metal molecules rotational cooling

          // CO
          lognhat = logCO[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH2[i], my_rates->LCO.props.dimension,
              my_rates->LCO.props.parameters[0],
              my_rates->LCO.props.parameter_spacing[0],
              my_rates->LCO.props.parameters[1],
              my_rates->LCO.props.parameter_spacing[1],
              my_rates->LCO.props.parameters[2],
              my_rates->LCO.props.parameter_spacing[2],
              my_rates->LCO.props.data_size, my_rates->LCO.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH2[i], my_rates->LCO.props.dimension,
                my_rates->LCO.props.parameters[0],
                my_rates->LCO.props.parameter_spacing[0],
                my_rates->LCO.props.parameters[1],
                my_rates->LCO.props.parameter_spacing[1],
                my_rates->LCO.props.parameters[2],
                my_rates->LCO.props.parameter_spacing[2],
                my_rates->LCO.props.data_size, my_rates->LCO.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LCO[i] = (G - L) / dom * CO(i, idx_range.j, idx_range.k) / 28.e0;
          if (LCO[i] != LCO[i]) {
            LCO[i] = 0.e0;
          }
          edot[i] = edot[i] + LCO[i];

          // OH
          lognhat = logOH[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH2[i], my_rates->LOH.props.dimension,
              my_rates->LOH.props.parameters[0],
              my_rates->LOH.props.parameter_spacing[0],
              my_rates->LOH.props.parameters[1],
              my_rates->LOH.props.parameter_spacing[1],
              my_rates->LOH.props.parameters[2],
              my_rates->LOH.props.parameter_spacing[2],
              my_rates->LOH.props.data_size, my_rates->LOH.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH2[i], my_rates->LOH.props.dimension,
                my_rates->LOH.props.parameters[0],
                my_rates->LOH.props.parameter_spacing[0],
                my_rates->LOH.props.parameters[1],
                my_rates->LOH.props.parameter_spacing[1],
                my_rates->LOH.props.parameters[2],
                my_rates->LOH.props.parameter_spacing[2],
                my_rates->LOH.props.data_size, my_rates->LOH.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LOH[i] = (G - L) / dom * OH(i, idx_range.j, idx_range.k) / 17.e0;
          if (LOH[i] != LOH[i]) {
            LOH[i] = 0.e0;
          }
          edot[i] = edot[i] + LOH[i];

          // H2O
          lognhat = logH2O[i] - logdvdr[i];

          log_Linv = grackle::impl::fortran_wrapper::interpolate_3d_g(
              lognhat, logT[i], logH2[i], my_rates->LH2O.props.dimension,
              my_rates->LH2O.props.parameters[0],
              my_rates->LH2O.props.parameter_spacing[0],
              my_rates->LH2O.props.parameters[1],
              my_rates->LH2O.props.parameter_spacing[1],
              my_rates->LH2O.props.parameters[2],
              my_rates->LH2O.props.parameter_spacing[2],
              my_rates->LH2O.props.data_size, my_rates->LH2O.data);
          L = std::pow(1.e1, (-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1) {
            log_Ginv = grackle::impl::fortran_wrapper::interpolate_3d_g(
                lognhat, logTcmb[i], logH2[i], my_rates->LH2O.props.dimension,
                my_rates->LH2O.props.parameters[0],
                my_rates->LH2O.props.parameter_spacing[0],
                my_rates->LH2O.props.parameters[1],
                my_rates->LH2O.props.parameter_spacing[1],
                my_rates->LH2O.props.parameters[2],
                my_rates->LH2O.props.parameter_spacing[2],
                my_rates->LH2O.props.data_size, my_rates->LH2O.data);
            G = std::pow(1.e1, (-log_Ginv));
          } else {
            G = tiny8;
          }

          LH2O[i] = (G - L) / dom * H2O(i, idx_range.j, idx_range.k) / 18.e0;
          if (LH2O[i] != LH2O[i]) {
            LH2O[i] = 0.e0;
          }
          edot[i] = edot[i] + LH2O[i];
        }
      }
    }
  }

  // Add user-provided volumetric and/or specific heating terms

  if (my_chemistry->use_volumetric_heating_rate == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        edot[i] = edot[i] + Vheat(i, idx_range.j, idx_range.k) / coolunit /
                                std::pow(dom, 2);
      }
    }
  }

  if (my_chemistry->use_specific_heating_rate == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        edot[i] = edot[i] + Mheat(i, idx_range.j, idx_range.k) *
                                d(i, idx_range.j, idx_range.k) * mh_local_var /
                                coolunit / dom;
      }
    }
  }

  // Continuum opacity

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      if (tau_con[i] > 1.e0) {
        if (tau_con[i] < 1.e2) {
          edot[i] = edot[i] * std::pow(tau_con[i], (-2.e0));
        } else {
          edot[i] = 0.e0;
        }
      }
    }
  }

  // Set tgasold

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      cool1dmulti_buf.tgasold[i] = tgas[i];
    }
  }

  // Free memory
  grackle::impl::drop_InternalDustPropBuf(&internal_dust_prop_buf);
  grackle::impl::drop_GrainSpeciesCollection(&grain_kappa);
  grackle::impl::drop_GrainSpeciesCollection(&gas_grainsp_heatrate);

  return;
}