// See LICENSE file for license and copyright information

/// @file cool1d_multi_g-cpp.C
/// @brief Declares signature of cool1d_multi_g

// This file was initially generated automatically during conversion of the
// cool1d_multi_g function from FORTRAN to C++

#include <cstdio>
#include <vector>
#include <iostream>
#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "cool1d_multi_g-cpp.h"
//TODO: manually removed from signature (double* comp1, double* comp2)
void cool1d_multi_g(
  int imetal, int iter, //double* comp1, double* comp2, 
  double* edot,
  double* tgas, double* mmw, double* p2d, double* tdust, double* metallicity,
  double* dust2gas, double* rhoH, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields,
  photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
  IndexRange idx_range,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf
)
{

  // SOLVE RADIATIVE COOLING/HEATING EQUATIONS
  
  // written by: Yu Zhang, Peter Anninos and Tom Abel
  // date:
  // modified1: January, 1996 by Greg Bryan; adapted to KRONOS
  // modified2: October, 1996 by GB; moved to AMR
  // modified3: February, 2003 by Robert Harkness; iteration mask
  // modified6: September, 2009 by BDS to include cloudy cooling
  
  // PURPOSE:
  //   Solve the energy cooling equations.
  
  // INPUTS:
  //   is,ie   - start and end indicies of active region (zero-based!)
  
  // PARAMETERS:

  // -----------------------------------------------------------------------


  // Arguments

  // -- removed line (previously just declared arg types) -- 

  // -- removed line (previously just declared arg types) -- 
  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(my_fields->HeI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(my_fields->HeII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(my_fields->HeIII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(my_fields->HM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(my_fields->HDI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(my_fields->metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(my_fields->dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Vheat(my_fields->volumetric_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mheat(my_fields->specific_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Tfloor(my_fields->temperature_floor, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> photogamma(my_fields->RT_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> isrf_habing(my_fields->isrf_habing, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(my_fields->CI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(my_fields->CII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(my_fields->CO_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(my_fields->OI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(my_fields->OH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(my_fields->H2O_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  // -- removed line (previously just declared arg types) --

  // Cloudy cooling data

  int iZscale, mycmbTfloor;
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 

  // Parameters

  const double mh_local_var = mh_grflt;
  // approx. mean molecular weight of metals
  const double mu_metal = 16.;
  const int ti_max = 20;

  // Locals

  int i, ti, iradfield;
  double dom, qq, vibl, logtem0, logtem9, dlogtem, zr, hdlte1, hdlow1, gamma2, x, fudge, fH2, gphdl1, dom_inv, tau, ciefudge, coolunit, dbase1, tbase1, xbase1, nH2, nother, nSSh, nratio, nssh_he, nratio_he, fSShHI, fSShHeI, pe_eps, pe_X, grbeta, ih2cox, min_metallicity;
  int itd;
  double comp1, comp2; //TODO: manually added and removed from signature

  // Slice locals
 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 

  // Cooling/heating slice locals

  std::vector<double> gaHI(my_fields->grid_dimension[0]);
  std::vector<double> gaH2(my_fields->grid_dimension[0]);
  std::vector<double> gaHe(my_fields->grid_dimension[0]);
  std::vector<double> gaHp(my_fields->grid_dimension[0]);
  std::vector<double> gael(my_fields->grid_dimension[0]);
  std::vector<double> h2lte(my_fields->grid_dimension[0]);
  std::vector<double> galdl(my_fields->grid_dimension[0]);
  std::vector<double> gasgr(my_fields->grid_dimension[0]);
  std::vector<double> myisrf(my_fields->grid_dimension[0]);
  int iden, item, itab;
  std::vector<double> cieY06(my_fields->grid_dimension[0]);
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  // opacity table
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  double logdom;
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
  double log_a, log_L, log_G;
  // grain growth
  std::vector<double> sgSiM(my_fields->grid_dimension[0]);
  std::vector<double> sgFeM(my_fields->grid_dimension[0]);
  std::vector<double> sgMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> sgMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> sgFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> sgAC(my_fields->grid_dimension[0]);
  std::vector<double> sgSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> sgMgO(my_fields->grid_dimension[0]);
  std::vector<double> sgFeS(my_fields->grid_dimension[0]);
  std::vector<double> sgAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> sgreforg(my_fields->grid_dimension[0]);
  std::vector<double> sgvolorg(my_fields->grid_dimension[0]);
  std::vector<double> sgH2Oice(my_fields->grid_dimension[0]);
  std::vector<double> sgtot(my_fields->grid_dimension[0]);
  std::vector<double> alSiM(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFeM(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMg2SiO4(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMgSiO3(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFe3O4(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alAC(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alSiO2D(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMgO(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFeS(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alAl2O3(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alreforg(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alvolorg(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alH2Oice(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> altot(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> kpSiM(my_fields->grid_dimension[0]);
  std::vector<double> kpFeM(my_fields->grid_dimension[0]);
  std::vector<double> kpMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> kpMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> kpFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> kpAC(my_fields->grid_dimension[0]);
  std::vector<double> kpSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> kpMgO(my_fields->grid_dimension[0]);
  std::vector<double> kpFeS(my_fields->grid_dimension[0]);
  std::vector<double> kpAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> kpreforg(my_fields->grid_dimension[0]);
  std::vector<double> kpvolorg(my_fields->grid_dimension[0]);
  std::vector<double> kpH2Oice(my_fields->grid_dimension[0]);
  std::vector<double> kptot(my_fields->grid_dimension[0]);
  // grain temperature
  // -- removed line (previously just declared arg types) -- 
  // -- removed line (previously just declared arg types) -- 
  std::vector<double> gasSiM(my_fields->grid_dimension[0]);
  std::vector<double> gasFeM(my_fields->grid_dimension[0]);
  std::vector<double> gasMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> gasMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> gasFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> gasAC(my_fields->grid_dimension[0]);
  std::vector<double> gasSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> gasMgO(my_fields->grid_dimension[0]);
  std::vector<double> gasFeS(my_fields->grid_dimension[0]);
  std::vector<double> gasAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> gasreforg(my_fields->grid_dimension[0]);
  std::vector<double> gasvolorg(my_fields->grid_dimension[0]);
  std::vector<double> gasH2Oice(my_fields->grid_dimension[0]);
  // Iteration mask

  gr_mask_type anydust, interp;
  std::vector<gr_mask_type> itmask_tab(my_fields->grid_dimension[0]);
  // !#define CALCULATE_TGAS_SELF_CONSISTENTLY
  // #ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
  int iter_tgas;
  double tgas_err, tgas0;
  // #endif /* NOT important */
  //      debug
  double edotunit;
  int i_max;
  gr_float d_max;

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Set flag for dust-related options

  if ((my_chemistry->h2_on_dust > 0)  ||  (my_chemistry->dust_chemistry > 0)  || 
      (my_chemistry->dust_recombination_cooling > 0))  {
    anydust = MASK_TRUE;
  } else {
    anydust = MASK_FALSE;
  }

  // Set flag for needing interpolation variables

  if ((my_chemistry->primordial_chemistry > 0)  ||  (my_chemistry->dust_chemistry > 0))  {
    interp = MASK_TRUE;
  } else {
    interp = MASK_FALSE;
  }
  // Set log values of start and end of lookup tables

  logtem0 = std::log(my_chemistry->TemperatureStart);
  logtem9 = std::log(my_chemistry->TemperatureEnd);
  dlogtem= (std::log(my_chemistry->TemperatureEnd) - std::log(my_chemistry->TemperatureStart))/(double)(my_chemistry->NumberOfTemperatureBins-1 );

  // Set units

  dom      = internalu.urho*(std::pow(internalu.a_value,3))/mh_local_var;
  dom_inv  = 1./dom;
  tbase1   = internalu.tbase1;
  xbase1   = internalu.uxyz/(internalu.a_value*internalu.a_units);    // uxyz is [x]*a      = [x]*[a]*a'        '
  dbase1   = internalu.urho*std::pow((internalu.a_value*internalu.a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '
  coolunit = (std::pow(internalu.a_units,5) * std::pow(xbase1,2) * std::pow(mh_local_var,2)) / (std::pow(tbase1,3) * dbase1);
  zr       = 1./(internalu.a_value*internalu.a_units) - 1.;
  fudge    = 1.;
  iradfield = -1;

  // Set compton cooling coefficients (and temperature)

  // (*comp1) = my_rates->comp  * std::pow((1. + zr),4);
  // (*comp2) = 2.73 * (1. + zr);
  comp1 = my_rates->comp  * std::pow((1. + zr),4);
  comp2 = 2.73 * (1. + zr);

  // multiplicative factor for including/excluding H2 cooling
  ih2cox = (double)(my_chemistry->ih2co );

  // ignore metal chemistry/cooling below this metallicity
  min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;

  // Initialize edot

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      edot[i-1] = 0.;
    }
  }

  // Compute Pressure

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      p2d[i-1] = (my_chemistry->Gamma - 1.)*d(i-1,idx_range.jp1-1,idx_range.kp1-1)*e(i-1,idx_range.jp1-1,idx_range.kp1-1);
    }
  }

  // Compute Temperature

  // If no chemistry, use a tabulated mean molecular weight
  // and iterate to convergence.

  if (my_chemistry->primordial_chemistry == 0)  {

    // fh is H mass fraction in metal-free gas.

    if (imetal == 1)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          rhoH[i-1] = my_chemistry->HydrogenFractionByMass * (d(i-1,idx_range.jp1-1,idx_range.kp1-1) - metal(i-1,idx_range.jp1-1,idx_range.kp1-1));
        }
      }
    } else {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          rhoH[i-1] = my_chemistry->HydrogenFractionByMass * d(i-1,idx_range.jp1-1,idx_range.kp1-1);
        }
      }
    }

     FORTRAN_NAME(calc_temp1d_cloudy_g)(d.data(), metal.data(), e.data(), rhoH,
         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
         tgas, mmw, &dom, &zr,
         &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
         &my_chemistry->Gamma, &internalu.utem, &imetal,
         &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
         my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2],
         &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.mmw_data,
         itmask);

  } else {

    // Compute mean molecular weight (and temperature) directly

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        mmw[i-1] =
             (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))/4. +
             HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1) + de(i-1,idx_range.jp1-1,idx_range.kp1-1);
        rhoH[i-1] = HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1);
        cool1dmulti_buf.myde[i-1] = de(i-1,idx_range.jp1-1,idx_range.kp1-1);
      }
    }

    // (include molecular hydrogen, but ignore deuterium)

    if (my_chemistry->primordial_chemistry > 1)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          mmw[i-1] = mmw[i-1] +
               HM(i-1,idx_range.jp1-1,idx_range.kp1-1) + (H2I(i-1,idx_range.jp1-1,idx_range.kp1-1) + H2II(i-1,idx_range.jp1-1,idx_range.kp1-1))/2.;
          rhoH[i-1] = rhoH[i-1] + H2I(i-1,idx_range.jp1-1,idx_range.kp1-1) + H2II(i-1,idx_range.jp1-1,idx_range.kp1-1);
        }
      }
    }

    // Include metal species

    if (imetal == 1)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          mmw[i-1] = mmw[i-1] + metal(i-1,idx_range.jp1-1,idx_range.kp1-1)/mu_metal;
        }
      }
    }

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        tgas[i-1] = std::fmax(p2d[i-1]*internalu.utem/mmw[i-1], my_chemistry->TemperatureStart);
        mmw[i-1] = d(i-1,idx_range.jp1-1,idx_range.kp1-1) / mmw[i-1];
      }
    }

    // Correct temperature for gamma from H2

    if (my_chemistry->primordial_chemistry > 1)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          nH2 = 0.5*(H2I(i-1,idx_range.jp1-1,idx_range.kp1-1) + H2II(i-1,idx_range.jp1-1,idx_range.kp1-1));
          nother = (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) +
               HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))/4. +
               HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1) + de(i-1,idx_range.jp1-1,idx_range.kp1-1);

          iter_tgas = 0;
          tgas_err = huge8;
          while ((iter_tgas < 100)
                && (tgas_err > 1.e-3)) {
            tgas0 = tgas[i-1];
            if (nH2/nother > 1.0e-3)  {
              x = 6100./tgas[i-1]; // not quite self-consistent
              if (x > 10.)  {
                gamma2 = 0.5*5.;
              } else {
                gamma2 = 0.5*(5. + 2.*std::pow(x,2)
                  * std::exp(x)/std::pow((std::exp(x)-1),2));
              }
            } else {
              gamma2 = 2.5;
            }
            gamma2 = 1. + (nH2 + nother)/
                 (nH2*gamma2 + nother/(my_chemistry->Gamma-1.));
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
            tgas[i-1] = std::fmax((gamma2 - 1.)*mmw[i-1]*e(i-1,idx_range.jp1-1,idx_range.kp1-1)*
                 internalu.utem, my_chemistry->TemperatureStart);
            tgas_err = grackle::impl::dabs(tgas0 - tgas[i-1]) / tgas0;
            iter_tgas = iter_tgas + 1;
#else
            tgas[i-1] = tgas[i-1] * (gamma2 - 1.)/
                 (my_chemistry->Gamma - 1.);
            iter_tgas = 101;
#endif
          }
        }
      }
    }

  }

  // Skip if below the temperature floor

  if (my_chemistry->use_temperature_floor == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        if (tgas[i-1] <= my_chemistry->temperature_floor_scalar)  {
          edot[i-1] = tiny_fortran_val;
          itmask[i-1] = MASK_FALSE;
        }
      }
    }
  } else if (my_chemistry->use_temperature_floor == 2)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        if (tgas[i-1] <= Tfloor(i-1,idx_range.jp1-1,idx_range.kp1-1))  {
          edot[i-1] = tiny_fortran_val;
          itmask[i-1] = MASK_FALSE;
        }
      }
    }
  }

  // Calculate metallicity and H number density

  if (imetal == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        metallicity[i-1] = metal(i-1,idx_range.jp1-1,idx_range.kp1-1) / d(i-1,idx_range.jp1-1,idx_range.kp1-1) / my_chemistry->SolarMetalFractionByMass;
      }
    }
  } else {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        metallicity[i-1] = tiny_fortran_val;
      }
    }
  }

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      cool1dmulti_buf.mynh[i-1] = rhoH[i-1] * dom;
    }
  }

  // If this is the first time through, just set tgasold to tgas

  if (iter == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        cool1dmulti_buf.tgasold[i-1] = tgas[i-1];
      }
    }
  }

  // Compute log densities

  //_// PORT:       logdom = log10(dom)
  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
  //_// PORT:             logT(i)   = log10(tgas(i))
  //_// PORT:             if(icmbTfloor .eq. 1)
  //_// PORT:      &         logTcmb(i) = log10(comp2)
  //_// PORT:             logrho(i) = log10(d(i,j,k) * dom*mh)
      if( my_chemistry->primordial_chemistry > 0)  {
  //_// PORT:                logH(i)   = log10(HI(i,j,k) * dom)
  //_// PORT:                logH2(i)  = log10(HI(i,j,k) * dom)
      }
      if( my_chemistry->primordial_chemistry > 1 )  {
  //_// PORT:               logH2(i)  = log10((HI(i,j,k) + H2I(i,j,k) / 2.d0) * dom)
  //_// PORT:               logH2I(i) = log10(H2I(i,j,k) * dom /  2.d0)
      }
      if( my_chemistry->primordial_chemistry > 2)  {
  //_// PORT:               logHDI(i) = log10(HDI(i,j,k) * dom /  3.d0)
      }
      if( my_chemistry->metal_cooling == 1 )  {
        if( my_chemistry->metal_chemistry == 1 )  {
  //_// PORT:               logCI(i)  = log10(CI (i,j,k) * dom / 12.d0)
  //_// PORT:               logCII(i) = log10(CII(i,j,k) * dom / 12.d0)
  //_// PORT:               logOI(i)  = log10(OI (i,j,k) * dom / 16.d0)
  //_// PORT:               logCO(i)  = log10(CO (i,j,k) * dom / 28.d0)
  //_// PORT:               logOH(i)  = log10(OH (i,j,k) * dom / 17.d0)
  //_// PORT:               logH2O(i) = log10(H2O(i,j,k) * dom / 18.d0)
        }
      }

      // From Chiaki & Wise (2019), approximate dv/dr as 1/(3 * t_ff)
      logdvdr[i-1] = -8.79947961814e0 + 0.5e0 * logrho[i-1]; // km/s / cm
      lshield_con[i-1] =
         std::sqrt((my_chemistry->Gamma * pi_fortran_val * kboltz_grflt * tgas[i-1]) /
         (GravConst_grflt *  mmw[i-1]*mh_local_var * d(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom*mh_local_var));

    }
  }
     
  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {

      // Compute log temperature and truncate if above/below table max/min

      logTlininterp_buf.logtem[i-1] = std::log(0.5*(tgas[i-1]+cool1dmulti_buf.tgasold[i-1]));
      logTlininterp_buf.logtem[i-1] = std::fmax(logTlininterp_buf.logtem[i-1], logtem0);
      logTlininterp_buf.logtem[i-1] = std::fmin(logTlininterp_buf.logtem[i-1], logtem9);

    }
  }

  // Compute interpolation indices

  if (interp != MASK_FALSE)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {

        // Compute index into the table and precompute parts of linear interp

        logTlininterp_buf.indixe[i-1] = std::fmin(my_chemistry->NumberOfTemperatureBins-1,
                std::fmax(1,(long long)((logTlininterp_buf.logtem[i-1]-logtem0)/dlogtem )+1));
        logTlininterp_buf.t1[i-1] = (logtem0 + (logTlininterp_buf.indixe[i-1] - 1)*dlogtem);
        logTlininterp_buf.t2[i-1] = (logtem0 + (logTlininterp_buf.indixe[i-1]    )*dlogtem);
        logTlininterp_buf.tdef[i-1] = (logTlininterp_buf.logtem[i-1] - logTlininterp_buf.t1[i-1]) / (logTlininterp_buf.t2[i-1] - logTlininterp_buf.t1[i-1]);

      }
    }
  }

  // --- 6 species cooling ---

  if (my_chemistry->primordial_chemistry > 0)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {

        // Lookup cooling values and do a linear temperature in log(T)

        coolingheating_buf.ceHI[i-1] = my_rates->ceHI[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ceHI[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ceHI[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ceHeI[i-1] = my_rates->ceHeI[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ceHeI[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ceHeI[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ceHeII[i-1] = my_rates->ceHeII[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ceHeII[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ceHeII[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ciHI[i-1] = my_rates->ciHI[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ciHI[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ciHI[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ciHeI[i-1] = my_rates->ciHeI[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ciHeI[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ciHeI[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ciHeIS[i-1] = my_rates->ciHeIS[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ciHeIS[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ciHeIS[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.ciHeII[i-1] = my_rates->ciHeII[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->ciHeII[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->ciHeII[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.reHII[i-1] = my_rates->reHII[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->reHII[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->reHII[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.reHeII1[i-1]=my_rates->reHeII1[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->reHeII1[logTlininterp_buf.indixe[i-1]+1-1]-my_rates->reHeII1[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.reHeII2[i-1]=my_rates->reHeII2[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->reHeII2[logTlininterp_buf.indixe[i-1]+1-1]-my_rates->reHeII2[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.reHeIII[i-1]=my_rates->reHeIII[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->reHeIII[logTlininterp_buf.indixe[i-1]+1-1]-my_rates->reHeIII[logTlininterp_buf.indixe[i-1]-1]);
        coolingheating_buf.brem[i-1] = my_rates->brem[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
              *(my_rates->brem[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->brem[logTlininterp_buf.indixe[i-1]-1]);

      }
    }

    // Compute the cooling function

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        edot[i-1] = (

        // Collisional excitations

                  - coolingheating_buf.ceHI  [i-1]*HI  (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)              // ce of HI
                  - coolingheating_buf.ceHeI [i-1]*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*std::pow(de(i-1,idx_range.jp1-1,idx_range.kp1-1),2)*dom/4.  // ce of HeI
                  - coolingheating_buf.ceHeII[i-1]*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.         // ce of HeII

        // Collisional ionizations

                  - coolingheating_buf.ciHI  [i-1]*HI  (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)              // ci of HI
                  - coolingheating_buf.ciHeI [i-1]*HeI (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.         // ci of HeI
                  - coolingheating_buf.ciHeII[i-1]*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.         // ci of HeII
                  - coolingheating_buf.ciHeIS[i-1]*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*std::pow(de(i-1,idx_range.jp1-1,idx_range.kp1-1),2)*dom/4.  // ci of HeIS

        // Recombinations

                  - coolingheating_buf.reHII  [i-1]*HII  (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)           // re of HII
                  - coolingheating_buf.reHeII1[i-1]*HeII (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.      // re of HeII
                  - coolingheating_buf.reHeII2[i-1]*HeII (i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.      // re of HeII
                  - coolingheating_buf.reHeIII[i-1]*HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1)*de(i-1,idx_range.jp1-1,idx_range.kp1-1)/4.      // re of HeIII

        // Bremsstrahlung

                  - coolingheating_buf.brem[i-1]*(HII(i-1,idx_range.jp1-1,idx_range.kp1-1)+HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)/4. +
                HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1)) * de(i-1,idx_range.jp1-1,idx_range.kp1-1)

                  );
        Lpri[i-1] = edot[i-1];

        if (edot[i-1] != edot[i-1])  {
          OMP_PRAGMA_CRITICAL
          {
            eprintf("NaN in edot[1]:  %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
                    i,
                    idx_range.jp1,
                    idx_range.kp1,
                    edot [ i-1 ],
                    HI ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    HII ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    HeI ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    HeII ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    HeIII ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    de ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    d ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    tgas [ i-1 ],
                    p2d [ i-1 ]);
          }
        }
         
      }
    }

  }

  // --- H2 cooling ---

  if (my_chemistry->primordial_chemistry > 1)  {
    // Chiaki & Wise (2019) H2 cooling rate
    if (my_chemistry->h2_cooling_rate == 3)  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {

          lognhat = logH2I[i-1] - logdvdr[i-1];

           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH[i-1], my_rates->LH2.props.dimension,
            my_rates->LH2.props.parameters[0], &my_rates->LH2.props.parameter_spacing[0], my_rates->LH2.props.parameters[1], &my_rates->LH2.props.parameter_spacing[1], my_rates->LH2.props.parameters[2], &my_rates->LH2.props.parameter_spacing[2],
            &my_rates->LH2.props.data_size, my_rates->LH2.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH[i-1], my_rates->LH2.props.dimension,
              my_rates->LH2.props.parameters[0], &my_rates->LH2.props.parameter_spacing[0], my_rates->LH2.props.parameters[1], &my_rates->LH2.props.parameter_spacing[1], my_rates->LH2.props.parameters[2], &my_rates->LH2.props.parameter_spacing[2],
              &my_rates->LH2.props.data_size, my_rates->LH2.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }

          LH2[i-1] = ih2cox * (G - L) / dom * H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)/2.e0;
          if (LH2[i-1] != LH2[i-1]) { LH2[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LH2[i-1];

        }
      }

      // Glover & Abel (2008) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 2)  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          gaHI[i-1] = my_rates->GAHI[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GAHI[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GAHI[logTlininterp_buf.indixe[i-1]-1]);
          gaH2[i-1] = my_rates->GAH2[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GAH2[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GAH2[logTlininterp_buf.indixe[i-1]-1]);
          gaHe[i-1] = my_rates->GAHe[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GAHe[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GAHe[logTlininterp_buf.indixe[i-1]-1]);
          gaHp[i-1] = my_rates->GAHp[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GAHp[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GAHp[logTlininterp_buf.indixe[i-1]-1]);
          gael[i-1] = my_rates->GAel[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GAel[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GAel[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.gphdl[i-1] = my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]-1]);
          h2lte[i-1] = my_rates->H2LTE[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->H2LTE[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->H2LTE[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.cieco[i-1] = my_rates->cieco[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->cieco[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->cieco[logTlininterp_buf.indixe[i-1]-1]);
        }
      }

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1);
          nother = (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) +
               HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))/4. +
               HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1) + de(i-1,idx_range.jp1-1,idx_range.kp1-1);
          fH2 = nH2/(nH2 + nother);
  //_// PORT:             fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
  //_// PORT:      &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
  //_// PORT:      &           ((nH2 + nother)*dom) )
          fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */
          // Note that this optical depth approximation comes from
          // RA04.
          if (my_chemistry->h2_optical_depth_approximation==1)  {
            fudge = std::pow((0.76*d(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom/
                8.e9),(-0.45));
            fudge = std::fmin(fudge, 1.);
          } else {
            fudge = 1.;
          }
          galdl[i-1] = gaHI[i-1] * HI(i-1,idx_range.jp1-1,idx_range.kp1-1)
                   + gaH2[i-1] * H2I(i-1,idx_range.jp1-1,idx_range.kp1-1) / 2.
                   + gaHe[i-1] * HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) / 4.
                   + gaHp[i-1] * HII(i-1,idx_range.jp1-1,idx_range.kp1-1)
                   + gael[i-1] * de(i-1,idx_range.jp1-1,idx_range.kp1-1);
          // gphdl1 = gphdl(i)/dom
          gphdl1 = h2lte[i-1]/dom;
          edot[i-1] = edot[i-1] - ih2cox*fudge*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)*
               h2lte[i-1]/(1. + gphdl1/galdl[i-1]) / (2.*dom);

        }
      }

      // Galli & Palla (1998) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 1)  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          coolingheating_buf.gpldl[i-1] = my_rates->GP99LowDensityLimit[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GP99LowDensityLimit[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GP99LowDensityLimit[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.gphdl[i-1] = my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->GP99HighDensityLimit[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.cieco[i-1] = my_rates->cieco[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->cieco[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->cieco[logTlininterp_buf.indixe[i-1]-1]);
        }
      }

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {

#define NO_OPTICAL_DEPTH_FUDGE
#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1);
          nother = (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) +
               HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))/4. +
               HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1) + de(i-1,idx_range.jp1-1,idx_range.kp1-1);
          fH2 = nH2/(nH2 + nother);
  //_// PORT:             fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
  //_// PORT:      &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
  //_// PORT:      &           ((nH2 + nother)*dom) )
          fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */
          // Note that this optical depth approximation comes from
          // RA04.
          if (my_chemistry->h2_optical_depth_approximation==1)  {
            fudge = std::pow((0.76*d(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom/
                8.e9),(-0.45));
            fudge = std::fmin(fudge, 1.);
          } else {
            fudge = 1.;
          }
          gphdl1 = coolingheating_buf.gphdl[i-1]/(HI(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom);
          edot[i-1] = edot[i-1] - ih2cox*fudge*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)*
               coolingheating_buf.gphdl[i-1]/(1. + gphdl1/coolingheating_buf.gpldl[i-1]) / (2.*dom);

        }
      }

      // Lepp & Shull (1983) H2 cooling rate
    } else if (my_chemistry->h2_cooling_rate == 0)  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          coolingheating_buf.hyd01k[i-1] = my_rates->hyd01k[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->hyd01k[logTlininterp_buf.indixe[i-1]+1-1]-my_rates->hyd01k[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.h2k01[i-1] = my_rates->h2k01[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->h2k01[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->h2k01[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.vibh[i-1] = my_rates->vibh[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->vibh[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->vibh[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.roth[i-1] = my_rates->roth[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->roth[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->roth[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.rotl[i-1] = my_rates->rotl[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->rotl[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->rotl[logTlininterp_buf.indixe[i-1]-1]);
          coolingheating_buf.cieco[i-1] = my_rates->cieco[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->cieco[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->cieco[logTlininterp_buf.indixe[i-1]-1]);
        }
      }

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          qq   = 1.2*std::pow((HI(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom),0.77) +
                    std::pow((H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom/2.),0.77);
          vibl = (HI(i-1,idx_range.jp1-1,idx_range.kp1-1)*coolingheating_buf.hyd01k[i-1] +
                 H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)/2.*coolingheating_buf.h2k01[i-1])
                 *dom*8.18e-13;

#ifdef OPTICAL_DEPTH_FUDGE
          nH2 = 0.5*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1);
          nother = (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) +
               HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))/4. +
               HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1) + de(i-1,idx_range.jp1-1,idx_range.kp1-1);
          fH2 = nH2/(nH2 + nother);
  //_// PORT:             fudge = sqrt((40._DKIND * 10._DKIND**(4.8_DKIND * 
  //_// PORT:      &           sqrt(max(log10(tgas(i)),2._DKIND)-2._DKIND)) / fH2**2)/
  //_// PORT:      &           ((nH2 + nother)*dom) )
          fudge = std::fmin(fudge, 1.);
#endif /* OPTICAL_DEPTH_FUDGE */

          edot[i-1] = edot[i-1] - ih2cox*fudge*H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)*(
               coolingheating_buf.vibh[i-1]/(1.+coolingheating_buf.vibh[i-1]/std::fmax(   vibl,     tiny_fortran_val)) +
               coolingheating_buf.roth[i-1]/(1.+coolingheating_buf.roth[i-1]/std::fmax(qq*coolingheating_buf.rotl[i-1],tiny_fortran_val))
               )/2./dom;
        }
      }

    }

    // CIE
    // cooling from H2-H2 and He-H2 collisional induced emission comes
    //- with its own radiative transfer correction as discussed in
    //- Ripamonti & Abel 2003
    if (my_chemistry->cie_cooling==1)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if (itmask[i-1] != MASK_FALSE)  {
          // Only calculate if H2I(i) is a substantial fraction
          if (d(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom>1e10)  {
            ciefudge = 1.;
            tau = std::pow(((d(i-1,idx_range.jp1-1,idx_range.kp1-1)/2e16)*dom),2.8);  // 2e16 is in units of cm^-3
            tau = std::fmax(tau, 1.e-5);
            ciefudge = std::fmin((1.-std::exp(-tau))/tau,1.);
            // Matt's attempt at a second exponentialier cutoff
            tau = std::pow(((d(i-1,idx_range.jp1-1,idx_range.kp1-1)/2.e18)*dom),8.);  // 2e18 is in units of cm^-3
            tau = std::fmax(tau, 1.e-5);
            ciefudge = ciefudge*std::fmin((1.f-std::exp(-tau))/tau,1.);
            // ciefudge, which is applied to the continuum, is applied to edot
            edot[i-1] = ciefudge*(edot[i-1] -
                    H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)*(d(i-1,idx_range.jp1-1,idx_range.kp1-1)*coolingheating_buf.cieco[i-1]));
          }
        }
      }
      // CIE H2 cooling using Yoshida et al. (2006)
    } else if (my_chemistry->cie_cooling == 2)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if (itmask[i-1] != MASK_FALSE)  {
          cieY06[i-1] = my_rates->cieY06[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->cieY06[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->cieY06[logTlininterp_buf.indixe[i-1]-1]);
          LCIE[i-1] = - cieY06[i-1] * std::pow((H2I(i-1,idx_range.jp1-1,idx_range.kp1-1)/2.e0),2);
          edot[i-1] = edot[i-1] + LCIE[i-1];
        }
      }
    }

  }

  // --- Cooling from HD ---

  if (my_chemistry->primordial_chemistry > 2)  {

    // Chiaki & Wise (2019) HD cooling rate
    if (my_chemistry->hd_cooling_rate == 1 )  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {

          lognhat = logHDI[i-1] - logdvdr[i-1];

           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH[i-1], my_rates->LHD.props.dimension,
            my_rates->LHD.props.parameters[0], &my_rates->LHD.props.parameter_spacing[0], my_rates->LHD.props.parameters[1], &my_rates->LHD.props.parameter_spacing[1], my_rates->LHD.props.parameters[2], &my_rates->LHD.props.parameter_spacing[2],
            &my_rates->LHD.props.data_size, my_rates->LHD.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));

          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH[i-1], my_rates->LHD.props.dimension,
              my_rates->LHD.props.parameters[0], &my_rates->LHD.props.parameter_spacing[0], my_rates->LHD.props.parameters[1], &my_rates->LHD.props.parameter_spacing[1], my_rates->LHD.props.parameters[2], &my_rates->LHD.props.parameter_spacing[2],
              &my_rates->LHD.props.data_size, my_rates->LHD.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }

          LHD[i-1] = (G - L) / dom * HDI(i-1,idx_range.jp1-1,idx_range.kp1-1)/3.e0;
          if (LHD[i-1] != LHD[i-1]) { LHD[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LHD[i-1];

        }
      }

      // Coppola et al (2011) and Wrathmall, Gusdorf, & Flower (2007) HD cooling rate
    } else if (my_chemistry->hd_cooling_rate == 0)  {

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          // CMB cooling floor
          if (tgas[i-1] > comp2)  {
            coolingheating_buf.hdlte[i-1] = my_rates->HDlte[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
            *(my_rates->HDlte[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->HDlte[logTlininterp_buf.indixe[i-1]-1]);
            coolingheating_buf.hdlow[i-1] = my_rates->HDlow[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
            *(my_rates->HDlow[logTlininterp_buf.indixe[i-1]+1-1] - my_rates->HDlow[logTlininterp_buf.indixe[i-1]-1]);
          } else {
            coolingheating_buf.hdlte[i-1] = tiny_fortran_val;
            coolingheating_buf.hdlow[i-1] = tiny_fortran_val;
          }
        }
      }

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          // old (incorrect) way:
          //              hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2._DKIND)
          //              hdlow1 = max(hdlow(i), tiny)
          //              edot(i) = edot(i) - HDI(i,j,k)*
          //    .                     (hdlte1/(1._DKIND + hdlte1/hdlow1)/(2._DKIND*dom))
          // new (correct) way: (april 4, 2007)
          hdlte1 = coolingheating_buf.hdlte[i-1]/(HI(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom);
          hdlow1 = std::fmax(coolingheating_buf.hdlow[i-1], tiny_fortran_val);
          edot[i-1] = edot[i-1] - HDI(i-1,idx_range.jp1-1,idx_range.kp1-1)*
               (coolingheating_buf.hdlte[i-1]/(1. + hdlte1/hdlow1)) /
               (3.*dom);
        }
      }

    }
  }

  // Iteration mask for metal-rich cells
  if (imetal == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (metallicity[i-1] >= min_metallicity)  {
        itmask_metal[i-1] = itmask[i-1];
      } else {
        itmask_metal[i-1] = MASK_FALSE;
      }
    }
  } else {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      itmask_metal[i-1] = MASK_FALSE;
    }
  }

  // Compute grain size increment

  if ( (my_chemistry->use_dust_density_field > 0)  &&  (my_chemistry->dust_species > 0) )  {

     FORTRAN_NAME(calc_grain_size_increment_1d)(
              &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->grain_growth, itmask_metal,
             &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, d.data(),
             my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
             my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
             my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
             metal.data(), my_fields->local_ISM_metal_density,
             my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
             my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
             my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density,
             &my_rates->SN0_N,
             my_rates->SN0_fSiM, my_rates->SN0_fFeM, my_rates->SN0_fMg2SiO4, my_rates->SN0_fMgSiO3,
             my_rates->SN0_fFe3O4, my_rates->SN0_fAC, my_rates->SN0_fSiO2D, my_rates->SN0_fMgO,
             my_rates->SN0_fFeS, my_rates->SN0_fAl2O3,
             my_rates->SN0_freforg, my_rates->SN0_fvolorg, my_rates->SN0_fH2Oice,
             my_rates->SN0_r0SiM, my_rates->SN0_r0FeM, my_rates->SN0_r0Mg2SiO4, my_rates->SN0_r0MgSiO3,
             my_rates->SN0_r0Fe3O4, my_rates->SN0_r0AC, my_rates->SN0_r0SiO2D, my_rates->SN0_r0MgO,
             my_rates->SN0_r0FeS, my_rates->SN0_r0Al2O3,
             my_rates->SN0_r0reforg, my_rates->SN0_r0volorg, my_rates->SN0_r0H2Oice,
             my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td,
             my_rates->SN0_kpSiM, my_rates->SN0_kpFeM, my_rates->SN0_kpMg2SiO4, my_rates->SN0_kpMgSiO3,
             my_rates->SN0_kpFe3O4, my_rates->SN0_kpAC, my_rates->SN0_kpSiO2D, my_rates->SN0_kpMgO,
             my_rates->SN0_kpFeS, my_rates->SN0_kpAl2O3,
             my_rates->SN0_kpreforg, my_rates->SN0_kpvolorg, my_rates->SN0_kpH2Oice,
             sgSiM.data(), sgFeM.data(), sgMg2SiO4.data(), sgMgSiO3.data(), sgFe3O4.data(), sgAC.data(),
             sgSiO2D.data(), sgMgO.data(), sgFeS.data(), sgAl2O3.data(),
             sgreforg.data(), sgvolorg.data(), sgH2Oice.data(), sgtot.data(),
             alSiM.data(), alFeM.data(), alMg2SiO4.data(), alMgSiO3.data(), alFe3O4.data(), alAC.data(),
             alSiO2D.data(), alMgO.data(), alFeS.data(), alAl2O3.data(),
             alreforg.data(), alvolorg.data(), alH2Oice.data(), altot.data()
          );

  }

  // Calculate dust to gas ratio AND interstellar radiation field
  // -> an earlier version of this logic would store values @ indices
  //    where `itmask_metal(i) .ne. MASK_FALSE`
  // -> this was undesirable, b/c these quantities are required for
  //    photo-electric heating, which can occur when
  //    `itmask_metal(i) .eq. MASK_FALSE` (we can revisit this choice
  //    later). Moreover, in most cases, these calculations will be
  //    faster when there is no branching

  if ((anydust != MASK_FALSE)  ||  (my_chemistry->photoelectric_heating > 0))  {
    if (my_chemistry->use_dust_density_field > 0)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        // REMINDER: use of `itmask` over `itmask_metal` is
        //   currently required by Photo-electric heating
        if (itmask[i-1] != MASK_FALSE)  {
          // it may be faster to remove this branching
          dust2gas[i-1] = dust(i-1,idx_range.jp1-1,idx_range.kp1-1) / d(i-1,idx_range.jp1-1,idx_range.kp1-1);
        }
      }
    } else {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        dust2gas[i-1] = my_chemistry->local_dust_to_gas_ratio * metallicity[i-1];
      }
    }
  }

  if ((anydust != MASK_FALSE)  ||  (my_chemistry->photoelectric_heating > 1))  {
    if (my_chemistry->use_isrf_field > 0)  {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        myisrf[i-1] = isrf_habing(i-1,idx_range.jp1-1,idx_range.kp1-1);
      }
    } else {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        myisrf[i-1] = my_chemistry->interstellar_radiation_field;
      }
    }
  }

  // compute dust temperature and cooling due to dust
  if (anydust != MASK_FALSE)  {

     FORTRAN_NAME(calc_all_tdust_gasgr_1d_g)(&my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                 &my_chemistry->use_dust_density_field, &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &my_chemistry->local_dust_to_gas_ratio, &my_rates->gamma_isrf,
                 &comp2, my_rates->gas_grain, logTlininterp_buf.indixe, logTlininterp_buf.tdef, tgas, tdust,
                 metallicity, dust2gas, cool1dmulti_buf.mynh, cool1dmulti_buf.gasgr_tdust,
                 itmask_metal,
                 &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT,
                 my_rates->gr_Td, grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
                 grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust],
                 grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], my_rates->gas_grain2, &my_rates->gamma_isrf2,
                 &coolunit, gasgr.data(), myisrf.data(), sgSiM.data(), sgFeM.data(), sgMg2SiO4.data(),
                 sgMgSiO3.data(), sgFe3O4.data(), sgAC.data(), sgSiO2D.data(), sgMgO.data(), sgFeS.data(),
                 sgAl2O3.data(), sgreforg.data(), sgvolorg.data(), sgH2Oice.data(), sgtot.data(),
                 alSiM.data(), alFeM.data(), alMg2SiO4.data(), alMgSiO3.data(), alFe3O4.data(), alAC.data(),
                 alSiO2D.data(), alMgO.data(), alFeS.data(), alAl2O3.data(), alreforg.data(),
                 alvolorg.data(), alH2Oice.data(), altot.data(), kpSiM.data(), kpFeM.data(),
                 kpMg2SiO4.data(), kpMgSiO3.data(), kpFe3O4.data(), kpAC.data(), kpSiO2D.data(),
                 kpMgO.data(), kpFeS.data(), kpAl2O3.data(), kpreforg.data(), kpvolorg.data(),
                 kpH2Oice.data(), kptot.data(), gasSiM.data(), gasFeM.data(), gasMg2SiO4.data(),
                 gasMgSiO3.data(), gasFe3O4.data(), gasAC.data(), gasSiO2D.data(), gasMgO.data(),
                 gasFeS.data(), gasAl2O3.data(), gasreforg.data(), gasvolorg.data(),
                 gasH2Oice.data()
               );

  }

  // Calculate dust cooling rate
  if (anydust != MASK_FALSE)  {


    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask_metal[i-1] != MASK_FALSE )  {

        if (my_chemistry->dust_species == 0)  {

          Ldst[i-1] =         -
               gasgr[i-1] * (tgas[i-1] - tdust[i-1]) *
               dust2gas[i-1] * rhoH[i-1] * rhoH[i-1];

        } else {

          if (my_chemistry->use_multiple_dust_temperatures == 0)  {
            Ldst[i-1] = - gasgr[i-1] * (tgas[i-1] - tdust[i-1])
               * d(i-1,idx_range.jp1-1,idx_range.kp1-1) * rhoH[i-1];
          } else {

            if (my_chemistry->dust_species > 0)  {
              Ldst[i-1] = - (
                  gasMgSiO3  [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust]  [i-1])
                + gasAC      [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::AC_dust]      [i-1])
                ) * d(i-1,idx_range.jp1-1,idx_range.kp1-1) * rhoH[i-1];
            }

            if (my_chemistry->dust_species > 1)  {
              Ldst[i-1] = Ldst[i-1] - (
                  gasSiM     [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::SiM_dust]     [i-1])
                + gasFeM     [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::FeM_dust]     [i-1])
                + gasMg2SiO4 [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust] [i-1])
                + gasFe3O4   [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust]   [i-1])
                + gasSiO2D   [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust]   [i-1])
                + gasMgO     [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::MgO_dust]     [i-1])
                + gasFeS     [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::FeS_dust]     [i-1])
                + gasAl2O3   [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust]   [i-1])
                ) * d(i-1,idx_range.jp1-1,idx_range.kp1-1) * rhoH[i-1];
            }

            if (my_chemistry->dust_species > 2)  {
              Ldst[i-1] = Ldst[i-1] - (
                  gasreforg  [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust]  [i-1])
                + gasvolorg  [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust]  [i-1])
                + gasH2Oice  [i-1] * (tgas[i-1] - grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust]  [i-1])
                ) * d(i-1,idx_range.jp1-1,idx_range.kp1-1) * rhoH[i-1];
            }
          }


        }

        edot[i-1] = edot[i-1] + Ldst[i-1];

      }
    }

  }

  // Compute continuum opacity

  if ( my_chemistry->use_primordial_continuum_opacity == 1 )  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {

        // ! primordial continuum opacity !!
         FORTRAN_NAME(interpolate_2d_g)(
          &logrho[i-1], &logT[i-1], my_rates->alphap.props.dimension, my_rates->alphap.props.parameters[0], &my_rates->alphap.props.parameter_spacing[0],
          my_rates->alphap.props.parameters[1], &my_rates->alphap.props.parameter_spacing[1], &my_rates->alphap.props.data_size,
          my_rates->alphap.data, &log_a);
        alpha[i-1] = std::pow(1.e1,log_a);
      }
    }

  } else {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        alpha[i-1] = 0.f;
      }
    }

  }

  // Add dust opacity.
  // if (idspecies .eq. 0), dust opacity is overestimated at Td > 50 K
  // We better not include dust opacity.
  if ((anydust != MASK_FALSE) && (my_chemistry->dust_species > 0))  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask_metal[i-1] != MASK_FALSE )  {

        if (my_chemistry->use_multiple_dust_temperatures == 0)  {

          alphad[i-1] = kptot[i-1];

        } else {

          if (my_chemistry->dust_species > 0)  {
            alphad[i-1] = kpMgSiO3  [i-1]
                      + kpAC      [i-1];
          }
          if (my_chemistry->dust_species > 1)  {
            alphad[i-1] = alphad[i-1]
                      + kpSiM     [i-1]
                      + kpFeM     [i-1]
                      + kpMg2SiO4 [i-1]
                      + kpFe3O4   [i-1]
                      + kpSiO2D   [i-1]
                      + kpMgO     [i-1]
                      + kpFeS     [i-1]
                      + kpAl2O3   [i-1];
          }
          if (my_chemistry->dust_species > 2)  {
            alphad[i-1] = alphad[i-1]
                      + kpreforg  [i-1]
                      + kpvolorg  [i-1]
                      + kpH2Oice  [i-1];
          }
        }

        alpha[i-1] = alpha[i-1] + alphad[i-1] * d(i-1,idx_range.jp1-1,idx_range.kp1-1)*dom*mh_local_var;

      }
    }

  }

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      tau_con[i-1] = alpha[i-1] * lshield_con[i-1];
    }
  }

  // --- Compute (external) radiative heating terms ---
  // Photoionization heating

  if (my_chemistry->primordial_chemistry > 0)  {

    if (my_chemistry->self_shielding_method == 0)  // no shielding
    {
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          edot[i-1] = edot[i-1] + (double)(my_chemistry->ipiht )*(
                 my_uvb_rates.piHI  *HI  (i-1,idx_range.jp1-1,idx_range.kp1-1) // pi of HI
               + my_uvb_rates.piHeI *HeI (i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25 // pi of HeI
               + my_uvb_rates.piHeII*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25 // pi of HeII
               )/dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 1)  {
      
      // approximate self shielding using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // to shield HI, while leaving HeI and HeII optically thin

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if (itmask[i-1] != MASK_FALSE)  {
          if (my_uvb_rates.k24 < tiny8)  {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
             std::pow((my_uvb_rates.crsHI /2.49e-18),(-2./3.))*
             std::pow((tgas[i-1]/1.0e4),(0.17))*
             std::pow((my_uvb_rates.k24/tbase1/1.0e-12),(2./3.));
            nratio = (HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1))*dom/nSSh;
            fSShHI =
             0.98*std::pow((1.+
             std::pow(nratio,(1.64))),(-2.28)) +
             0.02*std::pow((1.+
             nratio),(-0.84));
          }

          edot[i-1] = edot[i-1] + (double)(my_chemistry->ipiht)*(
                 my_uvb_rates.piHI  *HI  (i-1,idx_range.jp1-1,idx_range.kp1-1)* fSShHI
               + my_uvb_rates.piHeI * HeI(i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25
               + my_uvb_rates.piHeII*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25
                )/dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 2)   {
      
      // Better self-shielding in HI using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // approximate self shielding in HeI and HeII

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          
          // HI self shielding ratio
          if (my_uvb_rates.k24 < tiny8)  {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
             std::pow((my_uvb_rates.crsHI/2.49e-18),(-2./3.))*
             std::pow((tgas[i-1]/1.0e4),(0.17))*
             std::pow((my_uvb_rates.k24/tbase1/1.0e-12),(2./3.));
            nratio = (HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1))*dom/nSSh;
            fSShHI =
             0.98*std::pow((1.+
              std::pow(nratio,(1.64))),(-2.28))+
             0.02*std::pow((1.+
              nratio),(-0.84));
          }
          
          // HeI self shielding ratio
          if (my_uvb_rates.k26 < tiny8)  {
            fSShHeI = 1.;
          } else {
            nssh_he = 6.73e-3 *
             std::pow((my_uvb_rates.crsHeI/ 2.49e-18),(-2./3.))*
             std::pow((tgas[i-1]/1.0e4),(0.17))*
             std::pow((my_uvb_rates.k26/tbase1/1.0e-12),(2./3.));
            nratio_he = 0.25*
             (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))*dom/nssh_he;
            fSShHeI =
             0.98*std::pow((1.+
              std::pow(nratio_he,(1.64))),(-2.28))+
             0.02*std::pow((1.+
              nratio_he),(-0.84));
          }

          edot[i-1] = edot[i-1] + (double)(my_chemistry->ipiht )*(
                 my_uvb_rates.piHI * HI(i-1,idx_range.jp1-1,idx_range.kp1-1)* fSShHI
               + my_uvb_rates.piHeI * HeI(i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25* fSShHeI
               + my_uvb_rates.piHeII*HeII(i-1,idx_range.jp1-1,idx_range.kp1-1)*0.25
                 )/dom;
        }
      }

    } else if (my_chemistry->self_shielding_method == 3)  {
      
      // shielding using Eq. 13 and 14 from
      // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
      // in HI and HeI, but ignoring HeII heating entirely

      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask[i-1] != MASK_FALSE )  {
          
          // HI self shielding ratio
          if (my_uvb_rates.k24 < tiny8)  {
            fSShHI = 1.;
          } else {
            nSSh = 6.73e-3 *
             std::pow((my_uvb_rates.crsHI /2.49e-18),(-2./3.))*
             std::pow((tgas[i-1]/1.0e4),(0.17))*
             std::pow((my_uvb_rates.k24/tbase1/1.0e-12),(2./3.));
            nratio = (HI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HII(i-1,idx_range.jp1-1,idx_range.kp1-1))*dom/nSSh;
            fSShHI =
             0.98*std::pow((1.+
              std::pow(nratio,(1.64))),(-2.28))+
             0.02*std::pow((1.+
              nratio),(-0.84));
          }
          
          // HeI self shielding ratio
          if (my_uvb_rates.k26 < tiny8)  {
            fSShHeI = 1.;
          } else {
            nssh_he = 6.73e-3 *
             std::pow((my_uvb_rates.crsHeI /2.49e-18),(-2./3.))*
             std::pow((tgas[i-1]/1.0e4),(0.17))*
             std::pow((my_uvb_rates.k26/tbase1/1.0e-12),(2./3.));
            nratio_he = 0.25*
             (HeI(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeII(i-1,idx_range.jp1-1,idx_range.kp1-1) + HeIII(i-1,idx_range.jp1-1,idx_range.kp1-1))*dom/nssh_he;
            fSShHeI =
             0.98*std::pow((1.+
              std::pow(nratio_he,(1.64))),(-2.28))+
             0.02*std::pow((1.+
              nratio_he),(-0.84));
          }

          edot[i-1] = edot[i-1] + (double)(my_chemistry->ipiht )*(
                 my_uvb_rates.piHI * HI (i-1,idx_range.jp1-1,idx_range.kp1-1)* fSShHI
              + my_uvb_rates.piHeI * HeI(i-1,idx_range.jp1-1,idx_range.kp1-1)* fSShHeI
              )/dom;
          
          // Ignoring HeII heating (HeII heating rate -> 0)
        }
      }

    }

  }

  // --- Cloudy primordial cooling and heating ---

  if (my_chemistry->primordial_chemistry == 0)  {

    iZscale = 0;
    mycmbTfloor = 0;
     FORTRAN_NAME(cool1d_cloudy_g)(d.data(), rhoH, metallicity,
         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
         logTlininterp_buf.logtem, edot, &comp2, &dom, &zr,
         &mycmbTfloor, &my_chemistry->UVbackground, &iZscale,
         &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
         my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2],
         &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data,
         itmask);

    // Calculate electron density from mean molecular weight

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {

        cool1dmulti_buf.myde[i-1] = 1 - mmw[i-1] * (3.0 * my_chemistry->HydrogenFractionByMass + 1.0) /
             4.0;
        if (imetal == 1)  {
          cool1dmulti_buf.myde[i-1] = cool1dmulti_buf.myde[i-1] - mmw[i-1] * metal(i-1,idx_range.jp1-1,idx_range.kp1-1) /
               (d(i-1,idx_range.jp1-1,idx_range.kp1-1) * mu_metal);
        }
        cool1dmulti_buf.myde[i-1] = d(i-1,idx_range.jp1-1,idx_range.kp1-1) * cool1dmulti_buf.myde[i-1] / mmw[i-1];
        cool1dmulti_buf.myde[i-1] = std::fmax(cool1dmulti_buf.myde[i-1], 0.);

      }
    }

  }

  // Photo-electric heating by UV-irradiated dust

  if (my_chemistry->photoelectric_heating == 1)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        if ( tgas[i-1] > 2.e4 )  {
          cool1dmulti_buf.gammaha_eff[i-1] = 0.;
        } else {
          cool1dmulti_buf.gammaha_eff[i-1] = my_rates->gammah;
        }
      }
    }

    // Use eqn. 1 of Wolfire et al. (1995)
  } else if (my_chemistry->photoelectric_heating == 2)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        if ( tgas[i-1] > 2.e4 )  {
          cool1dmulti_buf.gammaha_eff[i-1] = 0.;
        } else {
          // Assume constant epsilon = 0.05.
          cool1dmulti_buf.gammaha_eff[i-1] = my_rates->gammah * 0.05 * myisrf[i-1];
        }
      }
    }

    // Full calculation of epsilon (eqn. 2 of Wolfire 1995)
  } else if (my_chemistry->photoelectric_heating == 3)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        pe_X = myisrf[i-1] * dom_inv * std::sqrt(tgas[i-1]) / cool1dmulti_buf.myde[i-1];
        pe_eps =
             (4.9e-2 /
              (1. + std::pow((pe_X / 1925.),0.73))) +
             ((3.7e-2 * std::pow((tgas[i-1] / 1.e4),0.7)) /
              (1. + (pe_X / 5000.)));
        cool1dmulti_buf.gammaha_eff[i-1] = my_rates->gammah * pe_eps * myisrf[i-1];
      }
    }

  }

  if (my_chemistry->photoelectric_heating > 0)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        edot[i-1] = edot[i-1] + cool1dmulti_buf.gammaha_eff[i-1] * rhoH[i-1] *
             dom_inv * dust2gas[i-1] / my_chemistry->local_dust_to_gas_ratio;
      }
    }
  }

  // Electron recombination onto dust grains (eqn. 9 of Wolfire 1995)

  if ((my_chemistry->dust_chemistry > 0)  ||  (my_chemistry->dust_recombination_cooling > 0))  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        cool1dmulti_buf.regr[i-1] = my_rates->regr[logTlininterp_buf.indixe[i-1]-1] + logTlininterp_buf.tdef[i-1]
             *(my_rates->regr[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->regr[logTlininterp_buf.indixe[i-1]-1]);
      }
    }

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        grbeta = 0.74 / std::pow(tgas[i-1],0.068);
        edot[i-1] = edot[i-1] -
             cool1dmulti_buf.regr[i-1] * std::pow((myisrf[i-1]*dom_inv / cool1dmulti_buf.myde[i-1]),grbeta) *
             cool1dmulti_buf.myde[i-1] * rhoH[i-1] * dust2gas[i-1] / my_chemistry->local_dust_to_gas_ratio;
      }
    }

  }

  // Compton cooling or heating and X-ray compton heating

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if (itmask[i-1] != MASK_FALSE)  {

      edot[i-1] = edot[i-1]

      // Compton cooling or heating

           - comp1      * (tgas[i-1] - comp2)     * cool1dmulti_buf.myde[i-1]*dom_inv

      // X-ray compton heating

           - my_uvb_rates.comp_xray * (tgas[i-1] - my_uvb_rates.temp_xray) * cool1dmulti_buf.myde[i-1]*dom_inv;

    }
  }
 
  // Photoheating from radiative transfer

  if (my_chemistry->use_radiative_transfer == 1)  {
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if (itmask[i-1] != MASK_FALSE)  {
        edot[i-1] = edot[i-1] + (double)(my_chemistry->ipiht ) * photogamma(i-1,idx_range.jp1-1,idx_range.kp1-1)
                          / coolunit * HI(i-1,idx_range.jp1-1,idx_range.kp1-1) / dom;

        if (edot[i-1] != edot[i-1])  {
          OMP_PRAGMA_CRITICAL
          {
            eprintf("NaN in edot[2]:  %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
                    i,
                    idx_range.jp1,
                    idx_range.kp1,
                    edot [ i-1 ],
                    photogamma ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    HI ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    de ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    d ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    e ( i-1, idx_range.jp1-1, idx_range.kp1-1 ),
                    p2d [ i-1 ],
                    tgas [ i-1 ],
                    dom,
                    internalu.urho,
                    internalu.a_value,
                    mh_local_var);
          }
        }

      }
    }
  }

  // --- Cloudy metal cooling and heating ---

  if (my_chemistry->metal_cooling == 1)  {

    // Determine if the temperature is above the threshold to do tabulated cooling.
    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      itmask_tab[i-1] = itmask_metal[i-1];
      if ( itmask_tab[i-1] != MASK_FALSE )  {
        if (( my_chemistry->tabulated_cooling_minimum_temperature > 0.0e0 )  && 
            ( tgas[i-1] < my_chemistry->tabulated_cooling_minimum_temperature ))  {
          itmask_tab[i-1] = MASK_FALSE;
        }
      }
    }

    if (my_rates->cloudy_data_new == 1)  {

      iZscale = 1;
       FORTRAN_NAME(cool1d_cloudy_g)(d.data(), rhoH, metallicity,
           &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
           logTlininterp_buf.logtem, edot, &comp2, &dom, &zr,
           &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &iZscale,
           &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
           my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2],
           &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data,
           itmask_tab.data());

    } else {

       FORTRAN_NAME(cool1d_cloudy_old_tables_g)(
           d.data(), de.data(), rhoH, metallicity,
           &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
           logTlininterp_buf.logtem, edot, &comp2, &my_chemistry->primordial_chemistry, &dom, &zr,
           &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground,
           &my_chemistry->cloudy_electron_fraction_factor, &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
           my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
           &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data,
           itmask_tab.data());

    }

    if (my_chemistry->metal_chemistry == 1)  {

      // --- C/O fine-structure, metal molecular rotational cooling for low temperatures ---
      
      // C/O fine-structure cooling
      for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if ( itmask_metal[i-1] != MASK_FALSE )  {
      
          // CI
          lognhat = logCI[i-1] - logdvdr[i-1];

           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH[i-1], my_rates->LCI.props.dimension,
            my_rates->LCI.props.parameters[0], &my_rates->LCI.props.parameter_spacing[0], my_rates->LCI.props.parameters[1], &my_rates->LCI.props.parameter_spacing[1], my_rates->LCI.props.parameters[2], &my_rates->LCI.props.parameter_spacing[2],
            &my_rates->LCI.props.data_size, my_rates->LCI.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH[i-1], my_rates->LCI.props.dimension,
              my_rates->LCI.props.parameters[0], &my_rates->LCI.props.parameter_spacing[0], my_rates->LCI.props.parameters[1], &my_rates->LCI.props.parameter_spacing[1], my_rates->LCI.props.parameters[2], &my_rates->LCI.props.parameter_spacing[2],
              &my_rates->LCI.props.data_size, my_rates->LCI.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LCI[i-1] = (G - L) / dom * CI(i-1,idx_range.jp1-1,idx_range.kp1-1)/12.e0;
          if (LCI[i-1] != LCI[i-1]) { LCI[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LCI[i-1];
      
      
          // CII
          lognhat = logCII[i-1] - logdvdr[i-1];
      
           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH[i-1], my_rates->LCII.props.dimension,
            my_rates->LCII.props.parameters[0], &my_rates->LCII.props.parameter_spacing[0], my_rates->LCII.props.parameters[1], &my_rates->LCII.props.parameter_spacing[1], my_rates->LCII.props.parameters[2], &my_rates->LCII.props.parameter_spacing[2],
            &my_rates->LCII.props.data_size, my_rates->LCII.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH[i-1], my_rates->LCII.props.dimension,
              my_rates->LCII.props.parameters[0], &my_rates->LCII.props.parameter_spacing[0], my_rates->LCII.props.parameters[1], &my_rates->LCII.props.parameter_spacing[1], my_rates->LCII.props.parameters[2], &my_rates->LCII.props.parameter_spacing[2],
              &my_rates->LCII.props.data_size, my_rates->LCII.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LCII[i-1] = (G - L) / dom * CII(i-1,idx_range.jp1-1,idx_range.kp1-1)/12.e0;
          if (LCII[i-1] != LCII[i-1]) { LCII[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LCII[i-1];
      
      
          // OI
          lognhat = logOI[i-1] - logdvdr[i-1];
      
           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH[i-1], my_rates->LOI.props.dimension,
            my_rates->LOI.props.parameters[0], &my_rates->LOI.props.parameter_spacing[0], my_rates->LOI.props.parameters[1], &my_rates->LOI.props.parameter_spacing[1], my_rates->LOI.props.parameters[2], &my_rates->LOI.props.parameter_spacing[2],
            &my_rates->LOI.props.data_size, my_rates->LOI.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH[i-1], my_rates->LOI.props.dimension,
              my_rates->LOI.props.parameters[0], &my_rates->LOI.props.parameter_spacing[0], my_rates->LOI.props.parameters[1], &my_rates->LOI.props.parameter_spacing[1], my_rates->LOI.props.parameters[2], &my_rates->LOI.props.parameter_spacing[2],
              &my_rates->LOI.props.data_size, my_rates->LOI.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LOI[i-1] = (G - L) / dom * OI(i-1,idx_range.jp1-1,idx_range.kp1-1)/16.e0;
          if (LOI[i-1] != LOI[i-1]) { LOI[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LOI[i-1];
      
      
          // Metal molecules rotational cooling

          // CO
          lognhat = logCO[i-1] - logdvdr[i-1];
      
           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH2[i-1], my_rates->LCO.props.dimension,
            my_rates->LCO.props.parameters[0], &my_rates->LCO.props.parameter_spacing[0], my_rates->LCO.props.parameters[1], &my_rates->LCO.props.parameter_spacing[1], my_rates->LCO.props.parameters[2], &my_rates->LCO.props.parameter_spacing[2],
            &my_rates->LCO.props.data_size, my_rates->LCO.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH2[i-1], my_rates->LCO.props.dimension,
              my_rates->LCO.props.parameters[0], &my_rates->LCO.props.parameter_spacing[0], my_rates->LCO.props.parameters[1], &my_rates->LCO.props.parameter_spacing[1], my_rates->LCO.props.parameters[2], &my_rates->LCO.props.parameter_spacing[2],
              &my_rates->LCO.props.data_size, my_rates->LCO.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LCO[i-1] = (G - L) / dom * CO(i-1,idx_range.jp1-1,idx_range.kp1-1)/28.e0;
          if (LCO[i-1] != LCO[i-1]) { LCO[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LCO[i-1];
      

          // OH
          lognhat = logOH[i-1] - logdvdr[i-1];
      
           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH2[i-1], my_rates->LOH.props.dimension,
            my_rates->LOH.props.parameters[0], &my_rates->LOH.props.parameter_spacing[0], my_rates->LOH.props.parameters[1], &my_rates->LOH.props.parameter_spacing[1], my_rates->LOH.props.parameters[2], &my_rates->LOH.props.parameter_spacing[2],
            &my_rates->LOH.props.data_size, my_rates->LOH.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH2[i-1], my_rates->LOH.props.dimension,
              my_rates->LOH.props.parameters[0], &my_rates->LOH.props.parameter_spacing[0], my_rates->LOH.props.parameters[1], &my_rates->LOH.props.parameter_spacing[1], my_rates->LOH.props.parameters[2], &my_rates->LOH.props.parameter_spacing[2],
              &my_rates->LOH.props.data_size, my_rates->LOH.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LOH[i-1] = (G - L) / dom * OH(i-1,idx_range.jp1-1,idx_range.kp1-1)/17.e0;
          if (LOH[i-1] != LOH[i-1]) { LOH[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LOH[i-1];


          // H2O
          lognhat = logH2O[i-1] - logdvdr[i-1];
      
           FORTRAN_NAME(interpolate_3d_g)(
            &lognhat, &logT[i-1], &logH2[i-1], my_rates->LH2O.props.dimension,
            my_rates->LH2O.props.parameters[0], &my_rates->LH2O.props.parameter_spacing[0], my_rates->LH2O.props.parameters[1], &my_rates->LH2O.props.parameter_spacing[1], my_rates->LH2O.props.parameters[2], &my_rates->LH2O.props.parameter_spacing[2],
            &my_rates->LH2O.props.data_size, my_rates->LH2O.data, &log_Linv);
          L = std::pow(1.e1,(-log_Linv));
      
          if (my_chemistry->cmb_temperature_floor == 1)  {
             FORTRAN_NAME(interpolate_3d_g)(
              &lognhat, &logTcmb[i-1], &logH2[i-1], my_rates->LH2O.props.dimension,
              my_rates->LH2O.props.parameters[0], &my_rates->LH2O.props.parameter_spacing[0], my_rates->LH2O.props.parameters[1], &my_rates->LH2O.props.parameter_spacing[1], my_rates->LH2O.props.parameters[2], &my_rates->LH2O.props.parameter_spacing[2],
              &my_rates->LH2O.props.data_size, my_rates->LH2O.data, &log_Ginv);
            G = std::pow(1.e1,(-log_Ginv));
          } else {
            G = tiny8;
          }
      
          LH2O[i-1] = (G - L) / dom * H2O(i-1,idx_range.jp1-1,idx_range.kp1-1)/18.e0;
          if (LH2O[i-1] != LH2O[i-1]) { LH2O[i-1] = 0.e0; }
          edot[i-1] = edot[i-1] + LH2O[i-1];
      
        }
      }

    }

  }

  // Add user-provided volumetric and/or specific heating terms

  if (my_chemistry->use_volumetric_heating_rate == 1)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        edot[i-1] = edot[i-1] + Vheat(i-1,idx_range.jp1-1,idx_range.kp1-1) / coolunit / std::pow(dom,2);
      }
    }

  }

  if (my_chemistry->use_specific_heating_rate == 1)  {

    for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
      if ( itmask[i-1] != MASK_FALSE )  {
        edot[i-1] = edot[i-1] + Mheat(i-1,idx_range.jp1-1,idx_range.kp1-1) * d(i-1,idx_range.jp1-1,idx_range.kp1-1) * mh_local_var
            / coolunit / dom;
      }
    }

  }

  // Continuum opacity

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      if ( tau_con[i-1] > 1.e0 )  {
        if ( tau_con[i-1] < 1.e2 )  {
          edot[i-1] = edot[i-1] * std::pow(tau_con[i-1],(-2.e0));
        } else {
          edot[i-1] = 0.e0;
        }
      }
    }
  }

  // Set tgasold

  for (i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
    if ( itmask[i-1] != MASK_FALSE )  {
      cool1dmulti_buf.tgasold[i-1] = tgas[i-1];
    }
  }

  return;
}
