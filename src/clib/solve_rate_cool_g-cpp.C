// See LICENSE file for license and copyright information

/// @file solve_rate_cool_g-cpp.C
/// @brief Declares signature of solve_rate_cool_g

// This file was initially generated automatically during conversion of the
// solve_rate_cool_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

#include "solve_rate_cool_g-cpp.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int solve_rate_cool_g(
  int imetal, double dt, double utem, double uxyz, double urho,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  code_units* my_units, grackle_field_data* my_fields,
  photo_rate_storage* my_uvb_rates
)
{

  // Density, energy and velocity fields fields

  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Radiative transfer fields

  grackle::impl::View<gr_float***> kphHI(my_fields->RT_HI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Constants

#ifdef GRACKLE_FLOAT_4
  const gr_float tolerance = (gr_float)(1.0e-05);
#else
  const gr_float tolerance = (gr_float)(1.0e-10);
#endif

  const double mh_local_var = mh_grflt;
  const double pi_local_var = pi_fortran_val;

  // Locals

  int i, j, k, iter;
  int t, dj, dk;
  double ttmin, dom, energy, comp1, comp2;
  double coolunit, dbase1, tbase1, xbase1, chunit, uvel;
  double heq1, heq2, eqk221, eqk222, eqk131, eqk132, eqt1, eqt2, eqtdef, dheq, heq, dlogtem, dx_cgs, c_ljeans, min_metallicity;
  gr_float factor;

  // row temporaries

  double olddtit;
  std::vector<double> dtit(my_fields->grid_dimension[0]);
  std::vector<double> ttot(my_fields->grid_dimension[0]);
  std::vector<double> p2d(my_fields->grid_dimension[0]);
  std::vector<double> tgas(my_fields->grid_dimension[0]);
  std::vector<double> tdust(my_fields->grid_dimension[0]);
  std::vector<double> metallicity(my_fields->grid_dimension[0]);
  std::vector<double> dust2gas(my_fields->grid_dimension[0]);
  std::vector<double> rhoH(my_fields->grid_dimension[0]);
  std::vector<double> mmw(my_fields->grid_dimension[0]);
  std::vector<double> ddom(my_fields->grid_dimension[0]);

  // Rate equation row temporaries

  std::vector<double> HIp(my_fields->grid_dimension[0]);
  std::vector<double> HIIp(my_fields->grid_dimension[0]);
  std::vector<double> HeIp(my_fields->grid_dimension[0]);
  std::vector<double> HeIIp(my_fields->grid_dimension[0]);
  std::vector<double> HeIIIp(my_fields->grid_dimension[0]);
  std::vector<double> HMp(my_fields->grid_dimension[0]);
  std::vector<double> H2Ip(my_fields->grid_dimension[0]);
  std::vector<double> H2IIp(my_fields->grid_dimension[0]);
  std::vector<double> dep(my_fields->grid_dimension[0]);
  std::vector<double> dedot(my_fields->grid_dimension[0]);
  std::vector<double> HIdot(my_fields->grid_dimension[0]);
  std::vector<double> dedot_prev(my_fields->grid_dimension[0]);
  std::vector<double> DIp(my_fields->grid_dimension[0]);
  std::vector<double> DIIp(my_fields->grid_dimension[0]);
  std::vector<double> HDIp(my_fields->grid_dimension[0]);
  std::vector<double> HIdot_prev(my_fields->grid_dimension[0]);
  std::vector<double> k24shield(my_fields->grid_dimension[0]);
  std::vector<double> k25shield(my_fields->grid_dimension[0]);
  std::vector<double> k26shield(my_fields->grid_dimension[0]);
  std::vector<double> k28shield(my_fields->grid_dimension[0]);
  std::vector<double> k29shield(my_fields->grid_dimension[0]);
  std::vector<double> k30shield(my_fields->grid_dimension[0]);
  std::vector<double> k31shield(my_fields->grid_dimension[0]);
  std::vector<double> k1(my_fields->grid_dimension[0]);
  std::vector<double> k2(my_fields->grid_dimension[0]);
  std::vector<double> k3(my_fields->grid_dimension[0]);
  std::vector<double> k4(my_fields->grid_dimension[0]);
  std::vector<double> k5(my_fields->grid_dimension[0]);
  std::vector<double> k6(my_fields->grid_dimension[0]);
  std::vector<double> k7(my_fields->grid_dimension[0]);
  std::vector<double> k8(my_fields->grid_dimension[0]);
  std::vector<double> k9(my_fields->grid_dimension[0]);
  std::vector<double> k10(my_fields->grid_dimension[0]);
  std::vector<double> k11(my_fields->grid_dimension[0]);
  std::vector<double> k12(my_fields->grid_dimension[0]);
  std::vector<double> k13(my_fields->grid_dimension[0]);
  std::vector<double> k14(my_fields->grid_dimension[0]);
  std::vector<double> k15(my_fields->grid_dimension[0]);
  std::vector<double> k16(my_fields->grid_dimension[0]);
  std::vector<double> k17(my_fields->grid_dimension[0]);
  std::vector<double> k18(my_fields->grid_dimension[0]);
  std::vector<double> k19(my_fields->grid_dimension[0]);
  std::vector<double> k22(my_fields->grid_dimension[0]);
  std::vector<double> k50(my_fields->grid_dimension[0]);
  std::vector<double> k51(my_fields->grid_dimension[0]);
  std::vector<double> k52(my_fields->grid_dimension[0]);
  std::vector<double> k53(my_fields->grid_dimension[0]);
  std::vector<double> k54(my_fields->grid_dimension[0]);
  std::vector<double> k55(my_fields->grid_dimension[0]);
  std::vector<double> k56(my_fields->grid_dimension[0]);
  std::vector<double> k57(my_fields->grid_dimension[0]);
  std::vector<double> k58(my_fields->grid_dimension[0]);
  std::vector<double> k13dd(my_fields->grid_dimension[0] * 14);
  std::vector<double> h2dust(my_fields->grid_dimension[0]);
  std::vector<double> ncrn(my_fields->grid_dimension[0]);
  std::vector<double> ncrd1(my_fields->grid_dimension[0]);
  std::vector<double> ncrd2(my_fields->grid_dimension[0]);
  std::vector<double> DMp(my_fields->grid_dimension[0]);
  std::vector<double> HDIIp(my_fields->grid_dimension[0]);
  std::vector<double> HeHIIp(my_fields->grid_dimension[0]);
  std::vector<double> CIp(my_fields->grid_dimension[0]);
  std::vector<double> CIIp(my_fields->grid_dimension[0]);
  std::vector<double> COp(my_fields->grid_dimension[0]);
  std::vector<double> CO2p(my_fields->grid_dimension[0]);
  std::vector<double> OIp(my_fields->grid_dimension[0]);
  std::vector<double> OHp(my_fields->grid_dimension[0]);
  std::vector<double> H2Op(my_fields->grid_dimension[0]);
  std::vector<double> O2p(my_fields->grid_dimension[0]);
  std::vector<double> SiIp(my_fields->grid_dimension[0]);
  std::vector<double> SiOIp(my_fields->grid_dimension[0]);
  std::vector<double> SiO2Ip(my_fields->grid_dimension[0]);
  std::vector<double> CHp(my_fields->grid_dimension[0]);
  std::vector<double> CH2p(my_fields->grid_dimension[0]);
  std::vector<double> COIIp(my_fields->grid_dimension[0]);
  std::vector<double> OIIp(my_fields->grid_dimension[0]);
  std::vector<double> OHIIp(my_fields->grid_dimension[0]);
  std::vector<double> H2OIIp(my_fields->grid_dimension[0]);
  std::vector<double> H3OIIp(my_fields->grid_dimension[0]);
  std::vector<double> O2IIp(my_fields->grid_dimension[0]);
  std::vector<double> Mgp(my_fields->grid_dimension[0]);
  std::vector<double> Alp(my_fields->grid_dimension[0]);
  std::vector<double> Sp(my_fields->grid_dimension[0]);
  std::vector<double> Fep(my_fields->grid_dimension[0]);
  std::vector<gr_float> SiMp(my_fields->grid_dimension[0]);
  std::vector<gr_float> FeMp(my_fields->grid_dimension[0]);
  std::vector<gr_float> Mg2SiO4p(my_fields->grid_dimension[0]);
  std::vector<gr_float> MgSiO3p(my_fields->grid_dimension[0]);
  std::vector<gr_float> Fe3O4p(my_fields->grid_dimension[0]);
  std::vector<gr_float> ACp(my_fields->grid_dimension[0]);
  std::vector<gr_float> SiO2Dp(my_fields->grid_dimension[0]);
  std::vector<gr_float> MgOp(my_fields->grid_dimension[0]);
  std::vector<gr_float> FeSp(my_fields->grid_dimension[0]);
  std::vector<gr_float> Al2O3p(my_fields->grid_dimension[0]);
  std::vector<gr_float> reforgp(my_fields->grid_dimension[0]);
  std::vector<gr_float> volorgp(my_fields->grid_dimension[0]);
  std::vector<gr_float> H2Oicep(my_fields->grid_dimension[0]);

  std::vector<double> k125(my_fields->grid_dimension[0]);
  std::vector<double> k129(my_fields->grid_dimension[0]);
  std::vector<double> k130(my_fields->grid_dimension[0]);
  std::vector<double> k131(my_fields->grid_dimension[0]);
  std::vector<double> k132(my_fields->grid_dimension[0]);
  std::vector<double> k133(my_fields->grid_dimension[0]);
  std::vector<double> k134(my_fields->grid_dimension[0]);
  std::vector<double> k135(my_fields->grid_dimension[0]);
  std::vector<double> k136(my_fields->grid_dimension[0]);
  std::vector<double> k137(my_fields->grid_dimension[0]);
  std::vector<double> k148(my_fields->grid_dimension[0]);
  std::vector<double> k149(my_fields->grid_dimension[0]);
  std::vector<double> k150(my_fields->grid_dimension[0]);
  std::vector<double> k151(my_fields->grid_dimension[0]);
  std::vector<double> k152(my_fields->grid_dimension[0]);
  std::vector<double> k153(my_fields->grid_dimension[0]);
  std::vector<double> kz15(my_fields->grid_dimension[0]);
  std::vector<double> kz16(my_fields->grid_dimension[0]);
  std::vector<double> kz17(my_fields->grid_dimension[0]);
  std::vector<double> kz18(my_fields->grid_dimension[0]);
  std::vector<double> kz19(my_fields->grid_dimension[0]);
  std::vector<double> kz20(my_fields->grid_dimension[0]);
  std::vector<double> kz21(my_fields->grid_dimension[0]);
  std::vector<double> kz22(my_fields->grid_dimension[0]);
  std::vector<double> kz23(my_fields->grid_dimension[0]);
  std::vector<double> kz24(my_fields->grid_dimension[0]);
  std::vector<double> kz25(my_fields->grid_dimension[0]);
  std::vector<double> kz26(my_fields->grid_dimension[0]);
  std::vector<double> kz27(my_fields->grid_dimension[0]);
  std::vector<double> kz28(my_fields->grid_dimension[0]);
  std::vector<double> kz29(my_fields->grid_dimension[0]);
  std::vector<double> kz30(my_fields->grid_dimension[0]);
  std::vector<double> kz31(my_fields->grid_dimension[0]);
  std::vector<double> kz32(my_fields->grid_dimension[0]);
  std::vector<double> kz33(my_fields->grid_dimension[0]);
  std::vector<double> kz34(my_fields->grid_dimension[0]);
  std::vector<double> kz35(my_fields->grid_dimension[0]);
  std::vector<double> kz36(my_fields->grid_dimension[0]);
  std::vector<double> kz37(my_fields->grid_dimension[0]);
  std::vector<double> kz38(my_fields->grid_dimension[0]);
  std::vector<double> kz39(my_fields->grid_dimension[0]);
  std::vector<double> kz40(my_fields->grid_dimension[0]);
  std::vector<double> kz41(my_fields->grid_dimension[0]);
  std::vector<double> kz42(my_fields->grid_dimension[0]);
  std::vector<double> kz43(my_fields->grid_dimension[0]);
  std::vector<double> kz44(my_fields->grid_dimension[0]);
  std::vector<double> kz45(my_fields->grid_dimension[0]);
  std::vector<double> kz46(my_fields->grid_dimension[0]);
  std::vector<double> kz47(my_fields->grid_dimension[0]);
  std::vector<double> kz48(my_fields->grid_dimension[0]);
  std::vector<double> kz49(my_fields->grid_dimension[0]);
  std::vector<double> kz50(my_fields->grid_dimension[0]);
  std::vector<double> kz51(my_fields->grid_dimension[0]);
  std::vector<double> kz52(my_fields->grid_dimension[0]);
  std::vector<double> kz53(my_fields->grid_dimension[0]);
  std::vector<double> kz54(my_fields->grid_dimension[0]);
  std::vector<double> kdSiM(my_fields->grid_dimension[0]);
  std::vector<double> kdFeM(my_fields->grid_dimension[0]);
  std::vector<double> kdMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> kdMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> kdFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> kdAC(my_fields->grid_dimension[0]);
  std::vector<double> kdSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> kdMgO(my_fields->grid_dimension[0]);
  std::vector<double> kdFeS(my_fields->grid_dimension[0]);
  std::vector<double> kdAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> kdreforg(my_fields->grid_dimension[0]);
  std::vector<double> kdvolorg(my_fields->grid_dimension[0]);
  std::vector<double> kdH2Oice(my_fields->grid_dimension[0]);

  // Cooling/heating row locals

  std::vector<double> ceHI(my_fields->grid_dimension[0]);
  std::vector<double> ceHeI(my_fields->grid_dimension[0]);
  std::vector<double> ceHeII(my_fields->grid_dimension[0]);
  std::vector<double> ciHI(my_fields->grid_dimension[0]);
  std::vector<double> ciHeI(my_fields->grid_dimension[0]);
  std::vector<double> ciHeIS(my_fields->grid_dimension[0]);
  std::vector<double> ciHeII(my_fields->grid_dimension[0]);
  std::vector<double> reHII(my_fields->grid_dimension[0]);
  std::vector<double> reHeII1(my_fields->grid_dimension[0]);
  std::vector<double> reHeII2(my_fields->grid_dimension[0]);
  std::vector<double> reHeIII(my_fields->grid_dimension[0]);
  std::vector<double> brem(my_fields->grid_dimension[0]);
  std::vector<double> edot(my_fields->grid_dimension[0]);
  std::vector<double> hyd01k(my_fields->grid_dimension[0]);
  std::vector<double> h2k01(my_fields->grid_dimension[0]);
  std::vector<double> vibh(my_fields->grid_dimension[0]);
  std::vector<double> roth(my_fields->grid_dimension[0]);
  std::vector<double> rotl(my_fields->grid_dimension[0]);
  std::vector<double> gpldl(my_fields->grid_dimension[0]);
  std::vector<double> gphdl(my_fields->grid_dimension[0]);
  std::vector<double> hdlte(my_fields->grid_dimension[0]);
  std::vector<double> hdlow(my_fields->grid_dimension[0]);
  std::vector<double> cieco(my_fields->grid_dimension[0]);

  // Iteration mask

  gr_mask_type anydust;
  std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_tmp(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_nr(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);
  std::vector<int> imp_eng(my_fields->grid_dimension[0]);

  
  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Set error indicator (we will return this value)
  int ierr = GR_SUCCESS;

  // Set flag for dust-related options

  if ((my_chemistry->h2_on_dust > 0)  ||  (my_chemistry->dust_chemistry > 0))  {
    anydust = MASK_TRUE;
  } else {
    anydust = MASK_FALSE;
  }

  // ignore metal chemistry/cooling below this metallicity
  min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;
      
  // Set units

  dom      = urho*(std::pow(my_units->a_value,3))/mh_local_var;
  tbase1   = my_units->time_units;
  xbase1   = uxyz/(my_units->a_value*my_units->a_units);    // uxyz is [x]*a      = [x]*[a]*a'        '
  dbase1   = urho*std::pow((my_units->a_value*my_units->a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '
  coolunit = (std::pow(my_units->a_units,5) * std::pow(xbase1,2) * std::pow(mh_local_var,2)) / (std::pow(tbase1,3) * dbase1);
  uvel = (uxyz/my_units->a_value) / my_units->time_units;
  // chunit = (1.60218e-12_DKIND)/(2._DKIND*uvel*uvel*mh)   ! 1 eV per H2 formed
  chunit = (1.60218e-12)/(uvel*uvel*mh_local_var);            // 1 eV per REACTION (Feb 2020, Gen Chiaki)

  dx_cgs = my_fields->grid_dx * xbase1;
  c_ljeans = std::sqrt((my_chemistry->Gamma * pi_local_var * kboltz_grflt) /
       (GravConst_grflt * mh_local_var * dbase1));

  dlogtem = (std::log(my_chemistry->TemperatureEnd) - std::log(my_chemistry->TemperatureStart))/(double)(my_chemistry->NumberOfTemperatureBins-1 );

  // We better make consistent at first GC202002

  if (my_chemistry->primordial_chemistry > 0)  {

#define ABUNDANCE_CORRECTION
#ifdef ABUNDANCE_CORRECTION
     FORTRAN_NAME(make_consistent_g)(de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                         d.data(), &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->HydrogenFractionByMass, &my_chemistry->DeuteriumToHydrogenRatio,
                        &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &dom,
                        my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                        my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                        my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                        my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                        my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                        my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                        my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                        my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                        my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                        my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                        &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
                        my_fields->local_ISM_metal_density,
                        my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                        my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                        my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density,
                        &my_rates->SN0_N,
                        my_rates->SN0_XC, my_rates->SN0_XO, my_rates->SN0_XMg, my_rates->SN0_XAl, my_rates->SN0_XSi,
                        my_rates->SN0_XS, my_rates->SN0_XFe,
                        my_rates->SN0_fC, my_rates->SN0_fO, my_rates->SN0_fMg, my_rates->SN0_fAl, my_rates->SN0_fSi,
                        my_rates->SN0_fS, my_rates->SN0_fFe
                          );
#endif

  }

  // Convert densities from comoving to proper

  if (my_units->comoving_coordinates == 1)  {

    factor = (gr_float)(std::pow(my_units->a_value,(-3)) );

     FORTRAN_NAME(scale_fields_g)(
      &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry,
      &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->multi_metals, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation, &factor,
      &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2], &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
      d.data(), de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
      my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
      my_fields->metal_density, my_fields->dust_density,
      my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
      my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
      my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
      my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
      my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
      my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
      my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
      my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
      my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
      my_fields->local_ISM_metal_density, my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
      my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
      my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);

  }

#ifdef ABUNDANCE_CORRECTION
   FORTRAN_NAME(ceiling_species_g)(d.data(), de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                       my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                       &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                       &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->use_dust_density_field,
                      my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                      &my_chemistry->metal_abundances, &my_chemistry->metal_chemistry, &my_chemistry->dust_species, &my_chemistry->multi_metals,
                      &my_chemistry->grain_growth, &my_chemistry->dust_sublimation,
                      my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                      my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                      my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                      my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                      my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                      my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                      my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                      my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                      my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                      my_fields->local_ISM_metal_density,
                      my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                      my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                      my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);
#endif


  // Loop over zones, and do an entire i-column in one go
  dk = my_fields->grid_end[2] - my_fields->grid_start[2] + 1;
  dj = my_fields->grid_end[1] - my_fields->grid_start[1] + 1;

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: ! ierr is declared as shared and should be modified with atomic operation
  //_// PORT: !$omp parallel do schedule(runtime) private(
  //_// PORT: !$omp&   i, j, k, iter,
  //_// PORT: !$omp&   ttmin, energy, comp1, comp2,
  //_// PORT: !$omp&   heq1, heq2, eqk221, eqk222, eqk131, eqk132,
  //_// PORT: !$omp&   eqt1, eqt2, eqtdef, dheq, heq,
  //_// PORT: !$omp&   dtit, ttot, p2d, tgas,
  //_// PORT: !$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
  //_// PORT: !$omp&   ddom,
  //_// PORT: !$omp&   olddtit,
  //_// PORT: !$omp&   HIp, HIIp, HeIp, HeIIp, HeIIIp,
  //_// PORT: !$omp&   HMp, H2Ip, H2IIp,
  //_// PORT: !$omp&   dep, dedot,HIdot, dedot_prev,
  //_// PORT: !$omp&   DIp, DIIp, HDIp, HIdot_prev,
  //_// PORT: !$omp&   k24shield, k25shield, k26shield,
  //_// PORT: !$omp&   k28shield, k29shield, k30shield,
  //_// PORT: !$omp&   k31shield,
  //_// PORT: !$omp&   k1 , k2 , k3 , k4 , k5,
  //_// PORT: !$omp&   k6 , k7 , k8 , k9 , k10,
  //_// PORT: !$omp&   k11, k12, k13, k14, k15,
  //_// PORT: !$omp&   k16, k17, k18, k19, k22,
  //_// PORT: !$omp&   k50, k51, k52, k53, k54,
  //_// PORT: !$omp&   k55, k56, k57, k58, k13dd, h2dust,
  //_// PORT: !$omp&   ncrn, ncrd1, ncrd2,
  //_// PORT: !$omp&   ceHI, ceHeI, ceHeII,
  //_// PORT: !$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
  //_// PORT: !$omp&   reHII, reHeII1, reHeII2, reHeIII,
  //_// PORT: !$omp&   brem, edot,
  //_// PORT: !$omp&   hyd01k, h2k01, vibh, roth, rotl,
  //_// PORT: !$omp&   gpldl, gphdl, hdlte, hdlow, cieco,
  //_// PORT: !$omp&   itmask, itmask_metal )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    // TODO: move more relevant variable declarations to here to replace the
    //       OMP private-clause

    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf =
      grackle::impl::new_Cool1DMultiScratchBuf(my_fields->grid_dimension[0]);

    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (t = 0; t<=(dk * dj - 1); t++) {
      k = t/dj      + my_fields->grid_start[2]+1;
      j = grackle::impl::mod(t,dj) + my_fields->grid_start[1]+1;

      // tolerance = 1.0e-06_DKIND * dt

      // Initialize iteration mask to true for all cells.

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        itmask[i-1] = MASK_TRUE;
      }

      // If we are using coupled radiation with intermediate stepping,
      // set iteration mask to include only cells with radiation in the
      // intermediate coupled chemistry / energy step
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_coupled_rate_solver == 1  &&  my_chemistry->radiative_transfer_intermediate_step == 1)  {
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if (kphHI(i-1,j-1,k-1) > 0)  {
              itmask[i-1] = MASK_TRUE;
            } else {
              itmask[i-1] = MASK_FALSE;
            }
          }
        }

        // Normal rate solver, but don't double count cells with radiation
        if (my_chemistry->radiative_transfer_coupled_rate_solver == 1  &&  my_chemistry->radiative_transfer_intermediate_step == 0)  {
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if (kphHI(i-1,j-1,k-1) > 0)  {
              itmask[i-1] = MASK_FALSE;
            } else {
              itmask[i-1] = MASK_TRUE;
            }
          }
        }
      }

      // Set time elapsed to zero for each cell in 1D section

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        ttot[i-1] = 0.;
      }

      // A useful slice variable since we do this a lot

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        ddom[i-1] = d(i-1,j-1,k-1) * dom;
      }

      // ------------------ Loop over subcycles ----------------

      for (iter = 1; iter<=(my_chemistry->max_iterations); iter++) {

        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          if (itmask[i-1] != MASK_FALSE)  {
            dtit[i-1] = huge8;
          }
        }

        // Compute the cooling rate, tgas, tdust, and metallicity for this row

         FORTRAN_NAME(cool1d_multi_g)(
                  d.data(), e.data(), my_fields->x_velocity, my_fields->y_velocity, my_fields->z_velocity, de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                  &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                  &my_units->comoving_coordinates, &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->metal_cooling,
                  &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
                  &my_fields->grid_rank, &my_fields->grid_start[0], &my_fields->grid_end[0], &j, &k, &my_chemistry->ih2co, &my_chemistry->ipiht, &iter, &my_chemistry->photoelectric_heating,
                  &my_units->a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
                  &utem, &uxyz, &my_units->a_units, &urho, &my_units->time_units,
                  &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
                  my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
                  my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
                  my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
                  &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates->comp_xray, &my_uvb_rates->temp_xray,
                  &my_uvb_rates->piHI, &my_uvb_rates->piHeI, &my_uvb_rates->piHeII, &comp1, &comp2,
                  my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                  my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
                  hyd01k.data(), h2k01.data(), vibh.data(), roth.data(), rotl.data(),
                  my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, gpldl.data(), gphdl.data(),
                  my_rates->HDlte, my_rates->HDlow, hdlte.data(), hdlow.data(),
                  my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp, my_rates->GAel,
                  my_rates->H2LTE, my_rates->gas_grain,
                  ceHI.data(), ceHeI.data(), ceHeII.data(), ciHI.data(), ciHeI.data(), ciHeIS.data(), ciHeII.data(),
                  reHII.data(), reHeII1.data(), reHeII2.data(), reHeIII.data(), brem.data(),
                  logTlininterp_buf.indixe, logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.logtem, logTlininterp_buf.tdef, edot.data(),
                  tgas.data(), cool1dmulti_buf.tgasold, mmw.data(), p2d.data(), tdust.data(), metallicity.data(),
                  dust2gas.data(), rhoH.data(), cool1dmulti_buf.mynh, cool1dmulti_buf.myde,
                  cool1dmulti_buf.gammaha_eff, cool1dmulti_buf.gasgr_tdust, cool1dmulti_buf.regr,
                  &my_chemistry->self_shielding_method, &my_uvb_rates->crsHI, &my_uvb_rates->crsHeI, &my_uvb_rates->crsHeII,
                  &my_uvb_rates->k24, &my_uvb_rates->k26,
                  &my_chemistry->use_radiative_transfer, my_fields->RT_heating_rate,
                  &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, cieco.data(),
                  &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
                  &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
                  my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
                  &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
                  &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
                  my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
                  &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
                  &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, my_fields->volumetric_heating_rate, my_fields->specific_heating_rate,
                  &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, my_fields->temperature_floor,
                  &my_chemistry->use_isrf_field, my_fields->isrf_habing, itmask.data(), itmask_metal.data(),
                 &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature,
                 my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                 my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                 my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                 my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                 my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                 my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                 my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                 my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                 my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                 my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                 my_rates->cieY06,
                 my_rates->LH2.props.dimension, &my_rates->LH2.props.data_size,
                 my_rates->LH2.props.parameters[0], my_rates->LH2.props.parameters[1], my_rates->LH2.props.parameters[2],
                 &my_rates->LH2.props.parameter_spacing[0], &my_rates->LH2.props.parameter_spacing[1], &my_rates->LH2.props.parameter_spacing[2], my_rates->LH2.data,
                 my_rates->LHD.props.dimension, &my_rates->LHD.props.data_size,
                 my_rates->LHD.props.parameters[0], my_rates->LHD.props.parameters[1], my_rates->LHD.props.parameters[2],
                 &my_rates->LHD.props.parameter_spacing[0], &my_rates->LHD.props.parameter_spacing[1], &my_rates->LHD.props.parameter_spacing[2], my_rates->LHD.data,
                 my_rates->LCI.props.dimension, &my_rates->LCI.props.data_size,
                 my_rates->LCI.props.parameters[0], my_rates->LCI.props.parameters[1], my_rates->LCI.props.parameters[2],
                 &my_rates->LCI.props.parameter_spacing[0], &my_rates->LCI.props.parameter_spacing[1], &my_rates->LCI.props.parameter_spacing[2], my_rates->LCI.data,
                 my_rates->LCII.props.dimension, &my_rates->LCII.props.data_size,
                 my_rates->LCII.props.parameters[0], my_rates->LCII.props.parameters[1], my_rates->LCII.props.parameters[2],
                 &my_rates->LCII.props.parameter_spacing[0], &my_rates->LCII.props.parameter_spacing[1], &my_rates->LCII.props.parameter_spacing[2], my_rates->LCII.data,
                 my_rates->LOI.props.dimension, &my_rates->LOI.props.data_size,
                 my_rates->LOI.props.parameters[0], my_rates->LOI.props.parameters[1], my_rates->LOI.props.parameters[2],
                 &my_rates->LOI.props.parameter_spacing[0], &my_rates->LOI.props.parameter_spacing[1], &my_rates->LOI.props.parameter_spacing[2], my_rates->LOI.data,
                 my_rates->LCO.props.dimension, &my_rates->LCO.props.data_size,
                 my_rates->LCO.props.parameters[0], my_rates->LCO.props.parameters[1], my_rates->LCO.props.parameters[2],
                 &my_rates->LCO.props.parameter_spacing[0], &my_rates->LCO.props.parameter_spacing[1], &my_rates->LCO.props.parameter_spacing[2], my_rates->LCO.data,
                 my_rates->LOH.props.dimension, &my_rates->LOH.props.data_size,
                 my_rates->LOH.props.parameters[0], my_rates->LOH.props.parameters[1], my_rates->LOH.props.parameters[2],
                 &my_rates->LOH.props.parameter_spacing[0], &my_rates->LOH.props.parameter_spacing[1], &my_rates->LOH.props.parameter_spacing[2], my_rates->LOH.data,
                 my_rates->LH2O.props.dimension, &my_rates->LH2O.props.data_size,
                 my_rates->LH2O.props.parameters[0], my_rates->LH2O.props.parameters[1], my_rates->LH2O.props.parameters[2],
                 &my_rates->LH2O.props.parameter_spacing[0], &my_rates->LH2O.props.parameter_spacing[1], &my_rates->LH2O.props.parameter_spacing[2], my_rates->LH2O.data,
                 my_rates->alphap.props.dimension, &my_rates->alphap.props.data_size,
                 my_rates->alphap.props.parameters[0], my_rates->alphap.props.parameters[1], &my_rates->alphap.props.parameter_spacing[0], &my_rates->alphap.props.parameter_spacing[1],
                 my_rates->alphap.data,
                 &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
                 my_fields->local_ISM_metal_density,
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
                 grain_temperatures.SiM, grain_temperatures.FeM, grain_temperatures.Mg2SiO4, grain_temperatures.MgSiO3, grain_temperatures.Fe3O4,
                 grain_temperatures.AC, grain_temperatures.SiO2D, grain_temperatures.MgO, grain_temperatures.FeS, grain_temperatures.Al2O3,
                 grain_temperatures.reforg, grain_temperatures.volorg, grain_temperatures.H2Oice,
                 my_rates->gas_grain2, &my_rates->gamma_isrf2
            );

        if (my_chemistry->primordial_chemistry == 0)  {
          // This is some basic book-keeping to ensure that itmask_tmp has
          // sensible values when ispecies is 0
          // -> There's room for optimization: when ispecies is 0, there is
          //    need to ever touch this variable. But, for now we focus on
          //    correct behavior before implementing this optimization
          std::memcpy(itmask_tmp.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);

        } else {

          // Look-up rates as a function of temperature for 1D set of zones
          //  (maybe should add itmask to this call)

           FORTRAN_NAME(lookup_cool_rates1d_g)(&my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->NumberOfTemperatureBins, &j, &k,
                   &my_fields->grid_start[0], &my_fields->grid_end[0], &my_chemistry->three_body_rate,
                   &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &anydust,
                   &my_chemistry->H2_self_shielding, &my_chemistry->self_shielding_method,
                   tgas.data(), mmw.data(), d.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                   my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density,
                   tdust.data(), dust2gas.data(),
                   my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
                   my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
                   my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
                   my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54, my_rates->k55, my_rates->k56,
                   my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
                   my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
                   &my_uvb_rates->crsHI, &my_uvb_rates->crsHeI, &my_uvb_rates->crsHeII, &my_uvb_rates->piHI, &my_uvb_rates->piHeI,
                   k1.data(), k2.data(), k3.data(), k4.data(), k5.data(), k6.data(), k7.data(), k8.data(), k9.data(), k10.data(),
                   k11.data(), k12.data(), k13.data(), k14.data(), k15.data(), k16.data(), k17.data(), k18.data(),
                   k19.data(), k22.data(), &my_uvb_rates->k24, &my_uvb_rates->k25, &my_uvb_rates->k26, &my_uvb_rates->k28, &my_uvb_rates->k29, &my_uvb_rates->k30, &my_uvb_rates->k31,
                   k50.data(), k51.data(), k52.data(), k53.data(), k54.data(), k55.data(), k56.data(), k57.data(),
                   k58.data(), k13dd.data(), k24shield.data(), k25shield.data(), k26shield.data(),
                   k28shield.data(), k29shield.data(), k30shield.data(),
                   k31shield.data(), h2dust.data(), ncrn.data(), ncrd1.data(), ncrd2.data(),
                   logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.tdef, logTlininterp_buf.logtem, logTlininterp_buf.indixe,
                   &dom, &coolunit, &tbase1, &xbase1, &dx_cgs, &c_ljeans,
                   &my_chemistry->use_radiative_transfer, my_fields->RT_H2_dissociation_rate, my_fields->H2_self_shielding_length, itmask.data(),
                   itmask_metal.data(),
                   &my_chemistry->HydrogenFractionByMass, my_fields->metal_density,
                   my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density, &imetal, &my_chemistry->metal_chemistry, &my_chemistry->grain_growth,
                   my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                   my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                   my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                   my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                   my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                   my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                   my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                   my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                   my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                   my_rates->k125, my_rates->k129, my_rates->k130, my_rates->k131, my_rates->k132,
                   my_rates->k133, my_rates->k134, my_rates->k135, my_rates->k136, my_rates->k137,
                   my_rates->k148, my_rates->k149, my_rates->k150, my_rates->k151, my_rates->k152,
                   my_rates->k153,
                   my_rates->kz15, my_rates->kz16, my_rates->kz17, my_rates->kz18, my_rates->kz19,
                   my_rates->kz20, my_rates->kz21, my_rates->kz22, my_rates->kz23, my_rates->kz24,
                   my_rates->kz25, my_rates->kz26, my_rates->kz27, my_rates->kz28, my_rates->kz29,
                   my_rates->kz30, my_rates->kz31, my_rates->kz32, my_rates->kz33, my_rates->kz34,
                   my_rates->kz35, my_rates->kz36, my_rates->kz37, my_rates->kz38, my_rates->kz39,
                   my_rates->kz40, my_rates->kz41, my_rates->kz42, my_rates->kz43, my_rates->kz44,
                   my_rates->kz45, my_rates->kz46, my_rates->kz47, my_rates->kz48, my_rates->kz49,
                   my_rates->kz50, my_rates->kz51, my_rates->kz52, my_rates->kz53, my_rates->kz54,
                   k125.data(),  k129.data(),  k130.data(),  k131.data(),  k132.data(),
                   k133.data(),  k134.data(),  k135.data(),  k136.data(),  k137.data(),
                   k148.data(),  k149.data(),  k150.data(),  k151.data(),  k152.data(),
                   k153.data(),
                   kz15.data(),  kz16.data(),  kz17.data(),  kz18.data(),  kz19.data(),
                   kz20.data(),  kz21.data(),  kz22.data(),  kz23.data(),  kz24.data(),
                   kz25.data(),  kz26.data(),  kz27.data(),  kz28.data(),  kz29.data(),
                   kz30.data(),  kz31.data(),  kz32.data(),  kz33.data(),  kz34.data(),
                   kz35.data(),  kz36.data(),  kz37.data(),  kz38.data(),  kz39.data(),
                   kz40.data(),  kz41.data(),  kz42.data(),  kz43.data(),  kz44.data(),
                   kz45.data(),  kz46.data(),  kz47.data(),  kz48.data(),  kz49.data(),
                   kz50.data(),  kz51.data(),  kz52.data(),  kz53.data(),  kz54.data(),
                   &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
                   my_fields->local_ISM_metal_density,
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
                   my_rates->h2dustS, my_rates->h2dustC, rhoH.data(), my_rates->grain_growth_rate, &dt,
                   kdSiM.data(), kdFeM.data(), kdMg2SiO4.data(),
                   kdMgSiO3.data(), kdFe3O4.data(), kdAC.data(), kdSiO2D.data(), kdMgO.data(), kdFeS.data(),
                   kdAl2O3.data(), kdreforg.data(), kdvolorg.data(), kdH2Oice.data(),
                   grain_temperatures.SiM, grain_temperatures.FeM, grain_temperatures.Mg2SiO4, grain_temperatures.MgSiO3, grain_temperatures.Fe3O4,
                   grain_temperatures.AC, grain_temperatures.SiO2D, grain_temperatures.MgO, grain_temperatures.FeS, grain_temperatures.Al2O3,
                   grain_temperatures.reforg, grain_temperatures.volorg, grain_temperatures.H2Oice, &my_chemistry->radiative_transfer_use_H2_shielding,
                   &my_chemistry->H2_custom_shielding, my_fields->H2_custom_shielding_factor
          );

          // Compute dedot and HIdot, the rates of change of de and HI
          //   (should add itmask to this call)

           FORTRAN_NAME(rate_timestep_g)(
                         dedot.data(), HIdot.data(), &my_chemistry->primordial_chemistry, &anydust,
                         de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, d.data(),
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density,
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &j, &k,
                         k1.data(), k2.data(), k3.data(), k4.data(), k5.data(), k6.data(), k7.data(), k8.data(), k9.data(), k10.data(), k11.data(),
                         k12.data(), k13.data(), k14.data(), k15.data(), k16.data(), k17.data(), k18.data(), k19.data(), k22.data(),
                         &my_uvb_rates->k24, &my_uvb_rates->k25, &my_uvb_rates->k26, &my_uvb_rates->k27, &my_uvb_rates->k28, &my_uvb_rates->k29, &my_uvb_rates->k30,
                         k50.data(), k51.data(), k52.data(), k53.data(), k54.data(), k55.data(), k56.data(), k57.data(), k58.data(),
                         h2dust.data(), ncrn.data(), ncrd1.data(), ncrd2.data(), rhoH.data(),
                         k24shield.data(), k25shield.data(), k26shield.data(),
                         k28shield.data(), k29shield.data(), k30shield.data(), k31shield.data(),
                         &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
                         kphHI.data(), my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate,
                         itmask.data(), edot.data(), &chunit, &dom, my_fields->metal_density,
                        my_fields->HDI_density, &my_chemistry->metal_chemistry, my_fields->CI_density, my_fields->OI_density, my_fields->OH_density, my_fields->CO_density, my_fields->H2O_density,
                        &my_chemistry->radiative_transfer_HDI_dissociation, my_fields->RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, my_fields->RT_CI_ionization_rate, my_fields->RT_OI_ionization_rate,
                        &my_chemistry->radiative_transfer_metal_dissociation, my_fields->RT_CO_dissociation_rate, my_fields->RT_OH_dissociation_rate, my_fields->RT_H2O_dissociation_rate
                              );

          // move itmask temporary array
          // then split cells with low densities
          //    => Gauss-Seidel scheme
          //              and with high densities
          //    => Newton-Raphson scheme

          std::memcpy(itmask_tmp.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);
          std::memcpy(itmask_nr.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if ( itmask_tmp[i-1] != MASK_FALSE )  {

              if ( ( (imetal == 0)
                &&  (ddom[i-1] < 1.e8) )
               ||  ( (imetal == 1)
                &&  ( ( (metallicity[i-1] <= min_metallicity)
                    &&  (ddom[i-1] < 1.e8) )
                   ||  ( (metallicity[i-1] > min_metallicity)
                    &&  (ddom[i-1] < 1.e6) ) ) ) )  {
                itmask_nr[i-1] = MASK_FALSE;
              } else {
                itmask[i-1] = MASK_FALSE;
              }

            }
          }

          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if (itmask_nr[i-1] != MASK_FALSE)  {
              if ( (my_chemistry->with_radiative_cooling == 1)  &&  (my_chemistry->primordial_chemistry > 1)  && 
                   ((ddom[i-1] > 1.e7)
               && (tgas[i-1] > 1650.e0)) )  {
                imp_eng[i-1] = 1;
              } else {
                imp_eng[i-1] = 0;
              }
            }
          }

          // Find timestep that keeps relative chemical changes below 10%

          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if (itmask[i-1] != MASK_FALSE)  {
              // Bound from below to prevent numerical errors
               
              if (std::fabs(dedot[i-1]) < tiny8)
                         { dedot[i-1] = std::fmin(tiny_fortran_val,de(i-1,j-1,k-1)); }
              if (std::fabs(HIdot[i-1]) < tiny8)
                         { HIdot[i-1] = std::fmin(tiny_fortran_val,HI(i-1,j-1,k-1)); }

              // If the net rate is almost perfectly balanced then set
              //     it to zero (since it is zero to available precision)

              if (std::fmin(std::fabs(k1[i-1]* de(i-1,j-1,k-1)*HI(i-1,j-1,k-1)),
                      std::fabs(k2[i-1]*HII(i-1,j-1,k-1)*de(i-1,j-1,k-1)))/
                  std::fmax(std::fabs(dedot[i-1]),std::fabs(HIdot[i-1])) >
                   1.0e6)  {
                dedot[i-1] = tiny8;
                HIdot[i-1] = tiny8;
              }

              // If the iteration count is high then take the smaller of
              //   the calculated dedot and last time step's actual dedot.
              //   This is intended to get around the problem of a low
              //   electron or HI fraction which is in equilibrium with high
              //   individual terms (which all nearly cancel).

              if (iter > 50)  {
                dedot[i-1] = std::fmin(std::fabs(dedot[i-1]), std::fabs(dedot_prev[i-1]));
                HIdot[i-1] = std::fmin(std::fabs(HIdot[i-1]), std::fabs(HIdot_prev[i-1]));
              }

              // compute minimum rate timestep

              olddtit = dtit[i-1];
              dtit[i-1] = grackle::impl::fmin(std::fabs(0.1*de(i-1,j-1,k-1)/dedot[i-1]),
                   std::fabs(0.1*HI(i-1,j-1,k-1)/HIdot[i-1]),
                   dt-ttot[i-1], 0.5*dt);

              if (ddom[i-1] > 1.e8  && 
                   edot[i-1] > 0.       && 
                  my_chemistry->primordial_chemistry > 1)  {
                // Equilibrium value for H is:
                // H = (-1._DKIND / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
                // We now want this to change by 10% or less, but we're only
                // differentiating by dT.  We have de/dt.  We need dT/de.
                // T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
                // dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
                // dH_eq / dT (see above; we can calculate the derivative here)
                // dT / de = utem * (gamma - 1._DKIND) / N == tgas / p2d
                // de / dt = edot
                // Now we use our estimate of dT/de to get the estimated
                // difference in the equilibrium
                eqt2 = std::fmin(std::log(tgas[i-1]) + 0.1*dlogtem, logTlininterp_buf.t2[i-1]);
                eqtdef = (eqt2 - logTlininterp_buf.t1[i-1])/(logTlininterp_buf.t2[i-1] - logTlininterp_buf.t1[i-1]);
                eqk222 = my_rates->k22[logTlininterp_buf.indixe[i-1]-1] +
                  (my_rates->k22[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i-1]-1])*eqtdef;
                eqk132 = my_rates->k13[logTlininterp_buf.indixe[i-1]-1] +
                  (my_rates->k13[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i-1]-1])*eqtdef;
                heq2 = (-1. / (4.*eqk222)) * (eqk132-
                     std::sqrt(8.*eqk132*eqk222*
                     my_chemistry->HydrogenFractionByMass*d(i-1,j-1,k-1)+std::pow(eqk132,2.)));

                eqt1 = std::fmax(std::log(tgas[i-1]) - 0.1*dlogtem, logTlininterp_buf.t1[i-1]);
                eqtdef = (eqt1 - logTlininterp_buf.t1[i-1])/(logTlininterp_buf.t2[i-1] - logTlininterp_buf.t1[i-1]);
                eqk221 = my_rates->k22[logTlininterp_buf.indixe[i-1]-1] +
                  (my_rates->k22[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i-1]-1])*eqtdef;
                eqk131 = my_rates->k13[logTlininterp_buf.indixe[i-1]-1] +
                  (my_rates->k13[logTlininterp_buf.indixe[i-1]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i-1]-1])*eqtdef;
                heq1 = (-1. / (4.*eqk221)) * (eqk131-
                     std::sqrt(8.*eqk131*eqk221*
                     my_chemistry->HydrogenFractionByMass*d(i-1,j-1,k-1)+std::pow(eqk131,2.)));

                dheq = (std::fabs(heq2-heq1)/(std::exp(eqt2) - std::exp(eqt1)))
                     * (tgas[i-1]/p2d[i-1]) * edot[i-1];
                heq = (-1. / (4.*k22[i-1])) * (k13[i-1]-
                     std::sqrt(8.*k13[i-1]*k22[i-1]*
                     my_chemistry->HydrogenFractionByMass*d(i-1,j-1,k-1)+std::pow(k13[i-1],2.)));
                dtit[i-1] = std::fmin(dtit[i-1], 0.1*heq/dheq);
              }
              if (iter>10LL)  {
                dtit[i-1] = std::fmin(olddtit*1.5, dtit[i-1]);
              }

            } else if ((itmask_nr[i-1]!=MASK_FALSE) && 
                     (imp_eng[i-1]==0))  {
              dtit[i-1] = grackle::impl::fmin(std::fabs(0.1*e(i-1,j-1,k-1)/edot[i-1]*d(i-1,j-1,k-1)),
                   dt-ttot[i-1], 0.5*dt);

            } else if ((itmask_nr[i-1]!=MASK_FALSE) && 
                     (imp_eng[i-1]==1))  {
              dtit[i-1] = dt - ttot[i-1];

            } else {
              dtit[i-1] = dt;
            }
          }

        }

        // Compute maximum timestep for cooling/heating

        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          if (itmask[i-1] != MASK_FALSE)  {
            // Set energy per unit volume of this cell based in the pressure
            // (the gamma used here is the right one even for H2 since p2d
            //  is calculated with this gamma).

            energy = std::fmax(p2d[i-1]/(my_chemistry->Gamma-1.), tiny8);

            // If the temperature is at the bottom of the temperature look-up
            // table and edot < 0, then shut off the cooling.

            if (tgas[i-1] <= 1.01*my_chemistry->TemperatureStart  && 
                 edot[i-1] < 0.)
                 { edot[i-1] = tiny8; }
            if (std::fabs(edot[i-1]) < tiny8) { edot[i-1] = tiny8; }

            // Compute timestep for 10% change

            dtit[i-1] = grackle::impl::fmin((double)(std::fabs(0.1*
              energy/edot[i-1]) ),
              dt-ttot[i-1], dtit[i-1]);

            if (dtit[i-1] != dtit[i-1])  {
              OMP_PRAGMA_CRITICAL
              {
                eprintf("HUGE dtit ::  %g %g %g %g %g %g %g\n",
                        energy,
                        edot [ i-1 ],
                        dtit [ i-1 ],
                        dt,
                        ttot [ i-1 ],
                        std::fabs ( 0.1 * energy / edot [ i-1 ] ),
                        (double) ( std::fabs ( 0.1 * energy / edot [ i-1 ] ) ));
              }
            }

          }
        }

        // Update total and gas energy

        if (my_chemistry->with_radiative_cooling == 1)  {
          for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
            if (itmask[i-1] != MASK_FALSE)  {
              e(i-1,j-1,k-1)  = e(i-1,j-1,k-1) +
                      (gr_float)(edot[i-1]/d(i-1,j-1,k-1)*dtit[i-1] );

            }
          }
        }

        if (my_chemistry->primordial_chemistry > 0)  {

          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask)

           FORTRAN_NAME(step_rate_g)(de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, d.data(),
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, dtit.data(),
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &j, &k,
                         &my_chemistry->primordial_chemistry, &anydust,
                         k1.data(), k2.data(), k3.data(), k4.data(), k5.data(), k6.data(), k7.data(), k8.data(), k9.data(), k10.data(), k11.data(),
                         k12.data(), k13.data(), k14.data(), k15.data(), k16.data(), k17.data(), k18.data(), k19.data(), k22.data(),
                         &my_uvb_rates->k24, &my_uvb_rates->k25, &my_uvb_rates->k26, &my_uvb_rates->k27, &my_uvb_rates->k28, &my_uvb_rates->k29, &my_uvb_rates->k30,
                         k50.data(), k51.data(), k52.data(), k53.data(), k54.data(), k55.data(), k56.data(), k57.data(), k58.data(),
                         h2dust.data(), rhoH.data(),
                         k24shield.data(), k25shield.data(), k26shield.data(),
                         k28shield.data(), k29shield.data(), k30shield.data(), k31shield.data(),
                         HIp.data(), HIIp.data(), HeIp.data(), HeIIp.data(), HeIIIp.data(), dep.data(),
                         HMp.data(), H2Ip.data(), H2IIp.data(), DIp.data(), DIIp.data(), HDIp.data(),
                         dedot_prev.data(), HIdot_prev.data(),
                         &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
                         kphHI.data(), my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate,
                         itmask.data(), itmask_metal.data(),
                        my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density, &imetal, my_fields->metal_density,
                        &my_chemistry->metal_chemistry, &my_chemistry->dust_species, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation,
                        my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                        my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                        my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                        my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                        my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                        my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                        my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                        my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                        my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                        k125.data(), k129.data(), k130.data(), k131.data(), k132.data(),
                        k133.data(), k134.data(), k135.data(), k136.data(), k137.data(),
                        k148.data(), k149.data(), k150.data(), k151.data(), k152.data(),
                        k153.data(),
                        kz15.data(),  kz16.data(),  kz17.data(),  kz18.data(),  kz19.data(),
                        kz20.data(),  kz21.data(),  kz22.data(),  kz23.data(),  kz24.data(),
                        kz25.data(),  kz26.data(),  kz27.data(),  kz28.data(),  kz29.data(),
                        kz30.data(),  kz31.data(),  kz32.data(),  kz33.data(),  kz34.data(),
                        kz35.data(),  kz36.data(),  kz37.data(),  kz38.data(),  kz39.data(),
                        kz40.data(),  kz41.data(),  kz42.data(),  kz43.data(),  kz44.data(),
                        kz45.data(),  kz46.data(),  kz47.data(),  kz48.data(),  kz49.data(),
                        kz50.data(),  kz51.data(),  kz52.data(),  kz53.data(),  kz54.data(),
                        DMp.data(), HDIIp.data(), HeHIIp.data(),
                        CIp.data(), CIIp.data(), COp.data(), CO2p.data(),
                        OIp.data(), OHp.data(), H2Op.data(), O2p.data(),
                        SiIp.data(), SiOIp.data(), SiO2Ip.data(),
                        CHp.data(), CH2p.data(), COIIp.data(), OIIp.data(),
                        OHIIp.data(), H2OIIp.data(), H3OIIp.data(), O2IIp.data(),
                        Mgp.data(), Alp.data(), Sp.data(), Fep.data(),
                        SiMp.data(), FeMp.data(), Mg2SiO4p.data(), MgSiO3p.data(), Fe3O4p.data(),
                        ACp.data(), SiO2Dp.data(), MgOp.data(), FeSp.data(), Al2O3p.data(),
                        reforgp.data(), volorgp.data(), H2Oicep.data(),
                        kdSiM.data(), kdFeM.data(), kdMg2SiO4.data(), kdMgSiO3.data(), kdFe3O4.data(),
                        kdAC.data(), kdSiO2D.data(), kdMgO.data(), kdFeS.data(), kdAl2O3.data(),
                        kdreforg.data(), kdvolorg.data(), kdH2Oice.data(),
                        &my_chemistry->radiative_transfer_HDI_dissociation, my_fields->RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, my_fields->RT_CI_ionization_rate, my_fields->RT_OI_ionization_rate,
                        &my_chemistry->radiative_transfer_metal_dissociation, my_fields->RT_CO_dissociation_rate, my_fields->RT_OH_dissociation_rate, my_fields->RT_H2O_dissociation_rate
               );


          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask_nr)

           FORTRAN_NAME(step_rate_newton_raphson)(&my_chemistry->with_radiative_cooling, d.data(), e.data(), my_fields->x_velocity, my_fields->y_velocity, my_fields->z_velocity, de.data(), HI.data(),
                    HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                    &my_units->comoving_coordinates, &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->metal_cooling, &my_chemistry->h2_on_dust,
                    &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_fields->grid_start[0], &my_fields->grid_end[0], &my_chemistry->ih2co, &my_chemistry->ipiht,
                    &my_chemistry->dust_recombination_cooling, &my_chemistry->photoelectric_heating, &my_units->a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &utem,
                    &uxyz, &my_units->a_units, &urho, &my_units->time_units, &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
                    my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9,
                    my_rates->k10, my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
                    my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22, &my_uvb_rates->k24, &my_uvb_rates->k25, &my_uvb_rates->k26, &my_uvb_rates->k27, &my_uvb_rates->k28,
                    &my_uvb_rates->k29, &my_uvb_rates->k30, &my_uvb_rates->k31, my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54,
                    my_rates->k55, my_rates->k56, my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart,
                    &my_chemistry->DustTemperatureEnd, my_rates->h2dust, my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2, my_rates->ceHI,
                    my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI, my_rates->ciHeIS, my_rates->ciHeII,
                    my_rates->reHII, my_rates->reHeII1, my_rates->reHeII2, my_rates->reHeIII, my_rates->brem,
                    &my_rates->comp, &my_rates->gammah, &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf,
                    &my_uvb_rates->comp_xray, &my_uvb_rates->temp_xray, &my_uvb_rates->piHI, &my_uvb_rates->piHeI, &my_uvb_rates->piHeII, my_fields->HM_density,
                    my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density, my_rates->hyd01k,
                    my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl, my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit,
                    my_rates->HDlte, my_rates->HDlow, my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp,
                    my_rates->GAel, my_rates->H2LTE, my_rates->gas_grain, &my_chemistry->H2_self_shielding, &my_chemistry->self_shielding_method,
                    &my_uvb_rates->crsHI, &my_uvb_rates->crsHeI, &my_uvb_rates->crsHeII, &my_chemistry->use_radiative_transfer,
                    &my_chemistry->radiative_transfer_hydrogen_only, kphHI.data(), my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate, my_fields->RT_H2_dissociation_rate,
                    my_fields->RT_heating_rate, my_fields->H2_self_shielding_length, &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling,
                    &my_chemistry->three_body_rate, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, &my_chemistry->cmb_temperature_floor,
                    &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor, &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
                    my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
                    &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
                    &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension, my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1],
                    my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4], &my_rates->cloudy_metal.data_size,
                    my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new, &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate,
                    my_fields->volumetric_heating_rate, my_fields->specific_heating_rate, &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, my_fields->temperature_floor,
                    &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature, my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                    my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density, my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density, my_fields->SiI_density, my_fields->SiOI_density,
                    my_fields->SiO2I_density, my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density, my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density,
                    my_fields->O2II_density, my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density, my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density,
                    my_fields->Fe3O4_dust_density, my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density, my_fields->ref_org_dust_density,
                    my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density, my_rates->k125, my_rates->k129, my_rates->k130, my_rates->k131,
                    my_rates->k132, my_rates->k133, my_rates->k134, my_rates->k135, my_rates->k136, my_rates->k137, my_rates->k148,
                    my_rates->k149, my_rates->k150, my_rates->k151, my_rates->k152, my_rates->k153, my_rates->kz15, my_rates->kz16,
                    my_rates->kz17, my_rates->kz18, my_rates->kz19, my_rates->kz20, my_rates->kz21, my_rates->kz22, my_rates->kz23,
                    my_rates->kz24, my_rates->kz25, my_rates->kz26, my_rates->kz27, my_rates->kz28, my_rates->kz29, my_rates->kz30,
                    my_rates->kz31, my_rates->kz32, my_rates->kz33, my_rates->kz34, my_rates->kz35, my_rates->kz36, my_rates->kz37,
                    my_rates->kz38, my_rates->kz39, my_rates->kz40, my_rates->kz41, my_rates->kz42, my_rates->kz43, my_rates->kz44,
                    my_rates->kz45, my_rates->kz46, my_rates->kz47, my_rates->kz48, my_rates->kz49, my_rates->kz50, my_rates->kz51,
                    my_rates->kz52, my_rates->kz53, my_rates->kz54, my_rates->cieY06, my_rates->LH2.props.dimension, &my_rates->LH2.props.data_size,
                    my_rates->LH2.props.parameters[0], my_rates->LH2.props.parameters[1], my_rates->LH2.props.parameters[2], &my_rates->LH2.props.parameter_spacing[0], &my_rates->LH2.props.parameter_spacing[1], &my_rates->LH2.props.parameter_spacing[2],
                    my_rates->LH2.data, my_rates->LHD.props.dimension, &my_rates->LHD.props.data_size, my_rates->LHD.props.parameters[0], my_rates->LHD.props.parameters[1], my_rates->LHD.props.parameters[2],
                    &my_rates->LHD.props.parameter_spacing[0], &my_rates->LHD.props.parameter_spacing[1], &my_rates->LHD.props.parameter_spacing[2], my_rates->LHD.data, my_rates->LCI.props.dimension, &my_rates->LCI.props.data_size,
                    my_rates->LCI.props.parameters[0], my_rates->LCI.props.parameters[1], my_rates->LCI.props.parameters[2], &my_rates->LCI.props.parameter_spacing[0], &my_rates->LCI.props.parameter_spacing[1], &my_rates->LCI.props.parameter_spacing[2],
                    my_rates->LCI.data, my_rates->LCII.props.dimension, &my_rates->LCII.props.data_size, my_rates->LCII.props.parameters[0], my_rates->LCII.props.parameters[1], my_rates->LCII.props.parameters[2],
                    &my_rates->LCII.props.parameter_spacing[0], &my_rates->LCII.props.parameter_spacing[1], &my_rates->LCII.props.parameter_spacing[2], my_rates->LCII.data, my_rates->LOI.props.dimension,
                    &my_rates->LOI.props.data_size, my_rates->LOI.props.parameters[0], my_rates->LOI.props.parameters[1], my_rates->LOI.props.parameters[2], &my_rates->LOI.props.parameter_spacing[0], &my_rates->LOI.props.parameter_spacing[1],
                    &my_rates->LOI.props.parameter_spacing[2], my_rates->LOI.data, my_rates->LCO.props.dimension, &my_rates->LCO.props.data_size, my_rates->LCO.props.parameters[0], my_rates->LCO.props.parameters[1],
                    my_rates->LCO.props.parameters[2], &my_rates->LCO.props.parameter_spacing[0], &my_rates->LCO.props.parameter_spacing[1], &my_rates->LCO.props.parameter_spacing[2], my_rates->LCO.data, my_rates->LOH.props.dimension,
                    &my_rates->LOH.props.data_size, my_rates->LOH.props.parameters[0], my_rates->LOH.props.parameters[1], my_rates->LOH.props.parameters[2], &my_rates->LOH.props.parameter_spacing[0], &my_rates->LOH.props.parameter_spacing[1],
                    &my_rates->LOH.props.parameter_spacing[2], my_rates->LOH.data, my_rates->LH2O.props.dimension, &my_rates->LH2O.props.data_size, my_rates->LH2O.props.parameters[0], my_rates->LH2O.props.parameters[1],
                    my_rates->LH2O.props.parameters[2], &my_rates->LH2O.props.parameter_spacing[0], &my_rates->LH2O.props.parameter_spacing[1], &my_rates->LH2O.props.parameter_spacing[2], my_rates->LH2O.data,
                    my_rates->alphap.props.dimension, &my_rates->alphap.props.data_size, my_rates->alphap.props.parameters[0], my_rates->alphap.props.parameters[1],
                    &my_rates->alphap.props.parameter_spacing[0], &my_rates->alphap.props.parameter_spacing[1], my_rates->alphap.data, &my_chemistry->multi_metals,
                    &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation, my_fields->local_ISM_metal_density,
                    my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                    my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                    my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density, &my_rates->SN0_N,
                    my_rates->SN0_fSiM, my_rates->SN0_fFeM, my_rates->SN0_fMg2SiO4, my_rates->SN0_fMgSiO3,
                    my_rates->SN0_fFe3O4, my_rates->SN0_fAC, my_rates->SN0_fSiO2D, my_rates->SN0_fMgO,
                    my_rates->SN0_fFeS, my_rates->SN0_fAl2O3, my_rates->SN0_freforg, my_rates->SN0_fvolorg,
                    my_rates->SN0_fH2Oice, my_rates->SN0_r0SiM, my_rates->SN0_r0FeM, my_rates->SN0_r0Mg2SiO4,
                    my_rates->SN0_r0MgSiO3, my_rates->SN0_r0Fe3O4, my_rates->SN0_r0AC, my_rates->SN0_r0SiO2D,
                    my_rates->SN0_r0MgO, my_rates->SN0_r0FeS, my_rates->SN0_r0Al2O3, my_rates->SN0_r0reforg,
                    my_rates->SN0_r0volorg, my_rates->SN0_r0H2Oice, my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT,
                    my_rates->gr_Td, my_rates->SN0_kpSiM, my_rates->SN0_kpFeM, my_rates->SN0_kpMg2SiO4,
                    my_rates->SN0_kpMgSiO3, my_rates->SN0_kpFe3O4, my_rates->SN0_kpAC, my_rates->SN0_kpSiO2D,
                    my_rates->SN0_kpMgO, my_rates->SN0_kpFeS, my_rates->SN0_kpAl2O3, my_rates->SN0_kpreforg,
                    my_rates->SN0_kpvolorg, my_rates->SN0_kpH2Oice, my_rates->h2dustS, my_rates->h2dustC,
                    my_rates->gas_grain2, &my_rates->gamma_isrf2, my_rates->grain_growth_rate, &my_chemistry->radiative_transfer_HDI_dissociation,
                    my_fields->RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, my_fields->RT_CI_ionization_rate, my_fields->RT_OI_ionization_rate, &my_chemistry->radiative_transfer_metal_dissociation, my_fields->RT_CO_dissociation_rate,
                    my_fields->RT_OH_dissociation_rate, my_fields->RT_H2O_dissociation_rate, &my_chemistry->radiative_transfer_use_H2_shielding, &my_chemistry->use_isrf_field,
                    my_fields->isrf_habing, &my_chemistry->H2_custom_shielding, my_fields->H2_custom_shielding_factor,
                    &j, &k, &iter, &dom, &comp1,
                    &comp2, &coolunit, &tbase1, &xbase1, &chunit, &dx_cgs,
                    &c_ljeans, logTlininterp_buf.indixe, logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.logtem, logTlininterp_buf.tdef, dtit.data(),
                    p2d.data(), tgas.data(), cool1dmulti_buf.tgasold, tdust.data(), metallicity.data(), dust2gas.data(),
                    rhoH.data(), mmw.data(), cool1dmulti_buf.mynh, cool1dmulti_buf.myde, cool1dmulti_buf.gammaha_eff, cool1dmulti_buf.gasgr_tdust,
                    cool1dmulti_buf.regr, h2dust.data(), ncrn.data(), ncrd1.data(), ncrd2.data(), grain_temperatures.SiM, grain_temperatures.FeM,
                    grain_temperatures.Mg2SiO4, grain_temperatures.MgSiO3, grain_temperatures.Fe3O4, grain_temperatures.AC, grain_temperatures.SiO2D, grain_temperatures.MgO,
                    grain_temperatures.FeS, grain_temperatures.Al2O3, grain_temperatures.reforg, grain_temperatures.volorg, grain_temperatures.H2Oice, ceHI.data(),
                    ceHeI.data(), ceHeII.data(), ciHI.data(), ciHeI.data(), ciHeIS.data(), ciHeII.data(),
                    reHII.data(), reHeII1.data(), reHeII2.data(), reHeIII.data(), brem.data(), edot.data(),
                    hyd01k.data(), h2k01.data(), vibh.data(), roth.data(), rotl.data(), gpldl.data(), gphdl.data(),
                    hdlte.data(), hdlow.data(), cieco.data(), &anydust, itmask_nr.data(),
                    itmask_metal.data(), imp_eng.data()
                  );

        }

        // return itmask
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          itmask[i-1] = itmask_tmp[i-1];
        }

        // Add the timestep to the elapsed time for each cell and find
        //  minimum elapsed time step in this row

        ttmin = huge8;
        for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
          ttot[i-1] = std::fmin(ttot[i-1] + dtit[i-1], dt);
          if (std::fabs(dt-ttot[i-1]) <
               tolerance*dt) { itmask[i-1] = MASK_FALSE; }
          if (ttot[i-1]<ttmin) { ttmin = ttot[i-1]; }
        }

        // If all cells are done (on this slice), break out of subcycle loop
        if (std::fabs(dt-ttmin) < tolerance*dt) { break; }

      }  // subcycle iteration loop (for current row)

      // review number of iterations that were spent in the subcycle loop

      if (iter > my_chemistry->max_iterations)  {
        OMP_PRAGMA_CRITICAL
        {
          printf("inside if statement solve rate cool: %d %d\n",
                 my_fields->grid_start[0],
                 my_fields->grid_end[0]);
          eprintf("MULTI_COOL iter >  %d  at j,k = %d %d\n",
                  my_chemistry->max_iterations,
                  j,
                  k);
          printf("FATAL error (2) in MULTI_COOL\n");
          printf("( dt = %.17e ttmin = %.17e )", dt, ttmin);
          grackle::impl::print_contiguous_row_(
            dtit.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            ttot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            edot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            itmask.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );

          if (my_chemistry->exit_after_iterations_exceeded == 1)  {
            ierr = GR_FAIL;
          }
        }  // OMP_PRAGMA_CRITICAL
      }

      if (iter > my_chemistry->max_iterations/2) { // WARNING_MESSAGE
        OMP_PRAGMA_CRITICAL
        {
          eprintf("MULTI_COOL iter,j,k = %d %d %d\n", iter, j, k);
        }
      }

    }  // loop over j,k pairs

    // cleanup manually allocated temporaries
    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_Cool1DMultiScratchBuf(&cool1dmulti_buf);

  }  // OMP_PRAGMA("omp parallel")

  // If an error has been produced, return now.

  if (ierr != GR_SUCCESS)  {
    return ierr;
  }

  // Convert densities back to comoving from proper

  if (my_units->comoving_coordinates == 1)  {

    factor = (gr_float)(std::pow(my_units->a_value,3) );

     FORTRAN_NAME(scale_fields_g)(
      &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry,
      &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->multi_metals, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation, &factor,
      &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2], &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
      d.data(), de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
      my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
      my_fields->metal_density, my_fields->dust_density,
      my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
      my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
      my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
      my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
      my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
      my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
      my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
      my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
      my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
      my_fields->local_ISM_metal_density, my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
      my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
      my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);

  }

  if (my_chemistry->primordial_chemistry > 0)  {

    // Correct the species to ensure consistency (i.e. type conservation)

#ifdef ABUNDANCE_CORRECTION
     FORTRAN_NAME(make_consistent_g)(de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                         d.data(), &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->HydrogenFractionByMass, &my_chemistry->DeuteriumToHydrogenRatio,
                        &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &dom,
                        my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                        my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                        my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                        my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                        my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                        my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                        my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                        my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                        my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                        my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                        &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
                        my_fields->local_ISM_metal_density,
                        my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                        my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                        my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density,
                        &my_rates->SN0_N,
                        my_rates->SN0_XC, my_rates->SN0_XO, my_rates->SN0_XMg, my_rates->SN0_XAl, my_rates->SN0_XSi,
                        my_rates->SN0_XS, my_rates->SN0_XFe,
                        my_rates->SN0_fC, my_rates->SN0_fO, my_rates->SN0_fMg, my_rates->SN0_fAl, my_rates->SN0_fSi,
                        my_rates->SN0_fS, my_rates->SN0_fFe
                          );

     FORTRAN_NAME(ceiling_species_g)(d.data(), de.data(), HI.data(), HII.data(), my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                         &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->use_dust_density_field,
                        my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                        &my_chemistry->metal_abundances, &my_chemistry->metal_chemistry, &my_chemistry->dust_species, &my_chemistry->multi_metals,
                        &my_chemistry->grain_growth, &my_chemistry->dust_sublimation,
                        my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                        my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                        my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                        my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                        my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                        my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                        my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                        my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                        my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                        my_fields->local_ISM_metal_density,
                        my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                        my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                        my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);
#endif

  }

  return ierr;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
