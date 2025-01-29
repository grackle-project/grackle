// See LICENSE file for license and copyright information

/// @file step_rate_newton_raphson.hpp
/// @brief Defines the step_rate_newton_raphson function

// This file was initially generated automatically during conversion of the
// step_rate_newton_raphson function from FORTRAN to C++

#ifndef STEP_RATE_NEWTON_RAPHSON_HPP
#define STEP_RATE_NEWTON_RAPHSON_HPP

#include <vector>

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "utils-cpp.hpp"

namespace grackle::impl {

/// An alternative to step_rate_g for evolving the species rate equations that
/// employs the Newton-Raphson scheme rather than Gauss-Seidel
///
/// Whereas the Gauss-Seidel scheme always operator-splits evolution of the
/// energy equation, this function internally evolves the energy equation
///
/// @note
/// The values in the various buffers (p2d, tgas, tdust, metallicity, dust2gas,
/// rhoH, mmw, h2dust, edot, grain_temperatures, logTlininterp_buf,
/// cool1dmulti_buf, coolingheating_buf, chemheatrates_buf), hold undefined
/// values after this function call wherever `itmask_nr` indicates that this
/// function runs.
/// - this has **ALWAYS** been the case. Historically these buffers have been
///   reused when computing finite differences
inline void step_rate_newton_raphson(
  int imetal, IndexRange idx_range, int iter, double dom, double chunit,
  double dx_cgs, double c_ljeans, double* dtit, double* p2d, double* tgas,
  double* tdust, double* metallicity, double* dust2gas, double* rhoH,
  double* mmw, double* h2dust, double* edot, gr_mask_type anydust,
  gr_mask_type* itmask_nr, gr_mask_type* itmask_metal, int* imp_eng,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  InternalGrUnits internalu,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
)
{

  // Density, energy and velocity fields fields

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
  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> u(my_fields->x_velocity, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> v(my_fields->y_velocity, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> w(my_fields->z_velocity, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(my_fields->metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(my_fields->dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Vheat(my_fields->volumetric_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mheat(my_fields->specific_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Tfloor(my_fields->temperature_floor, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
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
  grackle::impl::View<gr_float***> metal_loc(my_fields->local_ISM_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C13(my_fields->ccsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C20(my_fields->ccsn20_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C25(my_fields->ccsn25_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C30(my_fields->ccsn30_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F13(my_fields->fsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F15(my_fields->fsn15_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F50(my_fields->fsn50_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F80(my_fields->fsn80_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P170(my_fields->pisn170_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P200(my_fields->pisn200_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_Y19(my_fields->y19_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Radiative transfer fields

  grackle::impl::View<gr_float***> kphHI(my_fields->RT_HI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeI(my_fields->RT_HeI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeII(my_fields->RT_HeII_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissH2I(my_fields->RT_H2_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> photogamma(my_fields->RT_heating_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> kdissHDI(my_fields->RT_HDI_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphCI(my_fields->RT_CI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphOI(my_fields->RT_OI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissCO(my_fields->RT_CO_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissOH(my_fields->RT_OH_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissH2O(my_fields->RT_H2O_dissociation_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // H2 self-shielding length-scale field

  grackle::impl::View<gr_float***> xH2shield(my_fields->H2_self_shielding_length, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Interstellar radiation field for dust heating

  grackle::impl::View<gr_float***> isrf_habing(my_fields->isrf_habing, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Custom H2 shielding factor

  grackle::impl::View<gr_float***> f_shield_custom(my_fields->H2_custom_shielding_factor, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // ierror local variable
  //   - this variable is only used internally by this subroutine for
  //     determining control-flow
  //   - it does NOT report whether or not this function succeeded (maybe
  //     we should change the name?)
  int ierror;
  // Local variable
  int itr, itr_time;
  int nsp, isp, jsp, id;
  double dspj, err, err_max;
  const int i_eng = 52;
  // There may be an argument for allocating the following at a higher
  // level function, but we will leave that for after transcription
  std::vector<double> dsp(i_eng);
  std::vector<double> dsp0(i_eng);
  std::vector<double> dsp1(i_eng);
  std::vector<double> dspdot(i_eng);
  std::vector<double> dspdot1(i_eng);
  std::vector<double> ddsp(i_eng);
  std::vector<double> der_data_(i_eng * i_eng);
  grackle::impl::View<double**> der(der_data_.data(), i_eng, i_eng);

  // (In the future, we may want to reconsider when/how we allocate
  // the following 3 variables)
  std::vector<int> idsp;
  std::vector<double> mtrx_data_;
  grackle::impl::View<double**> mtrx;
  std::vector<double> vec;
  // Another parameter
  const double eps = 1.e-4;

  // The following was extracted from another subroutine
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    const int ip1 = i+1;
    if (itmask_nr[ip1-1] != MASK_FALSE)  {

      // If density and temperature are low, update gas energy explicitly

      if (my_chemistry->with_radiative_cooling == 1)  {
        if (imp_eng[ip1-1] == 0)  {
          e(ip1-1,idx_range.jp1-1,idx_range.kp1-1)  = e(ip1-1,idx_range.jp1-1,idx_range.kp1-1) +
               (gr_float)(edot[ip1-1]/d(ip1-1,idx_range.jp1-1,idx_range.kp1-1)*dtit[ip1-1] );
        }
      }

      // initialize arrays
      if (my_chemistry->primordial_chemistry > 0) { nsp = 6; }
      if (my_chemistry->primordial_chemistry > 1) { nsp = nsp + 3; }
      if (my_chemistry->primordial_chemistry > 2) { nsp = nsp + 3; }
      if (my_chemistry->primordial_chemistry > 3) { nsp = nsp + 3; }
      if (itmask_metal[ip1-1] != MASK_FALSE)  {
        if (my_chemistry->metal_chemistry == 1)  {
          nsp = nsp + 19;
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0) { nsp = nsp + 1; }
            if (my_chemistry->dust_species > 1) { nsp = nsp + 3; }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0) { nsp = nsp + 2; }
          if (my_chemistry->dust_species > 1) { nsp = nsp + 8; }
          if (my_chemistry->dust_species > 2) { nsp = nsp + 3; }
        }
      }
      nsp = nsp + imp_eng[ip1-1];
      idsp.reserve(nsp);
      mtrx_data_.reserve(nsp * nsp);
      mtrx = grackle::impl::View<double**>(mtrx_data_.data(), nsp, nsp);
      vec.reserve(nsp);

      if ( my_chemistry->primordial_chemistry > 0 )  {
        dsp[ 1-1] = de(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 2-1] = HI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 3-1] = HII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 4-1] = HeI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 5-1] = HeII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 6-1] = HeIII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        dsp[ 7-1] = HM(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 8-1] = H2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[ 9-1] = H2II(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        dsp[10-1] = DI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[11-1] = DII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[12-1] = HDI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        dsp[13-1] = DM(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[14-1] = HDII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
        dsp[15-1] = HeHII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
      }
      if ( itmask_metal[ip1-1] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          dsp[16-1] = CI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[17-1] = CII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[18-1] = CO(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[19-1] = CO2(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[20-1] = OI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[21-1] = OH(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[22-1] = H2O(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[23-1] = O2(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[24-1] = SiI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[25-1] = SiOI(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[26-1] = SiO2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[27-1] = CH(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[28-1] = CH2(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[29-1] = COII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[30-1] = OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[31-1] = OHII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[32-1] = H2OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[33-1] = H3OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          dsp[34-1] = O2II(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0)  {
              dsp[35-1] = Mg(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            }
            if (my_chemistry->dust_species > 1)  {
              dsp[36-1] = Al(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
              dsp[37-1] = S(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
              dsp[38-1] = Fe(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            dsp[39-1] = MgSiO3(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[40-1] = AC(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          }
          if (my_chemistry->dust_species > 1)  {
            dsp[41-1] = SiM(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[42-1] = FeM(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[43-1] = Mg2SiO4(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[44-1] = Fe3O4(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[45-1] = SiO2D(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[46-1] = MgO(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[47-1] = FeS(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[48-1] = Al2O3(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          }
          if (my_chemistry->dust_species > 2)  {
            dsp[49-1] = reforg(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[50-1] = volorg(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
            dsp[51-1] = H2Oice(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          }
        }
      }
      dsp[i_eng-1] = e(ip1-1,idx_range.jp1-1,idx_range.kp1-1);

      id = 0;
      if (my_chemistry->primordial_chemistry > 0)  {
        for (isp = 1; isp<=(6); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 1)  {
        for (isp = 7; isp<=(9); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 2)  {
        for (isp = 10; isp<=(12); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 3)  {
        for (isp = 13; isp<=(15); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (itmask_metal[ip1-1] != MASK_FALSE)  {
        if (my_chemistry->metal_chemistry == 1)  {
          for (isp = 16; isp<=(34); isp++) {
            id = id + 1;
            idsp[id-1] = isp;
          }
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0)  {
              for (isp = 35; isp<=(35); isp++) {
                id = id + 1;
                idsp[id-1] = isp;
              }
            }
            if (my_chemistry->dust_species > 1)  {
              for (isp = 36; isp<=(38); isp++) {
                id = id + 1;
                idsp[id-1] = isp;
              }
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            for (isp = 39; isp<=(40); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
          if (my_chemistry->dust_species > 1)  {
            for (isp = 41; isp<=(48); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
          if (my_chemistry->dust_species > 2)  {
            for (isp = 49; isp<=(51); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
        }
      }
      if ( imp_eng[ip1-1] ==1 )  {
        id = id + 1;
        idsp[id-1] = i_eng;
      }

      // Save arrays at ttot(ip1)

      std::memcpy(dsp0.data(), dsp.data(), sizeof(double)*i_eng);
      std::memset(ddsp.data(), 0, sizeof(double)*i_eng);

      // Search for the timestep for which chemistry converges

      ierror=1;
      itr_time=0;
      while (ierror==1) {

        // If not converge, restore arrays at ttot(ip1)

        std::memcpy(dsp.data(), dsp0.data(), sizeof(double)*i_eng);
        std::memset(ddsp.data(), 0, sizeof(double)*i_eng);

        // Iteration to solve ODEs

        err_max=1.e2;
        itr=0;
        while (err_max>1.e-8) {
          if(itr>=20)  {
            ierror = 1;
            goto label_9996;
          }

           FORTRAN_NAME(lookup_cool_rates0d)(&dtit[ip1-1],
          &d(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &u(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &v(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &w(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          &nsp, dsp.data(), dspdot.data(), &my_chemistry->NumberOfTemperatureBins,
          &internalu.extfields_in_comoving, &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->metal_cooling,
          &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
          &my_chemistry->ih2co, &my_chemistry->ipiht, &iter, &my_chemistry->photoelectric_heating,
          &internalu.a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
          &internalu.utem, &internalu.uxyz, &internalu.a_units, &internalu.urho, &internalu.tbase1,
          &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
          my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
          my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
          my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
          &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates.comp_xray, &my_uvb_rates.temp_xray,
          &my_uvb_rates.piHI, &my_uvb_rates.piHeI, &my_uvb_rates.piHeII,
          &metal(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &dust(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
          &coolingheating_buf.hyd01k[ip1-1], &coolingheating_buf.h2k01[ip1-1], &coolingheating_buf.vibh[ip1-1], &coolingheating_buf.roth[ip1-1], &coolingheating_buf.rotl[ip1-1],
          my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, &coolingheating_buf.gpldl[ip1-1], &coolingheating_buf.gphdl[ip1-1],
          my_rates->HDlte, my_rates->HDlow, &coolingheating_buf.hdlte[ip1-1], &coolingheating_buf.hdlow[ip1-1],
          my_rates->GAHI, my_rates->GAH2, my_rates->GAHe,
          my_rates->GAHp, my_rates->GAel,
          my_rates->H2LTE, my_rates->gas_grain,
          &coolingheating_buf.ceHI[ip1-1], &coolingheating_buf.ceHeI[ip1-1], &coolingheating_buf.ceHeII[ip1-1], &coolingheating_buf.ciHI[ip1-1], &coolingheating_buf.ciHeI[ip1-1],
          &coolingheating_buf.ciHeIS[ip1-1], &coolingheating_buf.ciHeII[ip1-1],
          &coolingheating_buf.reHII[ip1-1], &coolingheating_buf.reHeII1[ip1-1], &coolingheating_buf.reHeII2[ip1-1], &coolingheating_buf.reHeIII[ip1-1], &coolingheating_buf.brem[ip1-1],
          &logTlininterp_buf.indixe[ip1-1], &logTlininterp_buf.t1[ip1-1], &logTlininterp_buf.t2[ip1-1], &logTlininterp_buf.logtem[ip1-1], &logTlininterp_buf.tdef[ip1-1], &edot[ip1-1],
          &tgas[ip1-1], &cool1dmulti_buf.tgasold[ip1-1], &mmw[ip1-1], &p2d[ip1-1], &tdust[ip1-1], &metallicity[ip1-1],
          &dust2gas[ip1-1], &rhoH[ip1-1], &cool1dmulti_buf.mynh[ip1-1], &cool1dmulti_buf.myde[ip1-1],
          &cool1dmulti_buf.gammaha_eff[ip1-1], &cool1dmulti_buf.gasgr_tdust[ip1-1], &cool1dmulti_buf.regr[ip1-1],
          &my_chemistry->self_shielding_method, &my_uvb_rates.crsHI, &my_uvb_rates.crsHeI, &my_uvb_rates.crsHeII,
          &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
          &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, &coolingheating_buf.cieco[ip1-1],
          &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
          &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
          my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
          &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
          &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
          my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
          &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
          &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, &Vheat(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &Mheat(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, &Tfloor(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          &my_chemistry->use_isrf_field, &isrf_habing(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          &my_chemistry->three_body_rate, &anydust, &my_chemistry->H2_self_shielding,
          my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
          my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
          my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
          &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30, &my_uvb_rates.k31,
          my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54, my_rates->k55, my_rates->k56,
          my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
          my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
          &h2dust[ip1-1], &chemheatrates_buf.n_cr_n[ip1-1], &chemheatrates_buf.n_cr_d1[ip1-1], &chemheatrates_buf.n_cr_d2[ip1-1],
          &dom, &internalu.coolunit, &internalu.tbase1, &internalu.xbase1, &dx_cgs, &c_ljeans,
          &kphHI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphHeI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphHeII(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissH2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
          &photogamma(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &xH2shield(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &chunit, &itmask_nr[ip1-1],
          &itmask_metal[ip1-1],
           &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature,
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
           &metal_loc(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C13(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C20(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
           &metal_C25(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C30(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F13(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
           &metal_F15(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F50(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F80(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
           &metal_P170(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_P200(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_Y19(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
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
           my_rates->h2dustS, my_rates->h2dustC, my_rates->grain_growth_rate,
           &grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][ip1-1],
           &grain_temperatures.data[OnlyGrainSpLUT::AC_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][ip1-1],
           &grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][ip1-1],
           my_rates->gas_grain2, &my_rates->gamma_isrf2,
           &imp_eng[ip1-1],
           &my_chemistry->radiative_transfer_HDI_dissociation, &kdissHDI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &my_chemistry->radiative_transfer_metal_ionization, &kphCI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphOI(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
           &my_chemistry->radiative_transfer_metal_dissociation, &kdissCO(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissOH(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissH2O(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
           &my_chemistry->radiative_transfer_use_H2_shielding, &my_chemistry->H2_custom_shielding, &f_shield_custom(ip1-1,idx_range.jp1-1,idx_range.kp1-1)
          );

          for (jsp = 1; jsp<=(nsp); jsp++) {
            dspj = eps * dsp[idsp[jsp-1]-1];
            for (isp = 1; isp<=(nsp); isp++) {
              if(isp == jsp)  {
                dsp1[idsp[isp-1]-1] = dsp[idsp[isp-1]-1] + dspj;
              } else {
                dsp1[idsp[isp-1]-1] = dsp[idsp[isp-1]-1];
              }
            }

             FORTRAN_NAME(lookup_cool_rates0d)(&dtit[ip1-1],
            &d(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &u(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &v(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &w(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            &nsp, dsp1.data(), dspdot1.data(), &my_chemistry->NumberOfTemperatureBins,
            &internalu.extfields_in_comoving, &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->metal_cooling,
            &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
            &my_chemistry->ih2co, &my_chemistry->ipiht, &iter, &my_chemistry->photoelectric_heating,
            &internalu.a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
            &internalu.utem, &internalu.uxyz, &internalu.a_units, &internalu.urho, &internalu.tbase1,
            &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
            my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
            my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
            my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
            &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates.comp_xray, &my_uvb_rates.temp_xray,
            &my_uvb_rates.piHI, &my_uvb_rates.piHeI, &my_uvb_rates.piHeII,
            &metal(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &dust(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
            &coolingheating_buf.hyd01k[ip1-1], &coolingheating_buf.h2k01[ip1-1], &coolingheating_buf.vibh[ip1-1], &coolingheating_buf.roth[ip1-1], &coolingheating_buf.rotl[ip1-1],
            my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, &coolingheating_buf.gpldl[ip1-1], &coolingheating_buf.gphdl[ip1-1],
            my_rates->HDlte, my_rates->HDlow, &coolingheating_buf.hdlte[ip1-1], &coolingheating_buf.hdlow[ip1-1],
            my_rates->GAHI, my_rates->GAH2, my_rates->GAHe,
            my_rates->GAHp, my_rates->GAel,
            my_rates->H2LTE, my_rates->gas_grain,
            &coolingheating_buf.ceHI[ip1-1], &coolingheating_buf.ceHeI[ip1-1], &coolingheating_buf.ceHeII[ip1-1], &coolingheating_buf.ciHI[ip1-1], &coolingheating_buf.ciHeI[ip1-1],
            &coolingheating_buf.ciHeIS[ip1-1], &coolingheating_buf.ciHeII[ip1-1],
            &coolingheating_buf.reHII[ip1-1], &coolingheating_buf.reHeII1[ip1-1], &coolingheating_buf.reHeII2[ip1-1], &coolingheating_buf.reHeIII[ip1-1], &coolingheating_buf.brem[ip1-1],
            &logTlininterp_buf.indixe[ip1-1], &logTlininterp_buf.t1[ip1-1], &logTlininterp_buf.t2[ip1-1], &logTlininterp_buf.logtem[ip1-1], &logTlininterp_buf.tdef[ip1-1], &edot[ip1-1],
            &tgas[ip1-1], &cool1dmulti_buf.tgasold[ip1-1], &mmw[ip1-1], &p2d[ip1-1], &tdust[ip1-1], &metallicity[ip1-1],
            &dust2gas[ip1-1], &rhoH[ip1-1], &cool1dmulti_buf.mynh[ip1-1], &cool1dmulti_buf.myde[ip1-1],
            &cool1dmulti_buf.gammaha_eff[ip1-1], &cool1dmulti_buf.gasgr_tdust[ip1-1], &cool1dmulti_buf.regr[ip1-1],
            &my_chemistry->self_shielding_method, &my_uvb_rates.crsHI, &my_uvb_rates.crsHeI, &my_uvb_rates.crsHeII,
            &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
            &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, &coolingheating_buf.cieco[ip1-1],
            &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
            &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
            my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
            &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
            &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
            my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
            &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
            &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, &Vheat(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &Mheat(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, &Tfloor(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            &my_chemistry->use_isrf_field, &isrf_habing(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            &my_chemistry->three_body_rate, &anydust, &my_chemistry->H2_self_shielding,
            my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
            my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
            my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
            &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30, &my_uvb_rates.k31,
            my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54, my_rates->k55, my_rates->k56,
            my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
            my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
            &h2dust[ip1-1], &chemheatrates_buf.n_cr_n[ip1-1], &chemheatrates_buf.n_cr_d1[ip1-1], &chemheatrates_buf.n_cr_d2[ip1-1],
            &dom, &internalu.coolunit, &internalu.tbase1, &internalu.xbase1, &dx_cgs, &c_ljeans,
            &kphHI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphHeI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphHeII(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissH2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
            &photogamma(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &xH2shield(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &chunit, &itmask_nr[ip1-1],
            &itmask_metal[ip1-1],
             &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature,
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
             &metal_loc(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C13(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C20(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
             &metal_C25(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_C30(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F13(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
             &metal_F15(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F50(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_F80(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
             &metal_P170(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_P200(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &metal_Y19(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
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
             my_rates->h2dustS, my_rates->h2dustC, my_rates->grain_growth_rate,
             &grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][ip1-1],
             &grain_temperatures.data[OnlyGrainSpLUT::AC_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][ip1-1],
             &grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][ip1-1], &grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][ip1-1],
             my_rates->gas_grain2, &my_rates->gamma_isrf2,
             &imp_eng[ip1-1],
             &my_chemistry->radiative_transfer_HDI_dissociation, &kdissHDI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &my_chemistry->radiative_transfer_metal_ionization, &kphCI(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kphOI(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
             &my_chemistry->radiative_transfer_metal_dissociation, &kdissCO(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissOH(ip1-1,idx_range.jp1-1,idx_range.kp1-1), &kdissH2O(ip1-1,idx_range.jp1-1,idx_range.kp1-1),
             &my_chemistry->radiative_transfer_use_H2_shielding, &my_chemistry->H2_custom_shielding, &f_shield_custom(ip1-1,idx_range.jp1-1,idx_range.kp1-1)
            );

            for (isp = 1; isp<=(nsp); isp++) {
              if ( (dsp[idsp[isp-1]-1]==0.e0)
               &&  (dspdot1[idsp[isp-1]-1]
               ==  dspdot[idsp[isp-1]-1]) )  {
                der(idsp[isp-1]-1,idsp[jsp-1]-1) = 0.e0;
              } else {
                der(idsp[isp-1]-1,idsp[jsp-1]-1) =
                   (dspdot1[idsp[isp-1]-1]
                   - dspdot[idsp[isp-1]-1]) / dspj;
              }
            }

          }

          for (isp = 1; isp<=(nsp); isp++) {
            for (jsp = 1; jsp<=(nsp); jsp++) {
              if(isp == jsp)  {
                mtrx(isp-1,jsp-1) = 1.e0 - dtit[ip1-1]
                   * der(idsp[isp-1]-1,idsp[jsp-1]-1);
              } else {
                mtrx(isp-1,jsp-1) =      - dtit[ip1-1]
                   * der(idsp[isp-1]-1,idsp[jsp-1]-1);
              }
            }
          }

          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = dspdot[idsp[isp-1]-1] * dtit[ip1-1]
                   - ddsp[idsp[isp-1]-1];
          }

          // to get more accuracy
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]/d(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          }

           FORTRAN_NAME(gaussj_g)(&nsp, mtrx.data(), vec.data(), &ierror);
          if(ierror == 1)  {
            goto label_9998;
          }

          // multiply with density again
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]*d(ip1-1,idx_range.jp1-1,idx_range.kp1-1);
          }

          for (isp = 1; isp<=(nsp); isp++) {
            ddsp[idsp[isp-1]-1] = ddsp[idsp[isp-1]-1] + vec[isp-1];
            dsp[idsp[isp-1]-1]  = dsp[idsp[isp-1]-1]  + vec[isp-1];
          }

          if (imp_eng[ip1-1] == 1)  {
            if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
              for (isp = 1; isp<=(nsp); isp++) {
                if ( (dsp[idsp[isp-1]-1] != dsp[idsp[isp-1]-1])
                 ||  (dsp[idsp[isp-1]-1] <= 0.) )  {
                  ierror = 1;
                  goto label_9997;
                }
              }
            }
          }

          err_max = 0.e0;
          for (isp = 1; isp<=(nsp); isp++) {
            if(dsp[idsp[isp-1]-1] > tiny8)  {
              err = grackle::impl::dabs(vec[isp-1] / dsp[idsp[isp-1]-1]);
            } else {
              err = 0.e0;
            }
            if(err > err_max)  {
              err_max = err;
            }
          }

          itr=itr+1;
        }

label_9998:
label_9997:
label_9996:


        // Check if the fractions are valid after an iteration

        if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
          for (isp = 1; isp<=(nsp); isp++) {
            if ( (dsp[idsp[isp-1]-1] != dsp[idsp[isp-1]-1])
             ||  (dsp[idsp[isp-1]-1] <= 0.) )  {
              ierror = 1;
            }
          }
        }
        if(ierror == 1)  {
          dtit[ip1-1] = 0.5e0*dtit[ip1-1];
        }

        itr_time=itr_time+1;
      }

      if ( my_chemistry->primordial_chemistry > 0 )  {
        de(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[ 1-1];
        HI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[ 2-1];
        HII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[ 3-1];
        HeI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[ 4-1];
        HeII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[ 5-1];
        HeIII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[ 6-1];
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        HM(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[ 7-1];
        H2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[ 8-1];
        H2II(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[ 9-1];
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        DI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[10-1];
        DII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[11-1];
        HDI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[12-1];
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        DM(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[13-1];
        HDII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[14-1];
        HeHII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[15-1];
      }
      if ( itmask_metal[ip1-1] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          CI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[16-1];
          CII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[17-1];
          CO(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[18-1];
          CO2(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[19-1];
          OI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[20-1];
          OH(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[21-1];
          H2O(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[22-1];
          O2(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[23-1];
          SiI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[24-1];
          SiOI(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[25-1];
          SiO2I(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[26-1];
          CH(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[27-1];
          CH2(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[28-1];
          COII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[29-1];
          OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[30-1];
          OHII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[31-1];
          H2OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[32-1];
          H3OII(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[33-1];
          O2II(ip1-1,idx_range.jp1-1,idx_range.kp1-1)    = dsp[34-1];
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
            if (my_chemistry->dust_species > 0)  {
              Mg(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[35-1];
            }
            if (my_chemistry->dust_species > 1)  {
              Al(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[36-1];
              S(ip1-1,idx_range.jp1-1,idx_range.kp1-1)       = dsp[37-1];
              Fe(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[38-1];
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
          if (my_chemistry->dust_species > 0)  {
            MgSiO3(ip1-1,idx_range.jp1-1,idx_range.kp1-1)  = dsp[39-1];
            AC(ip1-1,idx_range.jp1-1,idx_range.kp1-1)      = dsp[40-1];
          }
          if (my_chemistry->dust_species > 1)  {
            SiM(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[41-1];
            FeM(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[42-1];
            Mg2SiO4(ip1-1,idx_range.jp1-1,idx_range.kp1-1) = dsp[43-1];
            Fe3O4(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[44-1];
            SiO2D(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[45-1];
            MgO(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[46-1];
            FeS(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[47-1];
            Al2O3(ip1-1,idx_range.jp1-1,idx_range.kp1-1)   = dsp[48-1];
          }
          if (my_chemistry->dust_species > 2)  {
            reforg(ip1-1,idx_range.jp1-1,idx_range.kp1-1)  = dsp[49-1];
            volorg(ip1-1,idx_range.jp1-1,idx_range.kp1-1)  = dsp[50-1];
            H2Oice(ip1-1,idx_range.jp1-1,idx_range.kp1-1)  = dsp[51-1];
          }
        }
      }

      e(ip1-1,idx_range.jp1-1,idx_range.kp1-1)     = dsp[i_eng-1];

      idsp.clear();
      vec.clear();
      mtrx = grackle::impl::View<double**>();
      mtrx_data_.clear();

    }
  }
}


} // namespace grackle::impl


#endif /* STEP_RATE_NEWTON_RAPHSON_HPP */
