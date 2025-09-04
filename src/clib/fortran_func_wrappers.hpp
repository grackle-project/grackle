// See LICENSE file for license and copyright information

/// @file fortran_func_wrappers.hpp
/// @brief Declares wrappers around routines that haven't been transcribed yet
///     from Fortran.
///
/// THIS FILE WILL BE DELETED ONCE WE COMPLETE TRANSCRIPTION
///
/// This file exists to aid the transcription process
/// - This file holds C++ functions that wrap untranscribed Fortran routines.
///   The idea is that these wrapper functions take reduced argument lists
///   (that may includes structs) before converting extended arg lists used
///   within the fortran subroutine.
/// - ideally the reduced argument list used by the wrapper function should
///   roughly approximate the final argument list that the routines will use
///   after transcription is complete.

#ifndef FORTRAN_FUNC_WRAPPERS_HPP
#define FORTRAN_FUNC_WRAPPERS_HPP

#ifndef __cplusplus
#error "This file must be read by a c++ compiler"
#endif

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "LUT.hpp"
#include "utils-cpp.hpp"

// callers of these functions are generally expected to locally shorten the
// namespace name when they call these routines
namespace grackle::impl::fortran_wrapper {

inline void calc_all_tdust_gasgr_1d_g(
  double trad, double* tgas, double* tdust, double* metallicity,
  double* dust2gas, double* nh, double* gasgr_tdust, gr_mask_type* itmask_metal,
  double coolunit, double* gasgr, double* myisrf, double* kappa_tot,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate,
  grackle::impl::GrainSpeciesCollection grain_kappa,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf
) {

  FORTRAN_NAME(calc_all_tdust_gasgr_1d_g)(
    &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
    &my_chemistry->use_dust_density_field, &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &my_chemistry->local_dust_to_gas_ratio, &my_rates->gamma_isrf,
    &trad, my_rates->gas_grain, logTlininterp_buf.indixe, logTlininterp_buf.tdef, tgas, tdust,
    metallicity, dust2gas, nh, gasgr_tdust,
    itmask_metal,
    &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT,
    my_rates->gr_Td, grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
    grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust],
    grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], my_rates->gas_grain2, &my_rates->gamma_isrf2,
    &coolunit, gasgr, myisrf, internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.sigma_per_gas_mass_tot,
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Al2O3_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::ref_org_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.dyntab_kappa_tot, grain_kappa.data[OnlyGrainSpLUT::SiM_dust], grain_kappa.data[OnlyGrainSpLUT::FeM_dust],
    grain_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], grain_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], grain_kappa.data[OnlyGrainSpLUT::AC_dust], grain_kappa.data[OnlyGrainSpLUT::SiO2_dust],
    grain_kappa.data[OnlyGrainSpLUT::MgO_dust], grain_kappa.data[OnlyGrainSpLUT::FeS_dust], grain_kappa.data[OnlyGrainSpLUT::Al2O3_dust], grain_kappa.data[OnlyGrainSpLUT::ref_org_dust], grain_kappa.data[OnlyGrainSpLUT::vol_org_dust],
    grain_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], kappa_tot, gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust]
  );

}

/// Compute grain size increment
///
/// @note
/// The description could obviously be improved! My general sense is that we
/// are computing dust-properties in each zone. Among other things, this
/// computes size and precomputes the dust opacity table
///
/// @param[in] dom a standard quantity used throughout the codebase
/// @param[in] idx_range Specifies the current index-range
/// @param[in] itmask_metal Specifies the `idx_range`'s iteration-mask
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] my_rates holds assorted rate data.
/// @param[in] my_fields specifies the field data
/// @param[in,out] internal_dust_prop_buf Holds dust-specific information that
///     gets updated by this function
inline void calc_grain_size_increment_1d (
  double dom, IndexRange idx_range, gr_mask_type* itmask_metal,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf
) {

  FORTRAN_NAME(calc_grain_size_increment_1d)(
    &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->grain_growth, itmask_metal,
    &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &dom, my_fields->density,
    my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
    my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
    my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
    my_fields->metal_density, my_fields->local_ISM_metal_density,
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
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.sigma_per_gas_mass_tot,
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Al2O3_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::ref_org_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.dyntab_kappa_tot
  );
}

inline void calc_temp1d_cloudy_g(
  double* rhoH, IndexRange idx_range, double* tgas, double* mmw, double dom,
  double zr, int imetal, cloudy_data cloudy_primordial, gr_mask_type* itmask,
  chemistry_data* my_chemistry, grackle_field_data* my_fields,
  InternalGrUnits internalu
) {
  FORTRAN_NAME(calc_temp1d_cloudy_g)(
    my_fields->density, my_fields->metal_density, my_fields->internal_energy, rhoH,
    &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &idx_range.jp1, &idx_range.kp1,
    tgas, mmw, &dom, &zr,
    &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
    &my_chemistry->Gamma, &internalu.utem, &imetal,
    &cloudy_primordial.grid_rank, cloudy_primordial.grid_dimension,
    cloudy_primordial.grid_parameters[0], cloudy_primordial.grid_parameters[1], cloudy_primordial.grid_parameters[2],
    &cloudy_primordial.data_size, cloudy_primordial.mmw_data,
    itmask
  );

}

inline void ceiling_species_g(
  int imetal, chemistry_data* my_chemistry, grackle_field_data* my_fields
) {

   FORTRAN_NAME(ceiling_species_g)(my_fields->density, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
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

}

inline void cool1d_multi_g(
  int imetal, IndexRange idx_range, int iter, double* edot, double* tgas,
  double* mmw, double* p2d, double* tdust, double* metallicity,
  double* dust2gas, double* rhoH, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields,
  photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf
) {

  // the following 2 variables should not be arguments (they are only inside
  // of cool1d_multi_g to temporarily store a local value)
  //
  // We will fix this in the future
  double comp1, comp2;

  // TODO: we should really pass in a value for min_metallicity (this is
  // computed in a few independent places)

         FORTRAN_NAME(cool1d_multi_g)(
                  my_fields->density, my_fields->internal_energy, my_fields->x_velocity, my_fields->y_velocity, my_fields->z_velocity, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                  &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                  &internalu.extfields_in_comoving, &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->metal_cooling,
                  &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
                  &my_fields->grid_rank, &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &my_chemistry->ih2co, &my_chemistry->ipiht, &iter, &my_chemistry->photoelectric_heating,
                  &internalu.a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
                  &internalu.utem, &internalu.uxyz, &internalu.a_units, &internalu.urho, &internalu.tbase1,
                  &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
                  my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
                  my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
                  my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
                  &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates.comp_xray, &my_uvb_rates.temp_xray,
                  &my_uvb_rates.piHI, &my_uvb_rates.piHeI, &my_uvb_rates.piHeII, &comp1, &comp2,
                  my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                  my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
                  coolingheating_buf.hyd01k, coolingheating_buf.h2k01, coolingheating_buf.vibh, coolingheating_buf.roth, coolingheating_buf.rotl,
                  my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, coolingheating_buf.gpldl, coolingheating_buf.gphdl,
                  my_rates->HDlte, my_rates->HDlow, coolingheating_buf.hdlte, coolingheating_buf.hdlow,
                  my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp, my_rates->GAel,
                  my_rates->H2LTE, my_rates->gas_grain,
                  coolingheating_buf.ceHI, coolingheating_buf.ceHeI, coolingheating_buf.ceHeII, coolingheating_buf.ciHI, coolingheating_buf.ciHeI, coolingheating_buf.ciHeIS, coolingheating_buf.ciHeII,
                  coolingheating_buf.reHII, coolingheating_buf.reHeII1, coolingheating_buf.reHeII2, coolingheating_buf.reHeIII, coolingheating_buf.brem,
                  logTlininterp_buf.indixe, logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.logtem, logTlininterp_buf.tdef, edot,
                  tgas, cool1dmulti_buf.tgasold, mmw, p2d, tdust, metallicity,
                  dust2gas, rhoH, cool1dmulti_buf.mynh, cool1dmulti_buf.myde,
                  cool1dmulti_buf.gammaha_eff, cool1dmulti_buf.gasgr_tdust, cool1dmulti_buf.regr,
                  &my_chemistry->self_shielding_method, &my_uvb_rates.crsHI, &my_uvb_rates.crsHeI, &my_uvb_rates.crsHeII,
                  &my_uvb_rates.k24, &my_uvb_rates.k26,
                  &my_chemistry->use_radiative_transfer, my_fields->RT_heating_rate,
                  &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, coolingheating_buf.cieco,
                  &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
                  &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
                  my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
                  &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
                  &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
                  my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
                  &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
                  &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, my_fields->volumetric_heating_rate, my_fields->specific_heating_rate,
                  &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, my_fields->temperature_floor,
                  &my_chemistry->use_isrf_field, my_fields->isrf_habing, itmask, itmask_metal,
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
                 grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
                 grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust],
                 grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust], grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust],
                 my_rates->gas_grain2, &my_rates->gamma_isrf2
            );

}


/// Performs Gauss-Jordan elimination to solve the specified system of linear
/// equations and inverts the coefficient matrix.
///
/// In more detail, it solves the linear matrix `ax=b`, where `a` is the
/// square coefficient matrix, `b` is the right-hand side vector and `x`
/// is the solution vector
///
/// @param[in]     n The number of linear equations being solved
/// @param[in,out] coef_matrix An n by n array that initially specifies the
///    coefficient matrix. It's overwritten by the inverted matrix.
/// @param[in,out] vec An n element array that initially specifies the
///    right-hand side vector. It's overwritten by the solution vector.
///
/// @retval 0 indicates success
/// @retval 1 indicates that the matrix is singular
///
/// > [!important]
/// > This appears to be taken directly from a routine provided in
/// > "Numerical Recipes" and slightly adapted. The original text allows
/// > `b` to hold multiple right-hand side vectors (and then we compute
/// > solutions for each right-hand side vector at the same time). In
/// > contrast, only allow 1 vector).
/// >
/// > The original Numerical Recipes snippet can be found here
/// > https://phys.uri.edu/nigh/NumRec/bookfpdf/f2-1.pdf
/// >
/// > Code from Numerical Recipes CANNOT be merged into Grackle (the
/// > licensing is fundamentally incompatible)
///
/// @todo
/// We need to replace this for reasons highlighted up above. When we do that,
/// we should account that the only place that calls this routine only needs
/// the solution to the system of equations (i.e. it does not need the
/// inverted matrix). For that reason, we may want to rename this to something
/// a little more generic. We should also consider information highlighted
/// within https://github.com/grackle-project/grackle/issues/255
inline int gaussj_g(int n, double* coef_matrix, double* vector) {
  int ierr;
  FORTRAN_NAME(gaussj_g)(&n, coef_matrix, vector, &ierr);
  return ierr;
}

/// wrapper for 1d interpolation
inline double interpolate_1d_g(
  double input1,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 1 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {
  double value;
  FORTRAN_NAME(interpolate_1d_g)(
    &input1, gridDim, gridPar1, &dgridPar1, &dataSize, dataField, &value
  );
  return value;
}

/// wrapper for 2d interpolation
inline double interpolate_2d_g(
  double input1, double input2,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 2 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  gr_i64 dataSize, const double* dataField
) {
  double value;
  FORTRAN_NAME(interpolate_2d_g)(
    &input1, &input2,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2,
    &dataSize, dataField, &value
  );
  return value;
}

/// Helper function used to implement interpolate_3dz_g
inline double interpolate_2df3d_g(
  double input1, double input3,
  const gr_i64* GRIMPL_RESTRICT gridDim, // 3 elements
  const double* GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 index2,
  const double* GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize,
  const double* GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_2df3d_g)(
    &input1, &input3,
    gridDim,
    gridPar1, &dgridPar1, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// Similar to interpolate_3d_g except index2 is calculated ahead of time
/// because it corresponds to redshift, which will not change for the entire
/// grid during the current calculation.
inline double interpolate_3dz_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, gr_i64 index2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
  gr_i64 end_int
) {

  double value;
  FORTRAN_NAME(interpolate_3dz_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &end_int, &value
  );
  return value;

}

/// wraps 3d interpolator functions
inline double interpolate_3d_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_3d_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 4d interpolator functions
inline double interpolate_4d_g(
  double input1, double input2, double input3, double input4,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 4 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_4d_g)(
    &input1, &input2, &input3, &input4,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 5d interpolation routine
inline double interpolate_5d_g(
  double input1, double input2, double input3, double input4, double input5,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 5 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  const double * GRIMPL_RESTRICT gridPar5, double dgridPar5,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_5d_g)(
    &input1, &input2, &input3, &input4, &input5,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4, gridPar5, &dgridPar5,
    &dataSize, dataField, &value
  );
  return value;

}

// the following case was handcoded (so the argument order may shift when we
// actually transcribe the routine)
inline void make_consistent_g(
  int imetal, double dom, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields
){
  FORTRAN_NAME(make_consistent_g)(my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                         my_fields->density, &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
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
}

/// This routine calculates the electron and HI rates of change in order to
/// determine the maximum permitted timestep
inline void rate_timestep_g(
  double* dedot, double* HIdot, gr_mask_type anydust, IndexRange idx_range,
  double* h2dust, double* rhoH, gr_mask_type* itmask, double* edot,
  double chunit, double dom, chemistry_data* my_chemistry,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  grackle::impl::CollisionalRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
) {
  FORTRAN_NAME(rate_timestep_g)(
                         dedot, HIdot, &my_chemistry->primordial_chemistry, &anydust,
                         my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, my_fields->density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density,
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
                         kcr_buf.data[CollisionalRxnLUT::k1], kcr_buf.data[CollisionalRxnLUT::k2], kcr_buf.data[CollisionalRxnLUT::k3], kcr_buf.data[CollisionalRxnLUT::k4], kcr_buf.data[CollisionalRxnLUT::k5], kcr_buf.data[CollisionalRxnLUT::k6], kcr_buf.data[CollisionalRxnLUT::k7], kcr_buf.data[CollisionalRxnLUT::k8], kcr_buf.data[CollisionalRxnLUT::k9], kcr_buf.data[CollisionalRxnLUT::k10], kcr_buf.data[CollisionalRxnLUT::k11],
                         kcr_buf.data[CollisionalRxnLUT::k12], kcr_buf.data[CollisionalRxnLUT::k13], kcr_buf.data[CollisionalRxnLUT::k14], kcr_buf.data[CollisionalRxnLUT::k15], kcr_buf.data[CollisionalRxnLUT::k16], kcr_buf.data[CollisionalRxnLUT::k17], kcr_buf.data[CollisionalRxnLUT::k18], kcr_buf.data[CollisionalRxnLUT::k19], kcr_buf.data[CollisionalRxnLUT::k22],
                         &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30,
                         kcr_buf.data[CollisionalRxnLUT::k50], kcr_buf.data[CollisionalRxnLUT::k51], kcr_buf.data[CollisionalRxnLUT::k52], kcr_buf.data[CollisionalRxnLUT::k53], kcr_buf.data[CollisionalRxnLUT::k54], kcr_buf.data[CollisionalRxnLUT::k55], kcr_buf.data[CollisionalRxnLUT::k56], kcr_buf.data[CollisionalRxnLUT::k57], kcr_buf.data[CollisionalRxnLUT::k58],
                         h2dust, chemheatrates_buf.n_cr_n, chemheatrates_buf.n_cr_d1, chemheatrates_buf.n_cr_d2, rhoH,
                         kshield_buf.k24, kshield_buf.k25, kshield_buf.k26,
                         kshield_buf.k28, kshield_buf.k29, kshield_buf.k30, kshield_buf.k31,
                         &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
                         my_fields->RT_HI_ionization_rate, my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate,
                         itmask, edot, &chunit, &dom, my_fields->metal_density,
                        my_fields->HDI_density, &my_chemistry->metal_chemistry, my_fields->CI_density, my_fields->OI_density, my_fields->OH_density, my_fields->CO_density, my_fields->H2O_density,
                        &my_chemistry->radiative_transfer_HDI_dissociation, my_fields->RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, my_fields->RT_CI_ionization_rate, my_fields->RT_OI_ionization_rate,
                        &my_chemistry->radiative_transfer_metal_dissociation, my_fields->RT_CO_dissociation_rate, my_fields->RT_OH_dissociation_rate, my_fields->RT_H2O_dissociation_rate
                              );
}


inline void scale_fields_g(
  int imetal, gr_float factor, chemistry_data* my_chemistry,
  grackle_field_data* my_fields
) {

  FORTRAN_NAME(scale_fields_g)(
      &my_chemistry->primordial_chemistry, &imetal, &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry,
      &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->multi_metals, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation, &factor,
      &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2], &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
      my_fields->density, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
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

/// Uses one linearly implicit Gauss-Seidel sweep of a backward-Euler time
/// integrator to advance the rate equations by one (sub-)cycle (dtit).
inline void step_rate_g(
  double* dtit, IndexRange idx_range, gr_mask_type anydust, double* h2dust,
  double* rhoH, double* dedot_prev, double* HIdot_prev,
  gr_mask_type* itmask, gr_mask_type* itmask_metal, int imetal,
  chemistry_data* my_chemistry, grackle_field_data* my_fields,
  photo_rate_storage my_uvb_rates,
  grackle::impl::GrainSpeciesCollection grain_growth_rates,
  grackle::impl::SpeciesCollection species_tmpdens,
  grackle::impl::CollisionalRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf
) {
           FORTRAN_NAME(step_rate_g)(my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, my_fields->density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, dtit,
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
                         &my_chemistry->primordial_chemistry, &anydust,
                         kcr_buf.data[CollisionalRxnLUT::k1], kcr_buf.data[CollisionalRxnLUT::k2], kcr_buf.data[CollisionalRxnLUT::k3], kcr_buf.data[CollisionalRxnLUT::k4], kcr_buf.data[CollisionalRxnLUT::k5], kcr_buf.data[CollisionalRxnLUT::k6], kcr_buf.data[CollisionalRxnLUT::k7], kcr_buf.data[CollisionalRxnLUT::k8], kcr_buf.data[CollisionalRxnLUT::k9], kcr_buf.data[CollisionalRxnLUT::k10], kcr_buf.data[CollisionalRxnLUT::k11],
                         kcr_buf.data[CollisionalRxnLUT::k12], kcr_buf.data[CollisionalRxnLUT::k13], kcr_buf.data[CollisionalRxnLUT::k14], kcr_buf.data[CollisionalRxnLUT::k15], kcr_buf.data[CollisionalRxnLUT::k16], kcr_buf.data[CollisionalRxnLUT::k17], kcr_buf.data[CollisionalRxnLUT::k18], kcr_buf.data[CollisionalRxnLUT::k19], kcr_buf.data[CollisionalRxnLUT::k22],
                         &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30,
                         kcr_buf.data[CollisionalRxnLUT::k50], kcr_buf.data[CollisionalRxnLUT::k51], kcr_buf.data[CollisionalRxnLUT::k52], kcr_buf.data[CollisionalRxnLUT::k53], kcr_buf.data[CollisionalRxnLUT::k54], kcr_buf.data[CollisionalRxnLUT::k55], kcr_buf.data[CollisionalRxnLUT::k56], kcr_buf.data[CollisionalRxnLUT::k57], kcr_buf.data[CollisionalRxnLUT::k58],
                         h2dust, rhoH,
                         kshield_buf.k24, kshield_buf.k25, kshield_buf.k26,
                         kshield_buf.k28, kshield_buf.k29, kshield_buf.k30, kshield_buf.k31,
                         species_tmpdens.data[SpLUT::HI], species_tmpdens.data[SpLUT::HII], species_tmpdens.data[SpLUT::HeI], species_tmpdens.data[SpLUT::HeII], species_tmpdens.data[SpLUT::HeIII], species_tmpdens.data[SpLUT::e],
                         species_tmpdens.data[SpLUT::HM], species_tmpdens.data[SpLUT::H2I], species_tmpdens.data[SpLUT::H2II], species_tmpdens.data[SpLUT::DI], species_tmpdens.data[SpLUT::DII], species_tmpdens.data[SpLUT::HDI],
                         dedot_prev, HIdot_prev,
                         &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
                         my_fields->RT_HI_ionization_rate, my_fields->RT_HeI_ionization_rate, my_fields->RT_HeII_ionization_rate,
                         itmask, itmask_metal,
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
                        kcr_buf.data[CollisionalRxnLUT::k125], kcr_buf.data[CollisionalRxnLUT::k129], kcr_buf.data[CollisionalRxnLUT::k130], kcr_buf.data[CollisionalRxnLUT::k131], kcr_buf.data[CollisionalRxnLUT::k132],
                        kcr_buf.data[CollisionalRxnLUT::k133], kcr_buf.data[CollisionalRxnLUT::k134], kcr_buf.data[CollisionalRxnLUT::k135], kcr_buf.data[CollisionalRxnLUT::k136], kcr_buf.data[CollisionalRxnLUT::k137],
                        kcr_buf.data[CollisionalRxnLUT::k148], kcr_buf.data[CollisionalRxnLUT::k149], kcr_buf.data[CollisionalRxnLUT::k150], kcr_buf.data[CollisionalRxnLUT::k151], kcr_buf.data[CollisionalRxnLUT::k152],
                        kcr_buf.data[CollisionalRxnLUT::k153],
                        kcr_buf.data[CollisionalRxnLUT::kz15],  kcr_buf.data[CollisionalRxnLUT::kz16],  kcr_buf.data[CollisionalRxnLUT::kz17],  kcr_buf.data[CollisionalRxnLUT::kz18],  kcr_buf.data[CollisionalRxnLUT::kz19],
                        kcr_buf.data[CollisionalRxnLUT::kz20],  kcr_buf.data[CollisionalRxnLUT::kz21],  kcr_buf.data[CollisionalRxnLUT::kz22],  kcr_buf.data[CollisionalRxnLUT::kz23],  kcr_buf.data[CollisionalRxnLUT::kz24],
                        kcr_buf.data[CollisionalRxnLUT::kz25],  kcr_buf.data[CollisionalRxnLUT::kz26],  kcr_buf.data[CollisionalRxnLUT::kz27],  kcr_buf.data[CollisionalRxnLUT::kz28],  kcr_buf.data[CollisionalRxnLUT::kz29],
                        kcr_buf.data[CollisionalRxnLUT::kz30],  kcr_buf.data[CollisionalRxnLUT::kz31],  kcr_buf.data[CollisionalRxnLUT::kz32],  kcr_buf.data[CollisionalRxnLUT::kz33],  kcr_buf.data[CollisionalRxnLUT::kz34],
                        kcr_buf.data[CollisionalRxnLUT::kz35],  kcr_buf.data[CollisionalRxnLUT::kz36],  kcr_buf.data[CollisionalRxnLUT::kz37],  kcr_buf.data[CollisionalRxnLUT::kz38],  kcr_buf.data[CollisionalRxnLUT::kz39],
                        kcr_buf.data[CollisionalRxnLUT::kz40],  kcr_buf.data[CollisionalRxnLUT::kz41],  kcr_buf.data[CollisionalRxnLUT::kz42],  kcr_buf.data[CollisionalRxnLUT::kz43],  kcr_buf.data[CollisionalRxnLUT::kz44],
                        kcr_buf.data[CollisionalRxnLUT::kz45],  kcr_buf.data[CollisionalRxnLUT::kz46],  kcr_buf.data[CollisionalRxnLUT::kz47],  kcr_buf.data[CollisionalRxnLUT::kz48],  kcr_buf.data[CollisionalRxnLUT::kz49],
                        kcr_buf.data[CollisionalRxnLUT::kz50],  kcr_buf.data[CollisionalRxnLUT::kz51],  kcr_buf.data[CollisionalRxnLUT::kz52],  kcr_buf.data[CollisionalRxnLUT::kz53],  kcr_buf.data[CollisionalRxnLUT::kz54],
                        species_tmpdens.data[SpLUT::DM], species_tmpdens.data[SpLUT::HDII], species_tmpdens.data[SpLUT::HeHII],
                        species_tmpdens.data[SpLUT::CI], species_tmpdens.data[SpLUT::CII], species_tmpdens.data[SpLUT::CO], species_tmpdens.data[SpLUT::CO2],
                        species_tmpdens.data[SpLUT::OI], species_tmpdens.data[SpLUT::OH], species_tmpdens.data[SpLUT::H2O], species_tmpdens.data[SpLUT::O2],
                        species_tmpdens.data[SpLUT::SiI], species_tmpdens.data[SpLUT::SiOI], species_tmpdens.data[SpLUT::SiO2I],
                        species_tmpdens.data[SpLUT::CH], species_tmpdens.data[SpLUT::CH2], species_tmpdens.data[SpLUT::COII], species_tmpdens.data[SpLUT::OII],
                        species_tmpdens.data[SpLUT::OHII], species_tmpdens.data[SpLUT::H2OII], species_tmpdens.data[SpLUT::H3OII], species_tmpdens.data[SpLUT::O2II],
                        species_tmpdens.data[SpLUT::Mg], species_tmpdens.data[SpLUT::Al], species_tmpdens.data[SpLUT::S], species_tmpdens.data[SpLUT::Fe],
                        species_tmpdens.data[SpLUT::SiM_dust], species_tmpdens.data[SpLUT::FeM_dust], species_tmpdens.data[SpLUT::Mg2SiO4_dust], species_tmpdens.data[SpLUT::MgSiO3_dust], species_tmpdens.data[SpLUT::Fe3O4_dust],
                        species_tmpdens.data[SpLUT::AC_dust], species_tmpdens.data[SpLUT::SiO2_dust], species_tmpdens.data[SpLUT::MgO_dust], species_tmpdens.data[SpLUT::FeS_dust], species_tmpdens.data[SpLUT::Al2O3_dust],
                        species_tmpdens.data[SpLUT::ref_org_dust], species_tmpdens.data[SpLUT::vol_org_dust], species_tmpdens.data[SpLUT::H2O_ice_dust],
                        grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust], grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust], grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust], grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust],
                        grain_growth_rates.data[OnlyGrainSpLUT::AC_dust], grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust], grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust], grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust], grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust],
                        grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust], grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust], grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust],
                        &my_chemistry->radiative_transfer_HDI_dissociation, my_fields->RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, my_fields->RT_CI_ionization_rate, my_fields->RT_OI_ionization_rate,
                        &my_chemistry->radiative_transfer_metal_dissociation, my_fields->RT_CO_dissociation_rate, my_fields->RT_OH_dissociation_rate, my_fields->RT_H2O_dissociation_rate
               );

}


} // namespace grackle::impl::fortran_wrapper

#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

