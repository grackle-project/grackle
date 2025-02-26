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
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "LUT.hpp"

// callers of these functions are generally expected to locally shorten the
// namespace name when they call these routines
namespace grackle::impl::fortran_wrapper {

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

/// This routine uses the temperature to look up the chemical rates that are
/// tabulated in a log table as a function of temperature
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
  double* tdust, double* dust2gas, double* k13dd, double* h2dust,
  double dom, double dx_cgs, double c_ljeans, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, int imetal, gr_float* rhoH, double dt,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  InternalGrUnits internalu,
  grackle::impl::GrainSpeciesCollection grain_growth_rates,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::ColRecRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
) {

  FORTRAN_NAME(lookup_cool_rates1d_g)(&my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->NumberOfTemperatureBins, &idx_range.jp1, &idx_range.kp1,
                   &idx_range.i_start, &idx_range.i_end, &my_chemistry->three_body_rate,
                   &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->primordial_chemistry, &anydust,
                   &my_chemistry->H2_self_shielding, &my_chemistry->self_shielding_method,
                   tgas1d, mmw, my_fields->density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                   my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density,
                   tdust, dust2gas,
                   my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
                   my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
                   my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
                   my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54, my_rates->k55, my_rates->k56,
                   my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
                   my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
                   &my_uvb_rates.crsHI, &my_uvb_rates.crsHeI, &my_uvb_rates.crsHeII, &my_uvb_rates.piHI, &my_uvb_rates.piHeI,
                   kcr_buf.data[ColRecRxnLUT::k1], kcr_buf.data[ColRecRxnLUT::k2], kcr_buf.data[ColRecRxnLUT::k3], kcr_buf.data[ColRecRxnLUT::k4], kcr_buf.data[ColRecRxnLUT::k5], kcr_buf.data[ColRecRxnLUT::k6], kcr_buf.data[ColRecRxnLUT::k7], kcr_buf.data[ColRecRxnLUT::k8], kcr_buf.data[ColRecRxnLUT::k9], kcr_buf.data[ColRecRxnLUT::k10],
                   kcr_buf.data[ColRecRxnLUT::k11], kcr_buf.data[ColRecRxnLUT::k12], kcr_buf.data[ColRecRxnLUT::k13], kcr_buf.data[ColRecRxnLUT::k14], kcr_buf.data[ColRecRxnLUT::k15], kcr_buf.data[ColRecRxnLUT::k16], kcr_buf.data[ColRecRxnLUT::k17], kcr_buf.data[ColRecRxnLUT::k18],
                   kcr_buf.data[ColRecRxnLUT::k19], kcr_buf.data[ColRecRxnLUT::k22], &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30, &my_uvb_rates.k31,
                   kcr_buf.data[ColRecRxnLUT::k50], kcr_buf.data[ColRecRxnLUT::k51], kcr_buf.data[ColRecRxnLUT::k52], kcr_buf.data[ColRecRxnLUT::k53], kcr_buf.data[ColRecRxnLUT::k54], kcr_buf.data[ColRecRxnLUT::k55], kcr_buf.data[ColRecRxnLUT::k56], kcr_buf.data[ColRecRxnLUT::k57],
                   kcr_buf.data[ColRecRxnLUT::k58], k13dd, kshield_buf.k24, kshield_buf.k25, kshield_buf.k26,
                   kshield_buf.k28, kshield_buf.k29, kshield_buf.k30,
                   kshield_buf.k31, h2dust, chemheatrates_buf.n_cr_n, chemheatrates_buf.n_cr_d1, chemheatrates_buf.n_cr_d2,
                   logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.tdef, logTlininterp_buf.logtem, logTlininterp_buf.indixe,
                   &dom, &internalu.coolunit, &internalu.tbase1, &internalu.xbase1, &dx_cgs, &c_ljeans,
                   &my_chemistry->use_radiative_transfer, my_fields->RT_H2_dissociation_rate, my_fields->H2_self_shielding_length, itmask,
                   itmask_metal,
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
                   kcr_buf.data[ColRecRxnLUT::k125],  kcr_buf.data[ColRecRxnLUT::k129],  kcr_buf.data[ColRecRxnLUT::k130],  kcr_buf.data[ColRecRxnLUT::k131],  kcr_buf.data[ColRecRxnLUT::k132],
                   kcr_buf.data[ColRecRxnLUT::k133],  kcr_buf.data[ColRecRxnLUT::k134],  kcr_buf.data[ColRecRxnLUT::k135],  kcr_buf.data[ColRecRxnLUT::k136],  kcr_buf.data[ColRecRxnLUT::k137],
                   kcr_buf.data[ColRecRxnLUT::k148],  kcr_buf.data[ColRecRxnLUT::k149],  kcr_buf.data[ColRecRxnLUT::k150],  kcr_buf.data[ColRecRxnLUT::k151],  kcr_buf.data[ColRecRxnLUT::k152],
                   kcr_buf.data[ColRecRxnLUT::k153],
                   kcr_buf.data[ColRecRxnLUT::kz15],  kcr_buf.data[ColRecRxnLUT::kz16],  kcr_buf.data[ColRecRxnLUT::kz17],  kcr_buf.data[ColRecRxnLUT::kz18],  kcr_buf.data[ColRecRxnLUT::kz19],
                   kcr_buf.data[ColRecRxnLUT::kz20],  kcr_buf.data[ColRecRxnLUT::kz21],  kcr_buf.data[ColRecRxnLUT::kz22],  kcr_buf.data[ColRecRxnLUT::kz23],  kcr_buf.data[ColRecRxnLUT::kz24],
                   kcr_buf.data[ColRecRxnLUT::kz25],  kcr_buf.data[ColRecRxnLUT::kz26],  kcr_buf.data[ColRecRxnLUT::kz27],  kcr_buf.data[ColRecRxnLUT::kz28],  kcr_buf.data[ColRecRxnLUT::kz29],
                   kcr_buf.data[ColRecRxnLUT::kz30],  kcr_buf.data[ColRecRxnLUT::kz31],  kcr_buf.data[ColRecRxnLUT::kz32],  kcr_buf.data[ColRecRxnLUT::kz33],  kcr_buf.data[ColRecRxnLUT::kz34],
                   kcr_buf.data[ColRecRxnLUT::kz35],  kcr_buf.data[ColRecRxnLUT::kz36],  kcr_buf.data[ColRecRxnLUT::kz37],  kcr_buf.data[ColRecRxnLUT::kz38],  kcr_buf.data[ColRecRxnLUT::kz39],
                   kcr_buf.data[ColRecRxnLUT::kz40],  kcr_buf.data[ColRecRxnLUT::kz41],  kcr_buf.data[ColRecRxnLUT::kz42],  kcr_buf.data[ColRecRxnLUT::kz43],  kcr_buf.data[ColRecRxnLUT::kz44],
                   kcr_buf.data[ColRecRxnLUT::kz45],  kcr_buf.data[ColRecRxnLUT::kz46],  kcr_buf.data[ColRecRxnLUT::kz47],  kcr_buf.data[ColRecRxnLUT::kz48],  kcr_buf.data[ColRecRxnLUT::kz49],
                   kcr_buf.data[ColRecRxnLUT::kz50],  kcr_buf.data[ColRecRxnLUT::kz51],  kcr_buf.data[ColRecRxnLUT::kz52],  kcr_buf.data[ColRecRxnLUT::kz53],  kcr_buf.data[ColRecRxnLUT::kz54],
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
                   my_rates->h2dustS, my_rates->h2dustC, rhoH, my_rates->grain_growth_rate, &dt,
                   grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust], grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust], grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust],
                   grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust], grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust], grain_growth_rates.data[OnlyGrainSpLUT::AC_dust], grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust], grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust], grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust],
                   grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust], grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust], grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust], grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust],
                   grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
                   grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust],
                   grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust], grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], &my_chemistry->radiative_transfer_use_H2_shielding,
                   &my_chemistry->H2_custom_shielding, my_fields->H2_custom_shielding_factor
          );

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
  grackle::impl::ColRecRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
) {
  FORTRAN_NAME(rate_timestep_g)(
                         dedot, HIdot, &my_chemistry->primordial_chemistry, &anydust,
                         my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, my_fields->density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density,
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
                         kcr_buf.data[ColRecRxnLUT::k1], kcr_buf.data[ColRecRxnLUT::k2], kcr_buf.data[ColRecRxnLUT::k3], kcr_buf.data[ColRecRxnLUT::k4], kcr_buf.data[ColRecRxnLUT::k5], kcr_buf.data[ColRecRxnLUT::k6], kcr_buf.data[ColRecRxnLUT::k7], kcr_buf.data[ColRecRxnLUT::k8], kcr_buf.data[ColRecRxnLUT::k9], kcr_buf.data[ColRecRxnLUT::k10], kcr_buf.data[ColRecRxnLUT::k11],
                         kcr_buf.data[ColRecRxnLUT::k12], kcr_buf.data[ColRecRxnLUT::k13], kcr_buf.data[ColRecRxnLUT::k14], kcr_buf.data[ColRecRxnLUT::k15], kcr_buf.data[ColRecRxnLUT::k16], kcr_buf.data[ColRecRxnLUT::k17], kcr_buf.data[ColRecRxnLUT::k18], kcr_buf.data[ColRecRxnLUT::k19], kcr_buf.data[ColRecRxnLUT::k22],
                         &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30,
                         kcr_buf.data[ColRecRxnLUT::k50], kcr_buf.data[ColRecRxnLUT::k51], kcr_buf.data[ColRecRxnLUT::k52], kcr_buf.data[ColRecRxnLUT::k53], kcr_buf.data[ColRecRxnLUT::k54], kcr_buf.data[ColRecRxnLUT::k55], kcr_buf.data[ColRecRxnLUT::k56], kcr_buf.data[ColRecRxnLUT::k57], kcr_buf.data[ColRecRxnLUT::k58],
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
  grackle::impl::ColRecRxnRateCollection kcr_buf,
  grackle::impl::PhotoRxnRateCollection kshield_buf
) {
           FORTRAN_NAME(step_rate_g)(my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density, my_fields->density,
                         my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, dtit,
                         &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
                         &my_chemistry->primordial_chemistry, &anydust,
                         kcr_buf.data[ColRecRxnLUT::k1], kcr_buf.data[ColRecRxnLUT::k2], kcr_buf.data[ColRecRxnLUT::k3], kcr_buf.data[ColRecRxnLUT::k4], kcr_buf.data[ColRecRxnLUT::k5], kcr_buf.data[ColRecRxnLUT::k6], kcr_buf.data[ColRecRxnLUT::k7], kcr_buf.data[ColRecRxnLUT::k8], kcr_buf.data[ColRecRxnLUT::k9], kcr_buf.data[ColRecRxnLUT::k10], kcr_buf.data[ColRecRxnLUT::k11],
                         kcr_buf.data[ColRecRxnLUT::k12], kcr_buf.data[ColRecRxnLUT::k13], kcr_buf.data[ColRecRxnLUT::k14], kcr_buf.data[ColRecRxnLUT::k15], kcr_buf.data[ColRecRxnLUT::k16], kcr_buf.data[ColRecRxnLUT::k17], kcr_buf.data[ColRecRxnLUT::k18], kcr_buf.data[ColRecRxnLUT::k19], kcr_buf.data[ColRecRxnLUT::k22],
                         &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30,
                         kcr_buf.data[ColRecRxnLUT::k50], kcr_buf.data[ColRecRxnLUT::k51], kcr_buf.data[ColRecRxnLUT::k52], kcr_buf.data[ColRecRxnLUT::k53], kcr_buf.data[ColRecRxnLUT::k54], kcr_buf.data[ColRecRxnLUT::k55], kcr_buf.data[ColRecRxnLUT::k56], kcr_buf.data[ColRecRxnLUT::k57], kcr_buf.data[ColRecRxnLUT::k58],
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
                        kcr_buf.data[ColRecRxnLUT::k125], kcr_buf.data[ColRecRxnLUT::k129], kcr_buf.data[ColRecRxnLUT::k130], kcr_buf.data[ColRecRxnLUT::k131], kcr_buf.data[ColRecRxnLUT::k132],
                        kcr_buf.data[ColRecRxnLUT::k133], kcr_buf.data[ColRecRxnLUT::k134], kcr_buf.data[ColRecRxnLUT::k135], kcr_buf.data[ColRecRxnLUT::k136], kcr_buf.data[ColRecRxnLUT::k137],
                        kcr_buf.data[ColRecRxnLUT::k148], kcr_buf.data[ColRecRxnLUT::k149], kcr_buf.data[ColRecRxnLUT::k150], kcr_buf.data[ColRecRxnLUT::k151], kcr_buf.data[ColRecRxnLUT::k152],
                        kcr_buf.data[ColRecRxnLUT::k153],
                        kcr_buf.data[ColRecRxnLUT::kz15],  kcr_buf.data[ColRecRxnLUT::kz16],  kcr_buf.data[ColRecRxnLUT::kz17],  kcr_buf.data[ColRecRxnLUT::kz18],  kcr_buf.data[ColRecRxnLUT::kz19],
                        kcr_buf.data[ColRecRxnLUT::kz20],  kcr_buf.data[ColRecRxnLUT::kz21],  kcr_buf.data[ColRecRxnLUT::kz22],  kcr_buf.data[ColRecRxnLUT::kz23],  kcr_buf.data[ColRecRxnLUT::kz24],
                        kcr_buf.data[ColRecRxnLUT::kz25],  kcr_buf.data[ColRecRxnLUT::kz26],  kcr_buf.data[ColRecRxnLUT::kz27],  kcr_buf.data[ColRecRxnLUT::kz28],  kcr_buf.data[ColRecRxnLUT::kz29],
                        kcr_buf.data[ColRecRxnLUT::kz30],  kcr_buf.data[ColRecRxnLUT::kz31],  kcr_buf.data[ColRecRxnLUT::kz32],  kcr_buf.data[ColRecRxnLUT::kz33],  kcr_buf.data[ColRecRxnLUT::kz34],
                        kcr_buf.data[ColRecRxnLUT::kz35],  kcr_buf.data[ColRecRxnLUT::kz36],  kcr_buf.data[ColRecRxnLUT::kz37],  kcr_buf.data[ColRecRxnLUT::kz38],  kcr_buf.data[ColRecRxnLUT::kz39],
                        kcr_buf.data[ColRecRxnLUT::kz40],  kcr_buf.data[ColRecRxnLUT::kz41],  kcr_buf.data[ColRecRxnLUT::kz42],  kcr_buf.data[ColRecRxnLUT::kz43],  kcr_buf.data[ColRecRxnLUT::kz44],
                        kcr_buf.data[ColRecRxnLUT::kz45],  kcr_buf.data[ColRecRxnLUT::kz46],  kcr_buf.data[ColRecRxnLUT::kz47],  kcr_buf.data[ColRecRxnLUT::kz48],  kcr_buf.data[ColRecRxnLUT::kz49],
                        kcr_buf.data[ColRecRxnLUT::kz50],  kcr_buf.data[ColRecRxnLUT::kz51],  kcr_buf.data[ColRecRxnLUT::kz52],  kcr_buf.data[ColRecRxnLUT::kz53],  kcr_buf.data[ColRecRxnLUT::kz54],
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

