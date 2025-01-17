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

inline void wrapped_ceiling_species_g_(
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


// the following case was handcoded (so the argument order may shift when we
// actually transcribe the routine)
inline void wrapped_make_consistent_g_(
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


inline void wrapped_scale_fields_g_(
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


#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

