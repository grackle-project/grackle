/***********************************************************************
/
/ Grackle variable types
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_TYPES_H__
#define __GRACKLE_TYPES_H__
/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

#include "grackle_float.h"

// the following include-directive only exists because the definition of
// the code_units type got moved and we don't want to break anybody's
// code who expect that type to be defined in this file
#include "grackle_chemistry_data.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef GRACKLE_FLOAT_4
#define gr_float float
#endif

#ifdef GRACKLE_FLOAT_8
#define gr_float double
#endif

#if defined(GRACKLE_FLOAT_4) == defined(GRACKLE_FLOAT_8)
#error "Both GRACKLE_FLOAT_4 and GRACKLE_FLOAT_8 are defined. Only one can be defined."
#endif

typedef struct
{

  int grid_rank;
  int *grid_dimension;
  int *grid_start;
  int *grid_end;

  gr_float grid_dx;

  gr_float *density;
  gr_float *HI_density;
  gr_float *HII_density;
  gr_float *HM_density;
  gr_float *HeI_density;
  gr_float *HeII_density;
  gr_float *HeIII_density;
  gr_float *H2I_density;
  gr_float *H2II_density;
  gr_float *DI_density;
  gr_float *DII_density;
  gr_float *HDI_density;
  gr_float *DM_density;
  gr_float *HDII_density;
  gr_float *HeHII_density;
  gr_float *CI_density;
  gr_float *CII_density;
  gr_float *CO_density;
  gr_float *CO2_density;
  gr_float *OI_density;
  gr_float *OH_density;
  gr_float *H2O_density;
  gr_float *O2_density;
  gr_float *SiI_density;
  gr_float *SiOI_density;
  gr_float *SiO2I_density;
  gr_float *CH_density;
  gr_float *CH2_density;
  gr_float *COII_density;
  gr_float *OII_density;
  gr_float *OHII_density;
  gr_float *H2OII_density;
  gr_float *H3OII_density;
  gr_float *O2II_density;
  gr_float *Mg_density;
  gr_float *Al_density;
  gr_float *S_density;
  gr_float *Fe_density;
  gr_float *SiM_density;
  gr_float *FeM_density;
  gr_float *Mg2SiO4_density;
  gr_float *MgSiO3_density;
  gr_float *Fe3O4_density;
  gr_float *AC_density;
  gr_float *SiO2D_density;
  gr_float *MgO_density;
  gr_float *FeS_density;
  gr_float *Al2O3_density;
  gr_float *reforg_density;
  gr_float *volorg_density;
  gr_float *H2Oice_density;
  gr_float *e_density;
  gr_float *metal_density;
  gr_float *dust_density;
  gr_float *local_ISM_metal_density;
  gr_float *ccsn13_metal_density;
  gr_float *ccsn20_metal_density;
  gr_float *ccsn25_metal_density;
  gr_float *ccsn30_metal_density;
  gr_float *fsn13_metal_density;
  gr_float *fsn15_metal_density;
  gr_float *fsn50_metal_density;
  gr_float *fsn80_metal_density;
  gr_float *pisn170_metal_density;
  gr_float *pisn200_metal_density;
  gr_float *y19_metal_density;

  gr_float *internal_energy;
  gr_float *x_velocity;
  gr_float *y_velocity;
  gr_float *z_velocity;

  gr_float *volumetric_heating_rate;
  gr_float *specific_heating_rate;

  gr_float *temperature_floor;

  gr_float *RT_heating_rate;
  gr_float *RT_HI_ionization_rate;
  gr_float *RT_HeI_ionization_rate;
  gr_float *RT_HeII_ionization_rate;
  gr_float *RT_H2_dissociation_rate;

  gr_float *RT_HDI_dissociation_rate;
  gr_float *RT_CI_ionization_rate;
  gr_float *RT_OI_ionization_rate;
  gr_float *RT_CO_dissociation_rate;
  gr_float *RT_OH_dissociation_rate;
  gr_float *RT_H2O_dissociation_rate;

  gr_float *H2_self_shielding_length;
  gr_float *H2_custom_shielding_factor;

  gr_float *isrf_habing;

  // Temporary fields for primordial chemistry == 4
  gr_float *SiM_temperature;
  gr_float *FeM_temperature;
  gr_float *Mg2SiO4_temperature;
  gr_float *MgSiO3_temperature;
  gr_float *Fe3O4_temperature;
  gr_float *AC_temperature;
  gr_float *SiO2D_temperature;
  gr_float *MgO_temperature;
  gr_float *FeS_temperature;
  gr_float *Al2O3_temperature;
  gr_float *reforg_temperature;
  gr_float *volorg_temperature;
  gr_float *H2Oice_temperature;

} grackle_field_data;


typedef struct
{

  const char* version;
  const char* branch;
  const char* revision;

} grackle_version;

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __GRACKLE_TYPES_H__ */
