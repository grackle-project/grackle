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

// this should go before the header-guard
#if !defined(GRIMPL_PUBLIC_INCLUDE) && !defined(GRIMPL_COMPILING_CORE_LIB)
  #include "grackle_misc.h"
  GRIMPL_COMPTIME_WARNING(
    "You are using a deprecated header file; include the public \"grackle.h\" "
    "header file instead! In a future Grackle version, \"grackle_types.h\" "
    "may cease to exist (or contents may change in an incompatible manner)."
  );
#endif


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

#define GRIMPL_MAX_INJ_PATHWAYS 12

typedef struct
{

  int grid_rank;
  int *grid_dimension;
  int *grid_start;
  int *grid_end;

  gr_float grid_dx;

  gr_float *density;
  gr_float *internal_energy;
  gr_float *x_velocity;
  gr_float *y_velocity;
  gr_float *z_velocity;

  // metal_cooling = 1
  gr_float *metal_density;

  // use_dust_density_field = 1
  gr_float *dust_density;

  // primordial_chemistry = 1
  gr_float *e_density;
  gr_float *HI_density;
  gr_float *HII_density;
  gr_float *HeI_density;
  gr_float *HeII_density;
  gr_float *HeIII_density;

  // primordial_chemistry = 2
  gr_float *HM_density;
  gr_float *H2I_density;
  gr_float *H2II_density;

  // primordial_chemistry = 3
  gr_float *DI_density;
  gr_float *DII_density;
  gr_float *HDI_density;

  // primordial_chemistry = 4
  gr_float *DM_density;
  gr_float *HDII_density;
  gr_float *HeHII_density;

  // metal_chemistry = 1
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

  // dust_species = 1
  gr_float *Mg_density;

  // dust_species = 2
  gr_float *Al_density;
  gr_float *S_density;
  gr_float *Fe_density;

  // dust_species = 1
  gr_float *MgSiO3_dust_density; // enstatite
  gr_float *AC_dust_density; // amorphous carbon

  // dust_species = 2
  gr_float *SiM_dust_density; // metallic silicon
  gr_float *FeM_dust_density; // metallic iron
  gr_float *Mg2SiO4_dust_density; // forsterite
  gr_float *Fe3O4_dust_density; // magnetite
  gr_float *SiO2_dust_density; // silica
  gr_float *MgO_dust_density; // magnesia
  gr_float *FeS_dust_density; // troilite
  gr_float *Al2O3_dust_density; // alumina

  // dust_species = 3
  gr_float *ref_org_dust_density; // refractory organics
  gr_float *vol_org_dust_density; // volatile organics
  gr_float *H2O_ice_dust_density; // water ice

  // metal_chemistry = 1 && multi_metals = 1
  // -> the number and names of models are queried with the experimental
  //    ratequery interface
  // -> when metal_chemistry = 1 && multi_metals = 0, Grackle currently expects
  //    this to be empty and uses metal_density, instead (This may change in
  //    the near future as we start using hdf5 configuration files)
  gr_float* inject_pathway_metal_density[GRIMPL_MAX_INJ_PATHWAYS];

  // use_volumetric_heating_rate = 1
  gr_float *volumetric_heating_rate;
  // use_specific_heating_rate = 1
  gr_float *specific_heating_rate;
  // use_temperature_floor = 1
  gr_float *temperature_floor;

  // use_radiative_transfer = 1
  // primordial_chemistry = 1
  gr_float *RT_heating_rate;
  gr_float *RT_HI_ionization_rate;
  gr_float *RT_HeI_ionization_rate;
  gr_float *RT_HeII_ionization_rate;
  // primordial_chemistry = 2
  gr_float *RT_H2_dissociation_rate;

  // radiative_transfer_HDI_dissociation = 1
  gr_float *RT_HDI_dissociation_rate;
  // radiative_transfer_metal_ionization = 1
  gr_float *RT_CI_ionization_rate;
  gr_float *RT_OI_ionization_rate;
  // radiative_transfer_metal_dissociation = 1
  gr_float *RT_CO_dissociation_rate;
  gr_float *RT_OH_dissociation_rate;
  gr_float *RT_H2O_dissociation_rate;

  // H2_self_shielding = 2
  gr_float *H2_self_shielding_length;
  // H2_custom_shielding = 1
  gr_float *H2_custom_shielding_factor;

  // use_isrf_field = 1
  gr_float *isrf_habing;

  // use_multiple_dust_temperatures = 1
  gr_float *SiM_dust_temperature;
  gr_float *FeM_dust_temperature;
  gr_float *Mg2SiO4_dust_temperature;
  gr_float *MgSiO3_dust_temperature;
  gr_float *Fe3O4_dust_temperature;
  gr_float *AC_dust_temperature;
  gr_float *SiO2_dust_temperature;
  gr_float *MgO_dust_temperature;
  gr_float *FeS_dust_temperature;
  gr_float *Al2O3_dust_temperature;
  gr_float *ref_org_dust_temperature;
  gr_float *vol_org_dust_temperature;
  gr_float *H2O_ice_dust_temperature;

} grackle_field_data;


typedef struct
{

  const char* version;
  const char* branch;
  const char* revision;

} grackle_version;

#if !defined(GRIMPL_COMPILING_CORE_LIB)
// avoid leaking these macros used for defining implementation details out of
// the core library
#undef GRIMPL_MAX_INJ_PATHWAYS
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __GRACKLE_TYPES_H__ */
