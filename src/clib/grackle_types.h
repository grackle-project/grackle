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

#ifndef __GRACKLE_TYPES_H_
#define __GRACKLE_TYPES_H_
/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

#ifdef CONFIG_BFLOAT_4
#define gr_float float
#endif

#ifdef CONFIG_BFLOAT_8
#define gr_float double
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
//#ifdef GRACKLE_MD
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
//#endif
  gr_float *e_density;
  gr_float *metal_density;
  gr_float *dust_density;
//#ifdef GRACKLE_MD
  gr_float *metal_loc;
  gr_float *metal_C30;
  gr_float *metal_F13;
//#endif

  gr_float *internal_energy;
  gr_float *x_velocity;
  gr_float *y_velocity;
  gr_float *z_velocity;

  gr_float *volumetric_heating_rate;
  gr_float *specific_heating_rate;

  gr_float *RT_heating_rate;
  gr_float *RT_HI_ionization_rate;
  gr_float *RT_HeI_ionization_rate;
  gr_float *RT_HeII_ionization_rate;
  gr_float *RT_H2_dissociation_rate;

  gr_float *H2_self_shielding_length;

  gr_float *isrf_habing;

} grackle_field_data;

typedef struct
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;

} code_units;

#endif
