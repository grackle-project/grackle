/***********************************************************************
/
/ Calculate dust temperature field
/
/
/ Copyright (c) Grackle Development Team. All rights reserved.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature);

extern void FORTRAN_NAME(calc_tdust_3d_g)(
	gr_float *d, gr_float *de, gr_float *HI, gr_float *HII,
	gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
	int *in, int *jn, int *kn,
	int *nratec, int *iexpand,
	int *ispecies, int *idim,
	int *is, int *js, int *ks,
	int *ie, int *je, int *ke,
	double *aye, double *temstart, double *temend,
	double *fgr, double *gasgra,
        double *gamma_isrfa, double *isrf,
	double *utem, double *uxyz, double *uaye,
	double *urho, double *utim,
	gr_float *gas_temp, gr_float *dust_temp,
        int *iisrffield, gr_float* isrf_habing
      , int *imetal, int *idustfield, int *igrgr
      , double *z_solar, gr_float *metal, gr_float *dust
      , gr_float *SiM, gr_float *FeM, gr_float *Mg2SiO4, gr_float *MgSiO3, gr_float *Fe3O4
      , gr_float *AC, gr_float *SiO2D, gr_float *MgO, gr_float *FeS, gr_float *Al2O3
      , gr_float *reforg, gr_float *volorg, gr_float *H2Oice
      , int *immulti, int *impop3, int *idspecies, int *itdspecies, int *idsub
      , gr_float *metal_loc, gr_float *metal_C30, gr_float *metal_F13
      , double *loc_fFeM    , double *loc_fMg2SiO4, double *loc_fMgSiO3 , double *loc_fFeS
      , double *loc_freforg , double *loc_fvolorg , double *loc_fH2Oice
      , double *loc_r0FeM    , double *loc_r0Mg2SiO4, double *loc_r0MgSiO3 , double *loc_r0FeS
      , double *loc_r0reforg , double *loc_r0volorg , double *loc_r0H2Oice 
      , double *C30_fSiM, double *C30_fFeM, double *C30_fMg2SiO4, double *C30_fMgSiO3
      , double *C30_fAC, double *C30_fSiO2D, double *C30_fMgO, double *C30_fFeS, double *C30_fAl2O3
      , double *C30_r0SiM, double *C30_r0FeM, double *C30_r0Mg2SiO4, double *C30_r0MgSiO3
      , double *C30_r0AC, double *C30_r0SiO2D, double *C30_r0MgO, double *C30_r0FeS, double *C30_r0Al2O3
      , double *F13_fFeM, double *F13_fMg2SiO4, double *F13_fMgSiO3, double *F13_fFe3O4
      , double *F13_fAC , double *F13_fSiO2D , double *F13_fAl2O3
      , double *F13_r0FeM, double *F13_r0Mg2SiO4, double *F13_r0MgSiO3, double *F13_r0Fe3O4
      , double *F13_r0AC , double *F13_r0SiO2D , double *F13_r0Al2O3 
      , int *gr_N, int *gr_Size, double *gr_dT, double *gr_Td
      , double *loc_kpFeM    , double *loc_kpMg2SiO4, double *loc_kpMgSiO3 , double *loc_kpFeS
      , double *loc_kpreforg , double *loc_kpvolorg , double *loc_kpH2Oice 
      , double *C30_kpSiM, double *C30_kpFeM, double *C30_kpMg2SiO4, double *C30_kpMgSiO3
      , double *C30_kpAC, double *C30_kpSiO2D, double *C30_kpMgO, double *C30_kpFeS, double *C30_kpAl2O3
      , double *F13_kpFeM, double *F13_kpMg2SiO4, double *F13_kpMgSiO3, double *F13_kpFe3O4
      , double *F13_kpAC , double *F13_kpSiO2D , double *F13_kpAl2O3
      , double *gasgr2a, double *gamma_isrf2a
      , gr_float *SiM_temp, gr_float *FeM_temp, gr_float *Mg2SiO4_temp, gr_float *MgSiO3_temp, gr_float *Fe3O4_temp
      , gr_float *AC_temp, gr_float *SiO2D_temp, gr_float *MgO_temp, gr_float *FeS_temp, gr_float *Al2O3_temp
      , gr_float *reforg_temp, gr_float *volorg_temp, gr_float *H2Oice_temp );

int local_calculate_dust_temperature(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *dust_temperature
                                   , gr_float *SiM_temperature
                                   , gr_float *FeM_temperature
                                   , gr_float *Mg2SiO4_temperature
                                   , gr_float *MgSiO3_temperature
                                   , gr_float *Fe3O4_temperature
                                   , gr_float *AC_temperature
                                   , gr_float *SiO2D_temperature
                                   , gr_float *MgO_temperature
                                   , gr_float *FeS_temperature
                                   , gr_float *Al2O3_temperature
                                   , gr_float *reforg_temperature
                                   , gr_float *volorg_temperature
                                   , gr_float *H2Oice_temperature )
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  if (my_chemistry->dust_chemistry < 1 && my_chemistry->h2_on_dust < 1)
    return SUCCESS;

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }
  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Compute the size of the fields. */
 
  int i, dim, size = 1;
  for (dim = 0; dim < my_fields->grid_rank; dim++)
    size *= my_fields->grid_dimension[dim];

  gr_float *temperature;
  temperature = malloc(size * sizeof(gr_float));
  if (local_calculate_temperature(my_chemistry, my_rates, my_units,
                                  my_fields, temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_temperature.\n");
    return FAIL;
  }

  /* Call the appropriate FORTRAN routine to do the work. */

  FORTRAN_NAME(calc_tdust_3d_g)(
       my_fields->density,
       my_fields->e_density,
       my_fields->HI_density,
       my_fields->HII_density,
       my_fields->HeI_density,
       my_fields->HeII_density,
       my_fields->HeIII_density,
       my_fields->HM_density,
       my_fields->H2I_density,
       my_fields->H2II_density,
       my_fields->grid_dimension,
       my_fields->grid_dimension+1,
       my_fields->grid_dimension+2,
       &my_chemistry->NumberOfTemperatureBins,
       &my_units->comoving_coordinates,
       &my_chemistry->primordial_chemistry,
       &(my_fields->grid_rank),
       my_fields->grid_start,
       my_fields->grid_start+1,
       my_fields->grid_start+2,
       my_fields->grid_end,
       my_fields->grid_end+1,
       my_fields->grid_end+2,
       &my_units->a_value,
       &my_chemistry->TemperatureStart,
       &my_chemistry->TemperatureEnd,
       &my_chemistry->local_dust_to_gas_ratio,
       my_rates->gas_grain,
       &my_rates->gamma_isrf,
       &my_chemistry->interstellar_radiation_field,
       &temperature_units,
       &co_length_units,
       &my_units->a_units,
       &co_density_units,
       &my_units->time_units,
       temperature,
       dust_temperature,
       &my_chemistry->use_isrf_field,
       my_fields->isrf_habing
     ,&metal_field_present
     ,&my_chemistry->use_dust_density_field
     ,&my_chemistry->grain_growth
     ,&my_chemistry->SolarMetalFractionByMass
     , my_fields->metal_density
     , my_fields->dust_density
     , my_fields->SiM_density
     , my_fields->FeM_density
     , my_fields->Mg2SiO4_density
     , my_fields->MgSiO3_density
     , my_fields->Fe3O4_density
     , my_fields->AC_density
     , my_fields->SiO2D_density
     , my_fields->MgO_density
     , my_fields->FeS_density
     , my_fields->Al2O3_density
     , my_fields->reforg_density
     , my_fields->volorg_density
     , my_fields->H2Oice_density
     ,&my_chemistry->multi_metals
     ,&my_chemistry->metal_pop3
     ,&my_chemistry->dust_species
     ,&my_chemistry->dust_temperature_species
     ,&my_chemistry->dust_sublimation
     , my_fields->metal_loc
     , my_fields->metal_C30
     , my_fields->metal_F13
     ,&my_chemistry->loc_fFeM    
     ,&my_chemistry->loc_fMg2SiO4
     ,&my_chemistry->loc_fMgSiO3 
     ,&my_chemistry->loc_fFeS
     ,&my_chemistry->loc_freforg 
     ,&my_chemistry->loc_fvolorg 
     ,&my_chemistry->loc_fH2Oice
     , my_chemistry->loc_r0FeM    
     , my_chemistry->loc_r0Mg2SiO4
     , my_chemistry->loc_r0MgSiO3 
     , my_chemistry->loc_r0FeS
     , my_chemistry->loc_r0reforg 
     , my_chemistry->loc_r0volorg 
     , my_chemistry->loc_r0H2Oice 
     ,&my_chemistry->C30_fSiM    
     ,&my_chemistry->C30_fFeM    
     ,&my_chemistry->C30_fMg2SiO4
     ,&my_chemistry->C30_fMgSiO3 
     ,&my_chemistry->C30_fAC
     ,&my_chemistry->C30_fSiO2D
     ,&my_chemistry->C30_fMgO
     ,&my_chemistry->C30_fFeS
     ,&my_chemistry->C30_fAl2O3
     , my_chemistry->C30_r0SiM    
     , my_chemistry->C30_r0FeM    
     , my_chemistry->C30_r0Mg2SiO4
     , my_chemistry->C30_r0MgSiO3 
     , my_chemistry->C30_r0AC
     , my_chemistry->C30_r0SiO2D
     , my_chemistry->C30_r0MgO
     , my_chemistry->C30_r0FeS
     , my_chemistry->C30_r0Al2O3
     ,&my_chemistry->F13_fFeM    
     ,&my_chemistry->F13_fMg2SiO4
     ,&my_chemistry->F13_fMgSiO3 
     ,&my_chemistry->F13_fFe3O4
     ,&my_chemistry->F13_fAC
     ,&my_chemistry->F13_fSiO2D
     ,&my_chemistry->F13_fAl2O3  
     , my_chemistry->F13_r0FeM    
     , my_chemistry->F13_r0Mg2SiO4
     , my_chemistry->F13_r0MgSiO3 
     , my_chemistry->F13_r0Fe3O4
     , my_chemistry->F13_r0AC
     , my_chemistry->F13_r0SiO2D
     , my_chemistry->F13_r0Al2O3  
     , my_rates->gr_N
     ,&my_rates->gr_Size
     ,&my_rates->gr_dT
     , my_rates->gr_Td
     , my_rates->loc_kpFeM    
     , my_rates->loc_kpMg2SiO4
     , my_rates->loc_kpMgSiO3 
     , my_rates->loc_kpFeS
     , my_rates->loc_kpreforg 
     , my_rates->loc_kpvolorg 
     , my_rates->loc_kpH2Oice 
     , my_rates->C30_kpSiM    
     , my_rates->C30_kpFeM    
     , my_rates->C30_kpMg2SiO4 
     , my_rates->C30_kpMgSiO3
     , my_rates->C30_kpAC     
     , my_rates->C30_kpSiO2D
     , my_rates->C30_kpMgO    
     , my_rates->C30_kpFeS    
     , my_rates->C30_kpAl2O3  
     , my_rates->F13_kpFeM    
     , my_rates->F13_kpMg2SiO4
     , my_rates->F13_kpMgSiO3 
     , my_rates->F13_kpFe3O4
     , my_rates->F13_kpAC     
     , my_rates->F13_kpSiO2D  
     , my_rates->F13_kpAl2O3  
     , my_rates->gas_grain2
     ,&my_rates->gamma_isrf2
     , SiM_temperature
     , FeM_temperature
     , Mg2SiO4_temperature
     , MgSiO3_temperature
     , Fe3O4_temperature
     , AC_temperature
     , SiO2D_temperature
     , MgO_temperature
     , FeS_temperature
     , Al2O3_temperature
     , reforg_temperature
     , volorg_temperature
     , H2Oice_temperature
   );

  free(temperature);

  return SUCCESS;
}

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature
                             , gr_float *SiM_temperature
                             , gr_float *FeM_temperature
                             , gr_float *Mg2SiO4_temperature
                             , gr_float *MgSiO3_temperature
                             , gr_float *Fe3O4_temperature
                             , gr_float *AC_temperature
                             , gr_float *SiO2D_temperature
                             , gr_float *MgO_temperature
                             , gr_float *FeS_temperature
                             , gr_float *Al2O3_temperature
                             , gr_float *reforg_temperature
                             , gr_float *volorg_temperature
                             , gr_float *H2Oice_temperature )
{
  if (local_calculate_dust_temperature(
          grackle_data, &grackle_rates, my_units,
          my_fields, dust_temperature
                   , SiM_temperature
                   , FeM_temperature
                   , Mg2SiO4_temperature
                   , MgSiO3_temperature
                   , Fe3O4_temperature
                   , AC_temperature
                   , SiO2D_temperature
                   , MgO_temperature
                   , FeS_temperature
                   , Al2O3_temperature
                   , reforg_temperature
                   , volorg_temperature
                   , H2Oice_temperature ) == FAIL) {
    fprintf(stderr, "Error in local_calculate_dust_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
