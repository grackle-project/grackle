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

double get_temperature_units(code_units *my_units);

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
      , int *imetal, int *imchem, int *idustfield, int *igrgr
      , double *z_solar, gr_float *metal, gr_float *dust
      , gr_float *SiM, gr_float *FeM, gr_float *Mg2SiO4, gr_float *MgSiO3, gr_float *Fe3O4
      , gr_float *AC, gr_float *SiO2D, gr_float *MgO, gr_float *FeS, gr_float *Al2O3
      , gr_float *reforg, gr_float *volorg, gr_float *H2Oice
      , int *immulti, int *imabund, int *idspecies, int *itdmulti, int *idsub
      , gr_float *metal_loc, gr_float *metal_C13, gr_float *metal_C20, gr_float *metal_C25, gr_float *metal_C30
      , gr_float *metal_F13, gr_float *metal_F15, gr_float *metal_F50, gr_float *metal_F80
      , gr_float *metal_P170, gr_float *metal_P200, gr_float *metal_Y19
      , int *SN0_N
      , double *SN0_fSiM, double *SN0_fFeM, double *SN0_fMg2SiO4, double *SN0_fMgSiO3, double *SN0_fFe3O4
      , double *SN0_fAC, double *SN0_fSiO2D, double *SN0_fMgO, double *SN0_fFeS, double *SN0_fAl2O3
      , double *SN0_freforg , double *SN0_fvolorg , double *SN0_fH2Oice
      , double *SN0_r0SiM, double *SN0_r0FeM, double *SN0_r0Mg2SiO4, double *SN0_r0MgSiO3, double *SN0_r0Fe3O4
      , double *SN0_r0AC, double *SN0_r0SiO2D, double *SN0_r0MgO, double *SN0_r0FeS, double *SN0_r0Al2O3
      , double *SN0_r0reforg , double *SN0_r0volorg , double *SN0_r0H2Oice
      , int *gr_N, int *gr_Size, double *gr_dT, double *gr_Td
      , double *SN0_kpSiM, double *SN0_kpFeM, double *SN0_kpMg2SiO4, double *SN0_kpMgSiO3, double *SN0_kpFe3O4
      , double *SN0_kpAC, double *SN0_kpSiO2D, double *SN0_kpMgO, double *SN0_kpFeS, double *SN0_kpAl2O3
      , double *SN0_kpreforg , double *SN0_kpvolorg , double *SN0_kpH2Oice
      , double *gasgr2a, double *gamma_isrf2a
      , gr_float *SiM_temp, gr_float *FeM_temp, gr_float *Mg2SiO4_temp, gr_float *MgSiO3_temp, gr_float *Fe3O4_temp
      , gr_float *AC_temp, gr_float *SiO2D_temp, gr_float *MgO_temp, gr_float *FeS_temp, gr_float *Al2O3_temp
      , gr_float *reforg_temp, gr_float *volorg_temp, gr_float *H2Oice_temp );

int local_calculate_dust_temperature(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *dust_temperature)
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
  double temperature_units = get_temperature_units(my_units);

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
     ,&my_chemistry->metal_chemistry
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
     ,&my_chemistry->metal_abundances
     ,&my_chemistry->dust_species
     ,&my_chemistry->dust_temperature_multi
     ,&my_chemistry->dust_sublimation
     , my_fields->local_ISM_metal_density
     , my_fields->ccsn13_metal_density
     , my_fields->ccsn20_metal_density
     , my_fields->ccsn25_metal_density
     , my_fields->ccsn30_metal_density
     , my_fields->fsn13_metal_density
     , my_fields->fsn15_metal_density
     , my_fields->fsn50_metal_density
     , my_fields->fsn80_metal_density
     , my_fields->pisn170_metal_density
     , my_fields->pisn200_metal_density
     , my_fields->y19_metal_density
     ,&my_rates->SN0_N
     , my_rates->SN0_fSiM    
     , my_rates->SN0_fFeM    
     , my_rates->SN0_fMg2SiO4
     , my_rates->SN0_fMgSiO3 
     , my_rates->SN0_fFe3O4
     , my_rates->SN0_fAC
     , my_rates->SN0_fSiO2D
     , my_rates->SN0_fMgO
     , my_rates->SN0_fFeS
     , my_rates->SN0_fAl2O3  
     , my_rates->SN0_freforg 
     , my_rates->SN0_fvolorg 
     , my_rates->SN0_fH2Oice
     , my_rates->SN0_r0SiM    
     , my_rates->SN0_r0FeM    
     , my_rates->SN0_r0Mg2SiO4
     , my_rates->SN0_r0MgSiO3 
     , my_rates->SN0_r0Fe3O4
     , my_rates->SN0_r0AC
     , my_rates->SN0_r0SiO2D
     , my_rates->SN0_r0MgO
     , my_rates->SN0_r0FeS
     , my_rates->SN0_r0Al2O3  
     , my_rates->SN0_r0reforg 
     , my_rates->SN0_r0volorg 
     , my_rates->SN0_r0H2Oice
     , my_rates->gr_N
     ,&my_rates->gr_Size
     ,&my_rates->gr_dT
     , my_rates->gr_Td
     , my_rates->SN0_kpSiM    
     , my_rates->SN0_kpFeM    
     , my_rates->SN0_kpMg2SiO4
     , my_rates->SN0_kpMgSiO3 
     , my_rates->SN0_kpFe3O4
     , my_rates->SN0_kpAC
     , my_rates->SN0_kpSiO2D
     , my_rates->SN0_kpMgO
     , my_rates->SN0_kpFeS
     , my_rates->SN0_kpAl2O3  
     , my_rates->SN0_kpreforg 
     , my_rates->SN0_kpvolorg 
     , my_rates->SN0_kpH2Oice
     , my_rates->gas_grain2
     ,&my_rates->gamma_isrf2
     , my_fields->SiM_temperature
     , my_fields->FeM_temperature
     , my_fields->Mg2SiO4_temperature
     , my_fields->MgSiO3_temperature
     , my_fields->Fe3O4_temperature
     , my_fields->AC_temperature
     , my_fields->SiO2D_temperature
     , my_fields->MgO_temperature
     , my_fields->FeS_temperature
     , my_fields->Al2O3_temperature
     , my_fields->reforg_temperature
     , my_fields->volorg_temperature
     , my_fields->H2Oice_temperature
   );

  free(temperature);

  return SUCCESS;
}

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature)
{
  if (local_calculate_dust_temperature(
          grackle_data, &grackle_rates, my_units,
          my_fields, dust_temperature) == FAIL) {
    fprintf(stderr, "Error in local_calculate_dust_temperature.\n");
    return FAIL;
  }
  return SUCCESS;
}
