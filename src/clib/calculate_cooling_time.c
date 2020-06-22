/***********************************************************************
/
/ Calculate cooling time field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
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

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units);
 
extern void FORTRAN_NAME(cool_multi_time_g)(
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	gr_float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idustall,
        int *idustfield, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, 
        int *ih2co, int *ipiht, int *igammah,
	double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *z_solar, double *fgr,
	double *ceHIa, double *ceHeIa, double *ceHeIIa, 
        double *ciHIa, double *ciHeIa,
	double *ciHeISa, double *ciHeIIa, double *reHIIa, double *reHeII1a,
	double *reHeII2a, double *reHeIIIa, double *brema, double *compa, 
        double *gammaha, double *isrf, double *regra, double *gamma_isrfa,
        double *comp_xraya, double *comp_temp,
        double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II, 
        gr_float *DI, gr_float *DII, gr_float *HDI,
        gr_float *metal, gr_float *dust,
	double *hyd01ka, double *h2k01a, double *vibha, 
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra,
        int *iradshield, double *avgsighi, double *avgsighei, double *avgsigheii,
        double *k24, double *k26,
        int *iradtrans, double *photogamma,
	int *ih2optical, int *iciecool, double *ciecoa,
 	int *icmbTfloor, int *iClHeat, double *clEleFra,
        long long *priGridRank, long long *priGridDim,
 	double *priPar1, double *priPar2, double *priPar3, 
        double *priPar4, double *priPar5,
 	long long *priDataSize, double *priCooling,
        double *priHeating, double *priMMW,
        long long *metGridRank, long long *metGridDim,
 	double *metPar1, double *metPar2, double *metPar3, 
        double *metPar4, double *metPar5,
 	long long *metDataSize, double *metCooling,
        double *metHeating, int *clnew,
        int *iVheat, int *iMheat, gr_float *Vheat, gr_float *Mheat,
        int *iisrffield, gr_float* isrf_habing
#ifdef GRACKLE_MD
      , int *imchem, int *igrgr
      , gr_float *DM, gr_float *HDII, gr_float *HeHII
      , gr_float *CI, gr_float *CII, gr_float *CO, gr_float *CO2
      , gr_float *OI, gr_float *OH, gr_float *H2O, gr_float *O2
      , gr_float *SiI, gr_float *SiOI, gr_float *SiO2I
      , gr_float *CH, gr_float *CH2, gr_float *COII, gr_float *OII
      , gr_float *OHII, gr_float *H2OII, gr_float *H3OII, gr_float *O2II
      , gr_float *Mg, gr_float *Al, gr_float *S, gr_float *Fe
      , gr_float *SiM, gr_float *FeM, gr_float *Mg2SiO4, gr_float *MgSiO3, gr_float *Fe3O4
      , gr_float *AC, gr_float *SiO2D, gr_float *MgO, gr_float *FeS, gr_float *Al2O3
      , double *cieY06
      , int *LH2_N, int *LH2_Size
      , double *LH2_D, double *LH2_T, double *LH2_H
      , double *LH2_dD, double *LH2_dT, double *LH2_dH, double *LH2_L
      , int *LHD_N, int *LHD_Size
      , double *LHD_D, double *LHD_T, double *LHD_H
      , double *LHD_dD, double *LHD_dT, double *LHD_dH, double *LHD_L
      , int *LCI_N, int *LCI_Size
      , double *LCI_D, double *LCI_T, double *LCI_H
      , double *LCI_dD, double *LCI_dT, double *LCI_dH, double *LCI_L
      , int *LCII_N, int *LCII_Size
      , double *LCII_D, double *LCII_T, double *LCII_H
      , double *LCII_dD, double *LCII_dT, double *LCII_dH, double *LCII_L
      , int *LOI_N, int *LOI_Size
      , double *LOI_D, double *LOI_T, double *LOI_H
      , double *LOI_dD, double *LOI_dT, double *LOI_dH, double *LOI_L
      , int *LCO_N, int *LCO_Size
      , double *LCO_D, double *LCO_T, double *LCO_H
      , double *LCO_dD, double *LCO_dT, double *LCO_dH, double *LCO_L
      , int *LOH_N, int *LOH_Size
      , double *LOH_D, double *LOH_T, double *LOH_H
      , double *LOH_dD, double *LOH_dT, double *LOH_dH, double *LOH_L
      , int *LH2O_N, int *LH2O_Size
      , double *LH2O_D, double *LH2O_T, double *LH2O_H
      , double *LH2O_dD, double *LH2O_dT, double *LH2O_dH, double *LH2O_L
      , int *alphap_N, int *alphap_Size
      , double *alphap_D, double *alphap_T, double *alphap_dD, double *alphap_dT
      , double *alphap_Data
      , int *grain_N, int *grain_Size
      , double *grain_D, double *grain_T, double *grain_dD, double *grain_dT
      , double *Tgrain, double *Ograin, double *Lgrain
      , int *impop3
      , double *C30_XC , double *C30_XO, double *C30_XMg, double *C30_XAl
      , double *C30_XSi, double *C30_XS, double *C30_XFe
      , double *C30_fC , double *C30_fO, double *C30_fMg, double *C30_fAl
      , double *C30_fSi, double *C30_fS, double *C30_fFe
      , double *C30_fSiM, double *C30_fFeM, double *C30_fMg2SiO4, double *C30_fMgSiO3
      , double *C30_fFe3O4, double *C30_fAC, double *C30_fSiO2D, double *C30_fMgO
      , double *C30_fFeS, double *C30_fAl2O3
      , double *C30_r0SiM, double *C30_r0FeM, double *C30_r0Mg2SiO4, double *C30_r0MgSiO3
      , double *C30_r0Fe3O4, double *C30_r0AC, double *C30_r0SiO2D, double *C30_r0MgO
      , double *C30_r0FeS, double *C30_r0Al2O3
      , double *C30_a0SiM, double *C30_a0FeM, double *C30_a0Mg2SiO4, double *C30_a0MgSiO3
      , double *C30_a0Fe3O4, double *C30_a0AC, double *C30_a0SiO2D, double *C30_a0MgO
      , double *C30_a0FeS, double *C30_a0Al2O3
      , double *C30_v0SiM, double *C30_v0FeM, double *C30_v0Mg2SiO4, double *C30_v0MgSiO3
      , double *C30_v0Fe3O4, double *C30_v0AC, double *C30_v0SiO2D, double *C30_v0MgO
      , double *C30_v0FeS, double *C30_v0Al2O3
      , int *C30_N, int *C30_Size
      , double *C30_D, double *C30_T, double *C30_dD, double *C30_dT
      , double *C30_RSiM, double *C30_RFeM, double *C30_RMg2SiO4, double *C30_RMgSiO3
      , double *C30_RFe3O4, double *C30_RAC, double *C30_RSiO2D, double *C30_RMgO
      , double *C30_RFeS, double *C30_RAl2O3
      , double *C30_dRSiM, double *C30_dRFeM, double *C30_dRMg2SiO4, double *C30_dRMgSiO3
      , double *C30_dRFe3O4, double *C30_dRAC, double *C30_dRSiO2D, double *C30_dRMgO
      , double *C30_dRFeS, double *C30_dRAl2O3
      , double *OC30_SiM, double *OC30_FeM, double *OC30_Mg2SiO4, double *OC30_MgSiO3
      , double *OC30_Fe3O4, double *OC30_AC, double *OC30_SiO2D, double *OC30_MgO
      , double *OC30_FeS, double *OC30_Al2O3
      , double *LC30_SiM, double *LC30_FeM, double *LC30_Mg2SiO4, double *LC30_MgSiO3
      , double *LC30_Fe3O4, double *LC30_AC, double *LC30_SiO2D, double *LC30_MgO
      , double *LC30_FeS, double *LC30_Al2O3
      , double *F13_XC , double *F13_XO, double *F13_XMg, double *F13_XAl
      , double *F13_XSi, double *F13_XS, double *F13_XFe
      , double *F13_fC , double *F13_fO, double *F13_fMg, double *F13_fAl
      , double *F13_fSi, double *F13_fS, double *F13_fFe
      , double *F13_fSiM, double *F13_fFeM, double *F13_fMg2SiO4, double *F13_fMgSiO3
      , double *F13_fFe3O4, double *F13_fAC, double *F13_fSiO2D, double *F13_fMgO
      , double *F13_fFeS, double *F13_fAl2O3
      , double *F13_r0SiM, double *F13_r0FeM, double *F13_r0Mg2SiO4, double *F13_r0MgSiO3
      , double *F13_r0Fe3O4, double *F13_r0AC, double *F13_r0SiO2D, double *F13_r0MgO
      , double *F13_r0FeS, double *F13_r0Al2O3
      , double *F13_a0SiM, double *F13_a0FeM, double *F13_a0Mg2SiO4, double *F13_a0MgSiO3
      , double *F13_a0Fe3O4, double *F13_a0AC, double *F13_a0SiO2D, double *F13_a0MgO
      , double *F13_a0FeS, double *F13_a0Al2O3
      , double *F13_v0SiM, double *F13_v0FeM, double *F13_v0Mg2SiO4, double *F13_v0MgSiO3
      , double *F13_v0Fe3O4, double *F13_v0AC, double *F13_v0SiO2D, double *F13_v0MgO
      , double *F13_v0FeS, double *F13_v0Al2O3
      , int *F13_N, int *F13_Size
      , double *F13_D, double *F13_T, double *F13_dD, double *F13_dT
      , double *F13_RSiM, double *F13_RFeM, double *F13_RMg2SiO4, double *F13_RMgSiO3
      , double *F13_RFe3O4, double *F13_RAC, double *F13_RSiO2D, double *F13_RMgO
      , double *F13_RFeS, double *F13_RAl2O3
      , double *F13_dRSiM, double *F13_dRFeM, double *F13_dRMg2SiO4, double *F13_dRMgSiO3
      , double *F13_dRFe3O4, double *F13_dRAC, double *F13_dRSiO2D, double *F13_dRMgO
      , double *F13_dRFeS, double *F13_dRAl2O3
      , double *OF13_SiM, double *OF13_FeM, double *OF13_Mg2SiO4, double *OF13_MgSiO3
      , double *OF13_Fe3O4, double *OF13_AC, double *OF13_SiO2D, double *OF13_MgO
      , double *OF13_FeS, double *OF13_Al2O3
      , double *LF13_SiM, double *LF13_FeM, double *LF13_Mg2SiO4, double *LF13_MgSiO3
      , double *LF13_Fe3O4, double *LF13_AC, double *LF13_SiO2D, double *LF13_MgO
      , double *LF13_FeS, double *LF13_Al2O3
#endif
        );

int local_calculate_cooling_time(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units,
                                 grackle_field_data *my_fields,
                                 gr_float *cooling_time)
{
 
  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */
  photo_rate_storage my_uvb_rates;

  my_uvb_rates.k24 = my_uvb_rates.k25 = my_uvb_rates.k26 =
    my_uvb_rates.k27 = my_uvb_rates.k28 = my_uvb_rates.k29 =
    my_uvb_rates.k30 = my_uvb_rates.k31 = my_uvb_rates.piHI =
    my_uvb_rates.piHeI = my_uvb_rates.piHeII = my_uvb_rates.crsHI =
    my_uvb_rates.crsHeI = my_uvb_rates.crsHeII =
    my_uvb_rates.comp_xray = my_uvb_rates.temp_xray = 0.;

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_rates,
                                  &my_uvb_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return FAIL;
    }
  }
  else {
    my_uvb_rates.k24       = my_rates->k24;
    my_uvb_rates.k25       = my_rates->k25;
    my_uvb_rates.k26       = my_rates->k26;
    my_uvb_rates.k27       = my_rates->k27;
    my_uvb_rates.k28       = my_rates->k28;
    my_uvb_rates.k29       = my_rates->k29;
    my_uvb_rates.k30       = my_rates->k30;
    my_uvb_rates.k31       = my_rates->k31;
    my_uvb_rates.piHI      = my_rates->piHI;
    my_uvb_rates.piHeI     = my_rates->piHeI;
    my_uvb_rates.piHeII    = my_rates->piHeII;
    my_uvb_rates.crsHI     = my_rates->crsHI;
    my_uvb_rates.crsHeI    = my_rates->crsHeI;
    my_uvb_rates.crsHeII   = my_rates->crsHeII;
    my_uvb_rates.comp_xray = my_rates->comp_xray;
    my_uvb_rates.temp_xray = my_rates->temp_xray;
  }

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

  /* Calculate temperature units. */

  double temperature_units = mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Call the fortran routine to solve cooling equations. */

    FORTRAN_NAME(cool_multi_time_g)(
       my_fields->density,
       my_fields->internal_energy,
       my_fields->x_velocity,
       my_fields->y_velocity,
       my_fields->z_velocity,
       my_fields->e_density,
       my_fields->HI_density,
       my_fields->HII_density,
       my_fields->HeI_density,
       my_fields->HeII_density,
       my_fields->HeIII_density,
       cooling_time,
       my_fields->grid_dimension,
       my_fields->grid_dimension+1,
       my_fields->grid_dimension+2,
       &my_chemistry->NumberOfTemperatureBins,
       &my_units->comoving_coordinates,
       &my_chemistry->primordial_chemistry,
       &metal_field_present,
       &my_chemistry->metal_cooling,
       &my_chemistry->h2_on_dust,
       &my_chemistry->dust_chemistry,
       &my_chemistry->use_dust_density_field,
       &(my_fields->grid_rank),
       my_fields->grid_start,
       my_fields->grid_start+1,
       my_fields->grid_start+2,
       my_fields->grid_end,
       my_fields->grid_end+1,
       my_fields->grid_end+2,
       &my_chemistry->ih2co,
       &my_chemistry->ipiht,
       &my_chemistry->photoelectric_heating,
       &my_units->a_value,
       &my_chemistry->TemperatureStart,
       &my_chemistry->TemperatureEnd,
       &temperature_units,
       &co_length_units,
       &my_units->a_units,
       &co_density_units,
       &my_units->time_units,
       &my_chemistry->Gamma,
       &my_chemistry->HydrogenFractionByMass,
       &my_chemistry->SolarMetalFractionByMass,
       &my_chemistry->local_dust_to_gas_ratio,
       my_rates->ceHI,
       my_rates->ceHeI,
       my_rates->ceHeII,
       my_rates->ciHI,
       my_rates->ciHeI,
       my_rates->ciHeIS,
       my_rates->ciHeII,
       my_rates->reHII,
       my_rates->reHeII1,
       my_rates->reHeII2,
       my_rates->reHeIII,
       my_rates->brem,
       &my_rates->comp,
       &my_rates->gammah,
       &my_chemistry->interstellar_radiation_field,
       my_rates->regr,
       &my_rates->gamma_isrf,
       &my_uvb_rates.comp_xray,
       &my_uvb_rates.temp_xray,
       &my_uvb_rates.piHI,
       &my_uvb_rates.piHeI,
       &my_uvb_rates.piHeII,
       my_fields->HM_density,
       my_fields->H2I_density,
       my_fields->H2II_density,
       my_fields->DI_density,
       my_fields->DII_density,
       my_fields->HDI_density,
       my_fields->metal_density,
       my_fields->dust_density,
       my_rates->hyd01k,
       my_rates->h2k01,
       my_rates->vibh,
       my_rates->roth,
       my_rates->rotl,
       my_rates->GP99LowDensityLimit,
       my_rates->GP99HighDensityLimit,
       my_rates->HDlte,
       my_rates->HDlow,
       my_rates->GAHI,
       my_rates->GAH2,
       my_rates->GAHe,
       my_rates->GAHp,
       my_rates->GAel,
       my_rates->H2LTE,
       my_rates->gas_grain,
       &my_chemistry->self_shielding_method,
       &my_uvb_rates.crsHI,
       &my_uvb_rates.crsHeI,
       &my_uvb_rates.crsHeII,
       &my_uvb_rates.k24,
       &my_uvb_rates.k26,
       &my_chemistry->use_radiative_transfer,
       my_fields->RT_heating_rate,
       &my_chemistry->h2_optical_depth_approximation,
       &my_chemistry->cie_cooling, my_rates->cieco,
       &my_chemistry->cmb_temperature_floor,
       &my_chemistry->UVbackground,
       &my_chemistry->cloudy_electron_fraction_factor,
       &my_rates->cloudy_primordial.grid_rank,
       my_rates->cloudy_primordial.grid_dimension,
       my_rates->cloudy_primordial.grid_parameters[0],
       my_rates->cloudy_primordial.grid_parameters[1],
       my_rates->cloudy_primordial.grid_parameters[2],
       my_rates->cloudy_primordial.grid_parameters[3],
       my_rates->cloudy_primordial.grid_parameters[4],
       &my_rates->cloudy_primordial.data_size,
       my_rates->cloudy_primordial.cooling_data, 
       my_rates->cloudy_primordial.heating_data,
       my_rates->cloudy_primordial.mmw_data,
       &my_rates->cloudy_metal.grid_rank,
       my_rates->cloudy_metal.grid_dimension,
       my_rates->cloudy_metal.grid_parameters[0],
       my_rates->cloudy_metal.grid_parameters[1],
       my_rates->cloudy_metal.grid_parameters[2],
       my_rates->cloudy_metal.grid_parameters[3],
       my_rates->cloudy_metal.grid_parameters[4],
       &my_rates->cloudy_metal.data_size,
       my_rates->cloudy_metal.cooling_data,
       my_rates->cloudy_metal.heating_data,
       &my_rates->cloudy_data_new,
       &my_chemistry->use_volumetric_heating_rate,
       &my_chemistry->use_specific_heating_rate,
       my_fields->volumetric_heating_rate,
       my_fields->specific_heating_rate,
       &my_chemistry->use_isrf_field,
       my_fields->isrf_habing
#ifdef GRACKLE_MD
     ,&my_chemistry->metal_chemistry
     ,&my_chemistry->grain_growth
     , my_fields->DM_density
     , my_fields->HDII_density
     , my_fields->HeHII_density
     , my_fields->CI_density
     , my_fields->CII_density
     , my_fields->CO_density
     , my_fields->CO2_density
     , my_fields->OI_density
     , my_fields->OH_density
     , my_fields->H2O_density
     , my_fields->O2_density
     , my_fields->SiI_density
     , my_fields->SiOI_density
     , my_fields->SiO2I_density
     , my_fields->CH_density
     , my_fields->CH2_density
     , my_fields->COII_density
     , my_fields->OII_density
     , my_fields->OHII_density
     , my_fields->H2OII_density
     , my_fields->H3OII_density
     , my_fields->O2II_density
     , my_fields->Mg_density
     , my_fields->Al_density
     , my_fields->S_density
     , my_fields->Fe_density
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
     , my_rates->cieY06
     , my_rates->LH2_N
     ,&my_rates->LH2_Size
     , my_rates->LH2_D
     , my_rates->LH2_T
     , my_rates->LH2_H
     ,&my_rates->LH2_dD
     ,&my_rates->LH2_dT
     ,&my_rates->LH2_dH
     , my_rates->LH2_L
     , my_rates->LHD_N
     ,&my_rates->LHD_Size
     , my_rates->LHD_D
     , my_rates->LHD_T
     , my_rates->LHD_H
     ,&my_rates->LHD_dD
     ,&my_rates->LHD_dT
     ,&my_rates->LHD_dH
     , my_rates->LHD_L
     , my_rates->LCI_N
     ,&my_rates->LCI_Size
     , my_rates->LCI_D
     , my_rates->LCI_T
     , my_rates->LCI_H
     ,&my_rates->LCI_dD
     ,&my_rates->LCI_dT
     ,&my_rates->LCI_dH
     , my_rates->LCI_L
     , my_rates->LCII_N
     ,&my_rates->LCII_Size
     , my_rates->LCII_D
     , my_rates->LCII_T
     , my_rates->LCII_H
     ,&my_rates->LCII_dD
     ,&my_rates->LCII_dT
     ,&my_rates->LCII_dH
     , my_rates->LCII_L
     , my_rates->LOI_N
     ,&my_rates->LOI_Size
     , my_rates->LOI_D
     , my_rates->LOI_T
     , my_rates->LOI_H
     ,&my_rates->LOI_dD
     ,&my_rates->LOI_dT
     ,&my_rates->LOI_dH
     , my_rates->LOI_L
     , my_rates->LCO_N
     ,&my_rates->LCO_Size
     , my_rates->LCO_D
     , my_rates->LCO_T
     , my_rates->LCO_H
     ,&my_rates->LCO_dD
     ,&my_rates->LCO_dT
     ,&my_rates->LCO_dH
     , my_rates->LCO_L
     , my_rates->LOH_N
     ,&my_rates->LOH_Size
     , my_rates->LOH_D
     , my_rates->LOH_T
     , my_rates->LOH_H
     ,&my_rates->LOH_dD
     ,&my_rates->LOH_dT
     ,&my_rates->LOH_dH
     , my_rates->LOH_L
     , my_rates->LH2O_N
     ,&my_rates->LH2O_Size
     , my_rates->LH2O_D
     , my_rates->LH2O_T
     , my_rates->LH2O_H
     ,&my_rates->LH2O_dD
     ,&my_rates->LH2O_dT
     ,&my_rates->LH2O_dH
     , my_rates->LH2O_L
     , my_rates->alphap_N
     ,&my_rates->alphap_Size
     , my_rates->alphap_D
     , my_rates->alphap_T
     ,&my_rates->alphap_dD
     ,&my_rates->alphap_dT
     , my_rates->alphap_Data
     , my_rates->grain_N
     ,&my_rates->grain_Size
     , my_rates->grain_D
     , my_rates->grain_T
     ,&my_rates->grain_dD
     ,&my_rates->grain_dT
     , my_rates->Tgrain
     , my_rates->Ograin
     , my_rates->Lgrain
     ,&my_chemistry->metal_pop3
     ,&my_chemistry->C30_XC 
     ,&my_chemistry->C30_XO 
     ,&my_chemistry->C30_XMg
     ,&my_chemistry->C30_XAl
     ,&my_chemistry->C30_XSi
     ,&my_chemistry->C30_XS 
     ,&my_chemistry->C30_XFe
     ,&my_chemistry->C30_fC 
     ,&my_chemistry->C30_fO 
     ,&my_chemistry->C30_fMg
     ,&my_chemistry->C30_fAl
     ,&my_chemistry->C30_fSi
     ,&my_chemistry->C30_fS 
     ,&my_chemistry->C30_fFe
     ,&my_chemistry->C30_fSiM    
     ,&my_chemistry->C30_fFeM    
     ,&my_chemistry->C30_fMg2SiO4
     ,&my_chemistry->C30_fMgSiO3 
     ,&my_chemistry->C30_fFe3O4  
     ,&my_chemistry->C30_fAC     
     ,&my_chemistry->C30_fSiO2D  
     ,&my_chemistry->C30_fMgO    
     ,&my_chemistry->C30_fFeS    
     ,&my_chemistry->C30_fAl2O3  
     ,&my_rates->C30_r0SiM    
     ,&my_rates->C30_r0FeM    
     ,&my_rates->C30_r0Mg2SiO4
     ,&my_rates->C30_r0MgSiO3 
     ,&my_rates->C30_r0Fe3O4  
     ,&my_rates->C30_r0AC     
     ,&my_rates->C30_r0SiO2D  
     ,&my_rates->C30_r0MgO    
     ,&my_rates->C30_r0FeS    
     ,&my_rates->C30_r0Al2O3  
     ,&my_rates->C30_a0SiM    
     ,&my_rates->C30_a0FeM    
     ,&my_rates->C30_a0Mg2SiO4
     ,&my_rates->C30_a0MgSiO3 
     ,&my_rates->C30_a0Fe3O4  
     ,&my_rates->C30_a0AC     
     ,&my_rates->C30_a0SiO2D  
     ,&my_rates->C30_a0MgO    
     ,&my_rates->C30_a0FeS    
     ,&my_rates->C30_a0Al2O3  
     ,&my_rates->C30_v0SiM    
     ,&my_rates->C30_v0FeM    
     ,&my_rates->C30_v0Mg2SiO4
     ,&my_rates->C30_v0MgSiO3 
     ,&my_rates->C30_v0Fe3O4  
     ,&my_rates->C30_v0AC     
     ,&my_rates->C30_v0SiO2D  
     ,&my_rates->C30_v0MgO    
     ,&my_rates->C30_v0FeS    
     ,&my_rates->C30_v0Al2O3  
     , my_rates->C30_N
     ,&my_rates->C30_Size
     , my_rates->C30_D
     , my_rates->C30_T
     ,&my_rates->C30_dD
     ,&my_rates->C30_dT
     , my_rates->C30_RSiM    
     , my_rates->C30_RFeM    
     , my_rates->C30_RMg2SiO4
     , my_rates->C30_RMgSiO3 
     , my_rates->C30_RFe3O4  
     , my_rates->C30_RAC     
     , my_rates->C30_RSiO2D  
     , my_rates->C30_RMgO    
     , my_rates->C30_RFeS    
     , my_rates->C30_RAl2O3  
     ,&my_rates->C30_dRSiM    
     ,&my_rates->C30_dRFeM    
     ,&my_rates->C30_dRMg2SiO4
     ,&my_rates->C30_dRMgSiO3 
     ,&my_rates->C30_dRFe3O4  
     ,&my_rates->C30_dRAC     
     ,&my_rates->C30_dRSiO2D  
     ,&my_rates->C30_dRMgO    
     ,&my_rates->C30_dRFeS    
     ,&my_rates->C30_dRAl2O3  
     , my_rates->OC30_SiM    
     , my_rates->OC30_FeM    
     , my_rates->OC30_Mg2SiO4
     , my_rates->OC30_MgSiO3 
     , my_rates->OC30_Fe3O4  
     , my_rates->OC30_AC     
     , my_rates->OC30_SiO2D  
     , my_rates->OC30_MgO    
     , my_rates->OC30_FeS    
     , my_rates->OC30_Al2O3  
     , my_rates->LC30_SiM    
     , my_rates->LC30_FeM    
     , my_rates->LC30_Mg2SiO4
     , my_rates->LC30_MgSiO3 
     , my_rates->LC30_Fe3O4  
     , my_rates->LC30_AC     
     , my_rates->LC30_SiO2D  
     , my_rates->LC30_MgO    
     , my_rates->LC30_FeS    
     , my_rates->LC30_Al2O3  
     ,&my_chemistry->F13_XC 
     ,&my_chemistry->F13_XO 
     ,&my_chemistry->F13_XMg
     ,&my_chemistry->F13_XAl
     ,&my_chemistry->F13_XSi
     ,&my_chemistry->F13_XS 
     ,&my_chemistry->F13_XFe
     ,&my_chemistry->F13_fC 
     ,&my_chemistry->F13_fO 
     ,&my_chemistry->F13_fMg
     ,&my_chemistry->F13_fAl
     ,&my_chemistry->F13_fSi
     ,&my_chemistry->F13_fS 
     ,&my_chemistry->F13_fFe
     ,&my_chemistry->F13_fSiM    
     ,&my_chemistry->F13_fFeM    
     ,&my_chemistry->F13_fMg2SiO4
     ,&my_chemistry->F13_fMgSiO3 
     ,&my_chemistry->F13_fFe3O4  
     ,&my_chemistry->F13_fAC     
     ,&my_chemistry->F13_fSiO2D  
     ,&my_chemistry->F13_fMgO    
     ,&my_chemistry->F13_fFeS    
     ,&my_chemistry->F13_fAl2O3  
     ,&my_rates->F13_r0SiM    
     ,&my_rates->F13_r0FeM    
     ,&my_rates->F13_r0Mg2SiO4
     ,&my_rates->F13_r0MgSiO3 
     ,&my_rates->F13_r0Fe3O4  
     ,&my_rates->F13_r0AC     
     ,&my_rates->F13_r0SiO2D  
     ,&my_rates->F13_r0MgO    
     ,&my_rates->F13_r0FeS    
     ,&my_rates->F13_r0Al2O3  
     ,&my_rates->F13_a0SiM    
     ,&my_rates->F13_a0FeM    
     ,&my_rates->F13_a0Mg2SiO4
     ,&my_rates->F13_a0MgSiO3 
     ,&my_rates->F13_a0Fe3O4  
     ,&my_rates->F13_a0AC     
     ,&my_rates->F13_a0SiO2D  
     ,&my_rates->F13_a0MgO    
     ,&my_rates->F13_a0FeS    
     ,&my_rates->F13_a0Al2O3  
     ,&my_rates->F13_v0SiM    
     ,&my_rates->F13_v0FeM    
     ,&my_rates->F13_v0Mg2SiO4
     ,&my_rates->F13_v0MgSiO3 
     ,&my_rates->F13_v0Fe3O4  
     ,&my_rates->F13_v0AC     
     ,&my_rates->F13_v0SiO2D  
     ,&my_rates->F13_v0MgO    
     ,&my_rates->F13_v0FeS    
     ,&my_rates->F13_v0Al2O3  
     , my_rates->F13_N
     ,&my_rates->F13_Size
     , my_rates->F13_D
     , my_rates->F13_T
     ,&my_rates->F13_dD
     ,&my_rates->F13_dT
     , my_rates->F13_RSiM    
     , my_rates->F13_RFeM    
     , my_rates->F13_RMg2SiO4
     , my_rates->F13_RMgSiO3 
     , my_rates->F13_RFe3O4  
     , my_rates->F13_RAC     
     , my_rates->F13_RSiO2D  
     , my_rates->F13_RMgO    
     , my_rates->F13_RFeS    
     , my_rates->F13_RAl2O3  
     ,&my_rates->F13_dRSiM    
     ,&my_rates->F13_dRFeM    
     ,&my_rates->F13_dRMg2SiO4
     ,&my_rates->F13_dRMgSiO3 
     ,&my_rates->F13_dRFe3O4  
     ,&my_rates->F13_dRAC     
     ,&my_rates->F13_dRSiO2D  
     ,&my_rates->F13_dRMgO    
     ,&my_rates->F13_dRFeS    
     ,&my_rates->F13_dRAl2O3  
     , my_rates->OF13_SiM    
     , my_rates->OF13_FeM    
     , my_rates->OF13_Mg2SiO4
     , my_rates->OF13_MgSiO3 
     , my_rates->OF13_Fe3O4  
     , my_rates->OF13_AC     
     , my_rates->OF13_SiO2D  
     , my_rates->OF13_MgO    
     , my_rates->OF13_FeS    
     , my_rates->OF13_Al2O3  
     , my_rates->LF13_SiM    
     , my_rates->LF13_FeM    
     , my_rates->LF13_Mg2SiO4
     , my_rates->LF13_MgSiO3 
     , my_rates->LF13_Fe3O4  
     , my_rates->LF13_AC     
     , my_rates->LF13_SiO2D  
     , my_rates->LF13_MgO    
     , my_rates->LF13_FeS    
     , my_rates->LF13_Al2O3  
#endif
    );
 
  return SUCCESS;
}

int _calculate_cooling_time(chemistry_data *my_chemistry,
                            chemistry_data_storage *my_rates,
                            code_units *my_units,
                            int grid_rank, int *grid_dimension,
                            int *grid_start, int *grid_end,
                            gr_float *density, gr_float *internal_energy,
                            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                            gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                            gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                            gr_float *H2I_density, gr_float *H2II_density,
                            gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
//#ifdef GRACKLE_MD
                            gr_float *DM_density, gr_float *HDII_density, gr_float *HeHII_density,
                            gr_float *CI_density,  gr_float *CII_density,  gr_float *CO_density,  gr_float *CO2_density,
                            gr_float *OI_density,  gr_float *OH_density,  gr_float *H2O_density,  gr_float *O2_density,
                            gr_float *SiI_density,  gr_float *SiOI_density,  gr_float *SiO2I_density,
                            gr_float *CH_density,  gr_float *CH2_density,  gr_float *COII_density,  gr_float *OII_density,
                            gr_float *OHII_density,  gr_float *H2OII_density,  gr_float *H3OII_density,  gr_float *O2II_density,
                            gr_float *Mg_density,  gr_float *Al_density,  gr_float *S_density,  gr_float *Fe_density,
                            gr_float *SiM_density,  gr_float *FeM_density,  gr_float *Mg2SiO4_density,  gr_float *MgSiO3_density,  gr_float *Fe3O4_density,
                            gr_float *AC_density,  gr_float *SiO2D_density,  gr_float *MgO_density,  gr_float *FeS_density,  gr_float *Al2O3_density,
//#endif
                            gr_float *e_density, gr_float *metal_density, gr_float *dust_density,
                            gr_float *cooling_time, gr_float *RT_heating_rate,
                            gr_float *volumetric_heating_rate, gr_float *specific_heating_rate)
{

  grackle_field_data my_fields;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.x_velocity               = x_velocity;
  my_fields.y_velocity               = y_velocity;
  my_fields.z_velocity               = z_velocity;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
#ifdef GRACKLE_MD
  my_fields.DM_density               = DM_density;
  my_fields.HDII_density             = HDII_density;
  my_fields.HeHII_density            = HeHII_density;
  my_fields.CI_density               = CI_density;
  my_fields.CII_density              = CII_density;
  my_fields.CO_density               = CO_density;
  my_fields.CO2_density              = CO2_density;
  my_fields.OI_density               = OI_density;
  my_fields.OH_density               = OH_density;
  my_fields.H2O_density              = H2O_density;
  my_fields.O2_density               = O2_density;
  my_fields.SiI_density              = SiI_density;
  my_fields.SiOI_density             = SiOI_density;
  my_fields.SiO2I_density            = SiO2I_density;
  my_fields.CH_density               = CH_density;
  my_fields.CH2_density              = CH2_density;
  my_fields.COII_density             = COII_density;
  my_fields.OII_density              = OII_density;
  my_fields.OHII_density             = OHII_density;
  my_fields.H2OII_density            = H2OII_density;
  my_fields.H3OII_density            = H3OII_density;
  my_fields.O2II_density             = O2II_density;
  my_fields.Mg_density               = Mg_density;
  my_fields.Al_density               = Al_density;
  my_fields.S_density                = S_density;
  my_fields.Fe_density               = Fe_density;
  my_fields.SiM_density              = SiM_density;
  my_fields.FeM_density              = FeM_density;
  my_fields.Mg2SiO4_density          = Mg2SiO4_density;
  my_fields.MgSiO3_density           = MgSiO3_density;
  my_fields.Fe3O4_density            = Fe3O4_density;
  my_fields.AC_density               = AC_density;
  my_fields.SiO2D_density            = SiO2D_density;
  my_fields.MgO_density              = MgO_density;
  my_fields.FeS_density              = FeS_density;
  my_fields.Al2O3_density            = Al2O3_density;
#endif
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;
  my_fields.dust_density             = dust_density;
  my_fields.volumetric_heating_rate  = volumetric_heating_rate;
  my_fields.specific_heating_rate    = specific_heating_rate;
  my_fields.RT_heating_rate          = RT_heating_rate;

  if (local_calculate_cooling_time(my_chemistry, my_rates, my_units,
                                   &my_fields, cooling_time) == FAIL) {
    fprintf(stderr, "Error in local_calculate_cooling_time.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_cooling_time(code_units *my_units,
                           grackle_field_data *my_fields,
                           gr_float *cooling_time)
{
  if (local_calculate_cooling_time(grackle_data, &grackle_rates, my_units,
                                   my_fields, cooling_time) == FAIL) {
    fprintf(stderr, "Error in local_calculate_cooling_time.\n");
    return FAIL;
  }
  return SUCCESS;
}
