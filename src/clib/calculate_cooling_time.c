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
#include "utils.h"

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

double get_temperature_units(code_units *my_units);

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
        int *idustfield, int *idustrec, int *idim,
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
	int *ih2optical, int *iciecool, int *ih2cr, int *ihdcr, double *ciecoa,
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
        int *iTfloor, gr_float *Tfloor_scalar, gr_float *Tfloor,
        int *imchem, int *igrgr, int *ipcont, double *tmcool,
        gr_float *DM, gr_float *HDII, gr_float *HeHII,
        gr_float *CI, gr_float *CII, gr_float *CO, gr_float *CO2,
        gr_float *OI, gr_float *OH, gr_float *H2O, gr_float *O2,
        gr_float *SiI, gr_float *SiOI, gr_float *SiO2I,
        gr_float *CH, gr_float *CH2, gr_float *COII, gr_float *OII,
        gr_float *OHII, gr_float *H2OII, gr_float *H3OII, gr_float *O2II,
        gr_float *Mg, gr_float *Al, gr_float *S, gr_float *Fe,
        gr_float *SiM, gr_float *FeM, gr_float *Mg2SiO4, gr_float *MgSiO3, gr_float *Fe3O4,
        gr_float *AC, gr_float *SiO2D, gr_float *MgO, gr_float *FeS, gr_float *Al2O3,
        gr_float *reforg, gr_float *volorg, gr_float *H2Oice,
        double *cieY06,
        int *LH2_N, int *LH2_Size,
        double *LH2_D, double *LH2_T, double *LH2_H,
        double *LH2_dD, double *LH2_dT, double *LH2_dH, double *LH2_L,
        int *LHD_N, int *LHD_Size,
        double *LHD_D, double *LHD_T, double *LHD_H,
        double *LHD_dD, double *LHD_dT, double *LHD_dH, double *LHD_L,
        int *LCI_N, int *LCI_Size,
        double *LCI_D, double *LCI_T, double *LCI_H,
        double *LCI_dD, double *LCI_dT, double *LCI_dH, double *LCI_L,
        int *LCII_N, int *LCII_Size,
        double *LCII_D, double *LCII_T, double *LCII_H,
        double *LCII_dD, double *LCII_dT, double *LCII_dH, double *LCII_L,
        int *LOI_N, int *LOI_Size,
        double *LOI_D, double *LOI_T, double *LOI_H,
        double *LOI_dD, double *LOI_dT, double *LOI_dH, double *LOI_L,
        int *LCO_N, int *LCO_Size,
        double *LCO_D, double *LCO_T, double *LCO_H,
        double *LCO_dD, double *LCO_dT, double *LCO_dH, double *LCO_L,
        int *LOH_N, int *LOH_Size,
        double *LOH_D, double *LOH_T, double *LOH_H,
        double *LOH_dD, double *LOH_dT, double *LOH_dH, double *LOH_L,
        int *LH2O_N, int *LH2O_Size,
        double *LH2O_D, double *LH2O_T, double *LH2O_H,
        double *LH2O_dD, double *LH2O_dT, double *LH2O_dH, double *LH2O_L,
        int *alphap_N, int *alphap_Size,
        double *alphap_D, double *alphap_T, double *alphap_dD, double *alphap_dT,
        double *alphap_Data,
        int *immulti, int *imabund, int *idspecies, int *itdmulti, int *idsub,
        gr_float *metal_loc, gr_float *metal_C13, gr_float *metal_C20, gr_float *metal_C25, gr_float *metal_C30,
        gr_float *metal_F13, gr_float *metal_F15, gr_float *metal_F50, gr_float *metal_F80,
        gr_float *metal_P170, gr_float *metal_P200, gr_float *metal_Y19,
        int *SN0_N,
        double *SN0_fSiM, double *SN0_fFeM, double *SN0_fMg2SiO4, double *SN0_fMgSiO3, double *SN0_fFe3O4,
        double *SN0_fAC, double *SN0_fSiO2D, double *SN0_fMgO, double *SN0_fFeS, double *SN0_fAl2O3,
        double *SN0_freforg , double *SN0_fvolorg , double *SN0_fH2Oice,
        double *SN0_r0SiM, double *SN0_r0FeM, double *SN0_r0Mg2SiO4, double *SN0_r0MgSiO3, double *SN0_r0Fe3O4,
        double *SN0_r0AC, double *SN0_r0SiO2D, double *SN0_r0MgO, double *SN0_r0FeS, double *SN0_r0Al2O3,
        double *SN0_r0reforg , double *SN0_r0volorg , double *SN0_r0H2Oice,
        int *gr_N, int *gr_Size, double *gr_dT, double *gr_Td,
        double *SN0_kpSiM, double *SN0_kpFeM, double *SN0_kpMg2SiO4, double *SN0_kpMgSiO3, double *SN0_kpFe3O4,
        double *SN0_kpAC, double *SN0_kpSiO2D, double *SN0_kpMgO, double *SN0_kpFeS, double *SN0_kpAl2O3,
        double *SN0_kpreforg , double *SN0_kpvolorg , double *SN0_kpH2Oice,
        double *gasgr2a, double *gamma_isrf2a,
        int *iisrffield, gr_float* isrf_habing);

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

  /* Error checking for H2 shielding approximation */
  if (self_shielding_err_check(my_chemistry, my_fields,
                               "local_calculate_temperature") == FAIL) {
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

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
       my_fields->grid_dimension+0,
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
       &my_chemistry->dust_recombination_cooling,
       &(my_fields->grid_rank),
       my_fields->grid_start+0,
       my_fields->grid_start+1,
       my_fields->grid_start+2,
       my_fields->grid_end+0,
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
       &my_chemistry->cie_cooling,
       &my_chemistry->h2_cooling_rate,
       &my_chemistry->hd_cooling_rate,
       my_rates->cieco,
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
       &my_chemistry->use_temperature_floor,
       &my_chemistry->temperature_floor_scalar,
       my_fields->temperature_floor,
       &my_chemistry->metal_chemistry,
       &my_chemistry->grain_growth,
       &my_chemistry->use_primordial_continuum_opacity,
       &my_chemistry->tabulated_cooling_minimum_temperature,
       my_fields->DM_density,
       my_fields->HDII_density,
       my_fields->HeHII_density,
       my_fields->CI_density,
       my_fields->CII_density,
       my_fields->CO_density,
       my_fields->CO2_density,
       my_fields->OI_density,
       my_fields->OH_density,
       my_fields->H2O_density,
       my_fields->O2_density,
       my_fields->SiI_density,
       my_fields->SiOI_density,
       my_fields->SiO2I_density,
       my_fields->CH_density,
       my_fields->CH2_density,
       my_fields->COII_density,
       my_fields->OII_density,
       my_fields->OHII_density,
       my_fields->H2OII_density,
       my_fields->H3OII_density,
       my_fields->O2II_density,
       my_fields->Mg_density,
       my_fields->Al_density,
       my_fields->S_density,
       my_fields->Fe_density,
       my_fields->SiM_dust_density,
       my_fields->FeM_dust_density,
       my_fields->Mg2SiO4_dust_density,
       my_fields->MgSiO3_dust_density,
       my_fields->Fe3O4_dust_density,
       my_fields->AC_dust_density,
       my_fields->SiO2_dust_density,
       my_fields->MgO_dust_density,
       my_fields->FeS_dust_density,
       my_fields->Al2O3_dust_density,
       my_fields->ref_org_dust_density,
       my_fields->vol_org_dust_density,
       my_fields->H2O_ice_dust_density,
       my_rates->cieY06,
       my_rates->LH2_N,
       &my_rates->LH2_Size,
       my_rates->LH2_D,
       my_rates->LH2_T,
       my_rates->LH2_H,
       &my_rates->LH2_dD,
       &my_rates->LH2_dT,
       &my_rates->LH2_dH,
       my_rates->LH2_L,
       my_rates->LHD_N,
       &my_rates->LHD_Size,
       my_rates->LHD_D,
       my_rates->LHD_T,
       my_rates->LHD_H,
       &my_rates->LHD_dD,
       &my_rates->LHD_dT,
       &my_rates->LHD_dH,
       my_rates->LHD_L,
       my_rates->LCI_N,
       &my_rates->LCI_Size,
       my_rates->LCI_D,
       my_rates->LCI_T,
       my_rates->LCI_H,
       &my_rates->LCI_dD,
       &my_rates->LCI_dT,
       &my_rates->LCI_dH,
       my_rates->LCI_L,
       my_rates->LCII_N,
       &my_rates->LCII_Size,
       my_rates->LCII_D,
       my_rates->LCII_T,
       my_rates->LCII_H,
       &my_rates->LCII_dD,
       &my_rates->LCII_dT,
       &my_rates->LCII_dH,
       my_rates->LCII_L,
       my_rates->LOI_N,
       &my_rates->LOI_Size,
       my_rates->LOI_D,
       my_rates->LOI_T,
       my_rates->LOI_H,
       &my_rates->LOI_dD,
       &my_rates->LOI_dT,
       &my_rates->LOI_dH,
       my_rates->LOI_L,
       my_rates->LCO_N,
       &my_rates->LCO_Size,
       my_rates->LCO_D,
       my_rates->LCO_T,
       my_rates->LCO_H,
       &my_rates->LCO_dD,
       &my_rates->LCO_dT,
       &my_rates->LCO_dH,
       my_rates->LCO_L,
       my_rates->LOH_N,
       &my_rates->LOH_Size,
       my_rates->LOH_D,
       my_rates->LOH_T,
       my_rates->LOH_H,
       &my_rates->LOH_dD,
       &my_rates->LOH_dT,
       &my_rates->LOH_dH,
       my_rates->LOH_L,
       my_rates->LH2O_N,
       &my_rates->LH2O_Size,
       my_rates->LH2O_D,
       my_rates->LH2O_T,
       my_rates->LH2O_H,
       &my_rates->LH2O_dD,
       &my_rates->LH2O_dT,
       &my_rates->LH2O_dH,
       my_rates->LH2O_L,
       my_rates->alphap_N,
       &my_rates->alphap_Size,
       my_rates->alphap_D,
       my_rates->alphap_T,
       &my_rates->alphap_dD,
       &my_rates->alphap_dT,
       my_rates->alphap_Data,
       &my_chemistry->multi_metals,
       &my_chemistry->metal_abundances,
       &my_chemistry->dust_species,
       &my_chemistry->use_multiple_dust_temperatures,
       &my_chemistry->dust_sublimation,
       my_fields->local_ISM_metal_density,
       my_fields->ccsn13_metal_density,
       my_fields->ccsn20_metal_density,
       my_fields->ccsn25_metal_density,
       my_fields->ccsn30_metal_density,
       my_fields->fsn13_metal_density,
       my_fields->fsn15_metal_density,
       my_fields->fsn50_metal_density,
       my_fields->fsn80_metal_density,
       my_fields->pisn170_metal_density,
       my_fields->pisn200_metal_density,
       my_fields->y19_metal_density,
       &my_rates->SN0_N,
       my_rates->SN0_fSiM,
       my_rates->SN0_fFeM,
       my_rates->SN0_fMg2SiO4,
       my_rates->SN0_fMgSiO3,
       my_rates->SN0_fFe3O4,
       my_rates->SN0_fAC,
       my_rates->SN0_fSiO2D,
       my_rates->SN0_fMgO,
       my_rates->SN0_fFeS,
       my_rates->SN0_fAl2O3,
       my_rates->SN0_freforg,
       my_rates->SN0_fvolorg,
       my_rates->SN0_fH2Oice,
       my_rates->SN0_r0SiM,
       my_rates->SN0_r0FeM,
       my_rates->SN0_r0Mg2SiO4,
       my_rates->SN0_r0MgSiO3,
       my_rates->SN0_r0Fe3O4,
       my_rates->SN0_r0AC,
       my_rates->SN0_r0SiO2D,
       my_rates->SN0_r0MgO,
       my_rates->SN0_r0FeS,
       my_rates->SN0_r0Al2O3,
       my_rates->SN0_r0reforg,
       my_rates->SN0_r0volorg,
       my_rates->SN0_r0H2Oice,
       my_rates->gr_N,
       &my_rates->gr_Size,
       &my_rates->gr_dT,
       my_rates->gr_Td,
       my_rates->SN0_kpSiM,
       my_rates->SN0_kpFeM,
       my_rates->SN0_kpMg2SiO4,
       my_rates->SN0_kpMgSiO3,
       my_rates->SN0_kpFe3O4,
       my_rates->SN0_kpAC,
       my_rates->SN0_kpSiO2D,
       my_rates->SN0_kpMgO,
       my_rates->SN0_kpFeS,
       my_rates->SN0_kpAl2O3,
       my_rates->SN0_kpreforg,
       my_rates->SN0_kpvolorg,
       my_rates->SN0_kpH2Oice,
       my_rates->gas_grain2,
       &my_rates->gamma_isrf2,
       &my_chemistry->use_isrf_field,
       my_fields->isrf_habing

    );
 
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
