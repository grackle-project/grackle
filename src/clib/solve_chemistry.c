/***********************************************************************
/
/ Solve the chemistry and cooling
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

/* function prototypes */

double get_temperature_units(code_units *my_units);

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units);

extern void FORTRAN_NAME(solve_rate_cool_g)(
        int *icool,
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idustall,
        int *idustfield, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
        int *ih2co, int *ipiht, int *idustrec, int *igammah,
	double *dx, double *dt, double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *dtoh, double *z_solar, double *fgr,
	double *k1a, double *k2a, double *k3a, double *k4a, double *k5a,
	double *k6a, double *k7a, double *k8a, double *k9a, double *k10a,
	double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a,
	double *k15a, double *k16a, double *k17a, double *k18a, double *k19a,
        double *k22a,	double *k24, double *k25, double *k26, double *k27,
        double *k28, double *k29, double *k30, double *k31,
	double *k50a, double *k51a, double *k52a, double *k53a, double *k54a,
	double *k55a, double *k56a, double *k57a, double *k58a,
	int *ndratec, double *dtemstart, double *dtemend, double *h2dusta,
	double *ncrna, double *ncrd1a, double *ncrd2a,
	double *ceHIa, double *ceHeIa, double *ceHeIIa, double *ciHIa,
	double *ciHeIa, double *ciHeISa, double *ciHeIIa,
        double *reHIIa, double *reHeII1a, double *reHeII2a, double *reHeIIIa,
        double *brema, double *compa, double *gammaha, double *isrf,
        double *regra, double *gamma_isrfa, double *comp_xraya, double *comp_temp,
	double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
        gr_float *DI, gr_float *DII, gr_float *HDI,
        gr_float *metal, gr_float *dust,
	double *hyd01ka, double *h2k01a, double *vibha,
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra, int *iH2shield,
        int *iradshield, double *avgsighi, double *avgsighei, double *avgsigheii,
        int *iradtrans, int *iradcoupled, int *iradstep, int *irt_honly,
        gr_float *kphHI, gr_float *kphHeI, gr_float *kphHeII, gr_float *kdissH2I,
        gr_float *photogamma, gr_float *xH2shield,
	int *ierr,
	int *ih2optical, int *iciecool, int *ithreebody, double *ciecoa,
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
        int *imchem, int *igrgr
      , gr_float *DM, gr_float *HDII, gr_float *HeHII
      , gr_float *CI, gr_float *CII, gr_float *CO, gr_float *CO2
      , gr_float *OI, gr_float *OH, gr_float *H2O, gr_float *O2
      , gr_float *SiI, gr_float *SiOI, gr_float *SiO2I
      , gr_float *CH, gr_float *CH2, gr_float *COII, gr_float *OII
      , gr_float *OHII, gr_float *H2OII, gr_float *H3OII, gr_float *O2II
      , gr_float *Mg, gr_float *Al, gr_float *S, gr_float *Fe
      , gr_float *SiM, gr_float *FeM, gr_float *Mg2SiO4, gr_float *MgSiO3, gr_float *Fe3O4
      , gr_float *AC, gr_float *SiO2D, gr_float *MgO, gr_float *FeS, gr_float *Al2O3
      , gr_float *reforg, gr_float *volorg, gr_float *H2Oice
      , double *k125a, double *k129a, double *k130a, double *k131a, double *k132a
      , double *k133a, double *k134a, double *k135a, double *k136a, double *k137a
      , double *k148a, double *k149a, double *k150a, double *k151a, double *k152a
      , double *k153a
      , double *kz15a, double *kz16a, double *kz17a, double *kz18a, double *kz19a
      , double *kz20a, double *kz21a, double *kz22a, double *kz23a, double *kz24a
      , double *kz25a, double *kz26a, double *kz27a, double *kz28a, double *kz29a
      , double *kz30a, double *kz31a, double *kz32a, double *kz33a, double *kz34a
      , double *kz35a, double *kz36a, double *kz37a, double *kz38a, double *kz39a
      , double *kz40a, double *kz41a, double *kz42a, double *kz43a, double *kz44a
      , double *kz45a, double *kz46a, double *kz47a, double *kz48a, double *kz49a
      , double *kz50a, double *kz51a, double *kz52a, double *kz53a, double *kz54a
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
      , double *Hgrain, double *Tgrain, double *Ograin, double *Lgrain
      , int *immulti, int *imabund, int *idspecies, int *itdmulti, int *idsub
      , gr_float *metal_loc, gr_float *metal_C13, gr_float *metal_C20, gr_float *metal_C25, gr_float *metal_C30
      , gr_float *metal_F13, gr_float *metal_F15, gr_float *metal_F50, gr_float *metal_F80
      , gr_float *metal_P170, gr_float *metal_P200, gr_float *metal_Y19
      , int *SN0_N
      , double *SN0_XC , double *SN0_XO, double *SN0_XMg, double *SN0_XAl
      , double *SN0_XSi, double *SN0_XS, double *SN0_XFe
      , double *SN0_fC , double *SN0_fO, double *SN0_fMg, double *SN0_fAl
      , double *SN0_fSi, double *SN0_fS, double *SN0_fFe
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
      , double *h2dustSa, double *h2dustCa, double *gasgr2a, double *gamma_isrf2a, double *grogra
      , int *idissHDI, gr_float *kdissHDI, int *iionZ, gr_float *kphCI, gr_float *kphOI
      , int *idissZ, gr_float *kdissCO, gr_float *kdissOH, gr_float *kdissH2O, int *iuseH2shield,
        int *iisrffield, gr_float* isrf_habing, 
        int *iH2shieldcustom, gr_float* f_shield_custom,
        int *itmax, int *exititmax);

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value)
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
  if (my_chemistry->H2_self_shielding == 1 && my_fields->grid_rank != 3){
    fprintf(stderr, "Error in solve_chemistry: H2 self-shielding option 1 "
                    "will only work for 3D Cartesian grids. Use option 2 "
                    "to provide an array of shielding lengths with "
                    "H2_self_shielding_length or option 3 to use the "
                    "local Jeans length.");
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  FORTRAN_NAME(solve_rate_cool_g)(
    &my_chemistry->with_radiative_cooling,
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
    &my_chemistry->dust_recombination_cooling,
    &my_chemistry->photoelectric_heating,
    &(my_fields->grid_dx),
    &dt_value,
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
    &my_chemistry->DeuteriumToHydrogenRatio,
    &my_chemistry->SolarMetalFractionByMass,
    &my_chemistry->local_dust_to_gas_ratio,
    my_rates->k1,
    my_rates->k2,
    my_rates->k3,
    my_rates->k4,
    my_rates->k5,
    my_rates->k6,
    my_rates->k7,
    my_rates->k8,
    my_rates->k9,
    my_rates->k10,
    my_rates->k11,
    my_rates->k12,
    my_rates->k13,
    my_rates->k13dd,
    my_rates->k14,
    my_rates->k15,
    my_rates->k16,
    my_rates->k17,
    my_rates->k18,
    my_rates->k19,
    my_rates->k22,
    &my_uvb_rates.k24,
    &my_uvb_rates.k25,
    &my_uvb_rates.k26,
    &my_uvb_rates.k27,
    &my_uvb_rates.k28,
    &my_uvb_rates.k29,
    &my_uvb_rates.k30,
    &my_uvb_rates.k31,
    my_rates->k50,
    my_rates->k51,
    my_rates->k52,
    my_rates->k53,
    my_rates->k54,
    my_rates->k55,
    my_rates->k56,
    my_rates->k57,
    my_rates->k58,
    &my_chemistry->NumberOfDustTemperatureBins,
    &my_chemistry->DustTemperatureStart,
    &my_chemistry->DustTemperatureEnd,
    my_rates->h2dust,
    my_rates->n_cr_n,
    my_rates->n_cr_d1,
    my_rates->n_cr_d2,
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
    &my_chemistry->H2_self_shielding,
    &my_chemistry->self_shielding_method,
    &my_uvb_rates.crsHI,
    &my_uvb_rates.crsHeI,
    &my_uvb_rates.crsHeII,
    &my_chemistry->use_radiative_transfer,
    &my_chemistry->radiative_transfer_coupled_rate_solver,
    &my_chemistry->radiative_transfer_intermediate_step,
    &my_chemistry->radiative_transfer_hydrogen_only,
    my_fields->RT_HI_ionization_rate,
    my_fields->RT_HeI_ionization_rate,
    my_fields->RT_HeII_ionization_rate,
    my_fields->RT_H2_dissociation_rate,
    my_fields->RT_heating_rate,
    my_fields-> H2_self_shielding_length,
    &ierr,
    &my_chemistry->h2_optical_depth_approximation,
    &my_chemistry->cie_cooling,
    &my_chemistry->three_body_rate,
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
   &my_chemistry->metal_chemistry
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
  , my_fields->reforg_density
  , my_fields->volorg_density
  , my_fields->H2Oice_density
  , my_rates->k125
  , my_rates->k129
  , my_rates->k130
  , my_rates->k131
  , my_rates->k132
  , my_rates->k133
  , my_rates->k134
  , my_rates->k135
  , my_rates->k136
  , my_rates->k137
  , my_rates->k148
  , my_rates->k149
  , my_rates->k150
  , my_rates->k151
  , my_rates->k152
  , my_rates->k153
  , my_rates->kz15
  , my_rates->kz16
  , my_rates->kz17
  , my_rates->kz18
  , my_rates->kz19
  , my_rates->kz20
  , my_rates->kz21
  , my_rates->kz22
  , my_rates->kz23
  , my_rates->kz24
  , my_rates->kz25
  , my_rates->kz26
  , my_rates->kz27
  , my_rates->kz28
  , my_rates->kz29
  , my_rates->kz30
  , my_rates->kz31
  , my_rates->kz32
  , my_rates->kz33
  , my_rates->kz34
  , my_rates->kz35
  , my_rates->kz36
  , my_rates->kz37
  , my_rates->kz38
  , my_rates->kz39
  , my_rates->kz40
  , my_rates->kz41
  , my_rates->kz42
  , my_rates->kz43
  , my_rates->kz44
  , my_rates->kz45
  , my_rates->kz46
  , my_rates->kz47
  , my_rates->kz48
  , my_rates->kz49
  , my_rates->kz50
  , my_rates->kz51
  , my_rates->kz52
  , my_rates->kz53
  , my_rates->kz54
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
  , my_rates->Hgrain
  , my_rates->Tgrain
  , my_rates->Ograin
  , my_rates->Lgrain
  ,&my_chemistry->multi_metals
  ,&my_chemistry->metal_abundances
  ,&my_chemistry->dust_species
  ,&my_chemistry->dust_temperature_multi
  ,&my_chemistry->dust_sublimation
  , my_fields->metal_loc
  , my_fields->metal_C13
  , my_fields->metal_C20
  , my_fields->metal_C25
  , my_fields->metal_C30
  , my_fields->metal_F13
  , my_fields->metal_F15
  , my_fields->metal_F50
  , my_fields->metal_F80
  , my_fields->metal_P170
  , my_fields->metal_P200
  , my_fields->metal_Y19
  ,&my_rates->SN0_N
  , my_rates->SN0_XC 
  , my_rates->SN0_XO
  , my_rates->SN0_XMg
  , my_rates->SN0_XAl
  , my_rates->SN0_XSi
  , my_rates->SN0_XS
  , my_rates->SN0_XFe
  , my_rates->SN0_fC 
  , my_rates->SN0_fO
  , my_rates->SN0_fMg
  , my_rates->SN0_fAl
  , my_rates->SN0_fSi
  , my_rates->SN0_fS
  , my_rates->SN0_fFe
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
  , my_rates->h2dustS
  , my_rates->h2dustC
  , my_rates->gas_grain2
  ,&my_rates->gamma_isrf2
  , my_rates->grain_growth_rate
  ,&my_chemistry->radiative_transfer_HDI_diss
  , my_fields->RT_HDI_dissociation_rate
  ,&my_chemistry->radiative_transfer_metal_ion
  , my_fields->RT_CI_ionization_rate
  , my_fields->RT_OI_ionization_rate
  ,&my_chemistry->radiative_transfer_metal_diss
  , my_fields->RT_CO_dissociation_rate
  , my_fields->RT_OH_dissociation_rate
  , my_fields->RT_H2O_dissociation_rate
  ,&my_chemistry->radiative_transfer_use_H2_shielding,
    &my_chemistry->use_isrf_field,
    my_fields->isrf_habing,
    &my_chemistry->H2_custom_shielding,
    my_fields->H2_custom_shielding_factor,
    &my_chemistry->max_iterations,
    &my_chemistry->exit_after_iterations_exceeded);

  if (ierr == FAIL) {
    fprintf(stderr, "Error in solve_rate_cool_g.\n");
  }

  return ierr;

}

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value)
{
  if (local_solve_chemistry(grackle_data, &grackle_rates,
                            my_units, my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}
