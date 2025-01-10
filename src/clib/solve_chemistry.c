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
#include "grackle.h"
#include "grackle_macros.h"
#include "phys_constants.h"
#include "solve_rate_cool_g-cpp.h"
#include "utils.h"

/* function prototypes */

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
	int *ih2optical, int *iciecool, int *ithreebody, int *ih2cr, int *ihdcr,
        double *ciecoa,
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
        double *k125a, double *k129a, double *k130a, double *k131a, double *k132a,
        double *k133a, double *k134a, double *k135a, double *k136a, double *k137a,
        double *k148a, double *k149a, double *k150a, double *k151a, double *k152a,
        double *k153a,
        double *kz15a, double *kz16a, double *kz17a, double *kz18a, double *kz19a,
        double *kz20a, double *kz21a, double *kz22a, double *kz23a, double *kz24a,
        double *kz25a, double *kz26a, double *kz27a, double *kz28a, double *kz29a,
        double *kz30a, double *kz31a, double *kz32a, double *kz33a, double *kz34a,
        double *kz35a, double *kz36a, double *kz37a, double *kz38a, double *kz39a,
        double *kz40a, double *kz41a, double *kz42a, double *kz43a, double *kz44a,
        double *kz45a, double *kz46a, double *kz47a, double *kz48a, double *kz49a,
        double *kz50a, double *kz51a, double *kz52a, double *kz53a, double *kz54a,
        double *cieY06,
        long long *LH2_N, long long *LH2_Size,
        double *LH2_D, double *LH2_T, double *LH2_H,
        double *LH2_dD, double *LH2_dT, double *LH2_dH, double *LH2_L,
        long long *LHD_N, long long *LHD_Size,
        double *LHD_D, double *LHD_T, double *LHD_H,
        double *LHD_dD, double *LHD_dT, double *LHD_dH, double *LHD_L,
        long long *LCI_N, long long *LCI_Size,
        double *LCI_D, double *LCI_T, double *LCI_H,
        double *LCI_dD, double *LCI_dT, double *LCI_dH, double *LCI_L,
        long long *LCII_N, long long *LCII_Size,
        double *LCII_D, double *LCII_T, double *LCII_H,
        double *LCII_dD, double *LCII_dT, double *LCII_dH, double *LCII_L,
        long long *LOI_N, long long *LOI_Size,
        double *LOI_D, double *LOI_T, double *LOI_H,
        double *LOI_dD, double *LOI_dT, double *LOI_dH, double *LOI_L,
        long long *LCO_N, long long *LCO_Size,
        double *LCO_D, double *LCO_T, double *LCO_H,
        double *LCO_dD, double *LCO_dT, double *LCO_dH, double *LCO_L,
        long long *LOH_N, long long *LOH_Size,
        double *LOH_D, double *LOH_T, double *LOH_H,
        double *LOH_dD, double *LOH_dT, double *LOH_dH, double *LOH_L,
        long long *LH2O_N, long long *LH2O_Size,
        double *LH2O_D, double *LH2O_T, double *LH2O_H,
        double *LH2O_dD, double *LH2O_dT, double *LH2O_dH, double *LH2O_L,
        long long *alphap_N, long long *alphap_Size,
        double *alphap_D, double *alphap_T, double *alphap_dD, double *alphap_dT,
        double *alphap_Data,
        int *immulti, int *imabund, int *idspecies, int *itdmulti, int *idsub,
        gr_float *metal_loc, gr_float *metal_C13, gr_float *metal_C20, gr_float *metal_C25, gr_float *metal_C30,
        gr_float *metal_F13, gr_float *metal_F15, gr_float *metal_F50, gr_float *metal_F80,
        gr_float *metal_P170, gr_float *metal_P200, gr_float *metal_Y19,
        int *SN0_N,
        double *SN0_XC , double *SN0_XO, double *SN0_XMg, double *SN0_XAl,
        double *SN0_XSi, double *SN0_XS, double *SN0_XFe,
        double *SN0_fC , double *SN0_fO, double *SN0_fMg, double *SN0_fAl,
        double *SN0_fSi, double *SN0_fS, double *SN0_fFe,
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
        double *h2dustSa, double *h2dustCa, double *gasgr2a, double *gamma_isrf2a, double *grogra,
        int *idissHDI, gr_float *kdissHDI, int *iionZ, gr_float *kphCI, gr_float *kphOI,
        int *idissZ, gr_float *kdissCO, gr_float *kdissOH, gr_float *kdissH2O, int *iuseH2shield,
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
  if (self_shielding_err_check(my_chemistry, my_fields,
                               "local_solve_chemistry") == FAIL) {
    return FAIL;
  }

  /* Calculate temperature units. */

  double temperature_units = get_temperature_units(my_units);

  /* Call the fortran routine to solve cooling equations. */

  int ierr = solve_rate_cool_g(
    &metal_field_present, &dt_value, &temperature_units, &co_length_units,
    &co_density_units, my_chemistry, my_rates, my_units, my_fields,
    &my_uvb_rates
  );

  if (ierr == GR_FAIL) {
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
