/// @file fortran_func_decls.h
/// @brief Header file defining forward declaring function signatures of
///        Fortran subroutines have not been ported to C++ yet.

#ifndef FORTRAN_FN_DECLARATIONS_HPP
#define FORTRAN_FN_DECLARATIONS_HPP

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "grackle_macros.h" // FORTRAN_NAME
#include "grackle.h"        // gr_float
#include "phys_constants.h" // physical constants
#include <stdint.h>         // int32_t
typedef int32_t gr_mask_type;
typedef long long gr_i64;

#define MASK_TRUE 1
#define MASK_FALSE 0

void FORTRAN_NAME(calc_all_tdust_gasgr_1d_g)(
  int* in, int* jn, int* kn, int* nratec, int* idustfield, int* is, int* ie,
  int* j, int* k, double* fgr, double* gamma_isrfa, double* trad,
  double* gasgra, long long* indixe, double* tdef, double* tgas, double* tdust,
  double* metallicity, double* dust2gas, double* nh, double* gasgr_tdust,
  gr_mask_type* itmask_metal, int* idspecies, int* itdmulti, int* gr_N,
  int* gr_Size, double* gr_dT, double* gr_Td, double* tSiM, double* tFeM,
  double* tMg2SiO4, double* tMgSiO3, double* tFe3O4, double* tAC,
  double* tSiO2D, double* tMgO, double* tFeS, double* tAl2O3, double* treforg,
  double* tvolorg, double* tH2Oice, double* gasgr2a, double* gamma_isrf2a,
  double* coolunit, double* gasgr, double* myisrf, double* sgSiM, double* sgFeM,
  double* sgMg2SiO4, double* sgMgSiO3, double* sgFe3O4, double* sgAC,
  double* sgSiO2D, double* sgMgO, double* sgFeS, double* sgAl2O3,
  double* sgreforg, double* sgvolorg, double* sgH2Oice, double* sgtot,
  double* alSiM_data_ptr, double* alFeM_data_ptr, double* alMg2SiO4_data_ptr,
  double* alMgSiO3_data_ptr, double* alFe3O4_data_ptr, double* alAC_data_ptr,
  double* alSiO2D_data_ptr, double* alMgO_data_ptr, double* alFeS_data_ptr,
  double* alAl2O3_data_ptr, double* alreforg_data_ptr,
  double* alvolorg_data_ptr, double* alH2Oice_data_ptr, double* altot_data_ptr,
  double* kpSiM, double* kpFeM, double* kpMg2SiO4, double* kpMgSiO3,
  double* kpFe3O4, double* kpAC, double* kpSiO2D, double* kpMgO, double* kpFeS,
  double* kpAl2O3, double* kpreforg, double* kpvolorg, double* kpH2Oice,
  double* kptot, double* gasSiM, double* gasFeM, double* gasMg2SiO4,
  double* gasMgSiO3, double* gasFe3O4, double* gasAC, double* gasSiO2D,
  double* gasMgO, double* gasFeS, double* gasAl2O3, double* gasreforg,
  double* gasvolorg, double* gasH2Oice
);

void FORTRAN_NAME(calc_grain_size_increment_1d)(
  int* immulti, int* imabund, int* idspecies, int* igrgr, gr_mask_type* itmask,
  int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k, double* dom,
  gr_float* d_data_ptr, gr_float* SiM_data_ptr, gr_float* FeM_data_ptr,
  gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, gr_float* metal_data_ptr,
  gr_float* metal_loc_data_ptr, gr_float* metal_C13_data_ptr,
  gr_float* metal_C20_data_ptr, gr_float* metal_C25_data_ptr,
  gr_float* metal_C30_data_ptr, gr_float* metal_F13_data_ptr,
  gr_float* metal_F15_data_ptr, gr_float* metal_F50_data_ptr,
  gr_float* metal_F80_data_ptr, gr_float* metal_P170_data_ptr,
  gr_float* metal_P200_data_ptr, gr_float* metal_Y19_data_ptr, int* SN0_N,
  double* SN0_fSiM, double* SN0_fFeM, double* SN0_fMg2SiO4, double* SN0_fMgSiO3,
  double* SN0_fFe3O4, double* SN0_fAC, double* SN0_fSiO2D, double* SN0_fMgO,
  double* SN0_fFeS, double* SN0_fAl2O3, double* SN0_freforg,
  double* SN0_fvolorg, double* SN0_fH2Oice, double* SN0_r0SiM_data_ptr,
  double* SN0_r0FeM_data_ptr, double* SN0_r0Mg2SiO4_data_ptr,
  double* SN0_r0MgSiO3_data_ptr, double* SN0_r0Fe3O4_data_ptr,
  double* SN0_r0AC_data_ptr, double* SN0_r0SiO2D_data_ptr,
  double* SN0_r0MgO_data_ptr, double* SN0_r0FeS_data_ptr,
  double* SN0_r0Al2O3_data_ptr, double* SN0_r0reforg_data_ptr,
  double* SN0_r0volorg_data_ptr, double* SN0_r0H2Oice_data_ptr, int* gr_N,
  int* gr_Size, double* gr_dT, double* gr_Td, double* SN0_kpSiM_data_ptr,
  double* SN0_kpFeM_data_ptr, double* SN0_kpMg2SiO4_data_ptr,
  double* SN0_kpMgSiO3_data_ptr, double* SN0_kpFe3O4_data_ptr,
  double* SN0_kpAC_data_ptr, double* SN0_kpSiO2D_data_ptr,
  double* SN0_kpMgO_data_ptr, double* SN0_kpFeS_data_ptr,
  double* SN0_kpAl2O3_data_ptr, double* SN0_kpreforg_data_ptr,
  double* SN0_kpvolorg_data_ptr, double* SN0_kpH2Oice_data_ptr, double* sgSiM,
  double* sgFeM, double* sgMg2SiO4, double* sgMgSiO3, double* sgFe3O4,
  double* sgAC, double* sgSiO2D, double* sgMgO, double* sgFeS, double* sgAl2O3,
  double* sgreforg, double* sgvolorg, double* sgH2Oice, double* sgtot,
  double* alSiM_data_ptr, double* alFeM_data_ptr, double* alMg2SiO4_data_ptr,
  double* alMgSiO3_data_ptr, double* alFe3O4_data_ptr, double* alAC_data_ptr,
  double* alSiO2D_data_ptr, double* alMgO_data_ptr, double* alFeS_data_ptr,
  double* alAl2O3_data_ptr, double* alreforg_data_ptr,
  double* alvolorg_data_ptr, double* alH2Oice_data_ptr, double* altot_data_ptr
);

void FORTRAN_NAME(calc_grain_size_increment_species_1d)(
  int* igrgr, gr_mask_type* itmask, int* SN0_N, int* in, int* jn, int* kn,
  int* is, int* ie, int* j, int* k, double* dom, gr_float* d_data_ptr, int* nSN,
  gr_float* dsp_data_ptr, gr_float* SN_metal_data_ptr, double* SN_fsp,
  double* SN_r0sp_data_ptr, double* ssp, double* sgsp, double* alsp_data_ptr,
  int* gr_N, int* gr_Size, double* gr_dT, double* gr_Td,
  double* SN_kp0sp_data_ptr
);

void FORTRAN_NAME(solve_cubic_equation)(
  double* a, double* b, double* c, double* root
);

void FORTRAN_NAME(calc_tdust_1d_g)(
  double* tdust, double* tgas, double* nh, double* gasgr, double* gamma_isrfa,
  double* isrf, gr_mask_type* itmask, double* trad, int* in, int* is, int* ie,
  int* j, int* k, int* gr_N, int* gr_Size, double* gr_dT, double* gr_Td,
  gr_float* alsp_data_ptr, double* kgr, int* idspecies
);

void FORTRAN_NAME(calc_kappa_gr_g)(
  double* tdust, double* kgr, gr_mask_type* itmask, int* in, int* is, int* ie,
  double* t_subl, int* gr_N, int* gr_Size, double* gr_dT, double* gr_Td,
  gr_float* logalsp_data_ptr, int* idspecies
);

void FORTRAN_NAME(calc_gr_balance_g)(
  double* tdust, double* tgas, double* kgr, double* trad4, double* gasgr,
  double* gamma_isrf, double* nh, gr_mask_type* itmask, double* sol, int* in,
  int* is, int* ie
);

void FORTRAN_NAME(calc_temp1d_cloudy_g)(
  gr_float* d_data_ptr, gr_float* metal_data_ptr, gr_float* e_data_ptr,
  double* rhoH, int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
  double* tgas, double* mmw, double* dom, double* zr, double* temstart,
  double* temend, double* gamma, double* utem, int* imetal,
  long long* clGridRank, long long* clGridDim, double* clPar1, double* clPar2,
  double* clPar3, long long* clDataSize, double* clMMW, gr_mask_type* itmask
);

void FORTRAN_NAME(cool1d_cloudy_g)(
  gr_float* d_data_ptr, double* rhoH, double* metallicity, int* in, int* jn,
  int* kn, int* is, int* ie, int* j, int* k, double* logtem, double* edot,
  double* comp2, double* dom, double* zr, int* icmbTfloor, int* iClHeat,
  int* iZscale, long long* clGridRank, long long* clGridDim, double* clPar1,
  double* clPar2, double* clPar3, long long* clDataSize, double* clCooling,
  double* clHeating, gr_mask_type* itmask
);

void FORTRAN_NAME(cool1d_cloudy_old_tables_g)(
  gr_float* d_data_ptr, gr_float* de_data_ptr, double* rhoH,
  double* metallicity, int* in, int* jn, int* kn, int* is, int* ie, int* j,
  int* k, double* logtem, double* edot, double* comp2, int* ispecies,
  double* dom, double* zr, int* icmbTfloor, int* iClHeat, double* clEleFra,
  long long* clGridRank, long long* clGridDim, double* clPar1, double* clPar2,
  double* clPar3, double* clPar4, double* clPar5, long long* clDataSize,
  double* clCooling, double* clHeating, gr_mask_type* itmask
);

void FORTRAN_NAME(cool1d_multi_g)(
  gr_float* d_data_ptr, gr_float* e_data_ptr, gr_float* u_data_ptr,
  gr_float* v_data_ptr, gr_float* w_data_ptr, gr_float* de_data_ptr,
  gr_float* HI_data_ptr, gr_float* HII_data_ptr, gr_float* HeI_data_ptr,
  gr_float* HeII_data_ptr, gr_float* HeIII_data_ptr, int* in, int* jn, int* kn,
  int* nratec, int* iexpand, int* ispecies, int* imetal, int* imcool,
  int* idust, int* idustall, int* idustfield, int* idustrec, int* idim, int* is,
  int* ie, int* j, int* k, int* ih2co, int* ipiht, int* iter, int* igammah,
  double* aye, double* temstart, double* temend, double* z_solar, double* fgr,
  double* utem, double* uxyz, double* uaye, double* urho, double* utim,
  double* gamma, double* fh, double* ceHIa, double* ceHeIa, double* ceHeIIa,
  double* ciHIa, double* ciHeIa, double* ciHeISa, double* ciHeIIa,
  double* reHIIa, double* reHeII1a, double* reHeII2a, double* reHeIIIa,
  double* brema, double* compa, double* gammaha, double* isrf, double* regra,
  double* gamma_isrfa, double* comp_xraya, double* comp_temp, double* piHI,
  double* piHeI, double* piHeII, double* comp1, double* comp2,
  gr_float* HM_data_ptr, gr_float* H2I_data_ptr, gr_float* H2II_data_ptr,
  gr_float* DI_data_ptr, gr_float* DII_data_ptr, gr_float* HDI_data_ptr,
  gr_float* metal_data_ptr, gr_float* dust_data_ptr, double* hyd01ka,
  double* h2k01a, double* vibha, double* rotha, double* rotla, double* hyd01k,
  double* h2k01, double* vibh, double* roth, double* rotl, double* gpldla,
  double* gphdla, double* gpldl, double* gphdl, double* hdltea, double* hdlowa,
  double* hdlte, double* hdlow, double* gaHIa, double* gaH2a, double* gaHea,
  double* gaHpa, double* gaela, double* h2ltea, double* gasgra, double* ceHI,
  double* ceHeI, double* ceHeII, double* ciHI, double* ciHeI, double* ciHeIS,
  double* ciHeII, double* reHII, double* reHeII1, double* reHeII2,
  double* reHeIII, double* brem, long long* indixe, double* t1, double* t2,
  double* logtem, double* tdef, double* edot, double* tgas, double* tgasold,
  double* mmw, double* p2d, double* tdust, double* metallicity,
  double* dust2gas, double* rhoH, double* mynh, double* myde,
  double* gammaha_eff, double* gasgr_tdust, double* regr, int* iradshield,
  double* avgsighi, double* avgsighei, double* avgsigheii, double* k24,
  double* k26, int* iradtrans, gr_float* photogamma_data_ptr, int* ih2optical,
  int* iciecool, int* ih2cr, int* ihdcr, double* ciecoa, double* cieco,
  int* icmbTfloor, int* iClHeat, double* clEleFra, long long* priGridRank,
  long long* priGridDim, double* priPar1, double* priPar2, double* priPar3,
  double* priPar4, double* priPar5, long long* priDataSize, double* priCooling,
  double* priHeating, double* priMMW, long long* metGridRank,
  long long* metGridDim, double* metPar1, double* metPar2, double* metPar3,
  double* metPar4, double* metPar5, long long* metDataSize, double* metCooling,
  double* metHeating, int* clnew, int* iVheat, int* iMheat,
  gr_float* Vheat_data_ptr, gr_float* Mheat_data_ptr, int* iTfloor,
  double* Tfloor_scalar, gr_float* Tfloor_data_ptr, int* iisrffield,
  gr_float* isrf_habing_data_ptr, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, int* imchem, int* igrgr, int* ipcont,
  double* tmcool, gr_float* DM_data_ptr, gr_float* HDII_data_ptr,
  gr_float* HeHII_data_ptr, gr_float* CI_data_ptr, gr_float* CII_data_ptr,
  gr_float* CO_data_ptr, gr_float* CO2_data_ptr, gr_float* OI_data_ptr,
  gr_float* OH_data_ptr, gr_float* H2O_data_ptr, gr_float* O2_data_ptr,
  gr_float* SiI_data_ptr, gr_float* SiOI_data_ptr, gr_float* SiO2I_data_ptr,
  gr_float* CH_data_ptr, gr_float* CH2_data_ptr, gr_float* COII_data_ptr,
  gr_float* OII_data_ptr, gr_float* OHII_data_ptr, gr_float* H2OII_data_ptr,
  gr_float* H3OII_data_ptr, gr_float* O2II_data_ptr, gr_float* Mg_data_ptr,
  gr_float* Al_data_ptr, gr_float* S_data_ptr, gr_float* Fe_data_ptr,
  gr_float* SiM_data_ptr, gr_float* FeM_data_ptr, gr_float* Mg2SiO4_data_ptr,
  gr_float* MgSiO3_data_ptr, gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr,
  gr_float* SiO2D_data_ptr, gr_float* MgO_data_ptr, gr_float* FeS_data_ptr,
  gr_float* Al2O3_data_ptr, gr_float* reforg_data_ptr,
  gr_float* volorg_data_ptr, gr_float* H2Oice_data_ptr, double* cieY06a,
  long long* LH2_N, long long* LH2_Size, double* LH2_D, double* LH2_T,
  double* LH2_H, double* LH2_dD, double* LH2_dT, double* LH2_dH, double* LH2_L,
  long long* LHD_N, long long* LHD_Size, double* LHD_D, double* LHD_T,
  double* LHD_H, double* LHD_dD, double* LHD_dT, double* LHD_dH, double* LHD_L,
  long long* LCI_N, long long* LCI_Size, double* LCI_D, double* LCI_T,
  double* LCI_H, double* LCI_dD, double* LCI_dT, double* LCI_dH, double* LCI_L,
  long long* LCII_N, long long* LCII_Size, double* LCII_D, double* LCII_T,
  double* LCII_H, double* LCII_dD, double* LCII_dT, double* LCII_dH,
  double* LCII_L, long long* LOI_N, long long* LOI_Size, double* LOI_D,
  double* LOI_T, double* LOI_H, double* LOI_dD, double* LOI_dT, double* LOI_dH,
  double* LOI_L, long long* LCO_N, long long* LCO_Size, double* LCO_D,
  double* LCO_T, double* LCO_H, double* LCO_dD, double* LCO_dT, double* LCO_dH,
  double* LCO_L, long long* LOH_N, long long* LOH_Size, double* LOH_D,
  double* LOH_T, double* LOH_H, double* LOH_dD, double* LOH_dT, double* LOH_dH,
  double* LOH_L, long long* LH2O_N, long long* LH2O_Size, double* LH2O_D,
  double* LH2O_T, double* LH2O_H, double* LH2O_dD, double* LH2O_dT,
  double* LH2O_dH, double* LH2O_L, long long* alphap_N, long long* alphap_Size,
  double* alphap_D, double* alphap_T, double* alphap_dD, double* alphap_dT,
  double* alphap_Data, int* immulti, int* imabund, int* idspecies,
  int* itdmulti, int* idsub, gr_float* metal_loc_data_ptr,
  gr_float* metal_C13_data_ptr, gr_float* metal_C20_data_ptr,
  gr_float* metal_C25_data_ptr, gr_float* metal_C30_data_ptr,
  gr_float* metal_F13_data_ptr, gr_float* metal_F15_data_ptr,
  gr_float* metal_F50_data_ptr, gr_float* metal_F80_data_ptr,
  gr_float* metal_P170_data_ptr, gr_float* metal_P200_data_ptr,
  gr_float* metal_Y19_data_ptr, int* SN0_N, double* SN0_fSiM, double* SN0_fFeM,
  double* SN0_fMg2SiO4, double* SN0_fMgSiO3, double* SN0_fFe3O4,
  double* SN0_fAC, double* SN0_fSiO2D, double* SN0_fMgO, double* SN0_fFeS,
  double* SN0_fAl2O3, double* SN0_freforg, double* SN0_fvolorg,
  double* SN0_fH2Oice, double* SN0_r0SiM_data_ptr, double* SN0_r0FeM_data_ptr,
  double* SN0_r0Mg2SiO4_data_ptr, double* SN0_r0MgSiO3_data_ptr,
  double* SN0_r0Fe3O4_data_ptr, double* SN0_r0AC_data_ptr,
  double* SN0_r0SiO2D_data_ptr, double* SN0_r0MgO_data_ptr,
  double* SN0_r0FeS_data_ptr, double* SN0_r0Al2O3_data_ptr,
  double* SN0_r0reforg_data_ptr, double* SN0_r0volorg_data_ptr,
  double* SN0_r0H2Oice_data_ptr, int* gr_N, int* gr_Size, double* gr_dT,
  double* gr_Td, double* SN0_kpSiM_data_ptr, double* SN0_kpFeM_data_ptr,
  double* SN0_kpMg2SiO4_data_ptr, double* SN0_kpMgSiO3_data_ptr,
  double* SN0_kpFe3O4_data_ptr, double* SN0_kpAC_data_ptr,
  double* SN0_kpSiO2D_data_ptr, double* SN0_kpMgO_data_ptr,
  double* SN0_kpFeS_data_ptr, double* SN0_kpAl2O3_data_ptr,
  double* SN0_kpreforg_data_ptr, double* SN0_kpvolorg_data_ptr,
  double* SN0_kpH2Oice_data_ptr, double* tSiM, double* tFeM, double* tMg2SiO4,
  double* tMgSiO3, double* tFe3O4, double* tAC, double* tSiO2D, double* tMgO,
  double* tFeS, double* tAl2O3, double* treforg, double* tvolorg,
  double* tH2Oice, double* gasgr2a, double* gamma_isrf2a
);

void FORTRAN_NAME(gaussj_g)(
  int* n, double* a_data_ptr, double* b, int* ierr
);

// in the following interpolate functions, all of the arguments are const
// pointers other than value. (but I have just annotated the subset of
// arguments with const that are necessary to get the desired wrapper function
// signature)

void FORTRAN_NAME(interpolate_1d_g)(
  double* input1, const long long* gridDim, const double* gridPar1,
  double* dgridPar1, long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_2d_g)(
  double* input1, double* input2, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  double* dgridPar2, long long* dataSize, const double* dataField,
  double* value
);

void FORTRAN_NAME(interpolate_3d_g)(
  double* input1, double* input2, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  double* dgridPar2, const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_3dz_g)(
  double* input1, double* input2, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, const double* gridPar2,
  long long* index2, const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, long long* end_int,
  double* value
);

void FORTRAN_NAME(interpolate_2df3d_g)(
  double* input1, double* input3, const long long* gridDim,
  const double* gridPar1, double* dgridPar1, long long* index2,
  const double* gridPar3, double* dgridPar3,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_4d_g)(
  double* input1, double* input2, double* input3, double* input4,
  const long long* gridDim, const double* gridPar1, double* dgridPar1,
  const double* gridPar2, double* dgridPar2, const double* gridPar3,
  double* dgridPar3, const double* gridPar4, double* dgridPar4,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(interpolate_5d_g)(
  double* input1, double* input2, double* input3, double* input4,
  double* input5, const long long* gridDim, const double* gridPar1,
  double* dgridPar1, const double* gridPar2, double* dgridPar2,
  const double* gridPar3, double* dgridPar3, const double* gridPar4,
  double* dgridPar4, const double* gridPar5, double* dgridPar5,
  long long* dataSize, const double* dataField, double* value
);

void FORTRAN_NAME(scale_fields_g)(
  int* ispecies, int* imetal, int* idustfield, int* imchem, int* imabund,
  int* idspecies, int* immulti, int* igrgr, int* idsub, gr_float* factor,
  int* is, int* ie, int* js, int* je, int* ks, int* ke, int* in, int* jn,
  int* kn, gr_float* d_data_ptr, gr_float* de_data_ptr, gr_float* HI_data_ptr,
  gr_float* HII_data_ptr, gr_float* HeI_data_ptr, gr_float* HeII_data_ptr,
  gr_float* HeIII_data_ptr, gr_float* HM_data_ptr, gr_float* H2I_data_ptr,
  gr_float* H2II_data_ptr, gr_float* DI_data_ptr, gr_float* DII_data_ptr,
  gr_float* HDI_data_ptr, gr_float* DM_data_ptr, gr_float* HDII_data_ptr,
  gr_float* HeHII_data_ptr, gr_float* metal_data_ptr, gr_float* dust_data_ptr,
  gr_float* CI_data_ptr, gr_float* CII_data_ptr, gr_float* CO_data_ptr,
  gr_float* CO2_data_ptr, gr_float* OI_data_ptr, gr_float* OH_data_ptr,
  gr_float* H2O_data_ptr, gr_float* O2_data_ptr, gr_float* SiI_data_ptr,
  gr_float* SiOI_data_ptr, gr_float* SiO2I_data_ptr, gr_float* CH_data_ptr,
  gr_float* CH2_data_ptr, gr_float* COII_data_ptr, gr_float* OII_data_ptr,
  gr_float* OHII_data_ptr, gr_float* H2OII_data_ptr, gr_float* H3OII_data_ptr,
  gr_float* O2II_data_ptr, gr_float* Mg_data_ptr, gr_float* Al_data_ptr,
  gr_float* S_data_ptr, gr_float* Fe_data_ptr, gr_float* SiM_data_ptr,
  gr_float* FeM_data_ptr, gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, gr_float* metal_loc_data_ptr,
  gr_float* metal_C13_data_ptr, gr_float* metal_C20_data_ptr,
  gr_float* metal_C25_data_ptr, gr_float* metal_C30_data_ptr,
  gr_float* metal_F13_data_ptr, gr_float* metal_F15_data_ptr,
  gr_float* metal_F50_data_ptr, gr_float* metal_F80_data_ptr,
  gr_float* metal_P170_data_ptr, gr_float* metal_P200_data_ptr,
  gr_float* metal_Y19_data_ptr
);

void FORTRAN_NAME(ceiling_species_g)(
  gr_float* d_data_ptr, gr_float* de_data_ptr, gr_float* HI_data_ptr,
  gr_float* HII_data_ptr, gr_float* HeI_data_ptr, gr_float* HeII_data_ptr,
  gr_float* HeIII_data_ptr, gr_float* HM_data_ptr, gr_float* H2I_data_ptr,
  gr_float* H2II_data_ptr, gr_float* DI_data_ptr, gr_float* DII_data_ptr,
  gr_float* HDI_data_ptr, gr_float* metal_data_ptr, gr_float* dust_data_ptr,
  int* is, int* ie, int* js, int* je, int* ks, int* ke, int* in, int* jn,
  int* kn, int* ispecies, int* imetal, int* idustfield, gr_float* DM_data_ptr,
  gr_float* HDII_data_ptr, gr_float* HeHII_data_ptr, int* imabund, int* imchem,
  int* idspecies, int* immulti, int* igrgr, int* idsub, gr_float* CI_data_ptr,
  gr_float* CII_data_ptr, gr_float* CO_data_ptr, gr_float* CO2_data_ptr,
  gr_float* OI_data_ptr, gr_float* OH_data_ptr, gr_float* H2O_data_ptr,
  gr_float* O2_data_ptr, gr_float* SiI_data_ptr, gr_float* SiOI_data_ptr,
  gr_float* SiO2I_data_ptr, gr_float* CH_data_ptr, gr_float* CH2_data_ptr,
  gr_float* COII_data_ptr, gr_float* OII_data_ptr, gr_float* OHII_data_ptr,
  gr_float* H2OII_data_ptr, gr_float* H3OII_data_ptr, gr_float* O2II_data_ptr,
  gr_float* Mg_data_ptr, gr_float* Al_data_ptr, gr_float* S_data_ptr,
  gr_float* Fe_data_ptr, gr_float* SiM_data_ptr, gr_float* FeM_data_ptr,
  gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, gr_float* metal_loc_data_ptr,
  gr_float* metal_C13_data_ptr, gr_float* metal_C20_data_ptr,
  gr_float* metal_C25_data_ptr, gr_float* metal_C30_data_ptr,
  gr_float* metal_F13_data_ptr, gr_float* metal_F15_data_ptr,
  gr_float* metal_F50_data_ptr, gr_float* metal_F80_data_ptr,
  gr_float* metal_P170_data_ptr, gr_float* metal_P200_data_ptr,
  gr_float* metal_Y19_data_ptr
);

void FORTRAN_NAME(lookup_cool_rates1d_g)(
  double* temstart, double* temend, int* nratec, int* j, int* k, int* is,
  int* ie, int* ithreebody, int* in, int* jn, int* kn, int* ispecies,
  gr_mask_type* anydust, int* iH2shield, int* iradshield, double* tgas1d,
  double* mmw, gr_float* d_data_ptr, gr_float* HI_data_ptr,
  gr_float* HII_data_ptr, gr_float* HeI_data_ptr, gr_float* HeII_data_ptr,
  gr_float* HeIII_data_ptr, gr_float* HM_data_ptr, gr_float* H2I_data_ptr,
  gr_float* H2II_data_ptr, gr_float* DI_data_ptr, gr_float* DII_data_ptr,
  gr_float* HDI_data_ptr, double* tdust, double* dust2gas, double* k1a,
  double* k2a, double* k3a, double* k4a, double* k5a, double* k6a, double* k7a,
  double* k8a, double* k9a, double* k10a, double* k11a, double* k12a,
  double* k13a, double* k13dda_data_ptr, double* k14a, double* k15a,
  double* k16a, double* k17a, double* k18a, double* k19a, double* k22a,
  double* k50a, double* k51a, double* k52a, double* k53a, double* k54a,
  double* k55a, double* k56a, double* k57a, double* k58a, int* ndratec,
  double* dtemstart, double* dtemend, double* h2dusta_data_ptr, double* ncrna,
  double* ncrd1a, double* ncrd2a, double* avgsighi, double* avgsighei,
  double* avgsigheii, double* piHI, double* piHeI, double* k1, double* k2,
  double* k3, double* k4, double* k5, double* k6, double* k7, double* k8,
  double* k9, double* k10, double* k11, double* k12, double* k13, double* k14,
  double* k15, double* k16, double* k17, double* k18, double* k19, double* k22,
  double* k24, double* k25, double* k26, double* k28, double* k29, double* k30,
  double* k31, double* k50, double* k51, double* k52, double* k53, double* k54,
  double* k55, double* k56, double* k57, double* k58, double* k13dd_data_ptr,
  double* k24shield, double* k25shield, double* k26shield, double* k28shield,
  double* k29shield, double* k30shield, double* k31shield, double* h2dust,
  double* ncrn, double* ncrd1, double* ncrd2, double* t1, double* t2,
  double* tdef, double* logtem, long long* indixe, double* dom,
  double* coolunit, double* tbase1, double* uxyz, double* xbase1, double* dx_cgs,
  double* c_ljeans, int* iradtrans, gr_float* kdissH2I_data_ptr,
  gr_float* xH2shield_data_ptr, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, double* fh, gr_float* metal_data_ptr,
  gr_float* DM_data_ptr, gr_float* HDII_data_ptr, gr_float* HeHII_data_ptr,
  int* imetal, int* imchem, int* igrgr, gr_float* CI_data_ptr,
  gr_float* CII_data_ptr, gr_float* CO_data_ptr, gr_float* CO2_data_ptr,
  gr_float* OI_data_ptr, gr_float* OH_data_ptr, gr_float* H2O_data_ptr,
  gr_float* O2_data_ptr, gr_float* SiI_data_ptr, gr_float* SiOI_data_ptr,
  gr_float* SiO2I_data_ptr, gr_float* CH_data_ptr, gr_float* CH2_data_ptr,
  gr_float* COII_data_ptr, gr_float* OII_data_ptr, gr_float* OHII_data_ptr,
  gr_float* H2OII_data_ptr, gr_float* H3OII_data_ptr, gr_float* O2II_data_ptr,
  gr_float* Mg_data_ptr, gr_float* Al_data_ptr, gr_float* S_data_ptr,
  gr_float* Fe_data_ptr, gr_float* SiM_data_ptr, gr_float* FeM_data_ptr,
  gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, double* k125a, double* k129a, double* k130a,
  double* k131a, double* k132a, double* k133a, double* k134a, double* k135a,
  double* k136a, double* k137a, double* k148a, double* k149a, double* k150a,
  double* k151a, double* k152a, double* k153a, double* kz15a, double* kz16a,
  double* kz17a, double* kz18a, double* kz19a, double* kz20a, double* kz21a,
  double* kz22a, double* kz23a, double* kz24a, double* kz25a, double* kz26a,
  double* kz27a, double* kz28a, double* kz29a, double* kz30a, double* kz31a,
  double* kz32a, double* kz33a, double* kz34a, double* kz35a, double* kz36a,
  double* kz37a, double* kz38a, double* kz39a, double* kz40a, double* kz41a,
  double* kz42a, double* kz43a, double* kz44a, double* kz45a, double* kz46a,
  double* kz47a, double* kz48a, double* kz49a, double* kz50a, double* kz51a,
  double* kz52a, double* kz53a, double* kz54a, double* k125, double* k129,
  double* k130, double* k131, double* k132, double* k133, double* k134,
  double* k135, double* k136, double* k137, double* k148, double* k149,
  double* k150, double* k151, double* k152, double* k153, double* kz15,
  double* kz16, double* kz17, double* kz18, double* kz19, double* kz20,
  double* kz21, double* kz22, double* kz23, double* kz24, double* kz25,
  double* kz26, double* kz27, double* kz28, double* kz29, double* kz30,
  double* kz31, double* kz32, double* kz33, double* kz34, double* kz35,
  double* kz36, double* kz37, double* kz38, double* kz39, double* kz40,
  double* kz41, double* kz42, double* kz43, double* kz44, double* kz45,
  double* kz46, double* kz47, double* kz48, double* kz49, double* kz50,
  double* kz51, double* kz52, double* kz53, double* kz54, int* immulti,
  int* imabund, int* idspecies, int* itdmulti, int* idsub,
  gr_float* metal_loc_data_ptr, gr_float* metal_C13_data_ptr,
  gr_float* metal_C20_data_ptr, gr_float* metal_C25_data_ptr,
  gr_float* metal_C30_data_ptr, gr_float* metal_F13_data_ptr,
  gr_float* metal_F15_data_ptr, gr_float* metal_F50_data_ptr,
  gr_float* metal_F80_data_ptr, gr_float* metal_P170_data_ptr,
  gr_float* metal_P200_data_ptr, gr_float* metal_Y19_data_ptr, int* SN0_N,
  double* SN0_fSiM, double* SN0_fFeM, double* SN0_fMg2SiO4, double* SN0_fMgSiO3,
  double* SN0_fFe3O4, double* SN0_fAC, double* SN0_fSiO2D, double* SN0_fMgO,
  double* SN0_fFeS, double* SN0_fAl2O3, double* SN0_freforg,
  double* SN0_fvolorg, double* SN0_fH2Oice, double* SN0_r0SiM_data_ptr,
  double* SN0_r0FeM_data_ptr, double* SN0_r0Mg2SiO4_data_ptr,
  double* SN0_r0MgSiO3_data_ptr, double* SN0_r0Fe3O4_data_ptr,
  double* SN0_r0AC_data_ptr, double* SN0_r0SiO2D_data_ptr,
  double* SN0_r0MgO_data_ptr, double* SN0_r0FeS_data_ptr,
  double* SN0_r0Al2O3_data_ptr, double* SN0_r0reforg_data_ptr,
  double* SN0_r0volorg_data_ptr, double* SN0_r0H2Oice_data_ptr, int* gr_N,
  int* gr_Size, double* gr_dT, double* gr_Td, double* SN0_kpSiM_data_ptr,
  double* SN0_kpFeM_data_ptr, double* SN0_kpMg2SiO4_data_ptr,
  double* SN0_kpMgSiO3_data_ptr, double* SN0_kpFe3O4_data_ptr,
  double* SN0_kpAC_data_ptr, double* SN0_kpSiO2D_data_ptr,
  double* SN0_kpMgO_data_ptr, double* SN0_kpFeS_data_ptr,
  double* SN0_kpAl2O3_data_ptr, double* SN0_kpreforg_data_ptr,
  double* SN0_kpvolorg_data_ptr, double* SN0_kpH2Oice_data_ptr,
  double* h2dustSa, double* h2dustCa, gr_float* rhoH, double* grogra,
  double* dt, double* kdSiM, double* kdFeM, double* kdMg2SiO4, double* kdMgSiO3,
  double* kdFe3O4, double* kdAC, double* kdSiO2D, double* kdMgO, double* kdFeS,
  double* kdAl2O3, double* kdreforg, double* kdvolorg, double* kdH2Oice,
  double* tSiM, double* tFeM, double* tMg2SiO4, double* tMgSiO3, double* tFe3O4,
  double* tAC, double* tSiO2D, double* tMgO, double* tFeS, double* tAl2O3,
  double* treforg, double* tvolorg, double* tH2Oice, int* iuseH2shield,
  int* iH2shieldcustom, gr_float* f_shield_custom_data_ptr
);

void FORTRAN_NAME(rate_timestep_g)(
  double* dedot, double* HIdot, int* ispecies, gr_mask_type* anydust,
  gr_float* de_data_ptr, gr_float* HI_data_ptr, gr_float* HII_data_ptr,
  gr_float* HeI_data_ptr, gr_float* HeII_data_ptr, gr_float* HeIII_data_ptr,
  gr_float* d_data_ptr, gr_float* HM_data_ptr, gr_float* H2I_data_ptr,
  gr_float* H2II_data_ptr, int* in, int* jn, int* kn, int* is, int* ie, int* j,
  int* k, double* k1, double* k2, double* k3, double* k4, double* k5,
  double* k6, double* k7, double* k8, double* k9, double* k10, double* k11,
  double* k12, double* k13, double* k14, double* k15, double* k16, double* k17,
  double* k18, double* k19, double* k22, double* k24, double* k25, double* k26,
  double* k27, double* k28, double* k29, double* k30, double* k50, double* k51,
  double* k52, double* k53, double* k54, double* k55, double* k56, double* k57,
  double* k58, double* h2dust, double* ncrn, double* ncrd1, double* ncrd2,
  double* rhoH, double* k24shield, double* k25shield, double* k26shield,
  double* k28shield, double* k29shield, double* k30shield, double* k31shield,
  int* iradtrans, int* irt_honly, gr_float* kphHI_data_ptr,
  gr_float* kphHeI_data_ptr, gr_float* kphHeII_data_ptr, gr_mask_type* itmask,
  double* edot, double* chunit, double* dom, gr_float* metal_data_ptr,
  gr_float* HDI_data_ptr, int* imchem, gr_float* CI_data_ptr,
  gr_float* OI_data_ptr, gr_float* OH_data_ptr, gr_float* CO_data_ptr,
  gr_float* H2O_data_ptr, int* idissHDI, gr_float* kdissHDI_data_ptr,
  int* iionZ, gr_float* kphCI_data_ptr, gr_float* kphOI_data_ptr, int* idissZ,
  gr_float* kdissCO_data_ptr, gr_float* kdissOH_data_ptr,
  gr_float* kdissH2O_data_ptr
);

void FORTRAN_NAME(step_rate_g)(
  gr_float* de_data_ptr, gr_float* HI_data_ptr, gr_float* HII_data_ptr,
  gr_float* HeI_data_ptr, gr_float* HeII_data_ptr, gr_float* HeIII_data_ptr,
  gr_float* d_data_ptr, gr_float* HM_data_ptr, gr_float* H2I_data_ptr,
  gr_float* H2II_data_ptr, gr_float* DI_data_ptr, gr_float* DII_data_ptr,
  gr_float* HDI_data_ptr, double* dtit, int* in, int* jn, int* kn, int* is,
  int* ie, int* j, int* k, int* ispecies, gr_mask_type* anydust, double* k1,
  double* k2, double* k3, double* k4, double* k5, double* k6, double* k7,
  double* k8, double* k9, double* k10, double* k11, double* k12, double* k13,
  double* k14, double* k15, double* k16, double* k17, double* k18, double* k19,
  double* k22, double* k24, double* k25, double* k26, double* k27, double* k28,
  double* k29, double* k30, double* k50, double* k51, double* k52, double* k53,
  double* k54, double* k55, double* k56, double* k57, double* k58,
  double* h2dust, double* rhoH, double* k24shield, double* k25shield,
  double* k26shield, double* k28shield, double* k29shield, double* k30shield,
  double* k31shield, double* HIp, double* HIIp, double* HeIp, double* HeIIp,
  double* HeIIIp, double* dep, double* HMp, double* H2Ip, double* H2IIp,
  double* DIp, double* DIIp, double* HDIp, double* dedot_prev,
  double* HIdot_prev, int* iradtrans, int* irt_honly, gr_float* kphHI_data_ptr,
  gr_float* kphHeI_data_ptr, gr_float* kphHeII_data_ptr, gr_mask_type* itmask,
  gr_mask_type* itmask_metal, gr_float* DM_data_ptr, gr_float* HDII_data_ptr,
  gr_float* HeHII_data_ptr, int* imetal, gr_float* metal_data_ptr, int* imchem,
  int* idspecies, int* igrgr, int* idsub, gr_float* CI_data_ptr,
  gr_float* CII_data_ptr, gr_float* CO_data_ptr, gr_float* CO2_data_ptr,
  gr_float* OI_data_ptr, gr_float* OH_data_ptr, gr_float* H2O_data_ptr,
  gr_float* O2_data_ptr, gr_float* SiI_data_ptr, gr_float* SiOI_data_ptr,
  gr_float* SiO2I_data_ptr, gr_float* CH_data_ptr, gr_float* CH2_data_ptr,
  gr_float* COII_data_ptr, gr_float* OII_data_ptr, gr_float* OHII_data_ptr,
  gr_float* H2OII_data_ptr, gr_float* H3OII_data_ptr, gr_float* O2II_data_ptr,
  gr_float* Mg_data_ptr, gr_float* Al_data_ptr, gr_float* S_data_ptr,
  gr_float* Fe_data_ptr, gr_float* SiM_data_ptr, gr_float* FeM_data_ptr,
  gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, double* k125, double* k129, double* k130,
  double* k131, double* k132, double* k133, double* k134, double* k135,
  double* k136, double* k137, double* k148, double* k149, double* k150,
  double* k151, double* k152, double* k153, double* kz15, double* kz16,
  double* kz17, double* kz18, double* kz19, double* kz20, double* kz21,
  double* kz22, double* kz23, double* kz24, double* kz25, double* kz26,
  double* kz27, double* kz28, double* kz29, double* kz30, double* kz31,
  double* kz32, double* kz33, double* kz34, double* kz35, double* kz36,
  double* kz37, double* kz38, double* kz39, double* kz40, double* kz41,
  double* kz42, double* kz43, double* kz44, double* kz45, double* kz46,
  double* kz47, double* kz48, double* kz49, double* kz50, double* kz51,
  double* kz52, double* kz53, double* kz54, double* DMp, double* HDIIp,
  double* HeHIIp, double* CIp, double* CIIp, double* COp, double* CO2p,
  double* OIp, double* OHp, double* H2Op, double* O2p, double* SiIp,
  double* SiOIp, double* SiO2Ip, double* CHp, double* CH2p, double* COIIp,
  double* OIIp, double* OHIIp, double* H2OIIp, double* H3OIIp, double* O2IIp,
  double* Mgp, double* Alp, double* Sp, double* Fep, gr_float* SiMp,
  gr_float* FeMp, gr_float* Mg2SiO4p, gr_float* MgSiO3p, gr_float* Fe3O4p,
  gr_float* ACp, gr_float* SiO2Dp, gr_float* MgOp, gr_float* FeSp,
  gr_float* Al2O3p, gr_float* reforgp, gr_float* volorgp, gr_float* H2Oicep,
  double* kdSiM, double* kdFeM, double* kdMg2SiO4, double* kdMgSiO3,
  double* kdFe3O4, double* kdAC, double* kdSiO2D, double* kdMgO, double* kdFeS,
  double* kdAl2O3, double* kdreforg, double* kdvolorg, double* kdH2Oice,
  int* idissHDI, gr_float* kdissHDI_data_ptr, int* iionZ,
  gr_float* kphCI_data_ptr, gr_float* kphOI_data_ptr, int* idissZ,
  gr_float* kdissCO_data_ptr, gr_float* kdissOH_data_ptr,
  gr_float* kdissH2O_data_ptr
);

void FORTRAN_NAME(make_consistent_g)(
  gr_float* de_data_ptr, gr_float* HI_data_ptr, gr_float* HII_data_ptr,
  gr_float* HeI_data_ptr, gr_float* HeII_data_ptr, gr_float* HeIII_data_ptr,
  gr_float* HM_data_ptr, gr_float* H2I_data_ptr, gr_float* H2II_data_ptr,
  gr_float* DI_data_ptr, gr_float* DII_data_ptr, gr_float* HDI_data_ptr,
  gr_float* metal_data_ptr, gr_float* dust_data_ptr, gr_float* d_data_ptr,
  int* is, int* ie, int* js, int* je, int* ks, int* ke, int* in, int* jn,
  int* kn, int* ispecies, int* imetal, double* fh, double* dtoh,
  int* idustfield, int* imchem, int* igrgr, double* dom, gr_float* DM_data_ptr,
  gr_float* HDII_data_ptr, gr_float* HeHII_data_ptr, gr_float* CI_data_ptr,
  gr_float* CII_data_ptr, gr_float* CO_data_ptr, gr_float* CO2_data_ptr,
  gr_float* OI_data_ptr, gr_float* OH_data_ptr, gr_float* H2O_data_ptr,
  gr_float* O2_data_ptr, gr_float* SiI_data_ptr, gr_float* SiOI_data_ptr,
  gr_float* SiO2I_data_ptr, gr_float* CH_data_ptr, gr_float* CH2_data_ptr,
  gr_float* COII_data_ptr, gr_float* OII_data_ptr, gr_float* OHII_data_ptr,
  gr_float* H2OII_data_ptr, gr_float* H3OII_data_ptr, gr_float* O2II_data_ptr,
  gr_float* Mg_data_ptr, gr_float* Al_data_ptr, gr_float* S_data_ptr,
  gr_float* Fe_data_ptr, gr_float* SiM_data_ptr, gr_float* FeM_data_ptr,
  gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
  gr_float* Fe3O4_data_ptr, gr_float* AC_data_ptr, gr_float* SiO2D_data_ptr,
  gr_float* MgO_data_ptr, gr_float* FeS_data_ptr, gr_float* Al2O3_data_ptr,
  gr_float* reforg_data_ptr, gr_float* volorg_data_ptr,
  gr_float* H2Oice_data_ptr, int* immulti, int* imabund, int* idspecies,
  int* itdmulti, int* idsub, gr_float* metal_loc_data_ptr,
  gr_float* metal_C13_data_ptr, gr_float* metal_C20_data_ptr,
  gr_float* metal_C25_data_ptr, gr_float* metal_C30_data_ptr,
  gr_float* metal_F13_data_ptr, gr_float* metal_F15_data_ptr,
  gr_float* metal_F50_data_ptr, gr_float* metal_F80_data_ptr,
  gr_float* metal_P170_data_ptr, gr_float* metal_P200_data_ptr,
  gr_float* metal_Y19_data_ptr, int* SN0_N, double* SN0_XC, double* SN0_XO,
  double* SN0_XMg, double* SN0_XAl, double* SN0_XSi, double* SN0_XS,
  double* SN0_XFe, double* SN0_fC, double* SN0_fO, double* SN0_fMg,
  double* SN0_fAl, double* SN0_fSi, double* SN0_fS, double* SN0_fFe
);


#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* FORTRAN_FN_DECLARATIONS_HPP */
