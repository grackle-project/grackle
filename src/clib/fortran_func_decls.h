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
  int* immulti, int* imabund, int* idspecies, int* igrgr,
  const gr_mask_type* itmask, int* in, int* jn, int* kn, int* is, int* ie,
  int* j, int* k, double* dom, gr_float* d_data_ptr, gr_float* SiM_data_ptr,
  gr_float* FeM_data_ptr, gr_float* Mg2SiO4_data_ptr, gr_float* MgSiO3_data_ptr,
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

void FORTRAN_NAME(cool1d_cloudy_old_tables_g)(
  gr_float* d_data_ptr, gr_float* de_data_ptr, double* rhoH,
  double* metallicity, int* in, int* jn, int* kn, int* is, int* ie, int* j,
  int* k, double* logtem, double* edot, double* comp2, int* ispecies,
  double* dom, double* zr, int* icmbTfloor, int* iClHeat, double* clEleFra,
  long long* clGridRank, long long* clGridDim, double* clPar1, double* clPar2,
  double* clPar3, double* clPar4, double* clPar5, long long* clDataSize,
  double* clCooling, double* clHeating, gr_mask_type* itmask
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
