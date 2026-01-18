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

void FORTRAN_NAME(calc_grain_size_increment_species_1d)(
  const int* igrgr, const gr_mask_type* itmask, const int* SN0_N, int* in, int* jn, int* kn,
  int* is, int* ie, int* j, int* k, double* dom, gr_float* d_data_ptr, int* nSN,
  const gr_float* dsp_data_ptr, gr_float* SN_metal_data_ptr, double* SN_fsp,
  double* SN_r0sp_data_ptr, double* ssp, double* sgsp, double* alsp_data_ptr,
  const int* gr_N, const int* gr_Size, const double* gr_dT, const double* gr_Td,
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
