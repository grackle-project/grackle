#include "grackle_macros.h"

// TODO: FORTRAN_NAME not needed for unit tests
extern void FORTRAN_NAME(cool1d_cloudy_g)(double* d, double* rhoH, double* metallicity,
                                          int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
                                          double* logtem, double* edot, double* comp2, double* dom, double* zr,
                                          int* icmbTfloor, int* iClHeat, int* iZscale,
                                          long long* clGridRank, long long*  clGridDim,
                                          double* clPar1, double* clPar2, double* clPar3,
                                          long long* clDataSize, double* clCooling, double* clHeating, int* itmask);

extern void FORTRAN_NAME(cool_rank1_interpolation)(int* get_heat, int* icmbTfloor,
                                                   double* log10tem_i, double* log10_tCMB,
                                                   long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                   double* clPar1, double* clCooling, double* clHeating,
                                                   double* dclPar_1,
                                                   double*  edot_met_i);

extern void FORTRAN_NAME(cool_rank2_interpolation)(int* get_heat, int* icmbTfloor,
                                                   double* log_n_h_i, double* log10tem_i, double* log10_tCMB,
                                                   long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                   double* clPar1, double* clPar2,
                                                   double* clCooling, double* clHeating,
                                                   double* dclPar_1, double* dclPar_2,
                                                   double* edot_met_i);

extern void FORTRAN_NAME(cool_rank3_interpolation)(int* get_heat, int* icmbTfloor,
                                                   double* log_n_h_i, double* zr, double* log10tem_i, double* log10_tCMB,
                                                   long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                   double* clPar1, double* clPar2, double* clPar3,
                                                   double* clCooling, double* clHeating,
                                                   double* dclPar_1, long long* zindex, double* dclPar_3,
                                                   double* edot_met_i);
