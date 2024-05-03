#include "grackle_macros.h"

// TODO: FORTRAN_NAME not needed for unit tests
extern void FORTRAN_NAME(cool1d_cloudy_g)(double* d, double* rhoH, double* metallicity,
                                          int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
                                          double* logtem, double* edot, double* comp2, double* dom, double* zr,
                                          int* icmbTfloor, int* iClHeat, int* iZscale,
                                          long long* clGridRank, long long*  clGridDim,
                                          double* clPar1, double* clPar2, double* clPar3,
                                          long long* clDataSize, double* clCooling, double* clHeating, int* itmask);
