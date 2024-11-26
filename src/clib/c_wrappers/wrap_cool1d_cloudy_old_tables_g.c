#include "grackle_macros.h"

extern void FORTRAN_NAME(cool1d_cloudy_old_tables_g)(double* d, double* de, double* rhoH, double* metallicity,
                                                     int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
                                                     double* logtem, double* edot, double* comp2, int* ispecies, 
                                                     double* dom, double* zr, int* icmbTfloor, int* iClHeat, double* clEleFra,
                                                     long long* clGridRank, long long*  clGridDim,
                                                     double* clPar1, double* clPar2, double* clPar3, double* clPar4, double* clPar5,
                                                     long long* clDataSize, double* clCooling, double* clHeating, int* itmask);

