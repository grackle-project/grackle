#include "grackle_macros.h"

extern void FORTRAN_NAME(calc_parameter_value_slopes)(double* dclPar, long long* clGridRank, long long* clGridDim,
                                                      double* clPar1, double* clPar2, double* clPar3,
                                                      double* clPar4, double* clPar5);

extern void FORTRAN_NAME(compute_redshift_dimension)(double* zr, int* clGridDim_2, double*  clPar2,
                                                     int* zindex, int* end_int, int* get_heat);
