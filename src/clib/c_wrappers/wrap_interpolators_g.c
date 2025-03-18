#include "grackle_macros.h"

// TODO: FORTRAN_NAME not needed for unit tests
extern void FORTRAN_NAME(interpolate_1d_g)(double *input1,
                                           long long *gridDim,
                                           double *gridPar1, double *dgridPar1,
                                           long long *dataSize, double *dataField,
                                           double *value);

extern void FORTRAN_NAME(interpolate_2d_g)(double *input1, double *input2,
                                           long long *gridDim,
                                           double *gridPar1, double *dgridPar1,
                                           double *gridPar2, double *dgridPar2,
                                           long long *dataSize, double *dataField,
                                           double *value);

extern void FORTRAN_NAME(interpolate_3d_g)(double *input1, double *input2, double *input3,
                                           long long *gridDim,
                                           double *gridPar1, double *dgridPar1,
                                           double *gridPar2, double *dgridPar2,
                                           double *gridPar3, double *dgridPar3,
                                           long long *dataSize, double *dataField,
                                           double *value);

extern void FORTRAN_NAME(interpolate_3dz_g)(double *input1, double *input2, double *input3,
                                            long long *gridDim,
                                            double *gridPar1, double *dgridPar1,
                                            double *gridPar2, long long *index2,
                                            double *gridPar3, double *dgridPar3,
                                            long long *dataSize, double *dataField,
                                            int *end_int,
                                            double *value);

extern void FORTRAN_NAME(interpolate_2df3d_g)(double *input1, double *input3,
                                              long long *gridDim,
                                              double *gridPar1, double *dgridPar1,
                                              long long *index2,
                                              double *gridPar3, double *dgridPar3,
                                              long long *dataSize, double *dataField,
                                              double *value);

extern void FORTRAN_NAME(interpolate_4d_g)(double *input1, double *input2, double *input3, double *input4,
                                           long long *gridDim,
                                           double *gridPar1, double *dgridPar1,
                                           double *gridPar2, double *dgridPar2,
                                           double *gridPar3, double *dgridPar3,
                                           double *gridPar4, double *dgridPar4,
                                           long long *dataSize, double *dataField,
                                           double *value);

extern void FORTRAN_NAME(interpolate_5d_g)(double *input1, double *input2, double *input3, double *input4, double *input5,
                                           long long *gridDim,
                                           double *gridPar1, double *dgridPar1,
                                           double *gridPar2, double *dgridPar2,
                                           double *gridPar3, double *dgridPar3,
                                           double *gridPar4, double *dgridPar4,
                                           double *gridPar5, double *dgridPar5,
                                           long long *dataSize, double *dataField,
                                           double *value);
