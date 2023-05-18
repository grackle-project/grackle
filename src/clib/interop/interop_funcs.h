/***********************************************************************
/
/ C interfaces of internal interop functions (and helper functions)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/
#ifndef __INTEROP_FUNCS_H
#define __INTEROP_FUNCS_H

#include <stdint.h> // int32_t

#include "grackle_macros.h"
#include "grackle_types.h"

typedef long long gr_int64;

void interpolate_1d_g(const double* input1_p,
                      const gr_int64* gridDim, // 1 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value);

void interpolate_2d_g(const double* input1_p, const double* input2_p,
                      const gr_int64* gridDim, // 2 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value);

void interpolate_3dz_g(const double* input1_p, const double* input2_p,
                       const double* input3_p,
                       const gr_int64* gridDim, // 3 elements
                       const double* gridPar1, const double* dgridPar1_p,
                       const double* gridPar2, const gr_int64* index2_p,
                       const double* gridPar3, const double* dgridPar3_p,
                       const gr_int64* dataSize, const double* dataField,
                       const gr_int64* end_int_p, double* value);

void interpolate_3d_g(const double* input1_p, const double* input2_p,
                      const double* input3_p,
                      const gr_int64* gridDim, // 3 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value);

void interpolate_4d_g(const double* input1_p, const double* input2_p,
                      const double* input3_p, const double* input4_p,
                      const gr_int64* gridDim, // 4 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const double* gridPar4, const double* dgridPar4_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value);

void interpolate_5d_g(const double* input1_p, const double* input2_p,
                      const double* input3_p, const double* input4_p,
                      const double* input5_p,
                      const gr_int64* gridDim, // 5 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const double* gridPar4, const double* dgridPar4_p,
                      const double* gridPar5, const double* dgridPar5_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value);

#endif /* __INTEROP_FUNCS_H */
