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

void interpolate_1d_g(double input1,
                      const gr_int64 * GR_RESTRICT gridDim, // 1 elements
                      const double * GR_RESTRICT gridPar1, double dgridPar1,
                      gr_int64 dataSize, const double * GR_RESTRICT dataField,
                      double * GR_RESTRICT value);

void interpolate_2d_g(double input1, double input2,
                      const gr_int64 * GR_RESTRICT gridDim, // 2 elements
                      const double * GR_RESTRICT gridPar1, double dgridPar1,
                      const double * GR_RESTRICT gridPar2, double dgridPar2,
                      gr_int64 dataSize, const double* dataField,
                      double * GR_RESTRICT value);

void interpolate_3dz_g(double input1, double input2, double input3,
                       const gr_int64 * GR_RESTRICT gridDim, // 3 elements
                       const double * GR_RESTRICT gridPar1, double dgridPar1,
                       const double * GR_RESTRICT gridPar2, gr_int64 index2,
                       const double * GR_RESTRICT gridPar3, double dgridPar3,
                       gr_int64 dataSize, const double * GR_RESTRICT dataField,
                       gr_int64 end_int, double * GR_RESTRICT value);

void interpolate_3d_g(double input1, double input2, double input3,
                      const gr_int64 * GR_RESTRICT gridDim, // 3 elements
                      const double * GR_RESTRICT gridPar1, double dgridPar1,
                      const double * GR_RESTRICT gridPar2, double dgridPar2,
                      const double * GR_RESTRICT gridPar3, double dgridPar3,
                      gr_int64 dataSize, const double * GR_RESTRICT dataField,
                      double* value);

void interpolate_4d_g(double input1, double input2, double input3,
                      double input4,
                      const gr_int64 * GR_RESTRICT gridDim, // 4 elements
                      const double * GR_RESTRICT gridPar1, double dgridPar1,
                      const double * GR_RESTRICT gridPar2, double dgridPar2,
                      const double * GR_RESTRICT gridPar3, double dgridPar3,
                      const double * GR_RESTRICT gridPar4, double dgridPar4,
                      gr_int64 dataSize, const double * GR_RESTRICT dataField,
                      double * GR_RESTRICT value);

void interpolate_5d_g(double input1, double input2, double input3,
                      double input4, double input5,
                      const gr_int64 * GR_RESTRICT gridDim, // 5 elements
                      const double * GR_RESTRICT gridPar1, double dgridPar1,
                      const double * GR_RESTRICT gridPar2, double dgridPar2,
                      const double * GR_RESTRICT gridPar3, double dgridPar3,
                      const double * GR_RESTRICT gridPar4, double dgridPar4,
                      const double * GR_RESTRICT gridPar5, double dgridPar5,
                      gr_int64 dataSize, const double * GR_RESTRICT dataField,
                      double * GR_RESTRICT value);

#endif /* __INTEROP_FUNCS_H */
