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
void FORTRAN_NAME(interpolate_1d_g)(
        const double* input1, const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_2d_g)(
        const double* input1, const double* input2,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const double* dgridPar2,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_3dz_g)(
        const double* input1, const double* input2, const double* input3,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const gr_int64* index2,
        const double* gridPar3, const double* dgridPar3,
        const gr_int64* dataSize, const double* dataField,
        const gr_int64* end_int,
        double* value);

// in the future, we will define GR_RESTRICT to be the restrict keyword
// introduced in C99 (or compiler-specific alternatives). For now, haven't done
// this to avoid headaches with the C++ unit-tests (restrict is not recognized
// by C++, but compiler-specific alternatives are)
//
// I don't actually think it would be problematic to use restrict within source
// files and skip it in header files (as long as assumptions about no-pointer
// aliasing) are enforced 
#ifndef GR_RESTRICT
#define GR_RESTRICT /* ... */
#endif /* GR_RESTRICT */


// helper function that retrieves index for redshift dimension (of cloudy
// tables) via bisection
// - the index is one-indexed
// - the names of variables have not been changed for backwards compatibility
//   (it may seem counter-intuitive that clGridDim[1] gives the length of
//    clPar2, but that's because in Fortran you would access clGridDim(2) )
// - NOTE: since we define this function in a header, we must declare it as
//   static inline (in C++ we could just declare it as inline)
static inline long long find_zindex(double zr, long long clGridRank,
                                    const long long* clGridDim,
                                    const double* clPar2){
  if (clGridRank > 2){
    long long zindex;
    if (zr <= clPar2[0]) {
      zindex = 1;
    } else if (zr >= clPar2[clGridDim[1]-2]) {
      zindex = clGridDim[1];
    } else if (zr >= clPar2[clGridDim[1]-3]) {
      zindex = clGridDim[1] - 2;
    } else {
      zindex = 1;
      long long zhighpt = clGridDim[1] - 2;
      while ((zhighpt - zindex) > 1) {
        long long zmidpt = (long long)((zhighpt + zindex) / 2);
        if (zr >= clPar2[zmidpt-1]){
          zindex = zmidpt;
        } else {
          zhighpt = zmidpt;
        }
      }
    }
    return zindex;
  } else {
    return 1;
  }
}


/// Calculate temperature and mean molecular weight for tabulated cooling.
/// @details This performs a calculation for a 1D slice from a 3D field.
///
/// @author Britton Smith
/// @author Matthew Abruzzo
/// @date May, 2015
/// @remark This was originally written by Britton Smith in May 2015 in
///     Fortran. It was transcribed to C by Matthew Abruzzo in April 2023.
///
/// @param[in]  d  3D density field
/// @param[in]  metal  3D metal density field
/// @param[in]  e  3D specific internal energy field
/// @param[in]  rhoH  (precomputed) total H mass density. This only holds
///     values for the 1D slice
/// @param[in]  in,jn,kn  dimensions of 3D fields (1D array hold ``in`` items)
/// @param[in]  is,ie  start and (inclusive) end indices of active region
///     (zero-based)
/// @param[in]  j,k  indices along other dimensions (one-based)
/// @param[out] tgas  1D array to store output temperature values
/// @param[out] mmw  1D array to store output mean molecular weight values
/// @param[in]  dom  unit conversion to proper number density in code units
/// @param[in]  zr  current redshift
/// @param[in]  temstart, temend  start and end of temperature range for rate
///     table
/// @param[in]  gamma  adiabatic index
/// @param[in]  utem  temperature units
/// @param[in]  imetal  flag if metal field is active (0 = no, 1 = yes)
/// @param[in]  clGridRank  rank of cloudy cooling data grid
/// @param[in]  clGridDim  array containing dimensions of cloudy data
/// @param[in]  clPar1, clPar2, clPar3  arrays containing cloudy grid parameter
///     values.
/// @param[in]  clDataSize  total size of flattened mmw data array
/// @param[in]  clMMW  cloudy mmw data
/// @param[in]  itmask  iteration mask
///
/// @note
/// None of the pointers are allowed to alias. Other than ``clPar2`` (when
/// ``clGridRank < 2``) or ``clPar3`` (when ``clGridRank < 3``), no parameters
/// should be passed ``NULL``
void calc_temp1d_cloudy_g(
        const gr_float* d, const gr_float* metal, const gr_float* e,
        const double* rhoH,
        int in, int jn, int kn, int is, int ie, int j, int k,
        GR_RESTRICT double* tgas, GR_RESTRICT double* mmw,
        double dom, double zr, double temstart, double temend, double gamma,
        double utem, int imetal,
        long long clGridRank,
        const long long* clGridDim,
        const double* clPar1,
        const double* clPar2,
        const double* clPar3,
        long long clDataSize,
        const double* clMMW,
        const int32_t* itmask);

/// Solve cloudy cooling by interpolating from data tables.
/// @details This performs a calculation for a 1D slice from a 3D field.
///
/// @author Britton Smith
/// @author Matthew Abruzzo
/// @date September, 2009
/// @remark This was originally written by Britton Smith in September 2009 in
///     Fortran. It was transcribed to C by Matthew Abruzzo in May 2023.
///
/// @param[in]     d  3D density field
/// @param[in]     rhoH  (precomputed) total H mass density. This only holds
///     values for the 1D slice
/// @param[in]     metallicity  (precomputed) metallicity. This only holds
///     values for the 1D slice
/// @param[in]     in,jn,kn  dimensions of 3D fields (1D array hold ``in``
///     items)
/// @param[in]     is,ie  start and (inclusive) end indices of active region
///     (zero-based)
/// @param[in]     j,k  indices along other dimensions (one-based)
/// @param[in]     logtem  (precomputed) natural log of temperature values.
///     This only holds values for the 1D slice
/// @param[in,out] edot the heating/cooling contributions computed in this
///     function are used to update this array (the net contributions are added
///     to prexisting values. This only holds values for the 1D slice
/// @param[in]     comp2 Temperature of the CMB at the redshift given by ``zr``
/// @param[in]     dom  unit conversion to proper number density in code units
/// @param[in]     zr  current redshift
/// @param[in]     icmbTfloor flag to include temperature floor from cmb
/// @param[in]     iClHeat flag to include cloudy heating
/// @param[in]     iZscale flag to scale cooling by metallicity
/// @param[in]     clGridRank  rank of cloudy cooling data grid
/// @param[in]     clGridDim  array containing dimensions of cloudy data
/// @param[in]     clPar1, clPar2, clPar3  arrays containing cloudy grid
///     parameter values.
/// @param[in]     clDataSize  total size of flattened cloudy data arrays
/// @param[in]     clCooling  cloudy cooling data
/// @param[in]     clHeating  cloudy heating data
/// @param[in]     itmask  iteration mask
void cool1d_cloudy_g(
        const gr_float* d, // 3D arrays
        const double* rhoH, const gr_float* metallicity, // 1D array
        int in, int jn, int kn, int is, int ie, int j, int k,
        const double* logtem, GR_RESTRICT double* edot, // 1D array
        double comp2, double dom, double zr,
        int icmbTfloor, int iClHeat, int iZscale,
        long long clGridRank,
        const long long* clGridDim,
        const double* clPar1, const double* clPar2, const double* clPar3,
        long long clDataSize, const double* clCooling, const double* clHeating,
        const int32_t* itmask);

#endif /* __INTEROP_FUNCS_H */
