/***********************************************************************
/
/ Calculate temperature field of a 1D slice from a cloudy table
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdint.h> // int32_t
#include <stdio.h>
#include <math.h> // log, log10
#include "../grackle_chemistry_data.h"
#include "../phys_constants.h"
#include "./interop_funcs.h"

void cool1d_cloudy_g(
        const gr_float* d, // 3D arrays
        const double* rhoH, const gr_float* metallicity, // 1D array
        int in, int jn, int kn, int is, int ie, int j, int k,
        const double* logtem, double * GR_RESTRICT edot, // 1D array
        double comp2, double dom, double zr,
        int icmbTfloor, int iClHeat, int iZscale,
        long long clGridRank,
        const long long* clGridDim,
        const double* clPar1, const double* clPar2, const double* clPar3,
        long long clDataSize, const double* clCooling, const double* clHeating,
        const int32_t* itmask)
{
  int get_heat = iClHeat;
  long long zindex;
  const double inv_log10 = 1.0 / log(10.0);
  const double log10_tCMB = log10(comp2);
  gr_int64 end_int = 0;

  // local slices
  double* log_n_h = malloc(sizeof(double)*in);
  double* log_cool = malloc(sizeof(double)*in);
  double* log_cool_cmb = malloc(sizeof(double)*in);
  double* log_heat = malloc(sizeof(double)*in);
  double* edot_met = malloc(sizeof(double)*in);
  double* log10tem = malloc(sizeof(double)*in);

  // Calculate parameter value slopes
  const double dclPar[3] = {
     ((clPar1[clGridDim[0]-1] - clPar1[0]) / (double)(clGridDim[0] - 1)),

     (clGridRank > 1)
     ? ((clPar2[clGridDim[1]-1] - clPar2[0]) / (double)(clGridDim[1] - 1)) : 0,

     (clGridRank > 2)
     ? ((clPar3[clGridDim[2]-1] - clPar3[0]) / (double)(clGridDim[2] - 1)) : 0
    };

  // do work
  for (int i = is + 1; i <= (ie + 1); i++) {
    if ( !itmask[i-1] ) { continue; }

    log10tem[i-1] = logtem[i-1] * inv_log10;
    log_n_h[i-1] = log10(rhoH[i-1] * dom); // Calculate proper log(n_H)

    // Calculate index for redshift dimension - intentionally kept 1-indexed
    // TODO: hoist this out of the loop
    zindex = find_zindex(zr, clGridRank, clGridDim, clPar2);
    if ((clGridRank > 2) && (zindex == clGridDim[1])){
      end_int = 1;
      get_heat = 0;
    }

    // Call interpolation functions to get heating and cooling
    if (clGridRank == 1) { // Interpolate over temperature.
      interpolate_1d_g(log10tem[i-1], clGridDim,
                       clPar1, dclPar[0],
                       clDataSize, clCooling, &(log_cool[i-1]));
      // pow returns NaN when a negative base is raised to a non-integer power
      edot_met[i-1] = -1.0 * pow(10.0,log_cool[i-1]);

      // Ignore CMB term if T >> T_CMB
      if ((icmbTfloor == 1) && ((log10tem[i-1] - log10_tCMB) < 2.0)) {
        interpolate_1d_g(log10_tCMB, clGridDim,
                         clPar1, dclPar[0],
                         clDataSize, clCooling, &(log_cool_cmb[i-1]));
        edot_met[i-1] += pow(10.0,log_cool_cmb[i-1]);
      }

      if (get_heat == 1){
        interpolate_1d_g(log10tem[i-1], clGridDim,
                         clPar1, dclPar[0],
                         clDataSize, clHeating, &(log_heat[i-1]));
        edot_met[i-1] += pow(10.0,log_heat[i-1]);
      }
    } else if ( clGridRank == 2) { // Interpolate over density & temperature.
      interpolate_2d_g(log_n_h[i-1], log10tem[i-1],
                       clGridDim,
                       clPar1, dclPar[0],
                       clPar2, dclPar[1],
                       clDataSize, clCooling, &(log_cool[i-1]));
      // pow returns NaN when a negative base is raised to a non-integer power
      edot_met[i-1] = -1.0 * pow(10.0,log_cool[i-1]);

      // Ignore CMB term if T >> T_CMB
      if ((icmbTfloor == 1) && ((log10tem[i-1] - log10_tCMB) < 2.0)) {
        interpolate_2d_g(log_n_h[i-1], log10_tCMB,
                         clGridDim,
                         clPar1, dclPar[0],
                         clPar2, dclPar[1],
                         clDataSize, clCooling, &(log_cool_cmb[i-1]));
        edot_met[i-1] += pow(10.0,log_cool_cmb[i-1]);
      }

      if (get_heat == 1){
        interpolate_2d_g(log_n_h[i-1], log10tem[i-1],
                         clGridDim,
                         clPar1, dclPar[0],
                         clPar2, dclPar[1],
                         clDataSize, clHeating, &(log_heat[i-1]));
        edot_met[i-1] += pow(10.0,log_heat[i-1]);
      }

    } else if ( clGridRank == 3) { // Interpolate over density, redshift,
                                   // & temperature.
      interpolate_3dz_g(log_n_h[i-1], zr, log10tem[i-1],
                        clGridDim,
                        clPar1, dclPar[0],
                        clPar2, zindex,
                        clPar3, dclPar[2],
                        clDataSize, clCooling,
                        end_int, &(log_cool[i-1]));
      // pow returns NaN when a negative base is raised to a non-integer power
      edot_met[i-1] = -1.0 * pow(10.0,log_cool[i-1]);

      // Ignore CMB term if T >> T_CMB
      if ((icmbTfloor == 1) && ((log10tem[i-1] - log10_tCMB) < 2.0)) {
        interpolate_3dz_g(log_n_h[i-1], zr, log10_tCMB,
                          clGridDim,
                          clPar1, dclPar[0],
                          clPar2, zindex,
                          clPar3, dclPar[2],
                          clDataSize, clCooling,
                          end_int, &(log_cool_cmb[i-1]));
        edot_met[i-1] += pow(10.0,log_cool_cmb[i-1]);
      }

      if (get_heat == 1){
        interpolate_3dz_g(log_n_h[i-1], zr, log10tem[i-1],
                          clGridDim,
                          clPar1, dclPar[0],
                          clPar2, zindex,
                          clPar3, dclPar[2],
                          clDataSize, clHeating,
                          end_int, &(log_heat[i-1]));
        edot_met[i-1] += pow(10.0,log_heat[i-1]);
      }

    } else {
      // no need to be within an openmp critical section
      printf("Maximum cooling data grid rank is 3!\n");
    }

    // scale cooling by metallicity
    if (iZscale == 1) {
      edot_met[i-1] = edot_met[i-1] * metallicity[i-1];
    }

    edot[i-1] += (edot_met[i-1] * rhoH[i-1] * rhoH[i-1]);
  }


  free(log_n_h);
  free(log_cool);
  free(log_cool_cmb);
  free(log_heat);
  free(edot_met);
  free(log10tem);
}
