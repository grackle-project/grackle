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

void calc_temp1d_cloudy_g(
        const gr_float* d, const gr_float* metal, // 3D arrays
        const gr_float* e, // 3D array
        const double* rhoH, // 1D array
        int in, int jn, int kn, int is, int ie, int j, int k,
        double * GR_RESTRICT tgas, double * GR_RESTRICT mmw, // 1D array
        double dom, double zr, double temstart, double temend, double gamma,
        double utem, int imetal,
        long long clGridRank,
        const long long* clGridDim,
        const double* clPar1,
        const double* clPar2,
        const double* clPar3,
        long long clDataSize,
        const double* clMMW,
        const int32_t* itmask)
{
  const double mu_metal = 16.0;
  const int ti_max = 20;

  const double inv_log10 = 1.0 / log(10.0);

  // Calculate parameter value slopes
  const double dclPar[3] = {
     ((clPar1[clGridDim[0]-1] - clPar1[0]) / (double)(clGridDim[0] - 1)),

     (clGridRank > 1)
     ? ((clPar2[clGridDim[1]-1] - clPar2[0]) / (double)(clGridDim[1] - 1)) : 0,

     (clGridRank > 2)
     ? ((clPar3[clGridDim[2]-1] - clPar3[0]) / (double)(clGridDim[2] - 1)) : 0
    };

  // Calculate index for redshift dimension - intentionally kept 1-indexed
  const long long zindex = find_zindex(zr, clGridRank, clGridDim, clPar2);
  const gr_int64 end_int = ((clGridRank > 2) && (zindex == clGridDim[1]));

  for (int i = is + 1; i <= (ie + 1); i++) {
    if ( !itmask[i-1] ) { continue; }

    const double log_n_h = log10(rhoH[i-1] * dom); // Calculate proper log(n_H)

    const long long ind_3D = (i-1) + in *( (j-1) + jn * (k-1));

    double munew = 1.0;
    double muold;
    for (int ti = 1; ti <= ti_max; ti++) {
      muold = munew;

      tgas[i-1] = max((gamma - 1.0) * e[ind_3D] * munew * utem, temstart);
      // the original version doesn't use log10*(tgas[i-1]) either
      const double log10tem = log(tgas[i-1]) * inv_log10;

      // Call interpolation functions to get mmw
      if (clGridRank == 1) { // Interpolate over temperature.
        interpolate_1d_g(log10tem, clGridDim,
                         clPar1, dclPar[0],
                         clDataSize, clMMW, &munew);
      } else if ( clGridRank == 2) { // Interpolate over density & temperature.
        interpolate_2d_g(log_n_h, log10tem,
                         clGridDim,
                         clPar1, dclPar[0],
                         clPar2, dclPar[1],
                         clDataSize, clMMW, &munew);
      } else if (clGridRank == 3) { // Interpolate over density, redshift,
                                    // & temperature.
        interpolate_3dz_g(log_n_h, zr, log10tem,
                          clGridDim,
                          clPar1, dclPar[0],
                          clPar2, zindex,
                          clPar3, dclPar[2],
                          clDataSize, clMMW,
                          end_int, &munew);
      } else {
        // no need to be within an openmp critical section
        printf("Maximum mmw data grid rank is 3!\n");
      }

      munew = 0.5 * (munew + muold);
      tgas[i-1] = tgas[i-1] * munew / muold;

      if (fabs((munew/muold) - 1.0) < 1.e-2) {
        muold = munew;

        if (imetal == 1){
          munew = d[ind_3D] / ((d[ind_3D] - metal[ind_3D])/ munew +
                               metal[ind_3D] / mu_metal);
          tgas[i-1] = tgas[i-1] * munew / muold;
        }

        mmw[i-1] = munew;
        goto completed_T_calc; // TODO: refactor to remove this goto statement!
      }
    }

    mmw[i-1] = munew;

    // no need to put this in an omp critical statement
    fprintf(stderr,
            "Mean molecular weight not converged! %#.16g, %#.16g, %#.16g\n",
            munew, muold, fabs((munew/muold) - 1.0));

  completed_T_calc:
    continue;

  }
}
