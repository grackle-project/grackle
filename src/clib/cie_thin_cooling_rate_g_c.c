/*******************************************************************************/
/*
 * ----------Function to Compute Optically Thin CIE Cooling Rate---------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 17/2/2021
 * Date of Last Modification: 17/2/2021
 * 
 * Purpose:
 *  Computes the optically thin cooling rate due to CIE as prescribed in 
 *  (Ripamonti and Abel 2003).
 * 
 * Inputs:
 *  T : Temperature of the gas (in Kelvin)
 *  
 * Outputs:
 *  cieRate : Optically thin CIE cooling rate.
 * 
 * Units:
 *  erg / s * cm^3 / gram (gram in H2 molecules) --> Units converted to erg/s cm^-3
 *                                                   by multiplying by rho * rho_mh
 *                                                   in calc_rates_g.c
 * 
 * This code is a refactoring of the original Grackle fortran function of the same name.
 */
/*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "cie_thin_cooling_rate_tables.h"

//* Define the function which will return the CIE thin cooling rate.
double cie_thin_cooling_rate_g_c(double T){

    // Maximal and minimal indices in the cie tables.
    int minInd = 0;
    int maxInd = 287;

    //* Compute rough extrapolations for extreme temperatures.
    // Low temperatures extrapolated with fourth power.
    if (T <= t_cie_c[minInd]) {
        return cie_table_c[maxInd]*pow(T/t_cie_c[maxInd], 4);
    }
    // High temperatures extrapolated with third power.
    if (T >= t_cie_c[maxInd]) {
        return cie_table_c[maxInd]*pow(T/t_cie_c[maxInd], 3);
    }

    //* Compute CIE cooling rate for moderate temperatures.
    for (int i = 0; i <= 200; i++) {
        // Calculate working index.
        int ind = (maxInd + minInd) / 2;

        // Update extreme indices.
        if (T >= t_cie_c[ind]) {
            minInd = ind;
        } else {
            maxInd = ind;
        }

        // Return rate if it has been calculated.
        if ( maxInd - minInd <= 1 ) {
            return (cie_table_c[maxInd] * ( T - t_cie_c[minInd] ) \
                   + cie_table_c[minInd] * ( t_cie_c[maxInd] - T)) \
                   / ( t_cie_c[maxInd] - t_cie_c[minInd] );
        }
    }

    //! The original fortran code seems to return 0 at the end. I am not sure why this happens right now but will put it here for now.
    return 0;
}