/*********************************************************************************************************************/
/*
 * ------------------------Function to Compute Density-Dependent H2 Dissociation Rate-----------------------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 16/2/2021
 * Date of Last Modification: 16/2/2021
 * 
 * Purpose:
 *  To generate the (density-dependent) rate coefficients for the collisional dissociation of H2 by HI,
 *  seven (temperature dependent) functions must first be computed. This code computes these seven
 *  functions.
 * 
 *  Within the code the is a flag 'idt'. If idt = 0, the code generates the functions corresponding to
 *  direct collisional dissociation:
 *  
 *          H2 + H --> 3H
 *  
 *  If idt = 1, the code generates functions corresponding to the destruction of H2 by dissociative
 *  tunneling:
 * 
 *          H2 + H --> (H2)* + H --> 3H      ; where (H2)* represents the quasi-bound state.
 * 
 *  The data for the computation of the functions is from Martin, Schwartz & Mandy (1996, ApJ 461, 265) 
 *  in both cases.
 *  
 * 
 * Inputs:
 *  Tbin_index : Index of the temperature bin within the calculation 
 *                       for which the rates are being computed
 *  T : Temperature of the gas (in Kelvin)
 *  idt : Flag which decides whether to calculate the direct collisional dissociation rate or the 
 *        dissociatve tunneling rate
 *  *my_chemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *my_rates : Pointer to the structure of type chemistry_data_storage
 *               which holds pointers to the rates.s
 * 
 * Outputs:
 *  f1 - f7 : Rates as described above
 * 
 * Units:
 * 
 * This code is a refactoring of the original Grackle fortran function of the same name.
 */
/*********************************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grackle_chemistry_data.h"
#include "grackle_macros.h"
#include "grackle_types.h"

//*Define the function which will calculate the rates.
int colh2diss_g_c(int Tbin_index, double T, int idt, chemistry_data *my_chemistry,  chemistry_data_storage *my_rates) {

    //*Define variables for each of the rates and give them a preliminary value.
    double f1 = tiny;
    double f2 = tiny;
    double f3 = tiny;
    double f4 = tiny;
    double f5 = 1.0;
    double f6 = 1.0;
    double f7 = 0.0;

    //*Define an array which will hold 21 fitting parameters necessary for calculating the rates.
    double fitParam[21];

    //*The gas temperature cannot be smaller than 500 Kelvin or larger than 1*e6 Kelvin.
    //Note that data and fits are only accurate for temperatures below 5*e5 Kelvin, but collisional dissociation by electrons 
    //dominates above this limit anyway so it is of minor significance.
    if (T <= 500.0) {
        T = 500.0;
    }
    if (T >= 1.0e6) {
        T = 1.0e6;
    }

    //*Set the values of the fitting parameters depending on the value of the idt flag.
    if (idt == 0) {
        fitParam[0]   =   -1.784239e2;
        fitParam[1]   =   -6.842243e1;
        fitParam[2]   =    4.320243e1;
        fitParam[3]   =   -4.633167e0;
        fitParam[4]   =    6.970086e1;
        fitParam[5]   =    4.087038e4;
        fitParam[6]   =   -2.370570e4;
        fitParam[7]   =    1.288953e2;
        fitParam[8]   =   -5.391334e1;
        fitParam[9]   =    5.315517e0;
        fitParam[10]  =   -1.973427e1;
        fitParam[11]  =    1.678095e4;
        fitParam[12]  =   -2.578611e4;
        fitParam[13]  =    1.482123e1;
        fitParam[14]  =   -4.890915e0;
        fitParam[15]  =    4.749030e-1;
        fitParam[16]  =   -1.338283e2;
        fitParam[17]  =   -1.164408e0;
        fitParam[18]  =    8.227443e-1;
        fitParam[19]  =    5.864073e-1;
        fitParam[20]  =   -2.056313e0;
    } else if (idt == 1) {
        fitParam[0]   =   -1.427664e+02;
        fitParam[1]   =    4.270741e+01;
        fitParam[2]   =   -2.027365e+00;
        fitParam[3]   =   -2.582097e-01;
        fitParam[4]   =    2.136094e+01;
        fitParam[5]   =    2.753531e+04;
        fitParam[6]   =   -2.146779e+04;
        fitParam[7]   =    6.034928e+01;
        fitParam[8]   =   -2.743096e+01;
        fitParam[9]   =    2.676150e+00;
        fitParam[10]  =   -1.128215e+01;
        fitParam[11]  =    1.425455e+04;
        fitParam[12]  =   -2.312520e+04;
        fitParam[13]  =    9.305564e+00;
        fitParam[14]  =   -2.464009e+00;
        fitParam[15]  =    1.985955e-01;
        fitParam[16]  =    7.430600e+02;
        fitParam[17]  =   -1.174242e+00;
        fitParam[18]  =    7.502286e-01;
        fitParam[19]  =    2.358848e-01;
        fitParam[20]  =    2.937507e+00;
    } else {
        //Print error message if value of idt is invalid and return failure.
        fprintf(stderr, 'idt has been set to an unknown value. Expected 0 or 1, received %d.\n', idt);
        return(FAIL);
    }

    //Define log10 of the temperature for convenience in the following calculations.
    double logT = log10(T);

    //*Calculate parameters needed to obtain the rates by using the fitting parameters.
    //High density limit.
    double a = fitParam[0] + fitParam[1]*logT + fitParam[2]*pow(logT, 2) + fitParam[3]*pow(logT, 3)
               + fitParam[4]*log10(1.0 + fitParam[5]/T);
    double a1 = fitParam[6]/T;
    //Low density limit.
    double b = fitParam[7] + fitParam[8]*logT + fitParam[9]*pow(logT, 2)
               + fitParam[10]*log10(1.0 + fitParam[11]/T); 
    double b1 = fitParam[12]/T;
    //Critical density.
    double c = fitParam[13] + fitParam[14]*logT + fitParam[15]*pow(logT, 2) + fitParam[16]/T;
    double c1 = fitParam[17] + c;
    double d = fitParam[18] + fitParam[19]*exp(-T/1850.0) + fitParam[20]*exp(-T/440.0);

    //*Calculate the rates from the parameters above.
    f1 = a;
    f2 = a - b;
    f3 = a1;
    f4 = a1 - b1;
    f5 = pow(10.0, c);
    f6 = pow(10.0, c1);
    f7 = d;
    
    //* Store the rates within k13dd in the position prescribed by the Tbin_index and rateIndex.
    //Get the number of temperature bins for which the computations are being computed.
    int noTempBins = my_chemistry->NumberOfTemperatureBins;
    //Store the rates appropriately.
    my_rates->k13dd[Tbin_index + noTempBins*(idt*7)] = f1;
    my_rates->k13dd[Tbin_index + noTempBins*(1 + idt*7)] = f2;
    my_rates->k13dd[Tbin_index + noTempBins*(2 + idt*7)] = f3;
    my_rates->k13dd[Tbin_index + noTempBins*(3 + idt*7)] = f4;
    my_rates->k13dd[Tbin_index + noTempBins*(4 + idt*7)] = f5;
    my_rates->k13dd[Tbin_index + noTempBins*(5 + idt*7)] = f6;
    my_rates->k13dd[Tbin_index + noTempBins*(6 + idt*7)] = f7;

    return(SUCCESS);
}