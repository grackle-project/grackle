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
 *  temperatureBinIndex : Index of the temperature bin within the calculation 
 *                       for which the rates are being computed
 *  T : Temperature of the gas (in Kelvin)
 *  idt : Flag which decides whether to calculate the direct collisional dissociation rate or the 
 *        dissociatve tunneling rate
 *  *userChemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *userRates : Pointer to the structure of type chemistry_data_storage
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
void colh2diss_g_c(int temperatureBinIndex, double T, int idt, chemistry_data *userChemistry,  chemistry_data_storage *userRates) {

    //*Define variables for each of the rates and give them a preliminary value.
    double f1 = tiny;
    double f2 = tiny;
    double f3 = tiny;
    double f4 = tiny;
    double f5 = 1.0;
    double f6 = 1.0;
    double f7 = 0.0;

    //*Define an array which will hold 21 fitting parameters necessary for calculating the rates.
    double fittingParams[21];

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
        fittingParams[0]   =   -1.784239e2;
        fittingParams[1]   =   -6.842243e1;
        fittingParams[2]   =    4.320243e1;
        fittingParams[3]   =   -4.633167e0;
        fittingParams[4]   =    6.970086e1;
        fittingParams[5]   =    4.087038e4;
        fittingParams[6]   =   -2.370570e4;
        fittingParams[7]   =    1.288953e2;
        fittingParams[8]   =   -5.391334e1;
        fittingParams[9]   =    5.315517e0;
        fittingParams[10]  =   -1.973427e1;
        fittingParams[11]  =    1.678095e4;
        fittingParams[12]  =   -2.578611e4;
        fittingParams[13]  =    1.482123e1;
        fittingParams[14]  =   -4.890915e0;
        fittingParams[15]  =    4.749030e-1;
        fittingParams[16]  =   -1.338283e2;
        fittingParams[17]  =   -1.164408e0;
        fittingParams[18]  =    8.227443e-1;
        fittingParams[19]  =    5.864073e-1;
        fittingParams[20]  =   -2.056313e0;
    } else if (idt == 1) {
        fittingParams[1]   =    4.270741e+01;
        fittingParams[2]   =   -2.027365e+00;
        fittingParams[3]   =   -2.582097e-01;
        fittingParams[4]   =    2.136094e+01;
        fittingParams[0]   =   -1.427664e+02;
        fittingParams[5]   =    2.753531e+04;
        fittingParams[6]   =   -2.146779e+04;
        fittingParams[7]   =    6.034928e+01;
        fittingParams[8]   =   -2.743096e+01;
        fittingParams[9]   =    2.676150e+00;
        fittingParams[10]  =   -1.128215e+01;
        fittingParams[11]  =    1.425455e+04;
        fittingParams[12]  =   -2.312520e+04;
        fittingParams[13]  =    9.305564e+00;
        fittingParams[14]  =   -2.464009e+00;
        fittingParams[15]  =    1.985955e-01;
        fittingParams[16]  =    7.430600e+02;
        fittingParams[17]  =   -1.174242e+00;
        fittingParams[18]  =    7.502286e-01;
        fittingParams[19]  =    2.358848e-01;
        fittingParams[20]  =    2.937507e+00;
    } else {
        //Print error message to terminal if value of idt is invalid and exit the program.
        printf("Incorrect value for idt encountered. \n");
        printf("Expected 0 or 1 -- received %d", idt);
        exit(0);
    }

    //Define log10 of the temperature for convenience in the following calculations.
    double logT = log10(T);

    //*Calculate parameters needed to obtain the rates by using the fitting parameters.
    //High density limit.
    double a = fittingParams[0] + fittingParams[1]*logT + fittingParams[2]*pow(logT, 2) \
               + fittingParams[3]*pow(logT, 3) + fittingParams[4]*log10(1.0 + fittingParams[5]/T);
    double a1 = fittingParams[6]/T;
    //Low density limit.
    double b = fittingParams[7] + fittingParams[8]*logT + fittingParams[9]*pow(logT, 2) \
               + fittingParams[10]*log10(1.0 + fittingParams[11]/T); 
    double b1 = fittingParams[12]/T;
    //Critical density.
    double c = fittingParams[13] + fittingParams[14]*logT + fittingParams[15]*pow(logT, 2) \
               + fittingParams[16]/T;
    double c1 = fittingParams[17] + c;
    double d = fittingParams[18] + fittingParams[19]*exp(-T/1850.0) + fittingParams[20]*exp(-T/440.0);

    //*Calculate the rates from the parameters above.
    f1 = a;
    f2 = a - b;
    f3 = a1;
    f4 = a1 - b1;
    f5 = pow(10.0, c);
    f6 = pow(10.0, c1);
    f7 = d;
    
    //* Store the rates within k13dd in the position prescribed by the temperatureBinIndex and rateIndex.
    //Get the number of temperature bins for which the computations are being computed.
    int noTempBins = userChemistry->NumberOfTemperatureBins;
    //Store the rates appropriately.
    userRates->k13dd[temperatureBinIndex + noTempBins*(idt*7)] = f1;
    userRates->k13dd[temperatureBinIndex + noTempBins*(1 + idt*7)] = f2;
    userRates->k13dd[temperatureBinIndex + noTempBins*(2 + idt*7)] = f3;
    userRates->k13dd[temperatureBinIndex + noTempBins*(3 + idt*7)] = f4;
    userRates->k13dd[temperatureBinIndex + noTempBins*(4 + idt*7)] = f5;
    userRates->k13dd[temperatureBinIndex + noTempBins*(5 + idt*7)] = f6;
    userRates->k13dd[temperatureBinIndex + noTempBins*(6 + idt*7)] = f7;
}