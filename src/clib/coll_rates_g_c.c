/*******************************************************************************/
/*
 * ----------Function to Calculate Multispecies Collisional Rates---------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 12/2/2021
 * Date of Last Modification: 16/2/2021
 * 
 * Purpose:
 *  Calculates multispecies collisional rates.
 * 
 * Inputs:
 *  temperatureBinIndex : Index of the temperature bin within the calculation 
 *                       for which the rates are being computed
 *  T : Temperature of the gas (in Kelvin)
 *  kUnit : Normalisation factor for the units
 *  *userChemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *userRates : Pointer to the structure of type chemistry_data_storage
 *               which holds pointers to the rates.
 * 
 * Outputs:
 *  k1 - k19 : Rates for the various chemical reactions shown below.
 * 
 * Units:
 *  cgs / kUnit
 * 
 * Rate Coefficients:
 *   --All rates are labelled as in Abel et al., 1996 (astro-ph/9608040)--
 *   @  k1  @     H + e --> H(+) + 2e
 *   @  k2  @     H(+) + e --> H + photon
 *   @  k3  @     He + e --> He(+) + 2e
 *   @  k4  @     He(+) + e --> He + photon
 *   @  k5  @     He(+) + e --> He(++) + 2e
 *   @  k6  @     He(++) + e --> He(+) + photon
 *   @  k7  @     H + e --> H(-) + photon
 *   @  k8  @     H + H(-) --> H2 + e
 *   @  k9  @     H + H(+) --> H2(+) + photon
 *   @  k10 @     H2(+) + H --> H2 + H(+)
 *   @  k11 @     H2 + H(+) --> H2(+) + H
 *   @  k12 @     H2 + e --> 2H + e
 *   @  k13 @     H2 + H --> 3H
 *   @  k14 @     H(-) + e --> H + 2e
 *   @  k15 @     H(-) + H --> 2H + e
 *   @  k16 @     H(-) + H(+) --> 2H
 *   @  k17 @     H(-) + H(+) --> H2(+) + e
 *   @  k18 @     H2(+) + e --> 2H
 *   @  k19 @     H2(+) + H(-) --> H2 + H
 * 
 * This code is a refactoring of the original Grackle fortran function of the same name.
 */
/*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grackle_chemistry_data.h"
#include "grackle_macros.h"
#include "grackle_types.h"

// Define the function which will calculate the collisional rates.
void coll_rates_g_c(int temperatureBinIndex, double T, double kUnit, chemistry_data *userChemistry,  chemistry_data_storage *userRates) {

    //*The formulae for the collisional rates rely on the temperature in various forms, these are calculated now for convenience.
    double logT = log(T);
    double T_ev = T/11605;
    double logT_ev = log(T_ev);


    //*Calculate k1, k3, k4 and k5. These depend on the temperature.
    userRates->k1[temperatureBinIndex] = exp( -32.71396786 \
                        + 13.536556*logT_ev \
                        - 5.73932875*pow(logT_ev, 2) \
                        + 1.56315498*pow(logT_ev, 3) \
                        - 0.2877056*pow(logT_ev, 4) \
                        + 0.0348255977*pow(logT_ev, 5) \
                        - 0.00263197617*pow(logT_ev, 6) \
                        + 0.000111954395*pow(logT_ev, 7) \
                        - 0.00000203914985*pow(logT_ev, 8)) / kUnit;

    if (T_ev > 0.8) {
        //k3.
        userRates->k3[temperatureBinIndex] = exp( -44.09864886561001 \
             + 23.91596563469*logT_ev \
             - 10.75323019821*pow(logT_ev, 2) \
             + 3.058038757198*pow(logT_ev, 3) \
             - 0.5685118909884001*pow(logT_ev, 4) \
             + 0.06795391233790001*pow(logT_ev, 5) \
             - 0.005009056101857001*pow(logT_ev, 6) \
             + 0.0002067236157507*pow(logT_ev, 7) \
             - 3.649161410833e-6*pow(logT_ev, 8)) / kUnit;

        //k4.
        userRates->k4[temperatureBinIndex] = (1.54e-9*(1.0 + 0.3 / exp(8.099328789667/T_ev)) \
             / (exp(40.49664394833662/T_ev)*pow(T_ev, 1.5)) \
             + 3.92e-13/pow(T_ev, 0.6353)) / kUnit;

        //k5.
        userRates->k5[temperatureBinIndex] = exp(-68.71040990212001 \
             + 43.93347632635*logT_ev \
             - 18.48066993568*pow(logT_ev, 2) \
             + 4.701626486759002*pow(logT_ev, 3) \
             - 0.7692466334492*pow(logT_ev, 4) \
             + 0.08113042097303*pow(logT_ev, 5) \
             - 0.005324020628287001*pow(logT_ev, 6) \
             + 0.0001975705312221*pow(logT_ev, 7) \
             - 3.165581065665e-6*pow(logT_ev, 8)) / kUnit;
    } else {
        userRates->k1[temperatureBinIndex] = fmax(tiny, userRates->k1[temperatureBinIndex]); 
        userRates->k3[temperatureBinIndex] = tiny;
        userRates->k4[temperatureBinIndex] = 3.92e-13/pow(T_ev, 0.6353) / kUnit;
        userRates->k5[temperatureBinIndex] = tiny;
    }

    //If the case B recombination rates have been requested, k4 must be recalculated.
    if (userChemistry->CaseBRecombination == 1) {
        userRates->k4[temperatureBinIndex] = 1.26e-14 * pow(5.7067e5/T, 0.75) / kUnit;
    }

    //*Calculate k2, which depends on temperature and has a case B form.
    if (userChemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            userRates->k2[temperatureBinIndex] = 4.881357e-6*pow(T, -1.5) \
                * pow((1.0 + 1.14813e2*pow(T, -0.407)), -2.242) / kUnit;
        } else {
            userRates->k2[temperatureBinIndex] = tiny;
        }  
    } else {
        if (T > 5500) {
            userRates->k2[temperatureBinIndex] = exp( -28.61303380689232 \
                - 0.7241125657826851*logT_ev \
                - 0.02026044731984691*pow(logT_ev, 2) \
                - 0.002380861877349834*pow(logT_ev, 3) \
                - 0.0003212605213188796*pow(logT_ev, 4) \
                - 0.00001421502914054107*pow(logT_ev, 5) \
                + 4.989108920299513e-6*pow(logT_ev, 6) \
                + 5.755614137575758e-7*pow(logT_ev, 7) \
                - 1.856767039775261e-8*pow(logT_ev, 8) \
                - 3.071135243196595e-9*pow(logT_ev, 9)) / kUnit;
        } else {
            userRates->k2[temperatureBinIndex] = userRates->k4[temperatureBinIndex];
        }
    }

    //*Calculate k6, which depends on temperature and has a case B form.
    if (userChemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            userRates->k6[temperatureBinIndex] = 7.8155e-5*pow(T, -1.5) \
                * pow((1.0 + 2.0189e2*pow(T, -0.407)), -2.242) / kUnit;
        } else {
            userRates->k6[temperatureBinIndex] = tiny;
        }
    } else {
        userRates->k6[temperatureBinIndex] = 3.36e-10/sqrt(T)/pow(T/1.0e3, 0.2) \
             / (1.0 + pow(T/1.0e6, 0.7)) / kUnit;
    }

    //*Calculate k7. 
    //Fit --> Stancil, Lepp & Dalgarno (1998, ApJ, 509, 1). Based on photodetachment cross-section --> Wishart (1979, MNRAS, 187, P59).
    userRates->k7[temperatureBinIndex] = 3.0e-16*pow(T/3.0e2, 0.95) * exp(-T/9.32e3) / kUnit;

    //*Calculate k8.
    //Fit based on experimental measurements --> Kreckel et al (2010, Science, 329, 69).
    userRates->k8[temperatureBinIndex] = 1.35e-9 * ( pow(T, 9.8493e-2) + 3.2852e-1*pow(T, 5.5610e-1) + 2.771e-7*pow(T, 2.1826) ) \
                / ( 1.0 + 6.191e-3*pow(T, 1.0461) + 8.9712e-11*pow(T, 3.0424) + 3.2576e-14*pow(T, 3.7741) ) \
                / kUnit;

    //*Calculate k9. This is dependent upon the temperature.
    //Fit --> Latif et al (2015, MNRAS, 446, 3163): valid for 1 < T < 32000 K.
    if (T < 30.0) {
        userRates->k9[temperatureBinIndex] = 2.10e-20*pow(T/30.0, -0.15) / kUnit;
    } else {
        //If temperature is exceeding 32000 K, it is instead fixed at 32000 K -- the behaviour at this temperature should not be important.
        double T_k9;
        if (T > 3.2e4) {
            T_k9 = 3.2e4;
        } else {
            T_k9 = T;
        }

        userRates->k9[temperatureBinIndex] = pow(10.0, -18.20  - 3.194*log10(T_k9) \
               + 1.786*pow(log10(T_k9), 2) - 0.2072*pow(log10(T_k9), 3)) \
               / kUnit;
    }

    //*Calculate k10, this has no temperature dependance. 
    userRates->k10[temperatureBinIndex] = 6.0e-10 / kUnit;

    //*Calculate k11, k12 and k13. These are dependent upon temperature and there are two schemes for the calculation of k11.
    if ( T_ev > 0.3) {
        //! NEED TO UNCOMMENT
        /*
        //k11 is calculated by using either Savin 2004 or Abel et al. 1997. The parameter to control this is within the chemistry_data struct.
        if (userChemistry->useSavin2004 == 1) {
            userRates->k11[temperatureBinIndex] = ( exp(-21237.15/T) * \
                (- 3.3232183e-07 \
                + 3.3735382e-07 * logT \
                - 1.4491368e-07 * pow(logT, 2) \
                + 3.4172805e-08 * pow(logT, 3) \
                - 4.7813720e-09 * pow(logT, 4) \
                + 3.9731542e-10 * pow(logT, 5) \
                - 1.8171411e-11 * pow(logT, 6) \
                + 3.5311932e-13 * pow(logT, 7))) / kUnit;
        } else if (userChemistry->useSavin2004 == 0) {
            userRates->k11[temperatureBinIndex] = exp( -24.24914687731536 \
                + 3.400824447095291*logT_ev \
                - 3.898003964650152*pow(logT_ev, 2) \
                + 2.045587822403071*pow(logT_ev, 3) \
                - 0.5416182856220388*pow(logT_ev, 4) \
                + 0.0841077503763412*pow(logT_ev, 5) \
                - 0.007879026154483455*pow(logT_ev, 6) \
                + 0.0004138398421504563*pow(logT_ev, 7) \
                - 9.36345888928611e-6*pow(logT_ev, 8)) / kUnit;
        } else {
            printf("Invalid value for flag 'useSavin2004' encountered. This flag must be set to either 1 or 0.");
            exit(0);
        }

        //k12 --> Trevisan & Tennyson (2002, Plasma Phys. Cont. Fus., 44, 1263).
        userRates->k12[temperatureBinIndex] = 4.4886e-9*pow(T, 0.109127) * exp(-101858.0/T) / kUnit;

        //k13.
        userRates->k13[temperatureBinIndex] = 1.0670825e-10*pow(T_ev, 2.012) \
             / ( exp(4.463/T_ev) * pow((1.0 + 0.2472*T_ev), 3.512) ) \
             / kUnit;       
    } else {
        userRates->k11[temperatureBinIndex] = tiny;
        userRates->k12[temperatureBinIndex] = tiny;
        userRates->k13[temperatureBinIndex] = tiny;
    */
    }

    //*Calculate k14 (dependant on temperature).
    if (T_ev > 0.04) {
        userRates->k14[temperatureBinIndex] = exp( -18.01849334273 \
             + 2.360852208681*logT_ev \
             - 0.2827443061704*pow(logT_ev, 2) \
             + 0.01623316639567*pow(logT_ev, 3) \
             - 0.03365012031362999*pow(logT_ev, 4) \
             + 0.01178329782711*pow(logT_ev, 5) \
             - 0.001656194699504*pow(logT_ev, 6) \
             + 0.0001068275202678*pow(logT_ev, 7) \
             - 2.631285809207e-6*pow(logT_ev, 8)) / kUnit;
    } else {
        userRates->k14[temperatureBinIndex] = tiny;
    }

    //*Calculate k15 (dependant on temperature).
    if (T_ev > 0.1) {
        userRates->k15[temperatureBinIndex] = exp( -20.37260896533324 \
             + 1.139449335841631*logT_ev \
             - 0.1421013521554148*pow(logT_ev, 2) \
             + 0.00846445538663*pow(logT_ev, 3) \
             - 0.0014327641212992*pow(logT_ev, 4) \
             + 0.0002012250284791*pow(logT_ev, 5) \
             + 0.0000866396324309*pow(logT_ev, 6) \
             - 0.00002585009680264*pow(logT_ev, 7) \
             + 2.4555011970392e-6*pow(logT_ev, 8) \
             - 8.06838246118e-8*pow(logT_ev, 9) ) / kUnit;
    } else {
        userRates->k15[temperatureBinIndex] = 2.56e-9*pow(T_ev, 1.78186) / kUnit;
    }

    //*Calculate k16.
    //Fit --> Croft et al (1999, MNRAS, 304, 327). Based on cross-section --> Fussen & Kubach (1986, J. Phys. B, 18 L31).
    userRates->k16[temperatureBinIndex] = 2.4e-6*(1.0 + T/2.0e4) / sqrt(T) / kUnit;

    //*Calculate k17 (dependant on temperature).
    if (T > 1.0e4) {
        userRates->k17[temperatureBinIndex] = 4.0e-4*pow(T, -1.4) * exp(-15100.0/T) / kUnit;
    } else {
        userRates->k17[temperatureBinIndex] = 1.0e-8*pow(T, -0.4) / kUnit;
    }

    //*Calculate k18 (dependant on temperature).
    if (T > 617.0) {
        userRates->k18[temperatureBinIndex] = 1.32e-6*pow(T, -0.76) / kUnit;
    } else {
        userRates->k18[temperatureBinIndex] = 1.0e-8 / kUnit;
    }

    //*Calculate k19.
    userRates->k19[temperatureBinIndex] = 5.e-7 * sqrt(100.0/T) / kUnit;

    //*Calculate k23.
    userRates->k23[temperatureBinIndex] = ( (8.125e-8/sqrt(T)) * exp(-52000.0/T) * (1.0 - exp(-6000.0/T)) ) / kUnit;
    userRates->k23[temperatureBinIndex] = max(tiny, userRates->k23[temperatureBinIndex]);

    //*------------------------------ END OF FUNCTION ------------------------------
}