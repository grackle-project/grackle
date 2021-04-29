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
 *  Tbin_index : Index of the temperature bin within the calculation 
 *                       for which the rates are being computed
 *  T : Temperature of the gas (in Kelvin)
 *  kUnit : Normalisation factor for the units
 *  *my_chemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *my_rates : Pointer to the structure of type chemistry_data_storage
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
int coll_rates_g_c(int Tbin_index, double T, double kUnit, chemistry_data *my_chemistry,  chemistry_data_storage *my_rates) {

    //*The formulae for the collisional rates rely on the temperature in various forms, these are calculated now for convenience.
    double logT = log(T);
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    

    //*Calculate k2, which depends on temperature and has a case B form.
    if (my_chemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            my_rates->k2[Tbin_index] = 4.881357e-6*pow(T, -1.5)
                * pow((1.0 + 1.14813e2*pow(T, -0.407)), -2.242) / kUnit;
        } else {
            my_rates->k2[Tbin_index] = tiny;
        }  
    } else {
        if (T > 5500) {
            my_rates->k2[Tbin_index] = exp( -28.61303380689232
                - 0.7241125657826851*logT_ev
                - 0.02026044731984691*pow(logT_ev, 2)
                - 0.002380861877349834*pow(logT_ev, 3)
                - 0.0003212605213188796*pow(logT_ev, 4)
                - 0.00001421502914054107*pow(logT_ev, 5)
                + 4.989108920299513e-6*pow(logT_ev, 6)
                + 5.755614137575758e-7*pow(logT_ev, 7)
                - 1.856767039775261e-8*pow(logT_ev, 8)
                - 3.071135243196595e-9*pow(logT_ev, 9)) / kUnit;
        } else {
            my_rates->k2[Tbin_index] = my_rates->k4[Tbin_index];
        }
    }    


    return(SUCCESS);
    //*------------------------------ END OF FUNCTION ------------------------------
}