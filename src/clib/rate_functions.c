/**********************************************************************
 * ---- File containing definitions of all rates ----
 *
 * Written by: Ewan Jones
 * Date of Creation: 26/4/2021
 * 
 * Purpose:
 *  Contains definitions for all rates. These will be called within 
 *  calc_rates_g_c.
**********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

// Calculation of k1.
double k1_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k1 = exp( -32.71396786375
                        + 13.53655609057*logT_ev
                        - 5.739328757388*pow(logT_ev, 2)
                        + 1.563154982022*pow(logT_ev, 3)
                        - 0.2877056004391*pow(logT_ev, 4)
                        + 0.03482559773736999*pow(logT_ev, 5)
                        - 0.00263197617559*pow(logT_ev, 6)
                        + 0.0001119543953861*pow(logT_ev, 7)
                        - 2.039149852002e-6*pow(logT_ev, 8)) / kUnit;
    if (T_ev <= 0.8){
        k1 = fmax(tiny, k1); 
    }
    return k1;
}

//Calculation of k2.
double k2_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    return 1;
}

//Calculation of k3.
double k3_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k3;
    if (T_ev > 0.8){
        k3 = exp( -44.09864886561001
                + 23.91596563469*logT_ev
                - 10.75323019821*pow(logT_ev, 2)
                + 3.058038757198*pow(logT_ev, 3)
                - 0.5685118909884001*pow(logT_ev, 4)
                + 0.06795391233790001*pow(logT_ev, 5)
                - 0.005009056101857001*pow(logT_ev, 6)
                + 0.0002067236157507*pow(logT_ev, 7)
                - 3.649161410833e-6*pow(logT_ev, 8)) / kUnit;
    } else {
        k3 = tiny;
    }
    return k3;
}   

//Calculation of k4.
double k4_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k4;
    //If case B recombination on.
    if (my_chemistry->CaseBRecombination == 1){
        k4 = 1.26e-14 * pow(5.7067e5/T, 0.75) / kUnit;
        return k4;
    }

    //If case B recombination off.
    if (T_ev > 0.8){
        k4 = (1.54e-9*(1.0 + 0.3 / exp(8.099328789667/T_ev))
             / (exp(40.49664394833662/T_ev)*pow(T_ev, 1.5))
             + 3.92e-13/pow(T_ev, 0.6353)) / kUnit;
    } else {
        k4 = 3.92e-13/pow(T_ev, 0.6353) / kUnit;
    }
    return k4;
}

//Calculation of k5.
double k5_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k5;
    if (T_ev > 0.8){
        k5 = exp(-68.71040990212001
                + 43.93347632635*logT_ev
                - 18.48066993568*pow(logT_ev, 2)
                + 4.701626486759002*pow(logT_ev, 3)
                - 0.7692466334492*pow(logT_ev, 4)
                + 0.08113042097303*pow(logT_ev, 5)
                - 0.005324020628287001*pow(logT_ev, 6)
                + 0.0001975705312221*pow(logT_ev, 7)
                - 3.165581065665e-6*pow(logT_ev, 8)) / kUnit;
    } else {
        k5 = tiny;
    }
    return k5;
}

//Calculation of k6.
double k6_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double k6;
    //Has case B recombination setting.
    if (my_chemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            k6 = 7.8155e-5*pow(T, -1.5)
                * pow((1.0 + 2.0189e2*pow(T, -0.407)), -2.242) / kUnit;
        } else {
            k6 = tiny;
        }
    } else {
        k6 = 3.36e-10/sqrt(T)/pow(T/1.0e3, 0.2)
             / (1.0 + pow(T/1.0e6, 0.7)) / kUnit;
    }
    return k6;
}

//Calculation of k7.
double k7_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    //Fit --> Stancil, Lepp & Dalgarno (1998, ApJ, 509, 1). Based on photodetachment cross-section --> Wishart (1979, MNRAS, 187, P59).
    return 3.0e-16*pow(T/3.0e2, 0.95) * exp(-T/9.32e3) / kUnit;
}

//Calculation of k8.
double k8_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    //Fit based on experimental measurements --> Kreckel et al (2010, Science, 329, 69).
    return 1.35e-9 * ( pow(T, 9.8493e-2) + 3.2852e-1*pow(T, 5.5610e-1) + 2.771e-7*pow(T, 2.1826) )
                / ( 1.0 + 6.191e-3*pow(T, 1.0461) + 8.9712e-11*pow(T, 3.0424) + 3.2576e-14*pow(T, 3.7741) )
                / kUnit;
}

//Calculation of k9.
double k9_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double k9;
    //Fit --> Latif et al (2015, MNRAS, 446, 3163): valid for 1 < T < 32000 K.
    if (T < 30.0) {
        k9 = 2.10e-20*pow(T/30.0, -0.15) / kUnit;
    } else {
        //If temperature is exceeding 32000 K, it is instead fixed at 32000 K -- the behaviour at this temperature should not be important.
        double T_k9 = fmin(T, 3.2e4);

        k9 = pow(10.0, -18.20  - 3.194*log10(T_k9)
                + 1.786*pow(log10(T_k9), 2) - 0.2072*pow(log10(T_k9), 3))
                / kUnit;
    }
    return k9;
}

//Calculation of k10.
double k10_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    return 6.0e-10 / kUnit;
}

//Calculation of k11.
double k11_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double logT = log(T);
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k11;
    if ( T_ev > 0.3) {
        //k11 is calculated by using either Savin 2004 or Abel et al. 1997. The parameter to control this is within the chemistry_data struct.
        if (my_chemistry->k11_rate == 1) {
            k11 = ( exp(-21237.15/T) *
                (- 3.3232183e-07
                + 3.3735382e-07 * logT
                - 1.4491368e-07 * pow(logT, 2)
                + 3.4172805e-08 * pow(logT, 3)
                - 4.7813720e-09 * pow(logT, 4)
                + 3.9731542e-10 * pow(logT, 5)
                - 1.8171411e-11 * pow(logT, 6)
                + 3.5311932e-13 * pow(logT, 7))) / kUnit;
        } else if (my_chemistry->k11_rate == 2) {
            k11 = exp( -24.24914687731536
                + 3.400824447095291*logT_ev
                - 3.898003964650152*pow(logT_ev, 2)
                + 2.045587822403071*pow(logT_ev, 3)
                - 0.5416182856220388*pow(logT_ev, 4)
                + 0.0841077503763412*pow(logT_ev, 5)
                - 0.007879026154483455*pow(logT_ev, 6)
                + 0.0004138398421504563*pow(logT_ev, 7)
                - 9.36345888928611e-6*pow(logT_ev, 8)) / kUnit;
        } else {
            fprintf(stderr, 'k11_rate flag set to unknown value. This must be either 1 \
                             or 2 but was set to %d.\n', my_chemistry->k11_rate);
            return(FAIL); //! NOT SURE HOW TO DEAL WITH THIS QUITE YET
        }
    } else {
        k11 = tiny;
    }
    return k11;
}

//Calculation of k12.
double k12_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    
    double k12;
    if ( T_ev > 0.3) {
        //k12 --> Trevisan & Tennyson (2002, Plasma Phys. Cont. Fus., 44, 1263).
        k12 = 4.4886e-9*pow(T, 0.109127) * exp(-101858.0/T) / kUnit;
    } else {
        k12 = tiny;
    }
    return k12;
}

//Calculation of k13.
double k13_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;

    double k13;
    if ( T_ev > 0.3) {
        k13 = 1.0670825e-10*pow(T_ev, 2.012)
             / ( exp(4.463/T_ev) * pow((1.0 + 0.2472*T_ev), 3.512) )
             / kUnit;       
    } else {
        k13 = tiny;
    }
    return k13;
}

//Calculation of k14.
double k14_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);

    double k14;
    if (T_ev > 0.04) {
        k14 = exp( -18.01849334273
             + 2.360852208681*logT_ev
             - 0.2827443061704*pow(logT_ev, 2)
             + 0.01623316639567*pow(logT_ev, 3)
             - 0.03365012031362999*pow(logT_ev, 4)
             + 0.01178329782711*pow(logT_ev, 5)
             - 0.001656194699504*pow(logT_ev, 6)
             + 0.0001068275202678*pow(logT_ev, 7)
             - 2.631285809207e-6*pow(logT_ev, 8)) / kUnit;
    } else {
        k14 = tiny;
    }
    return k14;
}

//Calculation of k15.
double k15_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double T_ev = T/11605.0;
    double logT_ev = log(T_ev);
    
    double k15;
    if (T_ev > 0.1) {
        k15 = exp( -20.37260896533324
             + 1.139449335841631*logT_ev
             - 0.1421013521554148*pow(logT_ev, 2)
             + 0.00846445538663*pow(logT_ev, 3)
             - 0.0014327641212992*pow(logT_ev, 4)
             + 0.0002012250284791*pow(logT_ev, 5)
             + 0.0000866396324309*pow(logT_ev, 6)
             - 0.00002585009680264*pow(logT_ev, 7)
             + 2.4555011970392e-6*pow(logT_ev, 8)
             - 8.06838246118e-8*pow(logT_ev, 9) ) / kUnit;
    } else {
        k15 = 2.56e-9*pow(T_ev, 1.78186) / kUnit;
    }
    return k15;
}

//Calculation of k16.
double k16_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    //Fit --> Croft et al (1999, MNRAS, 304, 327). Based on cross-section --> Fussen & Kubach (1986, J. Phys. B, 18 L31).
    return 2.4e-6*(1.0 + T/2.0e4) / sqrt(T) / kUnit;
}

//Calculation of k17.
double k17_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
     double k17;
    if (T > 1.0e4) {
        k17 = 4.0e-4*pow(T, -1.4) * exp(-15100.0/T) / kUnit;
    } else {
        k17 = 1.0e-8*pow(T, -0.4) / kUnit;
    }
    return k17;
}

//Calculation of k18.
double k18_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double k18;
    if (T > 617.0) {
        k18 = 1.32e-6*pow(T, -0.76) / kUnit;
    } else {
        k18 = 1.0e-8 / kUnit;
    }
    return k18;
}

//Calculation of k19.
double k19_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    return 5.e-7 * sqrt(100.0/T) / kUnit;
}

//Calculation of k23.
double k23_rate(double T, double kUnit, chemistry_data *my_chemistry)
{
    double k23;
    k23 = ( (8.125e-8/sqrt(T)) * exp(-52000.0/T) * (1.0 - exp(-6000.0/T)) ) / kUnit;
    k23 = max(tiny, k23);
    return k23;
}