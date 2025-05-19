/**********************************************************************
 * ---- File containing definitions of all rates ----
 *
 * Written by: Ewan Jones
 * Date of Creation: 26/4/2021
 * 
 * Purpose:
 *  Contains definitions for all rates. These will be called within 
 *  calc_rates_g_c.
 * 
 * Set "units" parameter = 1 for cgs unit outputs.
**********************************************************************/

//* Macro definitions for constants.
//Kelvin to eV conversion factor
#ifndef tevk
#define tevk 1.1605e4
#endif
//Comparison value
#ifndef dhuge
#define dhuge 1.0e30
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "cie_thin_cooling_rate_tables.h"
#include "phys_constants.h"


// Calculation of k1 (HI + e --> HII + 2e)
double k1_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
    double logT_ev = log(T_ev);

    double k1 = exp( -32.71396786375
                        + 13.53655609057*logT_ev
                        - 5.739328757388*pow(logT_ev, 2)
                        + 1.563154982022*pow(logT_ev, 3)
                        - 0.2877056004391*pow(logT_ev, 4)
                        + 0.03482559773736999*pow(logT_ev, 5)
                        - 0.00263197617559*pow(logT_ev, 6)
                        + 0.0001119543953861*pow(logT_ev, 7)
                        - 2.039149852002e-6*pow(logT_ev, 8)) / units;
    if (T_ev <= 0.8){
        k1 = fmax(tiny, k1); 
    }
    return k1;
}

//Calculation of k3 (HeI + e --> HeII + 2e)
double k3_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
    double logT_ev = log(T_ev);

    if (T_ev > 0.8){
        return exp( -44.09864886561001
                + 23.91596563469*logT_ev
                - 10.75323019821*pow(logT_ev, 2)
                + 3.058038757198*pow(logT_ev, 3)
                - 0.5685118909884001*pow(logT_ev, 4)
                + 0.06795391233790001*pow(logT_ev, 5)
                - 0.005009056101857001*pow(logT_ev, 6)
                + 0.0002067236157507*pow(logT_ev, 7)
                - 3.649161410833e-6*pow(logT_ev, 8)) / units;
    } else {
        return tiny;
    }
}   

//Calculation of k4 (HeII + e --> HeI + photon)
double k4_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
    double logT_ev = log(T_ev);

    double k4;
    //If case B recombination on.
    if (my_chemistry->CaseBRecombination == 1){
        return 1.26e-14 * pow(5.7067e5/T, 0.75) / units;
    }

    //If case B recombination off.
    if (T_ev > 0.8){
        return (1.54e-9*(1.0 + 0.3 / exp(8.099328789667/T_ev))
             / (exp(40.49664394833662/T_ev)*pow(T_ev, 1.5))
             + 3.92e-13/pow(T_ev, 0.6353)) / units;
    } else {
        return 3.92e-13/pow(T_ev, 0.6353) / units;
    }
}

//Calculation of k2 (HII + e --> HI + photon)
double k2_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            return 4.881357e-6*pow(T, -1.5) \
                * pow((1.0 + 1.14813e2*pow(T, -0.407)), -2.242) / units;
        } else {
            return tiny;
        }  
    } else {
        if (T > 5500) {
            //Convert temperature to appropriate form.
            double T_ev = T / tevk;
            double logT_ev = log(T_ev);

            return exp( -28.61303380689232 \
                - 0.7241125657826851*logT_ev \
                - 0.02026044731984691*pow(logT_ev, 2) \
                - 0.002380861877349834*pow(logT_ev, 3) \
                - 0.0003212605213188796*pow(logT_ev, 4) \
                - 0.00001421502914054107*pow(logT_ev, 5) \
                + 4.989108920299513e-6*pow(logT_ev, 6) \
                + 5.755614137575758e-7*pow(logT_ev, 7) \
                - 1.856767039775261e-8*pow(logT_ev, 8) \
                - 3.071135243196595e-9*pow(logT_ev, 9)) / units;
        } else {
            return k4_rate(T, units, my_chemistry);
        }
    }
}

//Calculation of k5 (HeII + e --> HeIII + 2e)
double k5_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
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
                - 3.165581065665e-6*pow(logT_ev, 8)) / units;
    } else {
        k5 = tiny;
    }
    return k5;
}

//Calculation of k6 (HeIII + e --> HeII + photon)
double k6_rate(double T, double units, chemistry_data *my_chemistry)
{
    double k6;
    //Has case B recombination setting.
    if (my_chemistry->CaseBRecombination == 1) {
        if (T < 1.0e9) {
            k6 = 7.8155e-5*pow(T, -1.5)
                * pow((1.0 + 2.0189e2*pow(T, -0.407)), -2.242) / units;
        } else {
            k6 = tiny;
        }
    } else {
        k6 = 3.36e-10/sqrt(T)/pow(T/1.0e3, 0.2)
             / (1.0 + pow(T/1.0e6, 0.7)) / units;
    }
    return k6;
}

//Calculation of k7 (HI + e --> HM + photon)
double k7_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Fit --> Stancil, Lepp & Dalgarno (1998, ApJ, 509, 1). Based on photodetachment cross-section --> Wishart (1979, MNRAS, 187, P59).
    return 3.0e-16*pow(T/3.0e2, 0.95) * exp(-T/9.32e3) / units;
}

//Calculation of k8 (HI + HM --> H2I* + e)
double k8_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Fit based on experimental measurements --> Kreckel et al (2010, Science, 329, 69).
    return 1.35e-9 * ( pow(T, 9.8493e-2) + 3.2852e-1*pow(T, 5.5610e-1) + 2.771e-7*pow(T, 2.1826) )
                / ( 1.0 + 6.191e-3*pow(T, 1.0461) + 8.9712e-11*pow(T, 3.0424) + 3.2576e-14*pow(T, 3.7741) )
                / units;
}

//Calculation of k9 (HI + HII --> H2II + photon)
double k9_rate(double T, double units, chemistry_data *my_chemistry)
{
    double k9;
    //Fit --> Latif et al (2015, MNRAS, 446, 3163): valid for 1 < T < 32000 K.
    if (T < 30.0) {
        k9 = 2.10e-20*pow(T/30.0, -0.15) / units;
    } else {
        //If temperature is exceeding 32000 K, it is instead fixed at 32000 K -- the behaviour at this temperature should not be important.
        double T_k9 = fmin(T, 3.2e4);

        k9 = pow(10.0, -18.20  - 3.194*log10(T_k9)
                + 1.786*pow(log10(T_k9), 2) - 0.2072*pow(log10(T_k9), 3))
                / units;
    }
    return k9;
}

//Calculation of k10 (H2II + HI --> H2I* + HII)
double k10_rate(double T, double units, chemistry_data *my_chemistry)
{
    return 6.0e-10 / units;
}

//Calculation of k11 (H2I + HII --> H2II + HI)
double k11_rate(double T, double units, chemistry_data *my_chemistry)
{
    double logT = log(T);
    double T_ev = T / 11605.0;
    double logT_ev = log(T_ev);

    double k11;
    if ( T_ev > 0.3) {
        //k11 is calculated by using either Savin 2004 or Abel et al. 1996. The parameter to control this is within the chemistry_data struct.
        if (my_chemistry->h2_charge_exchange_rate == 1) {
            k11 = ( exp(-21237.15/T) *
                (- 3.3232183e-07
                + 3.3735382e-07 * logT
                - 1.4491368e-07 * pow(logT, 2)
                + 3.4172805e-08 * pow(logT, 3)
                - 4.7813720e-09 * pow(logT, 4)
                + 3.9731542e-10 * pow(logT, 5)
                - 1.8171411e-11 * pow(logT, 6)
                + 3.5311932e-13 * pow(logT, 7))) / units;
        } else if (my_chemistry->h2_charge_exchange_rate == 2) {
            k11 = exp( -24.24914687731536
                + 3.400824447095291*logT_ev
                - 3.898003964650152*pow(logT_ev, 2)
                + 2.045587822403071*pow(logT_ev, 3)
                - 0.5416182856220388*pow(logT_ev, 4)
                + 0.0841077503763412*pow(logT_ev, 5)
                - 0.007879026154483455*pow(logT_ev, 6)
                + 0.0004138398421504563*pow(logT_ev, 7)
                - 9.36345888928611e-6*pow(logT_ev, 8)) / units;
        } else {
            fprintf(stderr, "k11_rate flag set to unknown value. This must be either 1 \
                             or 2 but was set to %d \n", my_chemistry->h2_charge_exchange_rate);
            exit(1);
        }
    } else {
        k11 = tiny;
    }
    return k11;
}

//Calculation of k12 (H2I + e --> 2HI + e)
double k12_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
    
    double k12;
    if ( T_ev > 0.3) {
        //k12 --> Trevisan & Tennyson (2002, Plasma Phys. Cont. Fus., 44, 1263).
        k12 = 4.4886e-9*pow(T, 0.109127) * exp(-101858.0/T) / units;
    } else {
        k12 = tiny;
    }
    return k12;
}

//Calculation of k13 (H2I + HI --> 3HI)
double k13_rate(double T, double units, chemistry_data *my_chemistry)
{   
    double T_ev = T / 11605.0;

    double k13;
    switch (my_chemistry->three_body_rate) {

        case 0:
            if ( T_ev > 0.3) {
                k13 = 1.0670825e-10*pow(T_ev, 2.012)
                    / ( exp(4.463/T_ev) * pow((1.0 + 0.2472*T_ev), 3.512) ); 
            } else {
                k13 = tiny * units;
            }
            break;

        case 1:
            //Inverse of PSS83 three-body rate
            k13 = (5.24e-7 / pow(T, 0.485)) * exp(-5.2e4 / T);
            break;

        case 2:
            //Inverse of CW83 rate
            k13 = 8.4e-11 * pow(T, 0.515) * exp(-5.2e4 / T);
            break;

        case 3:
            //Rate from JGC67, used by FH07 to construct their three-body rate
            k13 = (1.38e-4 / pow(T, 1.025)) * exp(-5.2e4 / T);
            break;

        case 4:
            //High density limiting rate from MSM96
            k13 = pow(1e1 ,(-178.4239 - 68.42243 * log10(T)
                            + 43.20243 * pow(log10(T), 2)
                            - 4.633167 * pow(log10(T), 3) 
                            + 69.70086 * log10(1.0 + 40870.38 / T)
                            - (23705.7 / T)));
            break;

        case 5:
            //From detailed balance from Forrey rate
            if (T <= 3000.0) {
                k13 = 2.4e-8 * exp(-5.2e4/T);
            } else {
                k13 = 2.2e-6 * pow(T, -0.565) * exp(-5.2e4/T); 
            }
            break;
        
        default:
            fprintf(stderr, "three_body_rate has been set to an unknown value: %d \n",
                my_chemistry->three_body_rate);
            exit(1);
    }
    return k13 / units;
}

//Calculates 7 k13dd coefficients for given idt (0 or 1).
//This is an internal function not meant to be accessed by the user -- use k13dd_rate.
void _k13dd_rate(double T, int idt, double units, double *k13dd_results, chemistry_data *my_chemistry)
{   
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
        fitParam[0]   =   -1.427664e+02;
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
        fprintf(stderr, "idt has been set to an unknown value. Expected 0 or 1, received %d \n", idt);
        exit(1);
    }

    //Define log10 of the temperature for convenience in the following calculations.
    double log10_T = log10(T);

    //*Calculate parameters needed to obtain the rates by using the fitting parameters.
    //High density limit.
    double a = fitParam[0] + fitParam[1]*log10_T + fitParam[2]*pow(log10_T, 2) \
               + fitParam[3]*pow(log10_T, 3) + fitParam[4]*log10(1.0 + fitParam[5]/T);
    double a1 = fitParam[6]/T;
    //Low density limit.
    double b = fitParam[7] + fitParam[8]*log10_T + fitParam[9]*pow(log10_T, 2) \
               + fitParam[10]*log10(1.0 + fitParam[11]/T); 
    double b1 = fitParam[12]/T;
    //Critical density.
    double c = fitParam[13] + fitParam[14]*log10_T + fitParam[15]*pow(log10_T, 2) \
               + fitParam[16]/T;
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
    
    //*Store the rates appropriately.
    k13dd_results[idt*7] = f1 - log10(units);
    k13dd_results[1 + idt*7] = f2;
    k13dd_results[2 + idt*7] = f3;
    k13dd_results[3 + idt*7] = f4;
    k13dd_results[4 + idt*7] = f5;
    k13dd_results[5 + idt*7] = f6;
    k13dd_results[6 + idt*7] = f7;

    //* Modify some of the rates if alternative scheme selected
    if (my_chemistry->three_body_rate){
        k13dd_results[idt*7]     = log10(1.12e-10 * exp(-7.035e4 / T) / units);
        k13dd_results[1 + idt*7] = log10(6.5e-7 / sqrt(T) * exp(-5.2e4 / T) *
                                         (1. - exp(-6.3e3 / T)) / units);
        k13dd_results[2 + idt*7] = pow(1.0e1, (4.0 - 0.416 * log10(T / 1.0e4) -
                                       0.327 * pow(log10(T / 1.0e4),2)));
    }
}

//Calculation of k13dd. k13dd_results is a pointer to an array of length 14 * sizeof(double).
void k13dd_rate(double T, double units, double *k13dd_results, chemistry_data *my_chemistry)
{
    for (int idt = 0; idt < 2; idt++){
        _k13dd_rate(T, idt, units, k13dd_results, my_chemistry);
    }
}

//Calculation of k14 (HM + e --> HI + 2e)
double k14_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
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
             - 2.631285809207e-6*pow(logT_ev, 8)) / units;
    } else {
        k14 = tiny;
    }
    return k14;
}

//Calculation of k15 (HM + HI --> 2HI + e)
double k15_rate(double T, double units, chemistry_data *my_chemistry)
{
    double T_ev = T / 11605.0;
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
             - 8.06838246118e-8*pow(logT_ev, 9) ) / units;
    } else {
        k15 = 2.56e-9*pow(T_ev, 1.78186) / units;
    }
    return k15;
}

//Calculation of k16 (HM + HI --> 2HI)
double k16_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Fit --> Croft et al (1999, MNRAS, 304, 327). Based on cross-section --> Fussen & Kubach (1986, J. Phys. B, 18 L31).
    return 2.4e-6*(1.0 + T/2.0e4) / sqrt(T) / units;
}

//Calculation of k17 (HM + HI --> H2I + e)
double k17_rate(double T, double units, chemistry_data *my_chemistry)
{
     double k17;
    if (T > 1.0e4) {
        k17 = 4.0e-4*pow(T, -1.4) * exp(-15100.0/T) / units;
    } else {
        k17 = 1.0e-8*pow(T, -0.4) / units;
    }
    return k17;
}

//Calculation of k18 (H2I + e --> 2HI)
double k18_rate(double T, double units, chemistry_data *my_chemistry)
{
    double k18;
    if (T > 617.0) {
        k18 = 1.32e-6*pow(T, -0.76) / units;
    } else {
        k18 = 1.0e-8 / units;
    }
    return k18;
}

//Calculation of k19 (H2I + HM --> H2I + HI)
double k19_rate(double T, double units, chemistry_data *my_chemistry)
{
    return 5.e-7 * sqrt(100.0/T) / units;
}

//Calculation of k20 (This is not currently used in the code)
double k20_rate(double T, double units, chemistry_data *my_chemistry)
{
    return tiny;
}

//Calculation of k21 (2HI + H2I --> H2I + H2I)
double k21_rate(double T, double units, chemistry_data *my_chemistry){
    return 2.8e-31 * pow(T, -0.6) / units;
}

//Calculation of k22 (2HI + HI --> H2I + HI)
double k22_rate(double T, double units, chemistry_data *my_chemistry)
{   
    double k22;
    switch (my_chemistry->three_body_rate) {

        case 0:
            if (T <= 300.0) {
                k22 = 1.3e-32 * pow(T/300.0, -0.38);
            } else {
                k22 = 1.3e-32 * pow(T/300.0, -1.0);
            }
            break;

        case 1:
            //PSS83 three-body rate
            k22 = 5.5e-29 / T;
            break;

        case 2:
            //CW83 three-body rate
            k22 = 8.8e-33;
            break;

        case 3:
            //FH07 three-body rate
            k22 = 1.44e-26 / pow(T, 1.54);
            break;

        case 4:
            //Rate from Glover (2008), derived by detailed balance from MSM96
            k22 = 7.7e-31 / pow(T, 0.464);
            break;
    
        case 5:
            //Rate from Forrey (2013)
            k22 = (6e-32 / pow(T, 0.25)) + (2e-31 / pow(T, 0.5));
            break;

        default:
            fprintf(stderr, "three_body_rate has been set to an unknown value: %d \n",
                my_chemistry->three_body_rate);
            exit(1);
    }
    return k22 / units;
}

//Calculation of k23.
double k23_rate(double T, double units, chemistry_data *my_chemistry)
{
    double k23;
    k23 = ( (8.125e-8/sqrt(T)) * exp(-52000.0/T) * (1.0 - exp(-6000.0/T)) ) / units;
    k23 = max(tiny, k23);
    return k23;
}

//Calculation of k50 (HII + DI --> HI + DII)
double k50_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // Fit taken from Savin (2002) which is valid for T < 2e5 K.
    // We extrapolate for higher temperatures.
    if (T <= 2.0e5) {
      return (2.0e-10 * pow(T, 0.402) * exp(-3.71e1/T)
              - 3.31e-17 * pow(T, 1.48)) / units;
    }
    else {
      return 2.5e-8 * pow(T/2.0e5, 0.402) / units;
    }
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 1.0e-9 * exp(-41.0 / T) / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}
        
//Calculation of k51 (HI + DII --> HII + DI)
double k51_rate(double T, double units, chemistry_data *my_chemistry)
{   
  if (my_chemistry->hd_reaction_rates == 0) {
      // Fit taken from Savin (2002) which is valid for T < 2e5 K.
      return (2.06e-10 * pow(T, 0.396) * exp(-3.30e1/T)
              + 2.03e-9 * pow(T, -0.332)) / units;
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 1.0e-9 / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of k52 (H2I + DII --> HDI + HII)
double k52_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // Fits from Galli & Palla (2002) to calculations by Gerlich (1982).
    // If T > 1e4 K use fixed value for k52 to avoid numerical issues with fitting function.
    // In this limit this reaction is not expected to be important anyway.
    if (T <= 1e4) {
      return 1.0e-9 * (0.417 + 0.846 * log10(T) - 0.137 * pow(log10(T), 2)) / units;
    } else {
      return 1.609e-9 / units;
    }
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 2.1e-9 / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of k53 (HDI + HII --> H2I + DII)
double k53_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // Fits from Galli & Palla (2002) to calculations by Gerlich (1982).
    return 1.1e-9 * exp(-4.88e2/T) / units;
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 1.0e-9 * exp(-464.0 / T) / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of k54 (H2I + DI --> HDI + HI)
double k54_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // Fit from Clark et al (2011), which is based on data in Mielke et al (2003).
    if (T <= 2.0e3) {
      return pow(1.0e1, (-5.64737e1 + 5.88886 * log10(T)
                         + 7.19692  * pow(log10(T), 2)
                         + 2.25069  * pow(log10(T), 3)
                         - 2.16903  * pow(log10(T), 4)
                         + 3.17887e-1 * pow(log10(T), 5)));
    } else {
      return 3.17e-10 * exp(-5.207e3 / T);
    }
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 7.5e-11 * exp(-3820.0 / T) / units;
    }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of k55 (HDI + HI --> H2I + DI)
double k55_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // Fit from Galli & Palla (2002), which is based on Shavitt (1959).
    // Fit has been modified at low temperature to avoid creating an
    // anomalously large rate coefficient -- as suggested by Ripamonti (2007)
    // and McGreer & Bryan (2008).
    if (T <= 2.0e2) {
      return 1.08e-22 / units;
    } else {
      return 5.25e-11 * exp(-4.43e3/T + 1.739e5/pow(T, 2)) / units;
    }
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 7.5e-11 * exp(-4240.0 / T) / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of k56 (DI + HM --> HDI + e)
double k56_rate(double T, double units, chemistry_data *my_chemistry)
{
  if (my_chemistry->hd_reaction_rates == 0) {
    // This is the same as DM + HI --> HDI + e
    // Measurements from Miller et al (2012) suggest there is no significant isotope effect
    // for this reaction.
    return k8_rate(T, units, my_chemistry);
  }
  // Stancil, Lepp & Dalgarno (1998)
  else if (my_chemistry->hd_reaction_rates == 1) {
    return 1.5e-9 * pow(T / 3.0e2, -0.1) / units;
  }
  else {
    fprintf(stderr, "hd_reaction_rates can only be 0 or 1.\n");
    exit(1);
  }
}
        
//Calculation of k57 (HI + HI --> HII + HI + e)
double k57_rate(double T, double units, chemistry_data *my_chemistry)
{
    // These rate coefficients are from Lenzuni, Chernoff & Salpeter (1991).
    // k57 value based on experimental cross-sections from Gealy & van Zyl (1987).
    if (T > 3.0e3) {
        return 1.2e-17  * pow(T, 1.2) * exp(-1.578e5 / T) / units;
    } else {
        return tiny;
    }
}

//Calculation of k58 (HI + HeI --> HII + HeI + e)
double k58_rate(double T, double units, chemistry_data *my_chemistry)
{
    // These rate coefficients are from Lenzuni, Chernoff & Salpeter (1991).
    // k58 value based on cross-sections from van Zyl, Le & Amme (1981).
    if (T > 3.0e3) {
        return 1.75e-17 * pow(T, 1.3) * exp(-1.578e5 / T) / units;
    } else {
        return tiny;
    }
}

//Calculation of h2dust (2H + grain --> H2 + grain)
double h2dust_rate(double T, double T_dust, double units, chemistry_data *my_chemistry)
{
    //Defined at the top of initialize_rates but pasted here for ease.
    double fgr = 0.009387;

    //Calculate convenient temperature parameters for this bin (in eV).
    double T_2 = T / 1.0e2;
    double T_dust_2 = T_dust / 1.0e2;

    //h2dust
    // The method used to calculate this is dependent upon what the user has selected
    // By default Omukai 2000 will be used.
    double h2dust;
    if (my_chemistry->h2_dust_rate == 1) {
        //k23 from Omukai (2000).

        h2dust = 6.0e-17 / fgr * pow(T / 300.0, 0.5) * 
                (pow(1.0 + exp(7.5e2 * ((1.0 / 75.0) - (1.0 / T_dust))), -1.0)) *
                (pow(1.0 + (4.0e-2 * pow(T + T_dust, 0.5))
                + (2.0e-3 * T) + (8.0e-6 * pow(T, 2.0)), -1.0));
        
    } else {
        //Equation 3.8 from Hollenbach & McKee (1979).

        h2dust = 3.0e-17 / fgr * pow(T_2, 0.5) / (1.0 + 0.4 * pow(T_2 + T_dust_2, 0.5)
                + 0.2 * T_2 + 8.0e-2 * pow(T_2, 2.0));

    }
    return h2dust / units;
}

//Calculation of n_cr_n.
double n_cr_n_rate(double T, double units, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.0e6 * pow(T, -0.5);
}

//Calculation of n_cr_d1.
double n_cr_d1_rate(double T, double units, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.6 * exp(-pow(400.0 / T, 2.0));
}

//Calculation of n_cr_d2.
double n_cr_d2_rate(double T, double units, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.4 * exp(-12000.0 / (T + 1200.0));
}

//Calculation of ceHI.
double ceHI_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates == 1){
        return 7.5e-19*exp( -fmin(log(dhuge), 118348.0 / T) )
                / ( 1.0 + sqrt(T / 1.0e5) ) / units;
    } else {
        return tiny;
    }
}

//Calculation of ceHeI.
double ceHeI_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates == 1){
        return 9.1e-27*exp(-fmin(log(dhuge), 13179.0/T))
                * pow(T, -0.1687) / ( 1.0 + sqrt(T/1.0e5) ) / units;
    } else {
        return tiny;
    }
}

//Calculation of ceHeII.
double ceHeII_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates == 1){
        return 5.54e-17*exp(-fmin(log(dhuge), 473638.0/T))
                * pow(T, -0.3970) / ( 1.0 + sqrt(T/1.0e5) ) / units;
    } else {
        return tiny;
    }
}

//Calculation of ciHeIS.
double ciHeIS_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 5.01e-27*pow(T, -0.1687) / ( 1.0 + sqrt(T/1.0e5) )
                * exp(-fmin(log(dhuge), 55338.0/T)) / units;
    } else {
        return tiny;
    }
}

//Calculation of ciHI.
double ciHI_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 2.18e-11 * k1_rate(T, 1, my_chemistry) / units;
    } else {
        return tiny;
    }
}

//Calculation of ciHeI.
double ciHeI_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 3.94e-11 * k3_rate(T, 1, my_chemistry) / units;
    } else {
        return tiny;
    }
}

//Calculation of ciHeII.
double ciHeII_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 8.72e-11 * k5_rate(T, 1, my_chemistry) / units; 
    } else {
        return tiny;
    }
}

//Calculation of reHII.
double reHII_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHI    = 2.0 * 157807.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 3.435e-30 * T * pow(lambdaHI, 1.970)
                    / pow( 1.0 + pow(lambdaHI/2.25, 0.376), 3.720)
                    / units;
        } else {
            return 1.778e-29 * T * pow(lambdaHI, 1.965)
                    / pow(1.0 + pow(lambdaHI/0.541, 0.502), 2.697)
                    / units; 
        }
    } else {
        return tiny;
    }
}

//Calculation of reHII.
double reHeII1_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHeII  = 2.0 * 285335.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 1.26e-14 * kboltz * T * pow(lambdaHeII, 0.75)
                            / units;
        } else {
            return 3e-14 * kboltz * T * pow(lambdaHeII, 0.654)
                    / units;
        }
    } else {
        return tiny;
    }
}

//Calculation of reHII2.
double reHeII2_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Dielectronic recombination (Cen, 1992).
    if (my_chemistry->recombination_cooling_rates == 1){
        return 1.24e-13 * pow(T, -1.5)
                * exp( -fmin(log(dhuge), 470000.0 / T) )
                * ( 1.0 + 0.3 * exp( -fmin(log(dhuge), 94000.0 / T) ) ) 
                / units;
    } else {
        return tiny;
    }
}

//Calculation of reHIII.
double reHeIII_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHeIII = 2.0 * 631515.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 8.0 * 3.435e-30 * T * pow(lambdaHeIII, 1.970)
                    / pow(1.0 + pow(lambdaHeIII / 2.25, 0.376), 3.720) 
                    / units;
        } else {
            return 8.0 * 1.778e-29 * T * pow(lambdaHeIII, 1.965)
                    / pow(1.0 + pow(lambdaHeIII / 0.541, 0.502), 2.697)
                    / units;
        }
    } else {
        return tiny;
    }
}

//Calculation of brem.
double brem_rate(double T, double units, chemistry_data *my_chemistry)
{
    if (my_chemistry->bremsstrahlung_cooling_rates == 1){
        return 1.43e-27 * sqrt(T)
                * ( 1.1 + 0.34 * exp( -pow(5.5 - log10(T), 2) / 3.0) )
                / units;
    } else {
        return tiny;
    }
}

//Calculation of vibh.
double vibh_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum;
    if (T > 1635.0) {
        par_dum = 1.0e-12 * sqrt(T) * exp(-1000.0 / T);
    } else {
        par_dum = 1.4e-13 * exp( (T / 125.0) - pow(T / 577.0, 2) );
    }

    return 1.1e-18 * exp( -fmin(log(dhuge), 6744.0 / T) ) / units;
}

//Calculation of hyd01k.
double hyd01k_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum;
    if (T > 1635.0) {
        par_dum = 1.0e-12 * sqrt(T) * exp(-1000.0 / T);
    } else {
        par_dum = 1.4e-13 * exp( (T / 125.0) - pow(T / 577.0, 2) );
    }

    return par_dum * exp( -fmin( log(dhuge), 8.152e-13 / (kboltz * T) ) )
            / units;
}
            
//Calculation of h2k01.
double h2k01_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum = 8.152e-13 * ( 4.2 / (kboltz * (T + 1190.0)) + 1.0 / (kboltz * T));

    return 1.45e-12 * sqrt(T) * exp(-fmin(log(dhuge), par_dum)) / units;
}

//Calculation of rotl.
double rotl_rate(double T, double units, chemistry_data *my_chemistry)
{
    double par_x = log10(T / 1.0e4); //Parameter used in the following calculation.

    if (T > 4031.0) {
        return 1.38e-22 * exp(-9243.0 / T) / units;
    } else {
        return pow(10.0, -22.9 - 0.553 * par_x - 1.148 * pow(par_x, 2)) / units;
    }
}

//Calculation of roth.
double roth_rate(double T, double units, chemistry_data *my_chemistry)
{
    double par_x = log10(T/1.0e4); //Parameter used in the following calculation.

    if(T > 1087.0) {
        return 3.9e-19 * exp(-6118.0 / T) / units;
    } else {
        return pow(10.0, -19.24 + 0.474*par_x - 1.247*pow(par_x, 2)) / units;
    }
}

//Calculation of GP99LowDensityLimit.
double GP99LowDensityLimit_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm = fmax(T, 13.0); //no cooling below 13 Kelvin
    tm = fmin(tm, 1.0e5); //fixes numerics
    double lt = log10(tm);

    return pow(10.0, -103.0 + 97.59*lt - 48.05*pow(lt, 2) + 10.8*pow(lt, 3)
            - 0.9032*pow(lt, 4)) / units;
}

//Calculatin of GP99HighDensityLimit.
double GP99HighDensityLimit_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm = fmax(T, 13.0); //no cooling below 13 Kelvin
    tm = fmin(tm, 1.0e5); //fixes numerics
    double t3 = tm / 1000.0;

    //Simplify formula for clarity.
    double HDLR = ( 9.5e-22*pow(t3, 3.76) ) / ( 1.0 + 0.12*pow(t3, 2.1) ) *
            exp( -pow((0.13 / t3), 3) ) + 3.0e-24 * exp(-0.51 / t3);
    double HDLV = 6.7e-19*exp(-5.86 / t3) + 1.6e-18*exp(-11.7 / t3);

    return (HDLR + HDLV) / units;
}

//Calculation of GAHI.
double GAHI_rate(double T, double units, chemistry_data *my_chemistry)
{
    /* -- Comments from original code -- 
    h2_h_cooling_rate == 1 :

        Revised low density H2 cooling rate due to H collisions, based on
        Lique (2015, MNRAS, 453, 810). This fit is accurate to within ~5% over the temperature
        range 100 < T < 5000 K. Lique (2015) doesn't present data above 5000 K, so at higher
        temperatures the rate has been calculated assuming that the de-excitation rate
        coefficients have the same values that they have at 5000 K. Lique also doesn't give
        rates for T < 100 K, but since we don't expect H2 cooling to be important there, it
        should be OK to just set the rate to zero.

    h2_h_cooling_rate == 2:

        Formula from Glover and Abel 2008.
    */

    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    //Calculate GAHI using user-specified method.
    if (my_chemistry->h2_h_cooling_rate == 1) {
        if (tm < 1e2) {
            return 0.0;
        } else {
            return pow(10.0, -24.07950609
                            + 4.54182810 * lt3
                            - 2.40206896 * pow(lt3, 2)
                            - 0.75355292 * pow(lt3, 3)
                            + 4.69258178 * pow(lt3, 4)
                            - 2.79573574 * pow(lt3, 5)
                            - 3.14766075 * pow(lt3, 6)
                            + 2.50751333 * pow(lt3, 7)) / units;
        }
    } else if (my_chemistry->h2_h_cooling_rate == 2) {
        if (tm < 1.0e2) {
            return pow(10.0, -16.818342
                            + 37.383713 * lt3
                            + 58.145166 * pow(lt3, 2)
                            + 48.656103 * pow(lt3, 3)
                            + 20.159831 * pow(lt3, 4)
                            + 3.8479610 * pow(lt3, 5)) / units;
        } else if (tm < 1.0e3) {
            return pow(10.0, -24.311209
                            + 3.5692468 * lt3
                            - 11.332860 * pow(lt3, 2)
                            - 27.850082 * pow(lt3, 3)
                            - 21.328264 * pow(lt3, 4)
                            - 4.2519023 * pow(lt3, 5)) / units;
        } else {
            return pow(10.0, -24.311209
                            + 4.6450521 * lt3
                            - 3.7209846 * pow(lt3, 2)
                            + 5.9369081 * pow(lt3, 3)
                            - 5.5108047 * pow(lt3, 4)
                            + 1.5538288 * pow(lt3, 5)) / units;
        }
    } else {
        fprintf(stderr, "h2_h_cooling_rate must be 1 or 2, it has been set \
                to %d \n", my_chemistry->h2_h_cooling_rate);
        exit(1);

    }
}

//Calculation of GAH2.
double GAH2_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    return pow(10.0, -23.962112
            + 2.09433740  * lt3
            - 0.77151436 * pow(lt3, 2)
            + 0.43693353 * pow(lt3, 3)
            - 0.14913216 * pow(lt3, 4)
            - 0.033638326 * pow(lt3, 5)) / units;
}

//Calculation of GAHe.
double GAHe_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    return pow(10.0, -23.689237
            + 2.1892372  * lt3
            - 0.81520438 * pow(lt3, 2)
            + 0.29036281 * pow(lt3, 3)
            - 0.16596184 * pow(lt3, 4)
            + 0.19191375 * pow(lt3, 5)) / units;
}

//Calculation of GAHp.
double GAHp_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    return pow(10.0, -22.089523
            + 1.5714711   * lt3
            + 0.015391166 * pow(lt3, 2)
            - 0.23619985  * pow(lt3, 3)
            - 0.51002221  * pow(lt3, 4)
            + 0.32168730  * pow(lt3, 5)) / units;
}

//Calculation of GAel.
double GAel_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    if (tm < 100.0) {
        return 0.0;
    } else if (tm < 500.0) {
        return pow(10.0, -21.928796
                + 16.815730 * lt3
                + 96.743155 * pow(lt3, 2)
                + 343.19180 * pow(lt3, 3)
                + 734.71651 * pow(lt3, 4)
                + 983.67576 * pow(lt3, 5)
                + 801.81247 * pow(lt3, 6)
                + 364.14446 * pow(lt3, 7)
                + 70.609154 * pow(lt3, 8)) / units;
    } else {
        return pow(10.0, -22.921189
                + 1.6802758  * lt3
                + 0.93310622 * pow(lt3, 2)
                + 4.0406627  * pow(lt3, 3)
                - 4.7274036  * pow(lt3, 4)
                - 8.8077017  * pow(lt3, 5)
                + 8.9167183  * pow(lt3, 6)
                + 6.4380698  * pow(lt3, 7)
                - 6.3701156  * pow(lt3, 8)) / units;
    }
}

//Calculation of H2LTE.
double H2LTE_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Constrain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 1.0e4);
    double lt3 = log10(tm / 1.0e3);

    if (tm < 1.0e2) {
        //Simple extrapolation as H2 cooling insignificant at these temperatures.
        return 7.0e-27 * pow(tm, 1.5) * exp(-512.0 / tm) / units;
    } else {
        return pow(10.0, -20.584225
                + 5.0194035 * lt3
                - 1.5738805 * pow(lt3, 2)
                - 4.7155769 * pow(lt3, 3)
                + 2.4714161 * pow(lt3, 4)
                + 5.4710750 * pow(lt3, 5)
                - 3.9467356 * pow(lt3, 6)
                - 2.2148338 * pow(lt3, 7)
                + 1.8161874 * pow(lt3, 8)) / units;
    }
}

//Calculation of HDlte.
double HDlte_rate(double T, double units, chemistry_data *my_chemistry)
{   
    //Fit from Coppola et al 2011. LTE (ergs/s) -> hdlte (ergs cm3/s)

    //Constain temperature.
    double tm  = fmax(T, 10.0);
    tm  = fmin(tm, 3.0e4);

    double HDlte = -55.5725 + 56.649 * log10(tm)
                    - 37.9102  * pow(log10(tm), 2)
                    + 12.698   * pow(log10(tm), 3)
                    - 2.02424  * pow(log10(tm), 4)
                    + 0.122393 * pow(log10(tm), 5);

    return pow(10.0, fmin(HDlte, 0.0)) / units;

}

//Calculation of HDlow.
double HDlow_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Rate based on (Wrathmall, Gusdorf & Flower, 2007) HD-H collisional excitation rates.
    //Constrain temperature.
    double tm = fmax(T, 1.0e1);
    tm = fmin(tm, 6.0e3);

    double HDlow = -23.175780  + 1.5035261 * log10(tm/1.0e3)
                    + 0.40871403  * pow(log10(tm/1.0e3), 2)
                    + 0.17849311  * pow(log10(tm/1.0e3), 3)
                    - 0.077291388 * pow(log10(tm/1.0e3), 4)
                    + 0.10031326 * pow(log10(tm/1.0e3), 5);
    
    return pow(10.0, HDlow) / units;
}

//Calculation of cie_thin_cooling_rate.
double cie_thin_cooling_rate(double T){

    /* 
    *  Comments found in the original fortran code:
    * 
    *     Compute optically thin cooling rate due to CIE cooling
    *     as discussed in Ripamonti and Abel 2003.
    * 
    *     input: temp is temperature in Kelvin
    * 
    *     output: total CIE cooling rate (including H2-H2 and He-H2)
    *             in erg per second times cm^3 per gram per (gram in H2 molecules).
    * 
    *      Hence to get cooling rate in erg/s cm^-3 one has to multiply by the total mass 
    *      density and the mass density in H2 molecules (*rho*rho_H2). This again is an 
    *      approximation in that it assumes large H2 fractions (greater than 0.5). Currently
    *      it is not clear what happens in the regime where H2 becomes mostly dissociated.
    */

    //* Compute rough extrapolations for extreme temperatures.
    // Low temperatures extrapolated with fourth power.
    if (T <= t_cie_c[0]) {
        return cie_table_c[0]*pow(T/t_cie_c[0], 4);
    }
    // High temperatures extrapolated with third power.
    if (T >= t_cie_c[287]) {
        return cie_table_c[287]*pow(T/t_cie_c[287], 3);
    }

    //* Compute CIE cooling rate for moderate temperatures.
    // Maximal and minimal indices in the cie tables.
    int minInd = 0;
    int maxInd = 287;
    int ind;
    for (int i = 0; i <= 200; i++) {
        // Calculate working index.
        ind = (maxInd + minInd) / 2;

        // Update extreme indices.
        if (T >= t_cie_c[ind]) {
            minInd = ind;
        } else {
            maxInd = ind;
        }

        // Return rate after convergence of indices.
        if ( (maxInd - minInd) <= 1 ) {
            return (cie_table_c[maxInd] * ( T - t_cie_c[minInd] ) 
                   + cie_table_c[minInd] * ( t_cie_c[maxInd] - T))
                   / ( t_cie_c[maxInd] - t_cie_c[minInd] );
        }
    }

    // the following logic was added as a quick way to address the compiler
    // warnings about a control flow that doesn't return get addressed.
    // - In the future, we may want to consider a way to handle failure more
    //   gracefully (if failure is even possible)
    // - for now, loudly exitting with a failure is better than quietly
    //   continuing with a garbage value
    fprintf(
      stderr,
      "INTERNAL ERROR: something went horribly wrong while computing the "
      "optically thin cooling rate due to CIE cooling. Aborting\n"
    );
    abort();
}

//Calculation of cieco.
double cieco_rate(double T, double units, chemistry_data *my_chemistry)
{
    double cierate = cie_thin_cooling_rate(T);

    return cierate * (mh/2.0) / units;
}

//Calculation of gas_grain.
double gasGrain_rate(double T, double units, chemistry_data *my_chemistry)
{
  double grain_coef;
  double fgr = 0.009387;
  if (my_chemistry->gas_grain_cooling_rate == 0) {
    //Calculate energy transfer from gas to dust grains (Equation 2.15, Hollenbach & McKee, 1989).
    //Normalize to the HM89 dust-to-gas ratio.
    grain_coef = 1.2e-31 * pow(1.0e3, -0.5) / fgr;
    return grain_coef * pow(T, 0.5) * (1.0 - 0.8 * exp(-75.0 / T)) / units;
  }
  else if (my_chemistry->gas_grain_cooling_rate == 1) {
    /*
      The rate depends on mass fraction and size distributino of grains.
      The heat transfer rate (Hollenbach & McKee 1989) is for MRN size
      distribution from 0.01 um to 0.25 um.
      The ISRF heating rate (Krumholz 2014) used here is for a uniform
      grain size distribution (a = 0.17 um), and optical depth (Omukai 2000)
      is for a MRN-like broken power-law.
      GC racalculated these rates for Omukai's dust model.
    */
    grain_coef = 2.57033e-32 * pow(1.033,-0.5) / fgr;
    double f_vel = 0.5 / sqrt(2.0) + 0.0833333 / sqrt(4.0);
    // Hollenbach & McKee (1989) considered the contribution of other species
    // than protons and charged grains, but we now consider only H2 and He
    // and neglect charged grains (Schneider et al. 2006).
    return grain_coef * f_vel * pow(T, 0.5) / units;
  }
  else {
    fprintf(stderr, "gas_grain_cooling_rate can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of gas_grain2.
double gasGrain2_rate(double T, double units, chemistry_data *my_chemistry)
{
    //Variables.
    double f_vel = 0.5 / sqrt(2.0) + 0.0833333 / sqrt(4.0);
    double vH_avg = sqrt(kboltz * T / 2.0 / pi / mh);

    return f_vel * 4.0 * vH_avg * 2.0 * kboltz * mh / units;
    //Later multiplied by sigma_gr / mass_gr for arbitrary size distribution.
}


//Calculation of regr.
double regr_rate(double T, double units, chemistry_data *my_chemistry)
{
    //(Equation 9, Wolfire et al., 1995)
    double grbeta = 0.74 / pow(T, 0.068);
    return  4.65e-30 * pow(T, 0.94 + 0.5 * grbeta) / units;
}

//The below rates are scalar -- they have no temperature dependence.

//Calculation of comp.
double comp_rate(double units, chemistry_data *my_chemistry)
{
    return 5.65e-36 / units;
}

//Calculation of gammah.
double gammah_rate(double units, chemistry_data *my_chemistry)
{
    //Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3).
    if (my_chemistry->photoelectric_heating <= 1) {
        return my_chemistry->photoelectric_heating_rate / units;
    } else { //photoelectric_heating set to 2 or 3
        //User to specify G_0, epsilon set to 0.05 or calculated directly.
        return 1.0e-24 / units;
    }
}

//Calculation of gamma_isrf.
double gamma_isrf_rate(double units, chemistry_data *my_chemistry)
{
  // Parameter definition.
  double fgr = 0.009387;
  if (my_chemistry->uniform_grain_isrf_heating_rate == 0) {
    // (Equation B15, Krumholz, 2014)
    // Don't normalize by coolunit since tdust calculation is done in CGS.
    return 3.9e-24 / mh / fgr;
  }
  else if (my_chemistry->uniform_grain_isrf_heating_rate == 1) {
    // For uniform grain size (Goldsmith 2001; Krumholz 2014)
    return 8.60892e-24 / (2.0 * mh) / fgr;
    // F_isrf sigma_gr / mass_gr
    // For MRN-like broken power low size distribution (Omukai 2000)
    // The factor 2 to cancel out the molecular mass of H2.
  }
  else {
    fprintf(stderr, "uniform_grain_isrf_heating_rate can only be 0 or 1.\n");
    exit(1);
  }
}

//Calculation of gamma_isrf2.
double gamma_isrf2_rate(double units, chemistry_data *my_chemistry)
{
    return 5.3e-3;
    //MW interstellar radiation field (Goldsmith 2001)
    //Later multiplied by total grain cross-section per unit dust mass
    //for arbitrary size distribution
}

//Calculation of grain growth rate.
double grain_growth_rate(double T, double units, chemistry_data *my_chemistry)
{
    double vH_avg = sqrt( kboltz * T / 2.0 / pi / mh);
    return 4.0 * vH_avg * mh / units;
    // Factor of 4 because gas-phase molecules are accreted
    // onto the entire surface area of grains.
}

//Calculation of H2 formation rate on S grain surfaces. (Cazaux & Tielens 2002)
double h2dust_S_rate(double T, double T_dust, double units, chemistry_data *my_chemistry)
{
    // Constants. 
    double E_HC_Silicate = 200.0, E_HP_Silicate = 650.0, E_S_Silicate = 3.0e4;

    // Sticking probability.
    double S_H = pow(1.0 + 0.4*pow((T + T_dust)/100.0, 0.5) + 0.2 * (T/100.0) +
                 0.08*pow(T/100.0,2.0), -1.0);

    // Calculate variables.
    double bHP_aPC = 1.0/4.0 * pow(1.0 + sqrt((E_HC_Silicate - E_S_Silicate) /
                     (E_HP_Silicate - E_S_Silicate)), 2.0) * exp(-E_S_Silicate / T_dust);
    double epsilon_H2 = pow(1.0 + bHP_aPC, -1.0);
    double vH_avg = sqrt( kboltz * T / 2.0 / pi / mh);

    return 0.5 * 4.0 * vH_avg * S_H * epsilon_H2 * mh / units;
    // Factor of 4 because gas-phase molecules are accreted
    // onto the entire surface area of grains.
}

//Calculation of H2 formation rate on C grain surfaces. (Cazaux & Tielens 2002)
double h2dust_C_rate(double T, double T_dust, double units, chemistry_data *my_chemistry)
{
    // Constants.
    double E_HC_AmCarbon = 250.0, E_HP_AmCarbon = 800.0, E_S_AmCarbon = 3.0e4;

    // Sticking probability.
    double S_H = pow(1.0 + 0.4*pow((T + T_dust)/100.0, 0.5) + 0.2 * (T/100.0) +
                 0.08*pow(T/100.0,2.0), -1.0);

    // Calculate variables.
    double bHP_aPC = 1.0/4.0 * pow(1.0 + sqrt((E_HC_AmCarbon - E_S_AmCarbon) /
                     (E_HP_AmCarbon - E_S_AmCarbon)), 2.0) * exp(-E_S_AmCarbon / T_dust);
    double epsilon_H2 = pow(1.0 + bHP_aPC, -1.0);
    double vH_avg = sqrt( kboltz * T / 2.0 / pi / mh);

    return 0.5 * 4.0 * vH_avg * S_H * epsilon_H2 * mh / units;
    // Factor of 4 because gas-phase molecules are accreted
    // onto the entire surface area of grains.
}
