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
#include "phys_constants.h"

// Calculation of k1.
double k1_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k2_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    return 1;
}

//Calculation of k3.
double k3_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k4_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k5_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k6_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k7_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Fit --> Stancil, Lepp & Dalgarno (1998, ApJ, 509, 1). Based on photodetachment cross-section --> Wishart (1979, MNRAS, 187, P59).
    return 3.0e-16*pow(T/3.0e2, 0.95) * exp(-T/9.32e3) / kUnit;
}

//Calculation of k8.
double k8_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Fit based on experimental measurements --> Kreckel et al (2010, Science, 329, 69).
    return 1.35e-9 * ( pow(T, 9.8493e-2) + 3.2852e-1*pow(T, 5.5610e-1) + 2.771e-7*pow(T, 2.1826) )
                / ( 1.0 + 6.191e-3*pow(T, 1.0461) + 8.9712e-11*pow(T, 3.0424) + 3.2576e-14*pow(T, 3.7741) )
                / kUnit;
}

//Calculation of k9.
double k9_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k10_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    return 6.0e-10 / kUnit;
}

//Calculation of k11.
double k11_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k12_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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

//Calculation of k13. //! Check I have handled this correctly.
double k13_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->three_body_rate == 0) { //! Added this if statement as this bit always calculated and was then overwritten.
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
    //! Following if statements are from 3-body rate
    } else if (my_chemistry->three_body_rate == 1) {
        //Inverse of PSS83 three-body rate
        return (5.24e-7 / pow(T, 0.485)) * exp(-5.2e4 / T);
    } else if (my_chemistry->three_body_rate == 2) {
        //Inverse of CW83 rate
        return 8.4e-11 * pow(T, 0.515) * exp(-5.2e4 / T);
    } else if (my_chemistry->three_body_rate == 3) {
        //Rate from JGC67, used by FH07 to construct their three-body rate
        return (1.38e-4 / pow(T, 1.025)) * exp(-5.2e4 / T);
    } else if (my_chemistry->three_body_rate == 4) {
        //High density limiting rate from MSM96
        return pow(1e1 ,(-178.4239 - 68.42243 * log10(T)
                        + 43.20243 * pow(log10(T), 2)
                        - 4.633167 * pow(log10(T), 3) 
                        + 69.70086 * log10(1.0 + 40870.38 / T)
                        - (23705.7 / T)));
    } else if (my_chemistry->three_body_rate == 5) {
        //From detailed balance from Forrey rate
        if (T <= 3000.0) {
            return 2.4e-8 * exp(-5.2e4/T);
        } else {
            return 2.2e-6 * pow(T, -0.565) * exp(-5.2e4/T); 
        }
    } else {
        fprintf(stderr, 'three_body_rate has been set to an unknown value: %d.\n',
            my_chemistry->three_body_rate);
        return(FAIL); //! NOT SURE IF I CAN RETURN FAIL FROM THIS
    }
}

//Calculation of k14.
double k14_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k15_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k16_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Fit --> Croft et al (1999, MNRAS, 304, 327). Based on cross-section --> Fussen & Kubach (1986, J. Phys. B, 18 L31).
    return 2.4e-6*(1.0 + T/2.0e4) / sqrt(T) / kUnit;
}

//Calculation of k17.
double k17_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k18_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
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
double k19_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    return 5.e-7 * sqrt(100.0/T) / kUnit;
}

//Calculation of k21.
double k21_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry){
    return 2.8e-31 * pow(T, -0.6);
}

//Calculation of k22. //! CHECK WITHIN FUNCTION.
double k22_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{   
    double k22;
    if (my_chemistry->three_body_rate == 0) {
        if (T <= 300.0) {
            return 1.3e-32 * pow(T/300.0, -0.38);
        } else {
            return 1.3e-32 * pow(T/300.0, -1.0);
        }
    } else if (my_chemistry->three_body_rate == 1) {
        //PSS83 three-body rate
        return 5.5e-29 / T;
    } else if (my_chemistry->three_body_rate == 2) {
        //CW83 three-body rate
        return 8.8e-33;
    } else if (my_chemistry->three_body_rate == 3) {
        //FH07 three-body rate
        return 1.44e-26 / pow(T, 1.54);
    } else if (my_chemistry->three_body_rate == 4) {
        //Rate from Glover (2008), derived by detailed balance from MSM96
        return 7.7e-31 / pow(T, 0.464);
    } else if (my_chemistry->three_body_rate == 5) {
        //Rate from Forrey (2013)
        return (6e-32 / pow(T, 0.25)) + (2e-31 / pow(T, 0.5));
    } else {
        fprintf(stderr, 'three_body_rate has been set to an unknown value: %d.\n',
            my_chemistry->three_body_rate);
        return(FAIL); //! NOT SURE IF I CAN RETURN FAIL FROM THIS
    }
}

//Calculation of k23.
double k23_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    double k23;
    k23 = ( (8.125e-8/sqrt(T)) * exp(-52000.0/T) * (1.0 - exp(-6000.0/T)) ) / kUnit;
    k23 = max(tiny, k23);
    return k23;
}

//Calculation of k50.
double k50_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Fit taken from Savin (2002) which is valid for T < 2e5 K.
    //We extrapolate for higher temperatures.
    if (T <= 2.0e5) {
        return (2.0e-10 * pow(T, 0.402) * exp(-3.71e1/T)
                    - 3.31e-17 * pow(T, 1.48)) / kUnit;
    } else {
        return 2.5e-8 * pow(T/2.0e5, 0.402) / kUnit;
    }
}
        
//Calculation of k51.
double k51_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // Fit taken from Savin (2002) which is valid for T < 2e5 K.
    return (2.06e-10 * pow(T, 0.396) * exp(-3.30e1/T)
                        + 2.03e-9 * pow(T, -0.332)) / kUnit;
}

//Calculation of k52.
double k52_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // Fits from Galli & Palla (2002) to calculations by Gerlich (1982).
    // If T > 1e4 K use fixed value for k52 to avoid numerical issues with fitting function.
    // In this limit this reaction is not expected to be important anyway.
    if (T <= 1e4) {
    return 1.0e-9 * (0.417 + 0.846 * log10(T) - 0.137 * pow(log10(T), 2)) / kUnit;
    } else {
        return 1.609e-9 / kUnit;
    }   
}

//Calculation of k53.
double k53_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // Fits from Galli & Palla (2002) to calculations by Gerlich (1982).
    return 1.1e-9 * exp(-4.88e2/T) / kUnit;
}

//Calculation of k54.
double k54_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
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

//Calculation of k55.
double k55_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // Fit from Galli & Palla (2002), which is based on Shavitt (1959).
    // Fit has been modified at low temperature to avoid creating an 
    // anomalously large rate coefficient -- as suggested by Ripamonti (2007)
    // and McGreer & Bryan (2008).
    if (T <= 2.0e2) {
        return 1.08e-22 / kUnit;
    } else {
        return 5.25e-11 * exp(-4.43e3/T + 1.739e5/pow(T, 2)) / kUnit;
    }
}

//Calculation of k56. //! Not sure how to deal with this. Maybe just not add a function and copy the other array in calc_rates?
double k56_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // This is the same as HM + HI. 
    // Measurements from Miller et al (2012) suggest there is no significant isotope effect
    // for this reaction.
    my_rates->k56[i] = my_rates->k8[i];
}
        
//Calculation of k57.
double k57_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // These rate coefficients are from Lenzuni, Chernoff & Salpeter (1991).
    // k57 value based on experimental cross-sections from Gealy & van Zyl (1987).
    if (T > 3.0e3) {
        return 1.2e-17  * pow(T, 1.2) * exp(-1.578e5 / T) / kUnit;
    } else {
        return tiny;
    }
}

//Calculation of k58.
double k58_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    // These rate coefficients are from Lenzuni, Chernoff & Salpeter (1991).
    // k58 value based on cross-sections from van Zyl, Le & Amme (1981).
    if (T > 3.0e3) {
        return 1.75e-17 * pow(T, 1.3) * exp(-1.578e5 / T) / kUnit;
    } else {
        return tiny;
    }
}

//Calculation of n_cr_n.
double n_cr_n_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.0e6 * pow(T, -0.5);
}

//Calculation of n_cr_d1.
double n_cr_d1_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.6 * exp(-pow(400.0 / T, 2.0));
}

//Calculation of n_cr_d2.
double n_cr_d2_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //H2 formation heating terms from Equation 23, Omuaki (2000).
    return 1.4 * exp(-12000.0 / (T + 1200.0));
}

//Calculation of ceHI.
double ceHI_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates == 1){
        return 7.5e-19*exp(-fmin(log(1.0e30), 118348.0/T))
                / ( 1.0 + sqrt(T/1.0e5) ) / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ceHeI.
double ceHeI_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates == 1){
        return 9.1e-27*exp(-fmin(log(1.0e30), 13179.0/T))
                * pow(T, -0.1687) / ( 1.0 + sqrt(T/1.0e5) ) / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ceHeII.
double ceHeII_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_excitation_rates = 1){
        return 5.54e-17*exp(-fmin(log(1.0e30), 473638.0/T))
                * pow(T, -0.3970) / ( 1.0 + sqrt(T/1.0e5) ) / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ciHeIS.
double ciHeIS_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 5.01e-27*pow(T, -0.1687) / ( 1.0 + sqrt(T/1.0e5) )
                * exp(-fmin(log(1.0e30), 55338.0/T)) / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ciHI.
double ciHI_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry) //! Relies on another rate
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 2.18e-11 * my_rates->k1[i] * kUnit / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ciHeI.
double ciHeI_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry) //! Relies on another rate
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 3.94e-11 * my_rates->k3[i] * kUnit / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of ciHeII.
double ciHeII_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry) //! Relies on another rate
{
    //Collisional ionization. Polynomial fit from Tom Abel.
    if (my_chemistry->collisional_ionisation_rates == 1){
        return 8.72e-11 * my_rates->k5[i] * kUnit / coolingUnits; 
    } else {
        return tiny;
    }
}

//Calculation of reHII.
double reHII_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHI    = 2.0 * 157807.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 3.435e-30 * T * pow(lambdaHI, 1.970)
                    / pow( 1.0 + pow(lambdaHI/2.25, 0.376), 3.720)
                    / coolingUnits;
        } else {
            return 1.778e-29 * T * pow(lambdaHI, 1.965)
                    / pow(1.0 + pow(lambdaHI/0.541, 0.502), 2.697)
                    / coolingUnits; 
        }
    } else {
        return tiny;
    }
}

//Calculation of reHII.
double reHeII1_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHI    = 2.0 * 157807.0 / T;
        double lambdaHeII  = 2.0 * 285335.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 1.26e-14 * kboltz * T * pow(lambdaHeII, 0.75)
                            / coolingUnits;
        } else {
            return 3e-14 * kboltz * T * pow(lambdaHeII, 0.654)
                    / coolingUnits;
        }
    } else {
        return tiny;
    }
}

//Calculation of reHII2.
double reHeII2_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry){
    //Dielectronic recombination (Cen, 1992).
    if (my_chemistry->recombination_cooling_rates == 1){
        return 1.24e-13 * pow(T, -1.5)
                * exp(-fmin(log(1.0e30), 470000.0/T))
                * ( 1.0 + 0.3*exp(-fmin(log(1.0e30), 94000.0/T))) 
                / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of reHIII.
double reHeIII_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->recombination_cooling_rates == 1){
        //Define parameters used in the calculations.
        double lambdaHeII  = 2.0 * 285335.0 / T;
        double lambdaHeIII = 2.0 * 631515.0 / T;

        //These depend on if the user has chosen recombination case A or B.
        if (my_chemistry->CaseBRecombination == 1) {
            return 8.0 * 3.435e-30 * T * pow(lambdaHeIII, 1.970)
                    / pow(1.0 + pow(lambdaHeIII/2.25, 0.376), 3.720) 
                    / coolingUnits;
        } else {
            return 8.0 * 1.778e-29 * T * pow(lambdaHeIII, 1.965)
                    / pow(1.0 + pow(lambdaHeIII/0.541, 0.502), 2.697)
                    / coolingUnits;
        }
    } else {
        return tiny;
    }
}

//Calculation of brem.
double brem_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    if (my_chemistry->bremsstrahlung_cooling_rates == 1){
        return 1.43e-27 * sqrt(T)
                * (1.1 + 0.34 * exp(-pow(5.5 - log10(T), 2) / 3.0))
                / coolingUnits;
    } else {
        return tiny;
    }
}

//Calculation of vibh.
double vibh_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum;
    if (T > 1635.0) {
        par_dum = 1.0e-12 * sqrt(T) * exp(-1000.0/T);
    } else {
        par_dum = 1.4e-13 * exp((T/125.0) - pow(T/577.0, 2));
    }

    return 1.1e-18 * exp(-fmin(log(1.0e30),6744.0/T)) / coolingUnits;
}

//Calculation of hyd01k.
double hyd01k_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum;
    if (T > 1635.0) {
        par_dum = 1.0e-12 * sqrt(T) * exp(-1000.0/T);
    } else {
        par_dum = 1.4e-13 * exp((T/125.0) - pow(T/577.0, 2));
    }

    return par_dum * exp(-fmin(log(1.0e30), 8.152e-13/(kboltz*T)))
            / coolingUnits;
}
            
//Calculation of h2k01.
double h2k01_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Dummy parameter used in the calculation.
    double par_dum = 8.152e-13 * ( 4.2 / (kboltz * (T + 1190.0)) + 1.0 / (kboltz * T));

    return 1.45e-12 * sqrt(T) * exp(-fmin(log(1.0e30), par_dum)) / coolingUnits;
}

//Calculation of rotl.
double rotl_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    double par_x = log10(T/1.0e4); //Parameter used in the following calculation.

    if (T > 4031.0) {
        return 1.38e-22 * exp(-9243.0/T) / coolingUnits;
    } else {
        return pow(10.0, -22.9 - 0.553*par_x - 1.148*pow(par_x, 2)) / coolingUnits;
    }
}

//Calculation of roth.
double roth_rate(double T, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    double par_x = log10(T/1.0e4); //Parameter used in the following calculation.

    if(T > 1087.0) {
        return 3.9e-19 * exp(-6118.0/T) / coolingUnits;
    } else {
        return pow(10.0, -19.24 + 0.474*par_x - 1.247*pow(par_x, 2)) / coolingUnits;
    }
}




            