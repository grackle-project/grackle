/*******************************************************************************/
/*
 * ----------Function to Calculate Multispecies Rate Lookup Table---------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 18/2/2021
 * Date of Last Modification: 18/2/2021
 * 
 * Purpose:
 *  Calculates lookup tables for the rate coefficients (and some cooling functions) as a function 
 *  of log(temperature).
 * 
 * Inputs:
 *  *userChemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *userRates : Pointer to the structure of type chemistry_data_storage
 *               which holds pointers to the rates.
 *  *userUnits : Pointer to the structure which defines the units for 
 *               the variables.
 *  *userCoLengthUnit : Pointer to the comoving length unit.
 *  *userCoDensityUnit : Pointer to the comoving density unit.
 *  *iOutput : Integer value which is modified when function finishes.
 * 
 * Outputs:* 
 *  
 * Units:
 *  Most rate coefficients are calculated using cgs units and then converted
 *  for convenience. Units and constants are included explicitly where the 
 *  code units are not used. Also be aware that comoving densities are NOT
 *  used in the rate equations.
 * 
 * Rate Coefficients:
 *   --All rates are labelled as in Abel et al., 1996 (astro-ph/9608040)--
 *   @  k1  @     HI + e --> HII + 2e
 *   @  k2  @     HI + e --> HI + photon
 *   @  k3  @     HeI + e --> HeII + 2e
 *   @  k4  @     HeII + e --> HeI + photon
 *   @  k5  @     HeII + e --> HeIII + 2e
 *   @  k6  @     HeIII + e --> HeII + photon
 *   @  k7  @     HI + e --> HM + photon
 *   @  k8  @     HI + HM --> H2I* + e
 *   @  k9  @     HI + HII --> H2II + photon
 *   @  k10 @     H2II + HI --> H2I* + HII
 *   @  k11 @     H2I + HII --> H2II + HI
 *   @  k12 @     H2I + e --> 2HI + e
 *   @  k13 @     H2I + HI --> 3HI
 *   @  k14 @     HM + e --> HI + 2e
 *   @  k15 @     HM + HI --> 2HI + e
 *   @  k16 @     HM + HI --> 2HI
 *   @  k17 @     HM + HI --> H2I + e
 *   @  k18 @     H2I + e --> 2HI
 *   @  k19 @     H2I + HM --> H2I + HI
 *   @  k20 @     Not Used
 *   @  k21 @     2HI + H2I --> H2I + H2I
 *   @  k22 @     2HI + HI --> H2I + HI
 *   @  k24 @     HI + p --> HII + e
 *   @  k25 @     HeIII + p --> HeII + e
 *   @  k26 @     HeI + p --> HeII + e
 *   @  k27 @     HM + p --> HI + e
 *   @  k28 @     H2II + p --> HI + HII
 *   @  k29 @     H2I + p --> H2II + e
 *   @  k30 @     H2II + p --> 2HII + e
 *   @  k31 @     H2I + p --> 2HI
 *   @  k50 @     HII + DI --> HI + DII
 *   @  k51 @     HI + DII --> HII + DI
 *   @  k52 @     H2I + DII --> HDI + HII
 *   @  k53 @     HDI + HII --> H2I + DII
 *   @  k54 @     H2I + DI --> HDI + HI
 *   @  k55 @     HDI + HI --> H2I + DI
 *   @  k56 @     DI + HM --> HDI + e
 *   @      @     DM + HI --> HDI + e  //This is included implicitly by multiplying k56 by two as they are assumed to have the same rate.
 *   @  k57 @     HI + HI --> HII + HI + e
 *   @  k58 @     HI + HeI --> HII + HeI + e
 *   
 *   @  h2dust @     2H + grain --> H2 + grain  //(Hollenback and McKee 1989, Equation 3.8) 
 * 
 * 
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
#include "phys_constants.h"

//* Prototypes for functions used in calc_rates_g_c.c
//Multispecies collisional rates
void coll_rates_g_c(int temperatureBinIndex, double T, double kUnit, chemistry_data *userChemistry, chemistry_data_storage *userRates);
//Optically thin CIE cooling rate
double cie_thin_cooling_rate_g_c(double T);
//Density-dependent H2  dissociation rate
void colh2diss_g_c(int temperatureBinIndex, double T, int idt, chemistry_data *userChemistry, chemistry_data_storage *userRates);

//* Definition of the function.
int calc_rates_g_c(chemistry_data *userChemistry, chemistry_data_storage *userRates, code_units *userUnits, double userCoLengthUnit, double userCoDensityUnit)
{

    //* Set various constant parameters that will be used throughout.
    double tevk = 1.1605e4; //Kelvin to eV conversion factor
    double dhuge = 1.0e30; //Comparison value
    // Cutoff values (eV)
    double e24  = 13.6;
    double e25  = 54.4;
    double e26  = 24.6;
    double e27  = 0.755;
    double e28a = 2.65;
    double e28b = 11.27;
    double e28c = 21.0;
    double e29a = 15.42;
    double e29b = 16.5;
    double e29c = 17.7;
    double e30a = 30.0;
    double e30b = 70.0;


    //* Set the flag for dust calculations.
    int anyDust;
    if ( userChemistry->h2_on_dust > 0 || userChemistry->dust_chemistry > 0) {
        anyDust = TRUE;
    } else {
        anyDust = FALSE;
    }

    //* Obtain the conversion factors which convert between code and physical units. We define a_value = 1 at z = zInit such that a = a_value * [a].
    double timeBase1 = userUnits->time_units;
    double lengthBase1 = userCoLengthUnit/(userUnits->a_value * userUnits->a_units);
    double densityBase1 = userCoDensityUnit*pow(userUnits->a_value * userUnits->a_units, 3);

    /*
    * 1) Set the dimensions of the non-radiative rate coefficients.
    * 
    *   --Following explanation taken from original fortran code--
    * 
    * 
    *   Note that we have included the units that convert density to 
    *   number density, so the rate equations should look like:
    * 
    *       d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
    * 
    *   where k~ is the dimenionless rate coefficients and d0-2~ are three
    *   dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
    *   dimensionless expansion coefficient.
    * 
    *   rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
    *   rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
    *   rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
    * 
    *   So, k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] ) [~]
    *   
    *   Reminder: the number densities here are normalized with [dens] which is not
    *   a constant (it has a factor a^3), so the number densities must be converted
    *   from comoving to proper.
    */ 
    double kUnit = (pow(userUnits->a_units, 3) * mh) / (densityBase1 * timeBase1);
    double kUnit_3Bdy = kUnit * (pow(userUnits->a_units, 3) * mh) / densityBase1;

    /*
    * 2) Set the dimensions of the cooling coefficients.
    * 
    *   --Following explanation taken from original fortran code--
    * 
    * 
    *   This equation has a rho because e is the specific energy, not energy/unit volume).
    * 
    *   delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
    *   [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
    *   delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
    * 
    *   So, L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
    *   but [e] = ([a]*[x])**2 / [time]**2, and [a] = 1 / (1 + zri)
    *   
    *   Therefore, [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
    * 
    * 
    *   Note: some of the coffiecients have only one power of n.  These do not have the /a^3
    *   factor, also they have units:
    *     
    *     [L1] = ([a]**2 * [x]**2 * mh) / [time]**3  =  [L] * [dens] * [a]**3 / mh
    *       
    *   This is done through the dom variable in cool.src (some have three powers of n and they
    *   are different by the reciprocal of the above factor multiplying [L]).
    * 
    */
    double coolingUnits = (pow(userUnits->a_value, 5) * pow(lengthBase1, 2) * pow(mh, 2)) \
                          / (densityBase1 * pow(timeBase1, 3));

    //* 3) Units for radiative transfer coefficients are 1/[time].

    //* 4) Calculate energy transfer from gas to dust grains (Equation 2.15, Hollenbach & McKee, 1989).
    //*    Normalize to the HM89 dust-to-gas ratio.
    double fgr = 0.009387;
    double grainCoeff = 1.2e-31 * pow(1.0e3, -0.5) / fgr;

    //* Compute log spacing in temperature.
    double TempStart = userChemistry->TemperatureStart;
    double logTempStart = log(TempStart);
    double dlogTemp = (log(userChemistry->TemperatureEnd) - logTempStart) \
                           / fabs(userChemistry->NumberOfTemperatureBins - 1);  

    //* Compute log spacing in dust temperature.
    double dust_TempStart = userChemistry->DustTemperatureStart;
    double dust_logTempStart = log(dust_TempStart);
    double dust_dlogTemp = (log(userChemistry->DustTemperatureEnd) - dust_logTempStart) \
                                / fabs(userChemistry->NumberOfDustTemperatureBins - 1);

    //* Compute rates for primordial chemistry.
    //All rates are calculated if primordial_chemistry > 0 as it is more convenient and there
    //is not much cost to doing it this way as this code only has to run once.
    if (userChemistry->primordial_chemistry > 0) {
        //Initialise constants to tiny.
        for (int i = 0; i < userChemistry->NumberOfTemperatureBins; i++) {
            userRates->k1[i] = tiny;
            userRates->k2[i] = tiny;
            userRates->k3[i] = tiny;
            userRates->k4[i] = tiny;
            userRates->k5[i] = tiny;
            userRates->k6[i] = tiny;
            userRates->k7[i] = tiny;
            userRates->k8[i] = tiny;
            userRates->k9[i] = tiny;
            userRates->k10[i] = tiny;
            userRates->k11[i] = tiny;
            userRates->k12[i] = tiny;
            userRates->k13[i] = tiny;
            userRates->k14[i] = tiny;
            userRates->k15[i] = tiny;
            userRates->k16[i] = tiny;
            userRates->k17[i] = tiny;
            userRates->k18[i] = tiny;
            userRates->k19[i] = tiny;
            userRates->k20[i] = tiny;
            userRates->k21[i] = tiny;
            userRates->k22[i] = tiny;
            userRates->k50[i] = tiny;
            userRates->k51[i] = tiny;
            userRates->k52[i] = tiny;
            userRates->k53[i] = tiny;
            userRates->k54[i] = tiny;
            userRates->k55[i] = tiny;
            userRates->k56[i] = tiny;
            userRates->k57[i] = tiny;
            userRates->k58[i] = tiny;

            //k13dd and h2dust were originally 2D arrays in the fortran code: I use 1D arrays here with the following mappings,
            //k13dd:
            //  {k13dd(0,0), k13dd(1,0), ... k13dd(NumberOfTemperatureBins, 0), k13dd(0,1), 2D(1,1)..., k13dd(NumberOfTemperatureBins, 13)}
            //h2dust:
            //  {h2dust(0,0), h2dust(1,0), ... h2dust(NumberOfTemperatureBins, 0), h2dust(0,1), 2D(1,1)..., h2dust(NumberOfTemperatureBins, NumberOfDustTemperatureBins)}
            //where k13dd(i,j) and h2dust(i,j) are the elements in the corresponding 2D arrays.
            for (int j = 0; j < 14; j++) {
                userRates->k13dd[i + userChemistry->NumberOfTemperatureBins*j] = tiny;
            }
            for (int j = 0; j < userChemistry->NumberOfDustTemperatureBins; j++) {
                userRates->h2dust[i + userChemistry->NumberOfTemperatureBins*j] = tiny;
            }

            userRates->n_cr_n[i]  = tiny;
            userRates->n_cr_d1[i] = tiny;
            userRates->n_cr_d2[i] = tiny;

            userRates->ceHI[i]    = tiny;
            userRates->ceHeI[i]   = tiny;
            userRates->ceHeII[i]  = tiny;
            userRates->ciHI[i]    = tiny;
            userRates->ciHeI[i]   = tiny;
            userRates->ciHeIS[i]  = tiny;
            userRates->ciHeII[i]  = tiny;
            userRates->reHII[i]   = tiny;
            userRates->reHeII1[i] = tiny;
            userRates->reHeII2[i] = tiny;
            userRates->reHeIII[i] = tiny;
            userRates->brem[i]    = tiny;
            userRates->comp       = tiny;

            userRates->hyd01k[i]  = tiny;
            userRates->h2k01[i]   = tiny;
            userRates->vibh[i]    = tiny;
            userRates->roth[i]    = tiny;
            userRates->rotl[i]    = tiny;
            userRates->HDlte[i]   = tiny;
            userRates->HDlow[i]   = tiny;
            userRates->cieco[i]   = tiny;

            userRates->GAHI[i]    = tiny;
            userRates->GAH2[i]    = tiny;
            userRates->GAHe[i]    = tiny;
            userRates->GAHp[i]    = tiny;
            userRates->GAel[i]    = tiny;
            userRates->H2LTE[i]   = tiny;
        } //End of parameter initialisation for loop.

    //* Populate tables for values between the specificed start and end temperatures.
    //* _____________________________________________________________________________

    //* 1) Rate Coefficients (excluding the external radiation field)
    //Loop over temperature bins.
    for (int i = 0; i < userChemistry->NumberOfTemperatureBins; i++) {
        //Calculate convenient temperature parameters for this bin (in eV).
        double logBinTemp = logTempStart + fabs(i)*dlogTemp;
        double binTemp = exp(logBinTemp);
        double binTemp2 = binTemp / 1.0e2;
        double binTemp_eV = binTemp / tevk;
        double logBinTemp_eV = log(binTemp_eV);


        //--------Calculate multispecies collissional rates--------

        coll_rates_g_c(i, binTemp, kUnit, userChemistry, userRates);


        //--------Calculate coefficients for density-dependent collisional H2 dissociation rate--------
        //
        // idt = 0 calculates coefficients for direct collisional dissociation idt = 1 
        // calculates coefficients for dissociative tunneling as prescribed in 
        // (Martin et al, 1996, ApJ, 461, 265). Code is ran with both values of idt
        // so that coefficients for both cases are calculated.
        // 
        // Results are stored within the k13dd array in chemistry_data_storage. The coefficients
        // corresponding to idt=0 come first, followed by those corresponding to idt=1. There are
        // seven coefficients for each case, and each coefficient has a value at each temperature bin.
        // Therefore the array is structured as follows:
        // k13dd = {coeff1(idt=0, Tbin1), coeff1(idt=0, Tbin2), ..., coeff1(idt=0, TbinFinal), coeff2(idt=0, Tbin1), ..., 
        //          coeff7(idt=0, TbinFinal), coeff1(idt=1, Tbin1), ..., coeff7(idt=1, TbinFinal)}
        int idt = 0;
        colh2diss_g_c(i, binTemp, idt, userChemistry, userRates);
        userRates->k13dd[i] = userRates->k13dd[i] - log10(kUnit); //Normalization
        idt = 1;
        colh2diss_g_c(i, binTemp, idt, userChemistry, userRates);
        userRates->k13dd[i + 7*userChemistry->NumberOfTemperatureBins] = \
            userRates->k13dd[i + 7*userChemistry->NumberOfTemperatureBins] - log10(kUnit); //Normalization


        //--------Calculate 3-body H2 rate--------

        // Calculated by the same method as done in the original code. First is the fit to 
        // A.E. Orel 1987, J.Chem.Phys., 87, 314, which is matched to the 1/T of 
        // Palla etal (1983) -- which is four times smaller than the Palla rate.
        
        //Varying threebody and corresponding collisional dissociation rates from Simon.
        if (userChemistry->three_body_rate == 0) {
            if (binTemp <= 300.0) {
                userRates->k22[i] = 1.3e-32 * pow(binTemp/300.0, -0.38);
            } else {
                userRates->k22[i] = 1.3e-32 * pow(binTemp/300.0, -1.0);
            }
        } else if (userChemistry->three_body_rate == 1) {
            //PSS83 three-body rate
            userRates->k22[i] = 5.5e-29 / binTemp;
            //Inverse of PSS83 three-body rate
            userRates->k13[i] = (5.24e-7 / pow(binTemp, 0.485)) * exp(-5.2e4 / binTemp);
        } else if (userChemistry->three_body_rate == 2) {
            //CW83 three-body rate
            userRates->k22[i] = 8.8e-33;
            //Inverse of CW83 rate
            userRates->k13[i] = 8.4e-11 * pow(binTemp, 0.515) * exp(-5.2e4 / binTemp);
        } else if (userChemistry->three_body_rate == 3) {
            //FH07 three-body rate
            userRates->k22[i] = 1.44e-26 / pow(binTemp, 1.54);
            //Rate from JGC67, used by FH07 to construct their three-body rate
            userRates->k13[i] = (1.38e-4 / pow(binTemp, 1.025)) * exp(-5.2e4 / binTemp);
        } else if (userChemistry->three_body_rate == 4) {
            //Rate from Glover (2008), derived by detailed balance from MSM96
            userRates->k22[i] = 7.7e-31 / pow(binTemp, 0.464);
            //High density limiting rate from MSM96
            userRates->k13[i] = pow(1e1 ,(-178.4239 - 68.42243 * log10(binTemp) \
                                + 43.20243 * pow(log10(binTemp), 2) \
                                - 4.633167 * pow(log10(binTemp), 3) \ 
                                + 69.70086 * log10(1.0 + 40870.38 / binTemp) \
                                - (23705.7 / binTemp)));
        } else if (userChemistry->three_body_rate == 5) {
            //Rate from Forrey (2013)
            userRates->k22[i] = (6e-32 / pow(binTemp, 0.25)) + (2e-31 / pow(binTemp, 0.5));
            //From detailed balance from Forrey rate
            if (binTemp <= 3000.0) {
                userRates->k13[i] = 2.4e-8 * exp(-5.2e4/binTemp);
            } else {
                userRates->k13[i] = 2.2e-6 * pow(binTemp, -0.565) * exp(-5.2e4/binTemp); 
            }
        } else {
            printf('chemistry_data->three_body_rate has been set to an unknown value. Stopping code.');
            exit(0);
        }
        //Normalisation
        userRates->k13[i] = userRates->k13[i] / kUnit;
        userRates->k22[i] = userRates->k22[i] / kUnit_3Bdy;
        userRates->k21[i] = 2.8e-31 * pow(binTemp, -0.6) / kUnit_3Bdy;


        //--------Deuterium Rates--------
        
        //k50
        //Fit taken from Savin (2002) which is valid for T < 2e5 K.
        //We extrapolate for higher temperatures.
        if (binTemp <= 2.0e5) {
            userRates->k50[i] = (2.0e-10 * pow(binTemp, 0.402) * exp(-3.71e1/binTemp) \
                        - 3.31e-17 * pow(binTemp, 1.48)) / kUnit;
        } else {
            userRates->k50[i] = 2.5e-8 * pow(binTemp/2.0e5, 0.402) / kUnit;
        }

        //k51
        // Fit taken from Savin (2002) which is valid for T < 2e5 K.
        userRates->k51[i] = (2.06e-10 * pow(binTemp, 0.396) * exp(-3.30e1/binTemp) \
                            + 2.03e-9 * pow(binTemp, -0.332)) / kUnit;

        //k52 and k53
        // Fits from Galli & Palla (2002) to calculations by Gerlich (1982).
        // If T > 1e4 K use fixed value for k52 to avoid numerical issues with fitting function.
        // In this limit this reaction is not expected to be important anyway.
        if (binTemp <= 1e4) {
            userRates->k52[i] = 1.0e-9 * (0.417 + 0.846 * log10(binTemp) \ 
                                - 0.137 * pow(log10(binTemp), 2)) / kUnit;
        } else {
            userRates->k52[i] = 1.609e-9 / kUnit;
        }   

        userRates->k53[i] = 1.1e-9 * exp(-4.8832/binTemp) / kUnit;

        //k54
        // Fit from Clark et al (2011), which is based on data in Mielke et al (2003).
        if (binTemp <= 2.0e3) {
          userRates->k54[i] = pow(1.0e1, (-5.64737e1 + 5.88886 * log10(binTemp) \
                              + 7.19692  * pow(log10(binTemp), 2) \
                              + 2.25069  * pow(log10(binTemp), 3) \
                              - 2.16903  * pow(log10(binTemp), 4) \
                              + 3.17887e-1 * pow(log10(binTemp), 5)));
        } else {
          userRates->k54[i] = 3.17e-10 * exp(-5.207e3 / binTemp);
        }
        
        //k55
        // Fit from Galli & Palla (2002), which is based on Shavitt (1959).
        // Fit has been modified at low temperature to avoid creating an 
        // anomalously large rate coefficient -- as suggested by Ripamonti (2007)
        // and McGreer & Bryan (2008).
        if (binTemp <= 2.0e2) {
            userRates->k55[i] = 1.08e-22 / kUnit;
        } else {
            userRates->k55[i] = 5.25e-11 * exp(-4.43e3/binTemp + 1.739e5/pow(binTemp, 2)) / kUnit;
        }

        //k56
        // This is the same as HM + HI. 
        // Measurements from Miller et al (2012) suggest there is no significant isotope effect
        // for this reaction.
        userRates->k56[i] = userRates->k8[i];


        //--------New H Ionization Rates--------

        //k57 and k58
        // These rate coefficients are from Lenzuni, Chernoff & Salpeter (1991).
        // k57 value based on experimental cross-sections from Gealy & van Zyl (1987).
        // k58 value based on cross-sections from van Zyl, Le & Amme (1981).
        if (binTemp > 3.0e3) {
            userRates->k57[i] = 1.2e-17  * pow(binTemp, 1.2) * exp(-1.578e5 / binTemp) / kUnit;
            userRates->k58[i] = 1.75e-17 * pow(binTemp, 1.3) * exp(-1.578e5 / binTemp) / kUnit;
        } else {
            userRates->k57[i] = tiny;
            userRates->k58[i] = tiny;
        }

        //H2 formation on dust grains requires loop over the dust temperature.
        for (int j = 0; j < userChemistry->NumberOfDustTemperatureBins; j++) {
            //Compute dust temperature for this bin.
            double dust_logBinTemp = dust_logTempStart + abs(j)*dust_dlogTemp;
            double dust_binTemp = exp(dust_logBinTemp);
            double dust_binTemp2 = dust_binTemp / 1.0e2;

            //h2dust
            // The method used to calculate this is dependent upon what the user has selected
            // By default Omukai 2000 will be used.
            if (userChemistry->useOmukai2000 == 1) {
                //k23 from Omukai (2000).
                userRates->h2dust[i + userChemistry->NumberOfTemperatureBins*j] = 6.0e-17 / fgr * pow(binTemp / 300.0, 0.5) * \
                                                        (pow(1.0 + exp(7.5e2 * ((1.0 / 75.0) - (1.0 / dust_binTemp))), -1.0)) * \
                                                        (pow(1.0 + (4.0e-2 * pow(binTemp + dust_binTemp, 0.5)) + 
                                                        (2.0e-3 *binTemp) + (8.0e-6 * pow(binTemp, 2.0)), -1.0)) / kUnit;
            } else {
                //Equation 3.8 from Hollenbach & McKee (1979).
                userRates->h2dust[i + userChemistry->NumberOfTemperatureBins*j] = 3.0e-17 / fgr * pow(binTemp2, 0.5) / \
                                                                    (1.0 + 0.4 * pow(binTemp2 + dust_binTemp2, 0.5) + \
                                                                    0.2 * dust_binTemp2 + 8.0e-2 * pow(dust_binTemp2, 2.0)) / kUnit;
            }
            
        } // End of loop over dust temperature bins.

        //H2 formation heating terms from Equation 23, Omuaki (2000).
        userRates->n_cr_n[i]  = 1.0e6 * pow(binTemp, -0.5);
        userRates->n_cr_d1[i] = 1.6 * exp(-pow(400.0 / binTemp, 2.0));
        userRates->n_cr_d2[i] = 1.4 * exp(-12000.0 / (binTemp + 1200.0));
        } //End of loop over temperature bins.

        //* 2) Cooling and heating Rates (excluding external radiation field)
        for (int i = 0; i < userChemistry->NumberOfTemperatureBins; i++) {
            //Bin temperature.
            double logBinTemp = log(TempStart) + abs(i)*dlogTemp;
            double binTemp = exp(logBinTemp);

            //* a) Collisional excitations (Black 1981; Cen 1992).
            if (userChemistry->crg_coolExi == 1) {
                userRates->ceHI[i] = 7.5e-19*exp(-fmin(log(dhuge), 118348.0/binTemp)) \
                                    / ( 1.0 + sqrt(binTemp/1.0e5) ) / coolingUnits;
                userRates->ceHeI[i] = 9.1e-27*exp(-fmin(log(dhuge), 13179.0/binTemp)) \
                                    * pow(binTemp, -0.1687) / ( 1.0 + sqrt(binTemp/1.0e5) ) \
                                    / coolingUnits;
                userRates->ceHeII[i] = 5.54e-17*exp(-fmin(log(dhuge), 473638.0/binTemp)) \
                                    * pow(binTemp, -0.3970) / ( 1.0 + sqrt(binTemp/1.0e5) ) \
                                    / coolingUnits;
            } else {
                userRates->ceHI[i] = tiny;
                userRates->ceHeI[i] = tiny;
                userRates->ceHeII[i] = tiny;
            }

            //* b) Collisional ionizations (Cen 1992, Abel 1996)
            if (userChemistry->crg_collIon == 1) {
                userRates->ciHeIS[i] = 5.01e-27*pow(binTemp, -0.1687) / ( 1.0 + sqrt(binTemp/1.0e5) ) \
                                    * exp(-fmin(log(dhuge), 55338.0/binTemp)) / coolingUnits;
                
                //Collisional ionizations. Polynomial fits from Tom Abel.
                userRates->ciHI[i] = 2.18e-11 * userRates->k1[i] * kUnit / coolingUnits;
                userRates->ciHeI[i] = 3.94e-11 * userRates->k3[i] * kUnit / coolingUnits;
                userRates->ciHeII[i] = 8.72e-11 * userRates->k5[i] * kUnit / coolingUnits; 
            } else {
                userRates->ciHeIS[i] = tiny;
                userRates->ciHI[i] = tiny;
                userRates->ciHeI[i] = tiny;
                userRates->ceHeII[i] = tiny;
            }

            //* c) Recombinations (Hui & Gnedin 1997, except where noted)
            if (userChemistry->crg_recomCool == 1) {
                //Define parameters used in the calculations.
                double lambdaHI    = 2.0 * 157807.0 / binTemp;
                double lambdaHeII  = 2.0 * 285335.0 / binTemp;
                double lambdaHeIII = 2.0 * 631515.0 / binTemp;

                //These depend on if the user has chosen recombination case A or B.
                if (userChemistry->CaseBRecombination == 1) {
                    userRates->reHII[i] = 3.435e-30 * binTemp * pow(lambdaHI, 1.970) \
                                        / pow( 1.0 + pow(lambdaHI/2.25, 0.376), 3.720) \
                                        / coolingUnits;
                    userRates->reHeII1[i] = 1.26e-14 * kboltz * binTemp * pow(lambdaHeII, 0.75) \
                                        / coolingUnits;
                    userRates->reHeIII[i] = 8.0 * 3.435e-30 * binTemp * pow(lambdaHeIII, 1.970) \
                                        / pow(1.0 + pow(lambdaHeIII/2.25, 0.376), 3.720) 
                                        / coolingUnits;
                } else {
                    userRates->reHII[i] = 1.778e-29 * binTemp * pow(lambdaHI, 1.965) \
                                        / pow(1.0 + pow(lambdaHI/0.541, 0.502), 2.697) \
                                        / coolingUnits; 
                    userRates->reHeII1[i] = 3e-14 * kboltz * binTemp * pow(lambdaHeII, 0.654) \ 
                                        / coolingUnits;
                    userRates->reHeIII[i] = 8.0 * 1.778e-29 * binTemp * pow(lambdaHeIII, 1.965) \
                                        / pow(1.0 + pow(lambdaHeIII/0.541, 0.502), 2.697) \ 
                                        / coolingUnits;
                }

                //Dielectronic recombination (Cen, 1992).
                userRates->reHeII2[i] = 1.24e-13 * pow(binTemp, -1.5) \
                                    * exp(-fmin(log(dhuge), 470000.0/binTemp)) \
                                    * ( 1.0 + 0.3*exp(-fmin(log(dhuge), 94000.0/binTemp))) \
                                    / coolingUnits;
            } else {
                userRates->reHII[i]   = tiny;
                userRates->reHeII1[i] = tiny;
                userRates->reHeII2[i] = tiny;
                userRates->reHeIII[i] = tiny;
            }

            //* d) Bremsstrahlung (Black, 1981 and Spitzer & Hart, 1979)
            if ( userChemistry->crg_bremCool == 1 ) {
                userRates->brem[i] = 1.43e-27 * sqrt(binTemp) \
                            * (1.1 + 0.34 * exp(-pow(5.5 - log10(binTemp), 2) / 3.0)) \ 
                            / coolingUnits;
            } else {
                userRates->brem[i] = tiny;
            }

            //* e) Molecular hydrogen cooling

            //------Calculate Lepp & Shull rates------
            double par_x = log10(binTemp/1.0e4); //Parameter used in the following calculations.
            double par_dum; //Dummy parameter used in the following calculations.
            if (binTemp > 1635.0) {
                par_dum = 1.0e-12 * sqrt(binTemp) * exp(-1000.0/binTemp);
            } else {
                par_dum = 1.4e-13 * exp((binTemp/125.0) - pow(binTemp/577.0, 2));
            }

            //First couple of rates.
            userRates->vibh[i] = 1.1e-18 * exp(-fmin(log(dhuge),6744.0/binTemp)) / coolingUnits;
            userRates->hyd01k[i] = par_dum * exp(-fmin(log(dhuge), 8.152e-13/(kboltz*binTemp))) \
                                / coolingUnits;

            //Update value of dummy parameter.
            par_dum = 8.152e-13 * ( 4.2 / (kboltz * (binTemp + 1190.0)) + 1.0 / (kboltz * binTemp));


            //Last few rates.
            userRates->h2k01[i] = 1.45e-12 * sqrt(binTemp) * exp(-fmin(log(dhuge), par_dum)) / coolingUnits;

            if (binTemp > 4031.0) {
            userRates->rotl[i] = 1.38e-22 * exp(-9243.0/binTemp) / coolingUnits;
            } else {
            userRates->rotl[i] = pow(10.0, -22.9 - 0.553*par_x - 1.148*pow(par_x, 2)) / coolingUnits;
            }

            if(binTemp > 1087.0) {
            userRates->roth[i] = 3.9e-19 * exp(-6118.0/binTemp) / coolingUnits;
            } else {
            userRates->roth[i] = pow(10.0, -19.24 + 0.474*par_x - 1.247*pow(par_x, 2)) / coolingUnits;
            }

            //------Calculate Galli & Palla (1999) Rates------ 
            //             (as fit by Tom Abel)

            //Constrain temperature.
            double tm;
            tm = fmax(binTemp, 13.0); //no cooling below 13 Kelvin
            tm = fmin(tm, 1.0e5); //fixes numerics
            double lt = log10(tm);

            //Low density limit (Galli & Palla).
            userRates->GP99LowDensityLimit[i] = pow(10.0, -103.0 + 97.59*lt - 48.05*pow(lt, 2) + 10.8*pow(lt, 3) \
                                                - 0.9032*pow(lt, 4)) / coolingUnits;

            //High density limit from HM79.
            double t3 = tm/1000.0;
            double HDLR = (9.5e-22*pow(t3, 3.76)) / (1.0 + 0.12*pow(t3, 2.1)) * \
                        exp(-pow((0.13/t3), 3)) + 3.0e-24 * exp(-0.51/t3);
            double HDLV = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
            userRates->GP99HighDensityLimit[i] = (HDLR + HDLV) / coolingUnits;

            //------Low density rates (Glover & Abel, 2008)------
            
            //Excitation by HI.
            //Constrain temperature.
            tm  = fmax(binTemp, 10.0);
            tm  = fmin(tm, 1.0e4);
            double lt3 = log10(tm / 1.0e3);

            //Calculate rate using method specified in chemistry_data struct.
            if (userChemistry->useLique2015 == 1) {
                if (tm < 1e2) {
                    userRates->GAHI[i] = 0.0;
                } else {
                    userRates->GAHI[i] = pow(10.0, -24.07950609 \
                                    + 4.54182810 * lt3 \
                                    - 2.40206896 * pow(lt3, 2) \
                                    - 0.75355292 * pow(lt3, 3) \
                                    + 4.69258178 * pow(lt3, 4) \
                                    - 2.79573574 * pow(lt3, 5) \
                                    - 3.14766075 * pow(lt3, 6) \
                                    + 2.50751333 * pow(lt3, 7)) / coolingUnits;
                }
            } else {
                if (tm < 1.0e2) {
                    userRates->GAHI[i] = pow(10.0, -16.818342 \
                                    + 37.383713 * lt3 \
                                    + 58.145166 * pow(lt3, 2) \
                                    + 48.656103 * pow(lt3, 3) \
                                    + 20.159831 * pow(lt3, 4) \
                                    + 3.8479610 * pow(lt3, 5)) / coolingUnits;
                } else if (tm < 1.0e3) {
                    userRates->GAHI[i] = pow(10.0, -24.311209 \
                                    + 3.5692468 * lt3 \
                                    - 11.332860 * pow(lt3, 2) \
                                    - 27.850082 * pow(lt3, 3) \
                                    - 21.328264 * pow(lt3, 4) \
                                    - 4.2519023 * pow(lt3, 5)) / coolingUnits;
                } else {
                    userRates->GAHI[i] = pow(10.0, -24.311209 \
                                    + 4.6450521 * lt3 \
                                    - 3.7209846 * pow(lt3, 2) \
                                    + 5.9369081 * pow(lt3, 3) \
                                    - 5.5108047 * pow(lt3, 4) \
                                    + 1.5538288 * pow(lt3, 5)) / coolingUnits;
                }
            }

            //Excitation by H2.
            userRates->GAH2[i] = pow(10.0, -23.962112 \
                                + 2.09433740  * lt3 \
                                - 0.77151436 * pow(lt3, 2) \
                                + 0.43693353 * pow(lt3, 3) \
                                - 0.14913216 * pow(lt3, 4) \
                                - 0.033638326 * pow(lt3, 5)) / coolingUnits;

            //Excitation by He.
            userRates->GAHe[i] = pow(10.0, -23.689237 \
                            + 2.1892372  * lt3 \
                            - 0.81520438 * pow(lt3, 2) \
                            + 0.29036281 * pow(lt3, 3) \
                            - 0.16596184 * pow(lt3, 4) \
                            + 0.19191375 * pow(lt3, 5)) / coolingUnits;

            //Excitation by HII.
            //Collisional excitation rates from Honvault et al (2011, Phys. Rev. Lett., 107, 023201) and 
            //Honvault et al (2012, Phys. Rev. Lett., 108, 109903).
            userRates->GAHp[i] = pow(10.0, -22.089523 \
                            + 1.5714711   * lt3 \
                            + 0.015391166 * pow(lt3, 2) \
                            - 0.23619985  * pow(lt3, 3) \
                            - 0.51002221  * pow(lt3, 4) \
                            + 0.32168730  * pow(lt3, 5)) / coolingUnits;

            //Excitation by electrons.
            //Based on data from  Yoon et al (2008, J. Phys. Chem. Ref. Data, 37, 913).
            if (tm < 100.0) {
                userRates->GAel[i] = 0.0;
            } else if (tm < 500.0) {
                userRates->GAel[i] = pow(10.0, -21.928796 \
                                + 16.815730 * lt3 \
                                + 96.743155 * pow(lt3, 2) \
                                + 343.19180 * pow(lt3, 3) \
                                + 734.71651 * pow(lt3, 4) \
                                + 983.67576 * pow(lt3, 5) \
                                + 801.81247 * pow(lt3, 6) \
                                + 364.14446 * pow(lt3, 7) \
                                + 70.609154 * pow(lt3, 8)) / coolingUnits;
            } else {
                userRates->GAel[i] = pow(10.0, -22.921189 \
                                + 1.6802758  * lt3 \
                                + 0.93310622 * pow(lt3, 2) \
                                + 4.0406627  * pow(lt3, 3) \
                                - 4.7274036  * pow(lt3, 4) \
                                - 8.8077017  * pow(lt3, 5) \
                                + 8.9167183  * pow(lt3, 6) \
                                + 6.4380698  * pow(lt3, 7) \
                                - 6.3701156  * pow(lt3, 8)) / coolingUnits;
            }

            //------New fit to LTE rate------ 
            //from Glover (2015, MNRAS, 451, 2082)
            if (tm < 1.0e2) {
                //Simple extrapolation as H2 cooling insignificant at these temperatures.
                userRates->H2LTE[i] = 7.0e-27 * pow(tm, 1.5) * exp(-512.0/tm) / coolingUnits;
            } else {
                userRates->H2LTE[i] = pow(10.0, -20.584225 \
                                + 5.0194035 * lt3 \
                                - 1.5738805 * pow(lt3, 2) \
                                - 4.7155769 * pow(lt3, 3) \
                                + 2.4714161 * pow(lt3, 4) \
                                + 5.4710750 * pow(lt3, 5) \
                                - 3.9467356 * pow(lt3, 6) \
                                - 2.2148338 * pow(lt3, 7) \
                                + 1.8161874 * pow(lt3, 8)) / coolingUnits;
            }

            //* f) HD cooling.
            //HD cooling function has units of ergs cm^3 / s

            //Fit from Coppola et al 2011. LTE (ergs/s) -> hdlte (ergs cm3/s)
            //Constain temperature.
            tm  = fmax(binTemp, 10.0);
            tm  = fmin(tm, 3.0e4);
            //Calculate rate.
            userRates->HDlte[i] = -55.5725 + 56.649 * log10(tm) \
                            - 37.9102  * pow(log10(tm), 2) \
                            + 12.698   * pow(log10(tm), 3) \
                            - 2.02424  * pow(log10(tm), 4) \
                            + 0.122393 * pow(log10(tm), 5);
            userRates->HDlte[i] = pow(10.0, fmin(userRates->HDlte[i], 0.0)) / coolingUnits;

            //Rate based on (Wrathmall, Gusdorf & Flower, 2007) HD-H collisional excitation rates.
            //Constrain temperature.
            tm = fmax(binTemp, 1.0e1);
            tm = fmin(tm, 6.0e3);
            //Calculate rate.
            userRates->HDlow[i] = -23.175780  + 1.5035261 * log10(tm/1.0e3) \
                            + 0.40871403  * pow(log10(tm/1.0e3), 2) \
                            + 0.17849311  * pow(log10(tm/1.0e3), 3) \
                            - 0.077291388 * pow(log10(tm/1.0e3), 4) \
                            + 0.10031326 * pow(log10(tm/1.0e3), 5);
            userRates->HDlow[i] = pow(10.0, userRates->HDlow[i]) / coolingUnits;

            //* g) CIE cooling (Ripamonti & Abel, 2003).

            //Get cooling rate.
            double cie_rate = cie_thin_cooling_rate_g_c(binTemp);
            
            //Below explanation originally written by Matt Turk:
            //  The above function returns the rate with units of ergs * s^-1 * cm^3 * gram^-1 *
            //  (gram in H2 molecules)^-1. To reproduce equation 5 in RA04, divide c_t_c_r by mh,
            // so to get erg/s cm^3 we multiply by mh. Divide by 2 for the H2 mass --> number.
            userRates->cieco[i] = cie_rate * (mh/2.0) / coolingUnits;

        } //End of loop over temperature bins.
    } //End of primordial_chemistry if statement.

    //* Dust-related Processes.
    if (anyDust == TRUE) {
        //Loop over temperature bins.
        for (int i = 0; i < userChemistry->NumberOfTemperatureBins; i++) {
            //Calculate bin temperature (eV).
            double logBinTemp = logTempStart + abs(i)*dlogTemp;
            double binTemp = exp(logBinTemp);

            //Energy transfer from gas to dust grains.
            //(Equation 2.15, Hollenbach & McKee, 1989)
            userRates->gas_grain[i] = grainCoeff * pow(binTemp, 0.5) * \
                                   (1.0 - 0.8 * exp(-75.0/binTemp)) / coolingUnits;

            //Electron recombination onto dust grains.
            //(Equation 9, Wolfire et al., 1995)

            double grbeta = pow(0.74/binTemp, 0.068); //! Power possibly in the incorrect place.
            userRates->regr[i] = 4.65e-30 * pow(binTemp, 0.94 + 0.5 * grbeta) \
                              / coolingUnits;
        } //End of loop over temperature bins.
    } //End of anyDust if statement.

    //* i) Compton cooling (Peebles 1971).
    userRates->comp = 5.65e-36 / coolingUnits;

    //* j) Photoelectric heating by UV-irradiated dust (Wolfire 1995).
    //Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3).
    if (userChemistry->photoelectric_heating <= 1) {
         userRates->gammah = userChemistry->photoelectric_heating_rate / coolingUnits;
    } else { //photoelectric_heating set to 2 or 3
        //User to specify G_0, epsilon set to 0.05 or calculated directly.
        userRates->gammah = 1.0e-24 / coolingUnits;
    }

    //Heating of dust by interstellar radiation field.
    //(Equation B15, Krumholz, 2014)
    //Don't normalize by coolunit since tdust calculation is done in CGS.
    userRates->gamma_isrf = 3.9e-24 / mh / fgr;
    

    return SUCCESS;
//End of function definition.
}
     