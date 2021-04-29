/*******************************************************************************/
/*
 * ----------Function to Calculate Multispecies Rate Lookup Table---------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 18/2/2021
 * 
 * Purpose:
 *  Calculates lookup tables for the rate coefficients (and some cooling functions) as a function 
 *  of log(temperature).
 * 
 * Inputs:
 *  *my_chemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *my_rates : Pointer to the structure of type chemistry_data_storage
 *               which holds pointers to the rates.
 *  *my_units : Pointer to the structure which defines the units for 
 *               the variables.
 *  *co_length_unit : Pointer to the comoving length unit.
 *  *co_density_unit : Pointer to the comoving density unit.
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
int coll_rates_g_c(int Tbin_index, double T, double kUnit, chemistry_data *my_chemistry,
                    chemistry_data_storage *my_rates);
//Optically thin CIE cooling rate
double cie_thin_cooling_rate_g_c(double T);
//Density-dependent H2  dissociation rate
int colh2diss_g_c(int Tbin_index, double T, int idt, chemistry_data *my_chemistry,
                   chemistry_data_storage *my_rates);

//* Definition of the function.
int calc_rates_g_c(chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
                   code_units *my_units, double co_length_unit, double co_density_unit){

    
    //* Compute rates for primordial chemistry.
    //All rates are calculated if primordial_chemistry > 0 as it is more convenient and there
    //is not much cost to doing it this way as this code only has to run once.
    if (my_chemistry->primordial_chemistry > 0) {
        //Initialise constants to tiny.
        for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
            my_rates->k2[i] = tiny;
            my_rates->k20[i] = tiny;
            my_rates->k21[i] = tiny;
            my_rates->k22[i] = tiny;
            my_rates->k50[i] = tiny;
            my_rates->k51[i] = tiny;
            my_rates->k52[i] = tiny;
            my_rates->k53[i] = tiny;
            my_rates->k54[i] = tiny;
            my_rates->k55[i] = tiny;
            my_rates->k56[i] = tiny;
            my_rates->k57[i] = tiny;
            my_rates->k58[i] = tiny;

            //k13dd and h2dust were originally 2D arrays in the fortran code: I use 1D arrays here with the following mappings,
            //k13dd:
            //  {k13dd(0,0), k13dd(1,0), ... k13dd(NumberOfTemperatureBins, 0), k13dd(0,1), 2D(1,1)..., k13dd(NumberOfTemperatureBins, 13)}
            //h2dust:
            //  {h2dust(0,0), h2dust(1,0), ... h2dust(NumberOfTemperatureBins, 0), h2dust(0,1), 2D(1,1)..., h2dust(NumberOfTemperatureBins, NumberOfDustTemperatureBins)}
            //where k13dd(i,j) and h2dust(i,j) are the elements in the corresponding 2D arrays.
            for (int j = 0; j < 14; j++) {
                my_rates->k13dd[i + my_chemistry->NumberOfTemperatureBins*j] = tiny;
            }
            for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++) {
                my_rates->h2dust[i + my_chemistry->NumberOfTemperatureBins*j] = tiny;
            }

            my_rates->n_cr_n[i]  = tiny;
            my_rates->n_cr_d1[i] = tiny;
            my_rates->n_cr_d2[i] = tiny;

            my_rates->ceHI[i]    = tiny;
            my_rates->ceHeI[i]   = tiny;
            my_rates->ceHeII[i]  = tiny;
            my_rates->ciHI[i]    = tiny;
            my_rates->ciHeI[i]   = tiny;
            my_rates->ciHeIS[i]  = tiny;
            my_rates->ciHeII[i]  = tiny;
            my_rates->reHII[i]   = tiny;
            my_rates->reHeII1[i] = tiny;
            my_rates->reHeII2[i] = tiny;
            my_rates->reHeIII[i] = tiny;
            my_rates->brem[i]    = tiny;
            my_rates->comp       = tiny;

            my_rates->hyd01k[i]  = tiny;
            my_rates->h2k01[i]   = tiny;
            my_rates->vibh[i]    = tiny;
            my_rates->roth[i]    = tiny;
            my_rates->rotl[i]    = tiny;
            my_rates->HDlte[i]   = tiny;
            my_rates->HDlow[i]   = tiny;
            my_rates->cieco[i]   = tiny;

            my_rates->GAHI[i]    = tiny;
            my_rates->GAH2[i]    = tiny;
            my_rates->GAHe[i]    = tiny;
            my_rates->GAHp[i]    = tiny;
            my_rates->GAel[i]    = tiny;
            my_rates->H2LTE[i]   = tiny;
        } //End of parameter initialisation for loop.

    //* Populate tables for values between the specificed start and end temperatures.
    //* _____________________________________________________________________________

    // Declare all variables used within loops.
    double logBinTemp, binTemp, binTemp2, binTemp_eV, logBinTemp_eV, dust_logBinTemp, dust_binTemp, dust_binTemp2,
           lambdaHI, lambdaHeII, lambdaHeIII, par_x, par_dum, tm, lt, t3, HDLR, HDLV, lt3, cierate;
    
    //* 1) Rate Coefficients (excluding the external radiation field)
    //Loop over temperature bins.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
        //Calculate convenient temperature parameters for this bin (in eV).
        logBinTemp = logTempStart + i*dlogTemp;
        binTemp = exp(logBinTemp);
        binTemp2 = binTemp / 1.0e2;
        binTemp_eV = binTemp / tevk;
        logBinTemp_eV = log(binTemp_eV);

        

        //H2 formation on dust grains requires loop over the dust temperature.
        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++) {
            //Compute dust temperature for this bin.
            dust_logBinTemp = dust_logTempStart + j*dust_dlogTemp;
            dust_binTemp = exp(dust_logBinTemp);
            dust_binTemp2 = dust_binTemp / 1.0e2;

            //h2dust
            // The method used to calculate this is dependent upon what the user has selected
            // By default Omukai 2000 will be used.
            if (my_chemistry->h2dust_rate == 1) {
                //k23 from Omukai (2000).
                my_rates->h2dust[i + my_chemistry->NumberOfTemperatureBins*j] = 6.0e-17 / fgr * pow(binTemp / 300.0, 0.5) *
                                                        (pow(1.0 + exp(7.5e2 * ((1.0 / 75.0) - (1.0 / dust_binTemp))), -1.0)) *
                                                        (pow(1.0 + (4.0e-2 * pow(binTemp + dust_binTemp, 0.5)) + 
                                                        (2.0e-3 * binTemp) + (8.0e-6 * pow(binTemp, 2.0)), -1.0)) / kUnit;
            } else {
                //Equation 3.8 from Hollenbach & McKee (1979).
                my_rates->h2dust[i + my_chemistry->NumberOfTemperatureBins*j] = 3.0e-17 / fgr * pow(binTemp2, 0.5) /
                                                                    (1.0 + 0.4 * pow(binTemp2 + dust_binTemp2, 0.5) +
                                                                    0.2 * dust_binTemp2 + 8.0e-2 * pow(dust_binTemp2, 2.0)) / kUnit;
            }
        } // End of loop over dust temperature bins.

        } //End of loop over temperature bins.





        for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
            


            //* e) Molecular hydrogen cooling

           
         

            //------Calculate Galli & Palla (1999) Rates------ 
            //             (as fit by Tom Abel)

            //Constrain temperature.
            tm = fmax(binTemp, 13.0); //no cooling below 13 Kelvin
            tm = fmin(tm, 1.0e5); //fixes numerics
            lt = log10(tm);

            //Low density limit (Galli & Palla).
            my_rates->GP99LowDensityLimit[i] = pow(10.0, -103.0 + 97.59*lt - 48.05*pow(lt, 2) + 10.8*pow(lt, 3)
                                                - 0.9032*pow(lt, 4)) / coolingUnits;

            //High density limit from HM79.
            t3 = tm/1000.0;
            HDLR = (9.5e-22*pow(t3, 3.76)) / (1.0 + 0.12*pow(t3, 2.1)) *
                        exp(-pow((0.13/t3), 3)) + 3.0e-24 * exp(-0.51/t3);
            HDLV = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
            my_rates->GP99HighDensityLimit[i] = (HDLR + HDLV) / coolingUnits;

            //------Low density rates (Glover & Abel, 2008)------
            
            //Excitation by HI.
            //Constrain temperature.
            tm  = fmax(binTemp, 10.0);
            tm  = fmin(tm, 1.0e4);
            lt3 = log10(tm / 1.0e3);

            //Calculate rate using method specified in chemistry_data struct.
            if (my_chemistry->h2_h_cooling_rate == 1) {
                if (tm < 1e2) {
                    my_rates->GAHI[i] = 0.0;
                } else {
                    my_rates->GAHI[i] = pow(10.0, -24.07950609
                                    + 4.54182810 * lt3
                                    - 2.40206896 * pow(lt3, 2)
                                    - 0.75355292 * pow(lt3, 3)
                                    + 4.69258178 * pow(lt3, 4)
                                    - 2.79573574 * pow(lt3, 5)
                                    - 3.14766075 * pow(lt3, 6)
                                    + 2.50751333 * pow(lt3, 7)) / coolingUnits;
                }
            } else {
                if (tm < 1.0e2) {
                    my_rates->GAHI[i] = pow(10.0, -16.818342
                                    + 37.383713 * lt3
                                    + 58.145166 * pow(lt3, 2)
                                    + 48.656103 * pow(lt3, 3)
                                    + 20.159831 * pow(lt3, 4)
                                    + 3.8479610 * pow(lt3, 5)) / coolingUnits;
                } else if (tm < 1.0e3) {
                    my_rates->GAHI[i] = pow(10.0, -24.311209
                                    + 3.5692468 * lt3
                                    - 11.332860 * pow(lt3, 2)
                                    - 27.850082 * pow(lt3, 3)
                                    - 21.328264 * pow(lt3, 4)
                                    - 4.2519023 * pow(lt3, 5)) / coolingUnits;
                } else {
                    my_rates->GAHI[i] = pow(10.0, -24.311209
                                    + 4.6450521 * lt3
                                    - 3.7209846 * pow(lt3, 2)
                                    + 5.9369081 * pow(lt3, 3)
                                    - 5.5108047 * pow(lt3, 4)
                                    + 1.5538288 * pow(lt3, 5)) / coolingUnits;
                }
            }

            //Excitation by H2.
            my_rates->GAH2[i] = pow(10.0, -23.962112
                                + 2.09433740  * lt3
                                - 0.77151436 * pow(lt3, 2)
                                + 0.43693353 * pow(lt3, 3)
                                - 0.14913216 * pow(lt3, 4)
                                - 0.033638326 * pow(lt3, 5)) / coolingUnits;

            //Excitation by He.
            my_rates->GAHe[i] = pow(10.0, -23.689237
                            + 2.1892372  * lt3
                            - 0.81520438 * pow(lt3, 2)
                            + 0.29036281 * pow(lt3, 3)
                            - 0.16596184 * pow(lt3, 4)
                            + 0.19191375 * pow(lt3, 5)) / coolingUnits;

            //Excitation by HII.
            //Collisional excitation rates from Honvault et al (2011, Phys. Rev. Lett., 107, 023201) and 
            //Honvault et al (2012, Phys. Rev. Lett., 108, 109903).
            my_rates->GAHp[i] = pow(10.0, -22.089523
                            + 1.5714711   * lt3
                            + 0.015391166 * pow(lt3, 2)
                            - 0.23619985  * pow(lt3, 3)
                            - 0.51002221  * pow(lt3, 4)
                            + 0.32168730  * pow(lt3, 5)) / coolingUnits;

            //Excitation by electrons.
            //Based on data from  Yoon et al (2008, J. Phys. Chem. Ref. Data, 37, 913).
            if (tm < 100.0) {
                my_rates->GAel[i] = 0.0;
            } else if (tm < 500.0) {
                my_rates->GAel[i] = pow(10.0, -21.928796
                                + 16.815730 * lt3
                                + 96.743155 * pow(lt3, 2)
                                + 343.19180 * pow(lt3, 3)
                                + 734.71651 * pow(lt3, 4)
                                + 983.67576 * pow(lt3, 5)
                                + 801.81247 * pow(lt3, 6)
                                + 364.14446 * pow(lt3, 7)
                                + 70.609154 * pow(lt3, 8)) / coolingUnits;
            } else {
                my_rates->GAel[i] = pow(10.0, -22.921189
                                + 1.6802758  * lt3
                                + 0.93310622 * pow(lt3, 2)
                                + 4.0406627  * pow(lt3, 3)
                                - 4.7274036  * pow(lt3, 4)
                                - 8.8077017  * pow(lt3, 5)
                                + 8.9167183  * pow(lt3, 6)
                                + 6.4380698  * pow(lt3, 7)
                                - 6.3701156  * pow(lt3, 8)) / coolingUnits;
            }

            //------New fit to LTE rate------ 
            //from Glover (2015, MNRAS, 451, 2082)
            if (tm < 1.0e2) {
                //Simple extrapolation as H2 cooling insignificant at these temperatures.
                my_rates->H2LTE[i] = 7.0e-27 * pow(tm, 1.5) * exp(-512.0/tm) / coolingUnits;
            } else {
                my_rates->H2LTE[i] = pow(10.0, -20.584225
                                + 5.0194035 * lt3
                                - 1.5738805 * pow(lt3, 2)
                                - 4.7155769 * pow(lt3, 3)
                                + 2.4714161 * pow(lt3, 4)
                                + 5.4710750 * pow(lt3, 5)
                                - 3.9467356 * pow(lt3, 6)
                                - 2.2148338 * pow(lt3, 7)
                                + 1.8161874 * pow(lt3, 8)) / coolingUnits;
            }




            //* f) HD cooling.
            //HD cooling function has units of ergs cm^3 / s

            //Fit from Coppola et al 2011. LTE (ergs/s) -> hdlte (ergs cm3/s)
            //Constain temperature.
            tm  = fmax(binTemp, 10.0);
            tm  = fmin(tm, 3.0e4);
            //Calculate rate.
            my_rates->HDlte[i] = -55.5725 + 56.649 * log10(tm)
                            - 37.9102  * pow(log10(tm), 2)
                            + 12.698   * pow(log10(tm), 3)
                            - 2.02424  * pow(log10(tm), 4)
                            + 0.122393 * pow(log10(tm), 5);
            my_rates->HDlte[i] = pow(10.0, fmin(my_rates->HDlte[i], 0.0)) / coolingUnits;

            //Rate based on (Wrathmall, Gusdorf & Flower, 2007) HD-H collisional excitation rates.
            //Constrain temperature.
            tm = fmax(binTemp, 1.0e1);
            tm = fmin(tm, 6.0e3);
            //Calculate rate.
            my_rates->HDlow[i] = -23.175780  + 1.5035261 * log10(tm/1.0e3)
                            + 0.40871403  * pow(log10(tm/1.0e3), 2)
                            + 0.17849311  * pow(log10(tm/1.0e3), 3)
                            - 0.077291388 * pow(log10(tm/1.0e3), 4)
                            + 0.10031326 * pow(log10(tm/1.0e3), 5);
            my_rates->HDlow[i] = pow(10.0, my_rates->HDlow[i]) / coolingUnits;

            //* g) CIE cooling (Ripamonti & Abel, 2003).
            
            //Below explanation originally written by Matt Turk:
            //  The above function returns the rate with units of ergs * s^-1 * cm^3 * gram^-1 *
            //  (gram in H2 molecules)^-1. To reproduce equation 5 in RA04, divide c_t_c_r by mh,
            // so to get erg/s cm^3 we multiply by mh. Divide by 2 for the H2 mass --> number.
            cierate = cie_thin_cooling_rate_g_c(binTemp);
            my_rates->cieco[i] =  cierate * (mh/2.0) / coolingUnits;

        } //End of loop over temperature bins.
    } //End of primordial_chemistry if statement.

    //* Dust-related Processes.
    if (anyDust == TRUE) {
        //Loop over temperature bins.
        for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
            //Calculate bin temperature (eV).
            double logBinTemp = logTempStart + abs(i)*dlogTemp;
            double binTemp = exp(logBinTemp);

            //Energy transfer from gas to dust grains.
            //(Equation 2.15, Hollenbach & McKee, 1989)
            my_rates->gas_grain[i] = grainCoeff * pow(binTemp, 0.5) *
                                   (1.0 - 0.8 * exp(-75.0/binTemp)) / coolingUnits;

            //Electron recombination onto dust grains.
            //(Equation 9, Wolfire et al., 1995)

            double grbeta = 0.74 / pow(binTemp, 0.068);
            my_rates->regr[i] = 4.65e-30 * pow(binTemp, 0.94 + 0.5 * grbeta)
                              / coolingUnits;
        } //End of loop over temperature bins.
    } //End of anyDust if statement.


    //* i) Compton cooling (Peebles 1971).
    my_rates->comp = 5.65e-36 / coolingUnits;

    //* j) Photoelectric heating by UV-irradiated dust (Wolfire 1995).
    //Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3).
    if (my_chemistry->photoelectric_heating <= 1) {
         my_rates->gammah = my_chemistry->photoelectric_heating_rate / coolingUnits;
    } else { //photoelectric_heating set to 2 or 3
        //User to specify G_0, epsilon set to 0.05 or calculated directly.
        my_rates->gammah = 1.0e-24 / coolingUnits;
    }

    //Heating of dust by interstellar radiation field.
    //(Equation B15, Krumholz, 2014)
    //Don't normalize by coolunit since tdust calculation is done in CGS.
    my_rates->gamma_isrf = 3.9e-24 / mh / fgr;

    /* //! Code that saves temperatures of all bins. Used only for testing. 
    FILE *fp;
    fp = fopen("temperatures.txt", "w");
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
        //Bin temperature.
        logBinTemp = log(TempStart) + abs(i)*dlogTemp;
        binTemp = exp(logBinTemp);
        fprintf(fp, "%E \n", binTemp);
    }
    fclose(fp);
    */     

    return SUCCESS;
//End of function definition.
}
     