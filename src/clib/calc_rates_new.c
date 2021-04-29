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
#include "rate_functions.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"

//Define the type of a generic rate function.
typedef double (*rate_function)(double, double, double, chemistry_data*);

//Define a function which is able to add a rate to the calculation.
int add_reaction_rate(double **rate_ptr, double logT_start, double d_logT, rate_function my_function, double kUnit, double coolingUnits, chemistry_data *my_chemistry)
{
    //Allocate memory for rate.
    rate_ptr = malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    double T, logT;
    //Calculate rate at each temperature.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++){
        //Set rate to tiny for safety.
        *rate_ptr[i] = tiny;

        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        //Calculate rate and store.
        *rate_ptr[i] = my_function(T, kUnit, coolingUnits, my_chemistry);
    }
    return SUCCESS;
}

//Definition of the calc_rates function.
int calc_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
                   code_units *my_units, double co_length_unit, double co_density_unit)
{
    //* Set various constant parameters that will be used throughout.
    //! Not using these anymore, their values are hardcoded.
    double tevk = 1.1605e4; //Kelvin to eV conversion factor
    double dhuge = 1.0e30; //Comparison value
    
    //* Set the flag for dust calculations.
    int anyDust;
    if ( my_chemistry->h2_on_dust > 0 || my_chemistry->dust_chemistry > 0) {
        anyDust = TRUE;
    } else {
        anyDust = FALSE;
    }

    //* Obtain the conversion factors which convert between code and physical units. We define a_value = 1 at z = zInit such that a = a_value * [a].
    double timeBase1 = my_units->time_units;
    double lengthBase1 = co_length_unit/(my_units->a_value * my_units->a_units);
    double densityBase1 = co_density_unit*pow(my_units->a_value * my_units->a_units, 3);

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
    double kUnit = (pow(my_units->a_units, 3) * mh) / (densityBase1 * timeBase1);
    double kUnit_3Bdy = kUnit * (pow(my_units->a_units, 3) * mh) / densityBase1;

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
    double coolingUnits = (pow(my_units->a_value, 5) * pow(lengthBase1, 2) * pow(mh, 2))
                          / (densityBase1 * pow(timeBase1, 3));

    //* 3) Units for radiative transfer coefficients are 1/[time].

    //* 4) Calculate energy transfer from gas to dust grains (Equation 2.15, Hollenbach & McKee, 1989).
    //*    Normalize to the HM89 dust-to-gas ratio.
    double fgr = 0.009387;
    double grainCoeff = 1.2e-31 * pow(1.0e3, -0.5) / fgr;

    //* Compute log spacing in temperature.
    double logT_start = log(my_chemistry->TemperatureStart);
    double d_logT = (log(my_chemistry->TemperatureEnd) - logT_start)
                           / (my_chemistry->NumberOfTemperatureBins - 1);

    //* Compute log spacing in dust temperature.
    double dust_logT_start = log(my_chemistry->DustTemperatureStart);
    double dust_d_logT = (log(my_chemistry->DustTemperatureEnd) - dust_logT_start)
                                / (my_chemistry->NumberOfDustTemperatureBins - 1);

    

    
    //* Compute rates for primordial chemistry.

    //* 1) Rate Coefficients (excluding the external radiation field)
    if (my_chemistry->primordial_chemistry > 0){ //! Check that I need this.
        //--------Calculate multispecies collissional rates--------
        add_reaction_rate(&my_rates->k1, logT_start, d_logT, k1_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k3, logT_start, d_logT, k3_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k4, logT_start, d_logT, k4_rate, kUnit, coolingUnits, my_chemistry);
        //Space for k2.
        add_reaction_rate(&my_rates->k5, logT_start, d_logT, k5_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k6, logT_start, d_logT, k6_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k7, logT_start, d_logT, k7_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k8, logT_start, d_logT, k8_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k9, logT_start, d_logT, k9_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k9, logT_start, d_logT, k9_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k10, logT_start, d_logT, k10_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k11, logT_start, d_logT, k11_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k12, logT_start, d_logT, k12_rate, kUnit, coolingUnits, my_chemistry);
        //k13 Calculated with 3-body H2 rate.
        add_reaction_rate(&my_rates->k14, logT_start, d_logT, k14_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k15, logT_start, d_logT, k15_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k16, logT_start, d_logT, k16_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k17, logT_start, d_logT, k17_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k18, logT_start, d_logT, k18_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k19, logT_start, d_logT, k19_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k23, logT_start, d_logT, k23_rate, kUnit, coolingUnits, my_chemistry);

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
        //! NEED TO FIGURE OUT HOW TO CONVERT colh2diss
        int idt = 0;
        colh2diss_g_c(i, binTemp, idt, my_chemistry, my_rates);
        my_rates->k13dd[i] -= log10(kUnit); //Normalization
        idt = 1;
        colh2diss_g_c(i, binTemp, idt, my_chemistry, my_rates);
        my_rates->k13dd[i + 7*my_chemistry->NumberOfTemperatureBins] -= log10(kUnit); //Normalization

        //--------Calculate 3-body H2 rate--------

        // Calculated by the same method as done in the original code. First is the fit to 
        // A.E. Orel 1987, J.Chem.Phys., 87, 314, which is matched to the 1/T of 
        // Palla etal (1983) -- which is four times smaller than the Palla rate.
        
        //Varying threebody and corresponding collisional dissociation rates from Simon.
        add_reaction_rate(&my_rates->k13, logT_start, d_logT, k13_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k21, logT_start, d_logT, k21_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k22, logT_start, d_logT, k22_rate, kUnit, coolingUnits, my_chemistry);
        //Normalisation //! Check that this is a sensible way of doing this.
        for (int i=0; i < my_chemistry->NumberOfDustTemperatureBins; i++) {
            my_rates->k13[i] = my_rates->k13[i] / kUnit;
            my_rates->k21[i] =  my_rates->k21[i] / kUnit_3Bdy;
            my_rates->k22[i] = my_rates->k22[i] / kUnit_3Bdy;
        }

        //--------Deuterium Rates--------
        add_reaction_rate(&my_rates->k50, logT_start, d_logT, k50_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k51, logT_start, d_logT, k51_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k52, logT_start, d_logT, k52_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k53, logT_start, d_logT, k53_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k54, logT_start, d_logT, k54_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k55, logT_start, d_logT, k55_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k56, logT_start, d_logT, k56_rate, kUnit, coolingUnits, my_chemistry);

        //--------New H Ionization Rates--------
        add_reaction_rate(&my_rates->k57, logT_start, d_logT, k57_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->k58, logT_start, d_logT, k58_rate, kUnit, coolingUnits, my_chemistry);

        //H2 formation on dust grains requires loop over the dust temperature.
        //! Left this for now.

        //H2 formation heating terms from Equation 23, Omuaki (2000).
        add_reaction_rate(&my_rates->n_cr_n, logT_start, d_logT, n_cr_n_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d1, logT_start, d_logT, n_cr_d1_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d2, logT_start, d_logT, n_cr_d2_rate, kUnit, coolingUnits, my_chemistry);

        //* 2) Cooling and heating Rates (excluding external radiation field)

        //* a) Collisional excitations (Black 1981; Cen 1992).
        add_reaction_rate(&my_rates->ceHI, logT_start, d_logT, ceHI_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeI, logT_start, d_logT, ceHeI_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeII, logT_start, d_logT, ceHeII_rate, kUnit, coolingUnits, my_chemistry);

        //* b) Collisional ionizations (Cen 1992, Abel 1996).
        add_reaction_rate(&my_rates->ciHeIS, logT_start, d_logT, ciHeIS_rate, kUnit, coolingUnits, my_chemistry);
        //Collisional ionizations. Polynomial fits from Tom Abel.
        add_reaction_rate(&my_rates->ciHI, logT_start, d_logT, ciHI_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeI, logT_start, d_logT, ciHeI_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeII, logT_start, d_logT, ciHeII_rate, kUnit, coolingUnits, my_chemistry);

        //* c) Recombinations (Hui & Gnedin 1997, except where noted).
        add_reaction_rate(&my_rates->reHII, logT_start, d_logT, reHII_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeII1, logT_start, d_logT, reHeII1_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeII2, logT_start, d_logT, reHeII2_rate, kUnit, coolingUnits, my_chemistry); //Dielectronic recombination (Cen, 1992).
        add_reaction_rate(&my_rates->reHeIII, logT_start, d_logT, reHeIII_rate, kUnit, coolingUnits, my_chemistry);

        //* d) Bremsstrahlung (Black, 1981 and Spitzer & Hart, 1979)
        add_reaction_rate(&my_rates->brem, logT_start, d_logT, brem_rate, kUnit, coolingUnits, my_chemistry);

        //* e) Molecular hydrogen cooling
        //------Calculate Lepp & Shull rates------
        add_reaction_rate(&my_rates->vibh, logT_start, d_logT, vibh_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->hyd01k, logT_start, d_logT, hyd01k_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->h2k01, logT_start, d_logT, h2k01_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->rotl, logT_start, d_logT, rotl_rate, kUnit, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->roth, logT_start, d_logT, roth_rate, kUnit, coolingUnits, my_chemistry);
         
        //------Calculate Galli & Palla (1999) Rates------ 
        //             (as fit by Tom Abel)
        
    }
}
