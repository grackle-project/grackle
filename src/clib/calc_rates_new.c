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
typedef double (*rate_function)(double, double, chemistry_data*);

//Define a function which is able to add a rate to the calculation.
int add_reaction_rate(double **rate_ptr, double logT_start, double d_logT, rate_function my_function, double kUnit, chemistry_data *my_chemistry)
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
        *rate_ptr[i] = my_function(T, kUnit, my_chemistry);
    }
    return SUCCESS;
}

//Definition of the calc_rates function.
int calc_rates(chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
                   code_units *my_units, double co_length_unit, double co_density_unit)
{
    //* Set various constant parameters that will be used throughout.
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
    //All rates are calculated if primordial_chemistry > 0 as it is more convenient and there
    //is not much cost to doing it this way as this code only has to run once.
    if (my_chemistry->primordial_chemistry > 0){
        add_reaction_rate(&my_rates->k1, logT_start, d_logT, k1_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k3, logT_start, d_logT, k3_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k4, logT_start, d_logT, k4_rate, kUnit, my_chemistry);
        //Space for k2.
        add_reaction_rate(&my_rates->k5, logT_start, d_logT, k5_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k6, logT_start, d_logT, k6_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k7, logT_start, d_logT, k7_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k8, logT_start, d_logT, k8_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k9, logT_start, d_logT, k9_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k9, logT_start, d_logT, k9_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k10, logT_start, d_logT, k10_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k11, logT_start, d_logT, k11_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k12, logT_start, d_logT, k12_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k13, logT_start, d_logT, k13_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k14, logT_start, d_logT, k14_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k15, logT_start, d_logT, k15_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k16, logT_start, d_logT, k16_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k17, logT_start, d_logT, k17_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k18, logT_start, d_logT, k18_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k19, logT_start, d_logT, k19_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k23, logT_start, d_logT, k23_rate, kUnit, my_chemistry);
    }
}
