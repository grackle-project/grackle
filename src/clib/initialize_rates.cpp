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
 *   @  k2  @     HII + e --> HI + photon
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
 *   @  h2dust @     2H + grain --> H2 + grain
 * 
 * 
 * 
 * This code is a refactoring of the original Grackle fortran function of the same name.
 */
/*******************************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "grackle.h"
#include "grackle_macros.h"
#include "grackle_rate_functions.h"
#include "initialize_dust_yields.hpp"  // initialize_dust_yields
#include "initialize_metal_chemistry_rates.hpp"  // initialize_metal_chemistry_rates
#include "initialize_rates.hpp"
#include "internal_types.hpp" // new_CollisionalRxnRateCollection
#include "LUT.hpp" // CollisionalRxnLUT
#include "opaque_storage.hpp" // gr_opaque_storage
#include "phys_constants.h"

// We define the function pointers inside an extern "C" block because they
// are used to describe the functions declared in grackle_rate_functions.h,
// and all of these functions have C linkage
extern "C" {

//Define the type of a scalar rate function.
typedef double (*scalar_rate_function)(double, chemistry_data*);

//Define the type of a generic rate function.
typedef double (*rate_function)(double, double, chemistry_data*);

} // extern "C"

namespace { // stuff inside an anonymous namespace are only locally visible

//Define a function which calculates start temperature and increments in log space.
void logT_spacing(double *logT_start, double *d_logT, chemistry_data *my_chemistry)
{
    *logT_start = log(my_chemistry->TemperatureStart);
    *d_logT = (log(my_chemistry->TemperatureEnd) - *logT_start)
                / (my_chemistry->NumberOfTemperatureBins - 1);
}

//Define a function which calculates start temperature and incremenets in log space for dust.
void logT_spacing_dust(double *logT_start_dust, double *d_logT_dust, chemistry_data *my_chemistry)
{
    *logT_start_dust = log(my_chemistry->DustTemperatureStart);
    *d_logT_dust = (log(my_chemistry->DustTemperatureEnd) - *logT_start_dust)
                                / (my_chemistry->NumberOfDustTemperatureBins - 1);
}

//Define a function which is able to add a scalar rate to the calculation.
int add_scalar_reaction_rate(double *rate_ptr, scalar_rate_function my_function, double units,
                        chemistry_data *my_chemistry)
{
    //As these are scalar they have no temperature-dependence.
    (*rate_ptr) = my_function(units, my_chemistry);

    return SUCCESS;
}

/// Initialize values in a pre-allocated rate_ptr
static int init_preallocated_rate(double *rate_ptr, rate_function my_function,
                                  double units, chemistry_data *my_chemistry)
{
  // Calculate temperature spacing.
  double logT_start, d_logT;
  logT_spacing(&logT_start, &d_logT, my_chemistry);

  //Calculate rate at each temperature.
  for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++){
    // Set rate to tiny for safety.
    rate_ptr[i] = tiny;

    // Calculate bin temperature.
    double logT = logT_start + i*d_logT;
    double T = exp(logT);

    // Calculate rate and store.
    rate_ptr[i] = my_function(T, units, my_chemistry);
  }

  return GR_SUCCESS;
}

/// Allocate and initialize a reaction rate
static int add_reaction_rate(double **rate_ptr, rate_function my_function,
                             double units,
                             chemistry_data *my_chemistry)
{
  *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins
                              * sizeof(double));
  return init_preallocated_rate(*rate_ptr, my_function, units, my_chemistry);
}

//Define a function which will calculate the k13dd rates.
int add_k13dd_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate array memory to the pointer.
    *rate_ptr = (double*)malloc(14 * my_chemistry->NumberOfTemperatureBins * sizeof(double));

    //Create a temporary array to retrieve the 14 coefficients calculated by the k13dd function.
    double temp_array[14];

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT;
    logT_spacing(&logT_start, &d_logT, my_chemistry);

    //Calculate k13dd for all temperatures. Store in results array.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++){
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        //Obtain 14 coefficients at this temperature and place them in results array.
        k13dd_rate(T, units, temp_array, my_chemistry);
        for (int j = 0; j < 14; j++){
            (*rate_ptr)[i + (my_chemistry->NumberOfTemperatureBins * j)] = temp_array[j];
        }
    }

    return SUCCESS;
}
    
//Define a function which will calculate h2dust rates.
int add_h2dust_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

// Define a function which will calculate h2dust_C rates.
int add_h2dust_C_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_C_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

// Define a function which will calculate h2dust_S rates.
int add_h2dust_S_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_S_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

} // anonymous namespace

//Definition of the initialise_rates function.
int grackle::impl::initialize_rates(
  chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
  code_units *my_units, double co_length_unit, double co_density_unit)
{ 
    // TODO: we REALLY need to do an error check that
    //    my_chemistry->NumberOfTemperatureBins >= 2
    //  is satisfied. We should give some thought about whether this should
    //  produce errors when primordial_chemistry == 0

    // TODO: also add error-checks for:
    //   TemperatureStart, TemperatureEnd, NumberOfDustTemperatureBins,
    //   DustTemperatureStart, DustTemperatureEnd

    // handle some allocations up front
    // TODO: we should really make a separate function that fully initializes
    //   the kcol_rate_tables (rather than what we do now and partially
    //   initialize the contents).
    if (my_chemistry->primordial_chemistry > 0) {
      // allocate storage for kcol_rate_tables
      my_rates->opaque_storage->kcol_rate_tables =
        (grackle::impl::CollisionalRxnRateCollection*) malloc
          (sizeof(grackle::impl::CollisionalRxnRateCollection));

      // allocate storage within kcol_rate_tables
      (*my_rates->opaque_storage->kcol_rate_tables) =
        grackle::impl::new_CollisionalRxnRateCollection
          (my_chemistry->NumberOfTemperatureBins);

      // set all of the entries within kcol_rate_tables to tiny
      // (this is very important for primordial_chemistry == 4 rates and
      // the metal chemistry rates)
      for (int i = 0; i < CollisionalRxnLUT::NUM_ENTRIES; i++) {
        double* ptr = my_rates->opaque_storage->kcol_rate_tables->data[i];
        for (int j = 0; j < my_chemistry->NumberOfTemperatureBins; j++) {
          ptr[j] = tiny;
        }
      }

      // make the my_rates->k[0-9]+ and my_rates->kz[0-9]+ alias the pointers
      // allocated within kcol_rate_tables
      // -> this logic is temporary! It exists to let us more gradually change
      //    other pieces of code so that it refers to kcol_rate_tables
      // -> (the logic will be removed before we try merging the commits
      //    introducing kcol_rate_tables into the newchem-cpp branch)
      grackle::impl::CollisionalRxnRateCollection* kcol_rate_tables =
        my_rates->opaque_storage->kcol_rate_tables;

      my_rates->k1 = kcol_rate_tables->data[CollisionalRxnLUT::k1];
      my_rates->k2 = kcol_rate_tables->data[CollisionalRxnLUT::k2];
      my_rates->k3 = kcol_rate_tables->data[CollisionalRxnLUT::k3];
      my_rates->k4 = kcol_rate_tables->data[CollisionalRxnLUT::k4];
      my_rates->k5 = kcol_rate_tables->data[CollisionalRxnLUT::k5];
      my_rates->k6 = kcol_rate_tables->data[CollisionalRxnLUT::k6];
      my_rates->k7 = kcol_rate_tables->data[CollisionalRxnLUT::k7];
      my_rates->k8 = kcol_rate_tables->data[CollisionalRxnLUT::k8];
      my_rates->k9 = kcol_rate_tables->data[CollisionalRxnLUT::k9];
      my_rates->k10 = kcol_rate_tables->data[CollisionalRxnLUT::k10];
      my_rates->k11 = kcol_rate_tables->data[CollisionalRxnLUT::k11];
      my_rates->k12 = kcol_rate_tables->data[CollisionalRxnLUT::k12];
      my_rates->k13 = kcol_rate_tables->data[CollisionalRxnLUT::k13];
      my_rates->k14 = kcol_rate_tables->data[CollisionalRxnLUT::k14];
      my_rates->k15 = kcol_rate_tables->data[CollisionalRxnLUT::k15];
      my_rates->k16 = kcol_rate_tables->data[CollisionalRxnLUT::k16];
      my_rates->k17 = kcol_rate_tables->data[CollisionalRxnLUT::k17];
      my_rates->k18 = kcol_rate_tables->data[CollisionalRxnLUT::k18];
      my_rates->k19 = kcol_rate_tables->data[CollisionalRxnLUT::k19];
      my_rates->k20 = kcol_rate_tables->data[CollisionalRxnLUT::k20];
      my_rates->k21 = kcol_rate_tables->data[CollisionalRxnLUT::k21];
      my_rates->k22 = kcol_rate_tables->data[CollisionalRxnLUT::k22];
      my_rates->k23 = kcol_rate_tables->data[CollisionalRxnLUT::k23];
      my_rates->k50 = kcol_rate_tables->data[CollisionalRxnLUT::k50];
      my_rates->k51 = kcol_rate_tables->data[CollisionalRxnLUT::k51];
      my_rates->k52 = kcol_rate_tables->data[CollisionalRxnLUT::k52];
      my_rates->k53 = kcol_rate_tables->data[CollisionalRxnLUT::k53];
      my_rates->k54 = kcol_rate_tables->data[CollisionalRxnLUT::k54];
      my_rates->k55 = kcol_rate_tables->data[CollisionalRxnLUT::k55];
      my_rates->k56 = kcol_rate_tables->data[CollisionalRxnLUT::k56];
      my_rates->k57 = kcol_rate_tables->data[CollisionalRxnLUT::k57];
      my_rates->k58 = kcol_rate_tables->data[CollisionalRxnLUT::k58];

      my_rates->k125 = kcol_rate_tables->data[CollisionalRxnLUT::k125];
      my_rates->k129 = kcol_rate_tables->data[CollisionalRxnLUT::k129];
      my_rates->k130 = kcol_rate_tables->data[CollisionalRxnLUT::k130];
      my_rates->k131 = kcol_rate_tables->data[CollisionalRxnLUT::k131];
      my_rates->k132 = kcol_rate_tables->data[CollisionalRxnLUT::k132];
      my_rates->k133 = kcol_rate_tables->data[CollisionalRxnLUT::k133];
      my_rates->k134 = kcol_rate_tables->data[CollisionalRxnLUT::k134];
      my_rates->k135 = kcol_rate_tables->data[CollisionalRxnLUT::k135];
      my_rates->k136 = kcol_rate_tables->data[CollisionalRxnLUT::k136];
      my_rates->k137 = kcol_rate_tables->data[CollisionalRxnLUT::k137];
      my_rates->k148 = kcol_rate_tables->data[CollisionalRxnLUT::k148];
      my_rates->k149 = kcol_rate_tables->data[CollisionalRxnLUT::k149];
      my_rates->k150 = kcol_rate_tables->data[CollisionalRxnLUT::k150];
      my_rates->k151 = kcol_rate_tables->data[CollisionalRxnLUT::k151];
      my_rates->k152 = kcol_rate_tables->data[CollisionalRxnLUT::k152];
      my_rates->k153 = kcol_rate_tables->data[CollisionalRxnLUT::k153];

      my_rates->kz15 = kcol_rate_tables->data[CollisionalRxnLUT::kz15];
      my_rates->kz16 = kcol_rate_tables->data[CollisionalRxnLUT::kz16];
      my_rates->kz17 = kcol_rate_tables->data[CollisionalRxnLUT::kz17];
      my_rates->kz18 = kcol_rate_tables->data[CollisionalRxnLUT::kz18];
      my_rates->kz19 = kcol_rate_tables->data[CollisionalRxnLUT::kz19];
      my_rates->kz20 = kcol_rate_tables->data[CollisionalRxnLUT::kz20];
      my_rates->kz21 = kcol_rate_tables->data[CollisionalRxnLUT::kz21];
      my_rates->kz22 = kcol_rate_tables->data[CollisionalRxnLUT::kz22];
      my_rates->kz23 = kcol_rate_tables->data[CollisionalRxnLUT::kz23];
      my_rates->kz24 = kcol_rate_tables->data[CollisionalRxnLUT::kz24];
      my_rates->kz25 = kcol_rate_tables->data[CollisionalRxnLUT::kz25];
      my_rates->kz26 = kcol_rate_tables->data[CollisionalRxnLUT::kz26];
      my_rates->kz27 = kcol_rate_tables->data[CollisionalRxnLUT::kz27];
      my_rates->kz28 = kcol_rate_tables->data[CollisionalRxnLUT::kz28];
      my_rates->kz29 = kcol_rate_tables->data[CollisionalRxnLUT::kz29];
      my_rates->kz30 = kcol_rate_tables->data[CollisionalRxnLUT::kz30];
      my_rates->kz31 = kcol_rate_tables->data[CollisionalRxnLUT::kz31];
      my_rates->kz32 = kcol_rate_tables->data[CollisionalRxnLUT::kz32];
      my_rates->kz33 = kcol_rate_tables->data[CollisionalRxnLUT::kz33];
      my_rates->kz34 = kcol_rate_tables->data[CollisionalRxnLUT::kz34];
      my_rates->kz35 = kcol_rate_tables->data[CollisionalRxnLUT::kz35];
      my_rates->kz36 = kcol_rate_tables->data[CollisionalRxnLUT::kz36];
      my_rates->kz37 = kcol_rate_tables->data[CollisionalRxnLUT::kz37];
      my_rates->kz38 = kcol_rate_tables->data[CollisionalRxnLUT::kz38];
      my_rates->kz39 = kcol_rate_tables->data[CollisionalRxnLUT::kz39];
      my_rates->kz40 = kcol_rate_tables->data[CollisionalRxnLUT::kz40];
      my_rates->kz41 = kcol_rate_tables->data[CollisionalRxnLUT::kz41];
      my_rates->kz42 = kcol_rate_tables->data[CollisionalRxnLUT::kz42];
      my_rates->kz43 = kcol_rate_tables->data[CollisionalRxnLUT::kz43];
      my_rates->kz44 = kcol_rate_tables->data[CollisionalRxnLUT::kz44];
      my_rates->kz45 = kcol_rate_tables->data[CollisionalRxnLUT::kz45];
      my_rates->kz46 = kcol_rate_tables->data[CollisionalRxnLUT::kz46];
      my_rates->kz47 = kcol_rate_tables->data[CollisionalRxnLUT::kz47];
      my_rates->kz48 = kcol_rate_tables->data[CollisionalRxnLUT::kz48];
      my_rates->kz49 = kcol_rate_tables->data[CollisionalRxnLUT::kz49];
      my_rates->kz50 = kcol_rate_tables->data[CollisionalRxnLUT::kz50];
      my_rates->kz51 = kcol_rate_tables->data[CollisionalRxnLUT::kz51];
      my_rates->kz52 = kcol_rate_tables->data[CollisionalRxnLUT::kz52];
      my_rates->kz53 = kcol_rate_tables->data[CollisionalRxnLUT::kz53];
      my_rates->kz54 = kcol_rate_tables->data[CollisionalRxnLUT::kz54];
    }

    int anyDust;
    if ( my_chemistry->h2_on_dust > 0 || my_chemistry->dust_chemistry > 0 || my_chemistry->dust_recombination_cooling > 0) {
        anyDust = TRUE;
    } else {
        anyDust = FALSE;
    }

    //* Obtain the conversion factors which convert between code and physical units. We define a_value = 1 at z = zInit such that a = a_value * [a].
    double timeBase1 = my_units->time_units;
    double lengthBase1 = co_length_unit / (my_units->a_value * my_units->a_units);
    double densityBase1 = co_density_unit * pow(my_units->a_value * my_units->a_units, 3);

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
    double coolingUnits = (pow(my_units->a_units, 5) * pow(lengthBase1, 2) * pow(mh, 2))
                          / (densityBase1 * pow(timeBase1, 3));

    // These always need to be allocated since we define other variables by them.
    my_rates->gr_N = (int*)calloc(2, sizeof(int));

    if (my_chemistry->use_primordial_continuum_opacity == 1) {
      initialize_primordial_opacity(my_chemistry, my_rates);
    }

    //* 3) Units for radiative transfer coefficients are 1/[time].

    //* Compute rates for primordial chemistry.

    //* 1) Rate Coefficients (excluding the external radiation field)
    if (my_chemistry->primordial_chemistry > 0){ 
        //--------Calculate multispecies collissional rates--------
        add_reaction_rate(&my_rates->k1, k1_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k3, k3_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k4, k4_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k2, k2_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k5, k5_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k6, k6_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k7, k7_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k8, k8_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k9, k9_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k10, k10_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k11, k11_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k12, k12_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k14, k14_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k15, k15_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k16, k16_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k17, k17_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k18, k18_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k19, k19_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k20, k20_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k23, k23_rate, kUnit, my_chemistry);

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
        add_k13dd_reaction_rate(&my_rates->k13dd, kUnit, my_chemistry);

        //--------Calculate 3-body H2 rate--------

        // Calculated by the same method as done in the original code. First is the fit to 
        // A.E. Orel 1987, J.Chem.Phys., 87, 314, which is matched to the 1/T of 
        // Palla etal (1983) -- which is four times smaller than the Palla rate.
        
        //Varying threebody and corresponding collisional dissociation rates from Simon.
        add_reaction_rate(&my_rates->k13, k13_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k21, k21_rate, kUnit_3Bdy, my_chemistry);
        add_reaction_rate(&my_rates->k22, k22_rate, kUnit_3Bdy, my_chemistry);
        
        //--------Deuterium Rates--------
        add_reaction_rate(&my_rates->k50, k50_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k51, k51_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k52, k52_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k53, k53_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k54, k54_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k55, k55_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k56, k56_rate, kUnit, my_chemistry);

        //--------New H Ionization Rates--------
        add_reaction_rate(&my_rates->k57, k57_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->k58, k58_rate, kUnit, my_chemistry);

        //H2 formation on dust grains requires loop over the dust temperature.
        add_h2dust_reaction_rate(&my_rates->h2dust, kUnit, my_chemistry);

        //H2 formation heating terms from Equation 23, Omuaki (2000).
        add_reaction_rate(&my_rates->n_cr_n, n_cr_n_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d1, n_cr_d1_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d2, n_cr_d2_rate, kUnit, my_chemistry);

        //* 2) Cooling and heating Rates (excluding external radiation field)

        //* a) Collisional excitations (Black 1981; Cen 1992).

        add_reaction_rate(&my_rates->ceHI, ceHI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeI, ceHeI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeII, ceHeII_rate, coolingUnits, my_chemistry);

        //* b) Collisional ionizations (Cen 1992, Abel 1996).

        add_reaction_rate(&my_rates->ciHeIS, ciHeIS_rate, coolingUnits, my_chemistry);
        //Collisional ionizations. Polynomial fits from Tom Abel.
        add_reaction_rate(&my_rates->ciHI, ciHI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeI, ciHeI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeII, ciHeII_rate, coolingUnits, my_chemistry);

        //* c) Recombinations (Hui & Gnedin 1997, except where noted).

        add_reaction_rate(&my_rates->reHII, reHII_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeII1, reHeII1_rate, coolingUnits, my_chemistry);
        //Dielectronic recombination (Cen, 1992).
        add_reaction_rate(&my_rates->reHeII2, reHeII2_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeIII, reHeIII_rate, coolingUnits, my_chemistry);

        //* d) Bremsstrahlung (Black, 1981 and Spitzer & Hart, 1979)
        add_reaction_rate(&my_rates->brem, brem_rate, coolingUnits, my_chemistry);

        //* e) Molecular hydrogen cooling

        //------Calculate Lepp & Shull rates------

        add_reaction_rate(&my_rates->vibh, vibh_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->hyd01k, hyd01k_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->h2k01, h2k01_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->rotl, rotl_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->roth, roth_rate, coolingUnits, my_chemistry);
         
        //------Calculate Galli & Palla (1999) Rates------ 
        //             (as fit by Tom Abel)
        
        //Low density limit (Galli & Palla).
        add_reaction_rate(&my_rates->GP99LowDensityLimit, GP99LowDensityLimit_rate, coolingUnits, my_chemistry);
        //High density limit from HM79.
        add_reaction_rate(&my_rates->GP99HighDensityLimit, GP99HighDensityLimit_rate, coolingUnits, my_chemistry);

        //------Low density rates (Glover & Abel, 2008)------

        //Excitation by HI.
        add_reaction_rate(&my_rates->GAHI, GAHI_rate, coolingUnits, my_chemistry);
        //Excitation by H2.
        add_reaction_rate(&my_rates->GAH2, GAH2_rate, coolingUnits, my_chemistry);
        //Excitation by He.
        add_reaction_rate(&my_rates->GAHe, GAHe_rate, coolingUnits, my_chemistry);
        //Excitation by HII.
        //Collisional excitation rates from Honvault et al (2011, Phys. Rev. Lett., 107, 023201) and 
        //Honvault et al (2012, Phys. Rev. Lett., 108, 109903).
        add_reaction_rate(&my_rates->GAHp, GAHp_rate, coolingUnits, my_chemistry);
        //Excitation by electrons.
        //Based on data from  Yoon et al (2008, J. Phys. Chem. Ref. Data, 37, 913).
        add_reaction_rate(&my_rates->GAel, GAel_rate, coolingUnits, my_chemistry);
        //------New fit to LTE rate------ 
        //from Glover (2015, MNRAS, 451, 2082)
        add_reaction_rate(&my_rates->H2LTE, H2LTE_rate, coolingUnits, my_chemistry);

        // Chiaki & Wise 2019 rates
        // Note: these are still defined in initialize_metal_chemistry_rates.c.
        // They should be moved someday.
        initialize_cooling_rate_H2(my_chemistry, my_rates, coolingUnits);
        initialize_cooling_rate_HD(my_chemistry, my_rates, coolingUnits);

        //* f) HD cooling.

        //HD cooling function has units of ergs cm^3 / s
        //Fit from Coppola et al 2011. LTE (ergs/s) -> hdlte (ergs cm3/s)
        add_reaction_rate(&my_rates->HDlte, HDlte_rate, coolingUnits, my_chemistry);
        //Rate based on (Wrathmall, Gusdorf & Flower, 2007) HD-H collisional excitation rates.
        add_reaction_rate(&my_rates->HDlow, HDlow_rate, coolingUnits, my_chemistry);

        //* g) CIE cooling (Ripamonti & Abel, 2003).

        //Below explanation originally written by Matt Turk:
        //  The above function returns the rate with units of ergs * s^-1 * cm^3 * gram^-1 *
        //  (gram in H2 molecules)^-1. To reproduce equation 5 in RA04, divide c_t_c_r by mh,
        //  so to get erg/s cm^3 we multiply by mh. Divide by 2 for the H2 mass --> number.
        add_reaction_rate(&my_rates->cieco, cieco_rate, coolingUnits, my_chemistry); 

    }//End of primordial chemistry if-statement.

    //* Dust-related processes.
    if (anyDust == TRUE) {

        //Energy transfer from gas to dust grains.
        //(Equation 2.15, Hollenbach & McKee, 1989)
        add_reaction_rate(&my_rates->gas_grain, gasGrain_rate, coolingUnits, my_chemistry); 

        //Electron recombination onto dust grains.
        //(Equation 9, Wolfire et al., 1995)
        add_reaction_rate(&my_rates->regr, regr_rate, coolingUnits, my_chemistry);

        //H2 formation on dust grains with C and S compositions
        add_h2dust_C_reaction_rate(&my_rates->h2dustC, kUnit, my_chemistry);
        add_h2dust_S_reaction_rate(&my_rates->h2dustS, kUnit, my_chemistry);

        //Heating of dust by interstellar radiation field, with an arbitrary grain size distribution
        add_scalar_reaction_rate(&my_rates->gamma_isrf2, gamma_isrf2_rate, coolingUnits, my_chemistry);

        //Gas-grain energy transfer, with an arbitrary grain size distribution
        add_reaction_rate(&my_rates->gas_grain2, gasGrain2_rate, coolingUnits, my_chemistry);

        //Grain growth rate
        add_reaction_rate(&my_rates->grain_growth_rate, grain_growth_rate, kUnit, my_chemistry);

    }//End of anyDust if-statement.

    //Below rates are scalar -- temperature independent.

    //* i) Compton cooling (Peebles 1971).
    add_scalar_reaction_rate(&my_rates->comp, comp_rate, coolingUnits, my_chemistry); 

    //* j) Photoelectric heating by UV-irradiated dust (Wolfire 1995).

    add_scalar_reaction_rate(&my_rates->gammah, gammah_rate, coolingUnits, my_chemistry); 
    //Heating of dust by interstellar radiation field.
    //(Equation B15, Krumholz, 2014)
    add_scalar_reaction_rate(&my_rates->gamma_isrf, gamma_isrf_rate, coolingUnits, my_chemistry); 

    /* Metal chemistry rates */
    if (grackle::impl::initialize_metal_chemistry_rates(my_chemistry, my_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in initialize_metal_chemistry_rates.\n");
      return FAIL;
    }
    /* Dust rates */
    if (grackle::impl::initialize_dust_yields(my_chemistry, my_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in initialize_dust_yields.\n");
      return FAIL;
    }

    //End of function definition.
    return SUCCESS;
}

