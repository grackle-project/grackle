//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement machinery for initializing the reaction and cooling rates related
/// to metal species
///
/// Table of rate coefficients for primordial_chemistry >= 4:
///
/// | rate-name | Reaction                          |
/// | --------- | --------------------------------- |
/// | k125      | HDII +  HI   ->  HII  +  HDI      |
/// | k129      | DI   +  HII  ->  HDII +  p        |
/// | k130      | DII  +  HI   ->  HDII +  p        |
/// | k131      | HDII +  e    ->  HI   +  DI       |
/// | k132      | DI   +  e    ->  DM   +  p        |
/// | k133      | DII  +  DM   ->  DI   +  DI       |
/// | k134      | HII  +  DM   ->  DI   +  HI       |
/// | k135      | HM   +  DI   ->  HI   +  DM       |
/// | k136      | DM   +  HI   ->  DI   +  HM       |
/// | k137      | DM   +  HI   ->  HDI  +  e        |
/// | -         | -                                 |
/// | k148      | HeI    +  HII  ->  HeHII  +  p    |
/// | k149      | HeI    +  HII  ->  HeHII  +  p    |
/// | k150      | HeI    +  H2II ->  HeHII  +  HI   |
/// | k151      | HeII   +  HI   ->  HeHII  +  p    |
/// | k152      | HeHII  +  HI   ->  HeI    +  H2II |
/// | k153      | HeHII  +  e    ->  HeI    +  HI   |
///
/// Table of rate coefficients for metal species:
///
/// | rate-name | Reactions                        |
/// | --------- | -------------------------------- |
/// | kz15      | HI     +  CH   ->  CI     +  H2I |
/// | kz16      | HI     +  CH2  ->  CH     +  H2I |
/// | kz17      | HI     +  OH   ->  H2I    +  OI  |
/// | kz18      | HI     +  H2O  ->  OH     +  H2I |
/// | kz19      | HI     +  O2   ->  OH     +  OI  |
/// | kz20      | CI     +  H2I  ->  CH     +  HI  |
/// | kz21      | OI     +  H2I  ->  OH     +  HI  |
/// | kz22      | HII    +  OI   ->  OII    +  HI  |
/// | kz23      | H2I    +  CH   ->  CH2    +  HI  |
/// | kz24      | H2I    +  OH   ->  H2O    +  HI  |
/// | kz25      | OH     +  OH   ->  H2O    +  OI  |
/// | kz26      | OH     +  CO   ->  CO2    +  HI  |
/// | kz27      | CI     +  HI   ->  CH     +  p   |
/// | kz28      | CI     +  OH   ->  CO     +  HI  |
/// | kz29      | CI     +  O2   ->  CO     +  OI  |
/// | kz30      | OI     +  HI   ->  OH     +  p   |
/// | kz31      | OI     +  OI   ->  O2     +  p   |
/// | kz32      | OI     +  CH   ->  CO     +  HI  |
/// | kz33      | OI     +  OH   ->  O2     +  HI  |
/// | kz34      | HII    +  OH   ->  OHII   +  HI  |
/// | kz35      | HII    +  H2O  ->  H2OII  +  HI  |
/// | kz36      | HII    +  O2   ->  O2II   +  HI  |
/// | kz37      | CII    +  OH   ->  COII   +  HI  |
/// | kz38      | CII    +  O2   ->  OII    +  CO  |
/// | kz39      | OII    +  HI   ->  HII    +  OI  |
/// | kz40      | OII    +  H2I  ->  OHII   +  HI  |
/// | kz41      | OHII   +  H2I  ->  H2OII  +  HI  |
/// | kz42      | H2OII  +  H2I  ->  H3OII  +  HI  |
/// | kz43      | COII   +  HI   ->  HII    +  CO  |
/// | kz44      | CII    +  e    ->  CI     +  p   |
/// | kz45      | OII    +  e    ->  OI     +  p   |
/// | kz46      | H2OII  +  e    ->  OH     +  HI  |
/// | kz47      | H2OII  +  e    ->  OI     +  H2I |
/// | kz48      | H3OII  +  e    ->  H2O    +  HI  |
/// | kz49      | H3OII  +  e    ->  OH     +  HI  |
/// | kz50      | O2II   +  e    ->  OI     +  OI  |
/// | kz51      | H2I    +  CI   ->  CH2    +  p   |
/// | kz52      | SiI    +  OH   ->  SiOI   +  HI  |
/// | kz53      | SiI    +  O2   ->  SiOI   +  OI  |
/// | kz54      | SiOI   +  OH   ->  SiO2I  +  HI  |
///
//===----------------------------------------------------------------------===//

#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "grackle.h"
#include "grackle_macros.h"
#include "interp_table_utils.h" // free_interp_grid_
#include "initialize_metal_chemistry_rates.hpp"  // forward declarations
#include "internal_units.h"  // InternalGrUnits
#include "phys_constants.h"
#include "grackle_rate_functions.h" // forward declarations of some funcs

#define tiny 1.0e-20
#define tevk 1.1605e+4

static int allocate_rates_metal(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

int grackle::impl::initialize_metal_chemistry_rates(
  chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
  code_units *my_units)
{

  /* TO-DO: k125 - k153 are primordial_chemistry=4.
     These should be moved to initialize_rates.c so this is only metal species. */
  if (my_chemistry->primordial_chemistry == 0)
    return SUCCESS;


  // temporarily construct the InternalGrUnits struct
  // -> the construction logic deduplicates a lot of logic that was
  //    previously copied and pasted across a lot of fortran files
  InternalGrUnits internalu = new_internalu_legacy_C_(my_units);
  const double kunit = internalu_calc_kunit_(internalu);
  const double coolunit = internalu.coolunit;

  // compute log spacing of the temperature table
  //
  // TODO: address the ubiquity of this variable (among transcribed routines).
  //       To ensure that every part of the code uses exactly the same value,
  //       we should either:
  //       1. cache this quantity
  //       2. have a function that we use everywhere to call it
  const double dlogtem = (
    (log(my_chemistry->TemperatureEnd) -
     log(my_chemistry->TemperatureStart)) /
    (double)(my_chemistry->NumberOfTemperatureBins-1)
  );


  // Allocate buffers to hold the rates
  allocate_rates_metal(my_chemistry, my_rates);

  // Initialize constants to tiny
  for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
        my_rates->cieY06[i] = tiny;

        my_rates->k125[i] = tiny;
        my_rates->k129[i] = tiny;
        my_rates->k130[i] = tiny;
        my_rates->k131[i] = tiny;
        my_rates->k132[i] = tiny;
        my_rates->k133[i] = tiny;
        my_rates->k134[i] = tiny;
        my_rates->k135[i] = tiny;
        my_rates->k136[i] = tiny;
        my_rates->k137[i] = tiny;
        my_rates->k148[i] = tiny;
        my_rates->k149[i] = tiny;
        my_rates->k150[i] = tiny;
        my_rates->k151[i] = tiny;
        my_rates->k152[i] = tiny;
        my_rates->k153[i] = tiny;

        my_rates->kz15[i] = tiny;
        my_rates->kz16[i] = tiny;
        my_rates->kz17[i] = tiny;
        my_rates->kz18[i] = tiny;
        my_rates->kz19[i] = tiny;
        my_rates->kz20[i] = tiny;
        my_rates->kz21[i] = tiny;
        my_rates->kz22[i] = tiny;
        my_rates->kz23[i] = tiny;
        my_rates->kz24[i] = tiny;
        my_rates->kz25[i] = tiny;
        my_rates->kz26[i] = tiny;
        my_rates->kz27[i] = tiny;
        my_rates->kz28[i] = tiny;
        my_rates->kz29[i] = tiny;
        my_rates->kz30[i] = tiny;
        my_rates->kz31[i] = tiny;
        my_rates->kz32[i] = tiny;
        my_rates->kz33[i] = tiny;
        my_rates->kz34[i] = tiny;
        my_rates->kz35[i] = tiny;
        my_rates->kz36[i] = tiny;
        my_rates->kz37[i] = tiny;
        my_rates->kz38[i] = tiny;
        my_rates->kz39[i] = tiny;
        my_rates->kz40[i] = tiny;
        my_rates->kz41[i] = tiny;
        my_rates->kz42[i] = tiny;
        my_rates->kz43[i] = tiny;
        my_rates->kz44[i] = tiny;
        my_rates->kz45[i] = tiny;
        my_rates->kz46[i] = tiny;
        my_rates->kz47[i] = tiny;
        my_rates->kz48[i] = tiny;
        my_rates->kz49[i] = tiny;
        my_rates->kz50[i] = tiny;
        my_rates->kz51[i] = tiny;
        my_rates->kz52[i] = tiny;
        my_rates->kz53[i] = tiny;
        my_rates->kz54[i] = tiny;
  }

  // Fill in tables of
  //   - CIE H2 cooling rates,
  //   - collisional rate coefficients for primordial_chemistry >= 4 species
  //   - collisional rate coefficients for metal species
  //
  // We do this for every temperature in the range spanned by
  // my_chemistry->TemperatureStart & my_chemistry->TemperatureEnd
  for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
    // Compute the current temperature
    // NOTE: an earlier version of this comment notes that temperature is in
    //       eV, but I think that's incorrect
    double logttt = log(my_chemistry->TemperatureStart) + (double)(i  )*dlogtem;
    double ttt = exp(logttt);
    double ttt300 = ttt / 300.0;

    // CIE H2 cooling rate from Yoshida et al. (2006)
    my_rates->cieY06[i] =
      pow(10.0,
          -116.6 + 96.34  * log10(ttt)
                 - 47.153 * pow(log10(ttt), 2)
                 + 10.744 * pow(log10(ttt), 3)
                 -  0.916 * pow(log10(ttt), 4) ) / coolunit;

        my_rates->k125[i] = 6.4e-10;
        my_rates->k129[i] = 3.9e-19 * pow(ttt300, 1.8) * exp(20.0/ttt);
        my_rates->k130[i] = 3.9e-19 * pow(ttt300, 1.8) * exp(20.0/ttt);
        my_rates->k131[i] = 3.4e-9 * pow(ttt300, -0.4);
        my_rates->k132[i] = 3.0e-16 * pow(ttt300, 0.95) * exp(-ttt/9320.0);
        my_rates->k133[i] = 5.7e-8 * pow(ttt300, -0.50);
        my_rates->k134[i] = 4.6e-8 * pow(ttt300, -0.50);
        my_rates->k135[i] = 4.6e-8 * pow(ttt300, -0.50);
        my_rates->k136[i] = 6.4e-9 * pow(ttt300, 0.41);
        my_rates->k137[i] = 1.5e-9 * pow(ttt300, -0.1);
  
        my_rates->k148[i] = 5.0e-21;
        if(ttt < 1000.0)
          my_rates->k149[i] = 7.60e-18 * pow(ttt, -0.50);
        else
          my_rates->k149[i] = 3.45e-16 * pow(ttt, -1.06);
        my_rates->k150[i] = 3.0e-10 * exp(-6717.0/ttt);
        if(ttt < 4000.0)
          my_rates->k151[i] = 1.6e-14 * pow(ttt, -0.33);
        else
          my_rates->k151[i] = 1.0e-15;
        my_rates->k152[i] = 9.1e-10;
        my_rates->k153[i] = 1.7e-7 * pow(ttt, -0.5);
  
        my_rates->kz15[i] = 4.98e-11;
        my_rates->kz16[i] = 2.70e-10;
        my_rates->kz17[i] = 7.00e-14 * pow(ttt300, 2.80) * exp(-1950.0/ttt);
        my_rates->kz18[i] = 6.83e-12 * pow(ttt300, 1.60) * exp(-9720.0/ttt);
        my_rates->kz19[i] = 3.30e-10 * exp(-8460.0/ttt);
        my_rates->kz20[i] = 6.64e-10 * exp(-11700.0/ttt);
        if(ttt < 1.0e7)
            my_rates->kz21[i] = 3.43e-13 * pow(ttt300, 2.67) * exp(-3160.0/ttt);
        else
            my_rates->kz21[i] = 3.43e-13 * pow(1.0e7/300.0, 2.67) * exp(-3160.0/1.0e7);
        // The rate comes from an experiment (297-3532 K).
        // We refrain to extrapolate it to high temperatures.
        my_rates->kz22[i] = 7.00e-10 * exp(-232.0/ttt);
        my_rates->kz23[i] = 2.38e-10 * exp(-1760.0/ttt);
        my_rates->kz24[i] = 1.55e-12 * pow(ttt300, 1.60) * exp(-1660.0/ttt);
        my_rates->kz25[i] = 1.65e-12 * pow(ttt300, 1.14) * exp(-50.0/ttt);
        my_rates->kz26[i] = 1.0e-13;
        my_rates->kz27[i] = 1.0e-17;
        my_rates->kz28[i] = 1.1e-10 * pow(ttt300, 0.5);
        my_rates->kz29[i] = 3.3e-11;
        my_rates->kz30[i] = 9.9e-19 * pow(ttt300, -0.38);
        my_rates->kz31[i] = 4.9e-20 * pow(ttt300, 1.58);
        my_rates->kz32[i] = 6.6e-11;
        my_rates->kz33[i] = 4.34e-11 * pow(ttt300, -0.5) * exp(-30.0/ttt);
        my_rates->kz34[i] = 2.1e-9;
        my_rates->kz35[i] = 6.9e-9;
        my_rates->kz36[i] = 2.0e-9;
        my_rates->kz37[i] = 7.7e-10;
        my_rates->kz38[i] = 6.2e-10;
        my_rates->kz39[i] = 6.8e-10;
        my_rates->kz40[i] = 1.7e-9;
        my_rates->kz41[i] = 1.01e-9;
        my_rates->kz42[i] = 8.3e-10;
        my_rates->kz43[i] = 7.5e-10;
        my_rates->kz44[i] = 4.4e-12 * pow(ttt300, -0.61);
        my_rates->kz45[i] = 3.4e-12 * pow(ttt300, -0.63);
        my_rates->kz46[i] = 1.6e-7 * pow(ttt300, -0.5);
        my_rates->kz47[i] = 2.0e-7 * pow(ttt300, -0.5);
        my_rates->kz48[i] = 3.5e-7 * pow(ttt300, -0.5);
        my_rates->kz49[i] = 6.5e-7 * pow(ttt300, -0.5);
        my_rates->kz50[i] = 1.95e-7 * pow(ttt300, -0.7);
        my_rates->kz51[i] = 1.0e-17;
        my_rates->kz52[i]  = 3.00e-11;
        my_rates->kz53[i]  = 1.30e-11 * exp(-111.0/ttt);
        my_rates->kz54[i]  = 2.00e-13;

//      printf("CHECK %13.5e %13.5e %13.5e %13.5e\n", ttt
//         , my_rates->k125[i]
//         , my_rates->k129[i]
//         , my_rates->k130[i]);
  }

  for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
        my_rates->k125[i] = fmax(my_rates->k125[i], tiny) / kunit;
        my_rates->k129[i] = fmax(my_rates->k129[i], tiny) / kunit;
        my_rates->k130[i] = fmax(my_rates->k130[i], tiny) / kunit;
        my_rates->k131[i] = fmax(my_rates->k131[i], tiny) / kunit;
        my_rates->k132[i] = fmax(my_rates->k132[i], tiny) / kunit;
        my_rates->k133[i] = fmax(my_rates->k133[i], tiny) / kunit;
        my_rates->k134[i] = fmax(my_rates->k134[i], tiny) / kunit;
        my_rates->k135[i] = fmax(my_rates->k135[i], tiny) / kunit;
        my_rates->k136[i] = fmax(my_rates->k136[i], tiny) / kunit;
        my_rates->k137[i] = fmax(my_rates->k137[i], tiny) / kunit;

        my_rates->k148[i] = fmax(my_rates->k148[i], tiny) / kunit;
        my_rates->k149[i] = fmax(my_rates->k149[i], tiny) / kunit;
        my_rates->k150[i] = fmax(my_rates->k150[i], tiny) / kunit;
        my_rates->k151[i] = fmax(my_rates->k151[i], tiny) / kunit;
        my_rates->k152[i] = fmax(my_rates->k152[i], tiny) / kunit;
        my_rates->k153[i] = fmax(my_rates->k153[i], tiny) / kunit;

        my_rates->kz15[i] = fmax(my_rates->kz15[i], tiny) / kunit;
        my_rates->kz16[i] = fmax(my_rates->kz16[i], tiny) / kunit;
        my_rates->kz17[i] = fmax(my_rates->kz17[i], tiny) / kunit;
        my_rates->kz18[i] = fmax(my_rates->kz18[i], tiny) / kunit;
        my_rates->kz19[i] = fmax(my_rates->kz19[i], tiny) / kunit;
        my_rates->kz20[i] = fmax(my_rates->kz20[i], tiny) / kunit;
        my_rates->kz21[i] = fmax(my_rates->kz21[i], tiny) / kunit;
        my_rates->kz22[i] = fmax(my_rates->kz22[i], tiny) / kunit;
        my_rates->kz23[i] = fmax(my_rates->kz23[i], tiny) / kunit;
        my_rates->kz24[i] = fmax(my_rates->kz24[i], tiny) / kunit;
        my_rates->kz25[i] = fmax(my_rates->kz25[i], tiny) / kunit;
        my_rates->kz26[i] = fmax(my_rates->kz26[i], tiny) / kunit;
        my_rates->kz27[i] = fmax(my_rates->kz27[i], tiny) / kunit;
        my_rates->kz28[i] = fmax(my_rates->kz28[i], tiny) / kunit;
        my_rates->kz29[i] = fmax(my_rates->kz29[i], tiny) / kunit;
        my_rates->kz30[i] = fmax(my_rates->kz30[i], tiny) / kunit;
        my_rates->kz31[i] = fmax(my_rates->kz31[i], tiny) / kunit;
        my_rates->kz32[i] = fmax(my_rates->kz32[i], tiny) / kunit;
        my_rates->kz33[i] = fmax(my_rates->kz33[i], tiny) / kunit;
        my_rates->kz34[i] = fmax(my_rates->kz34[i], tiny) / kunit;
        my_rates->kz35[i] = fmax(my_rates->kz35[i], tiny) / kunit;
        my_rates->kz36[i] = fmax(my_rates->kz36[i], tiny) / kunit;
        my_rates->kz37[i] = fmax(my_rates->kz37[i], tiny) / kunit;
        my_rates->kz38[i] = fmax(my_rates->kz38[i], tiny) / kunit;
        my_rates->kz39[i] = fmax(my_rates->kz39[i], tiny) / kunit;
        my_rates->kz40[i] = fmax(my_rates->kz40[i], tiny) / kunit;
        my_rates->kz41[i] = fmax(my_rates->kz41[i], tiny) / kunit;
        my_rates->kz42[i] = fmax(my_rates->kz42[i], tiny) / kunit;
        my_rates->kz43[i] = fmax(my_rates->kz43[i], tiny) / kunit;
        my_rates->kz44[i] = fmax(my_rates->kz44[i], tiny) / kunit;
        my_rates->kz45[i] = fmax(my_rates->kz45[i], tiny) / kunit;
        my_rates->kz46[i] = fmax(my_rates->kz46[i], tiny) / kunit;
        my_rates->kz47[i] = fmax(my_rates->kz47[i], tiny) / kunit;
        my_rates->kz48[i] = fmax(my_rates->kz48[i], tiny) / kunit;
        my_rates->kz49[i] = fmax(my_rates->kz49[i], tiny) / kunit;
        my_rates->kz50[i] = fmax(my_rates->kz50[i], tiny) / kunit;
        my_rates->kz51[i] = fmax(my_rates->kz51[i], tiny) / kunit;
        my_rates->kz52[i] = fmax(my_rates->kz52[i], tiny) / kunit;
        my_rates->kz53[i] = fmax(my_rates->kz53[i], tiny) / kunit;
        my_rates->kz54[i] = fmax(my_rates->kz54[i], tiny) / kunit;
  }

      initialize_cooling_rate_CI (my_chemistry, my_rates, coolunit);
      initialize_cooling_rate_CII(my_chemistry, my_rates, coolunit);
      initialize_cooling_rate_OI (my_chemistry, my_rates, coolunit);
      initialize_cooling_rate_CO (my_chemistry, my_rates, coolunit);
      initialize_cooling_rate_OH (my_chemistry, my_rates, coolunit);
      initialize_cooling_rate_H2O(my_chemistry, my_rates, coolunit);

  return SUCCESS;
}

int grackle::impl::free_metal_chemistry_rates(chemistry_data *my_chemistry,
                                              chemistry_data_storage *my_rates)
{

  /* TO-DO: k125 - k153 are primordial_chemistry=4.
     These should be moved to initialize_rates.c so this is only metal species. */
  if (my_chemistry->primordial_chemistry == 0)
    return SUCCESS;


  free_interp_grid_(&my_rates->LCI);
  free_interp_grid_(&my_rates->LCII);
  free_interp_grid_(&my_rates->LOI);

  free_interp_grid_(&my_rates->LCO);
  free_interp_grid_(&my_rates->LOH);
  free_interp_grid_(&my_rates->LH2O);

  GRACKLE_FREE(my_rates->k125);
  GRACKLE_FREE(my_rates->k129);
  GRACKLE_FREE(my_rates->k130);
  GRACKLE_FREE(my_rates->k131);
  GRACKLE_FREE(my_rates->k132);
  GRACKLE_FREE(my_rates->k133);
  GRACKLE_FREE(my_rates->k134);
  GRACKLE_FREE(my_rates->k135);
  GRACKLE_FREE(my_rates->k136);
  GRACKLE_FREE(my_rates->k137);
  GRACKLE_FREE(my_rates->k148);
  GRACKLE_FREE(my_rates->k149);
  GRACKLE_FREE(my_rates->k150);
  GRACKLE_FREE(my_rates->k151);
  GRACKLE_FREE(my_rates->k152);
  GRACKLE_FREE(my_rates->k153);

  GRACKLE_FREE(my_rates->kz15);
  GRACKLE_FREE(my_rates->kz16);
  GRACKLE_FREE(my_rates->kz17);
  GRACKLE_FREE(my_rates->kz18);
  GRACKLE_FREE(my_rates->kz19);
  GRACKLE_FREE(my_rates->kz20);
  GRACKLE_FREE(my_rates->kz21);
  GRACKLE_FREE(my_rates->kz22);
  GRACKLE_FREE(my_rates->kz23);
  GRACKLE_FREE(my_rates->kz24);
  GRACKLE_FREE(my_rates->kz25);
  GRACKLE_FREE(my_rates->kz26);
  GRACKLE_FREE(my_rates->kz27);
  GRACKLE_FREE(my_rates->kz28);
  GRACKLE_FREE(my_rates->kz29);
  GRACKLE_FREE(my_rates->kz30);
  GRACKLE_FREE(my_rates->kz31);
  GRACKLE_FREE(my_rates->kz32);
  GRACKLE_FREE(my_rates->kz33);
  GRACKLE_FREE(my_rates->kz34);
  GRACKLE_FREE(my_rates->kz35);
  GRACKLE_FREE(my_rates->kz36);
  GRACKLE_FREE(my_rates->kz37);
  GRACKLE_FREE(my_rates->kz38);
  GRACKLE_FREE(my_rates->kz39);
  GRACKLE_FREE(my_rates->kz40);
  GRACKLE_FREE(my_rates->kz41);
  GRACKLE_FREE(my_rates->kz42);
  GRACKLE_FREE(my_rates->kz43);
  GRACKLE_FREE(my_rates->kz44);
  GRACKLE_FREE(my_rates->kz45);
  GRACKLE_FREE(my_rates->kz46);
  GRACKLE_FREE(my_rates->kz47);
  GRACKLE_FREE(my_rates->kz48);
  GRACKLE_FREE(my_rates->kz49);
  GRACKLE_FREE(my_rates->kz50);
  GRACKLE_FREE(my_rates->kz51);
  GRACKLE_FREE(my_rates->kz52);
  GRACKLE_FREE(my_rates->kz53);
  GRACKLE_FREE(my_rates->kz54);

  GRACKLE_FREE(my_rates->cieY06 );

  return SUCCESS;
}

struct regular_range_{
  int count;
  double start;
  double step;
};

/// helper function to assist with setting up a generic interp_grid_props
static void setup_generic_grid_props_(gr_interp_grid_props* grid_props,
                                      int rank,
                                      const struct regular_range_* parameters)
{
  grid_props->rank = rank;

  long long data_size = 1ll ? (rank > 0) : 0ll;
  for (int i = 0; i < rank; i++) {
    const struct regular_range_ par_range = parameters[i];

    double* arr = (double*)malloc(par_range.count * sizeof(double));
    for(int j = 0; j < par_range.count; j++) {
      arr[j] = par_range.start + (double)j * par_range.step;
    }

    grid_props->dimension[i] = par_range.count;
    grid_props->parameters[i] = arr;
    grid_props->parameter_spacing[i] = par_range.step;

    data_size *= (long long)par_range.count;
  }
  grid_props->data_size = data_size;
}


/// helper function to assist with setting up a interp_grid for cooling
static void setup_cool_interp_grid_(gr_interp_grid* grid,
                                    int rank,
                                    const struct regular_range_* parameters,
                                    const double* data,
                                    double log_coolrate)
{
  setup_generic_grid_props_(&grid->props, rank, parameters);
  const long long data_size = grid->props.data_size;
  grid->data = (double*)malloc(data_size * sizeof(double));
  for(long long i = 0; i < data_size; i++) {
    grid->data[i] = data[i] + log_coolrate;
  }
}

// all of the following functions are declared extern "C" because (at the time
// of writing) 

extern "C" void initialize_cooling_rate_H2(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {16,  20.0, 1.0}, // log10(number-density like)
    {11,   1.6, 0.2}, // log10(temperature)
    {21, -10.0, 1.0}  // log10(H2 number density)
  };

  double L[] = 
    {41.47,  40.47,  39.47,  38.47,  37.47,  36.47,  35.47,  34.47,  33.47,  32.47,  31.47,  30.48,  29.58,  29.03,  28.91,  28.89,  28.89,  28.89,  28.89,  28.89,  28.89, 
     39.38,  38.38,  37.38,  36.38,  35.38,  34.38,  33.38,  32.38,  31.38,  30.38,  29.38,  28.39,  27.51,  27.03,  26.89,  26.86,  26.86,  26.86,  26.86,  26.86,  26.86, 
     37.96,  36.96,  35.96,  34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.98,  26.13,  25.64,  25.39,  25.33,  25.32,  25.32,  25.32,  25.32,  25.32, 
     36.94,  35.94,  34.94,  33.94,  32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.97,  25.16,  24.59,  24.25,  24.17,  24.16,  24.16,  24.16,  24.16,  24.16, 
     36.16,  35.16,  34.16,  33.16,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.17,  25.21,  24.43,  23.77,  23.44,  23.32,  23.29,  23.28,  23.28,  23.28,  23.28, 
     35.53,  34.53,  33.53,  32.53,  31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.54,  24.61,  23.80,  23.14,  22.76,  22.54,  22.48,  22.48,  22.48,  22.48,  22.48, 
     34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.95,  25.95,  24.97,  24.04,  23.15,  22.42,  21.91,  21.71,  21.64,  21.63,  21.63,  21.63,  21.63, 
     34.38,  33.38,  32.36,  31.36,  30.36,  29.36,  28.36,  27.36,  26.36,  25.36,  24.38,  23.29,  22.19,  21.40,  20.99,  20.73,  20.62,  20.61,  20.60,  20.60,  20.60, 
     33.62,  32.69,  31.69,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.69,  23.65,  22.17,  21.01,  20.44,  20.11,  19.81,  19.72,  19.71,  19.71,  19.71,  19.71, 
     33.03,  32.00,  31.00,  30.00,  29.00,  28.00,  27.00,  26.00,  25.00,  24.00,  22.76,  20.89,  19.98,  19.68,  19.38,  19.13,  19.08,  19.07,  19.07,  19.07,  19.07, 
     32.36,  31.37,  30.37,  29.37,  28.37,  27.37,  26.37,  25.37,  24.37,  23.33,  21.68,  19.76,  19.23,  19.09,  18.83,  18.66,  18.64,  18.64,  18.63,  18.63,  18.63, 
     41.47,  40.47,  39.47,  38.47,  37.47,  36.47,  35.47,  34.47,  33.47,  32.47,  31.47,  30.48,  29.58,  29.03,  28.91,  28.89,  28.89,  28.89,  28.89,  28.89,  28.89, 
     39.38,  38.38,  37.38,  36.38,  35.38,  34.38,  33.38,  32.38,  31.38,  30.38,  29.38,  28.39,  27.51,  27.03,  26.89,  26.86,  26.86,  26.86,  26.86,  26.86,  26.86, 
     37.96,  36.96,  35.96,  34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.98,  26.13,  25.64,  25.39,  25.33,  25.32,  25.32,  25.32,  25.32,  25.32, 
     36.94,  35.94,  34.94,  33.94,  32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.97,  25.16,  24.59,  24.25,  24.17,  24.16,  24.16,  24.16,  24.16,  24.16, 
     36.16,  35.16,  34.16,  33.16,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.17,  25.21,  24.43,  23.77,  23.44,  23.32,  23.29,  23.28,  23.28,  23.28,  23.28, 
     35.53,  34.53,  33.53,  32.53,  31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.54,  24.61,  23.80,  23.14,  22.76,  22.54,  22.48,  22.48,  22.48,  22.48,  22.48, 
     34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.95,  25.95,  24.97,  24.04,  23.15,  22.42,  21.91,  21.71,  21.64,  21.63,  21.63,  21.63,  21.63, 
     34.38,  33.38,  32.36,  31.36,  30.36,  29.36,  28.36,  27.36,  26.36,  25.36,  24.38,  23.29,  22.19,  21.40,  20.99,  20.73,  20.62,  20.61,  20.60,  20.60,  20.60, 
     33.62,  32.69,  31.69,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.69,  23.65,  22.17,  21.01,  20.44,  20.11,  19.81,  19.72,  19.71,  19.71,  19.71,  19.71, 
     33.03,  32.00,  31.00,  30.00,  29.00,  28.00,  27.00,  26.00,  25.00,  24.00,  22.76,  20.89,  19.98,  19.68,  19.38,  19.13,  19.08,  19.07,  19.07,  19.07,  19.07, 
     32.36,  31.37,  30.37,  29.37,  28.37,  27.37,  26.37,  25.37,  24.37,  23.33,  21.68,  19.76,  19.23,  19.09,  18.83,  18.66,  18.64,  18.64,  18.64,  18.63,  18.63, 
     41.47,  40.47,  39.47,  38.47,  37.47,  36.47,  35.47,  34.47,  33.47,  32.47,  31.48,  30.49,  29.58,  29.04,  28.91,  28.89,  28.89,  28.89,  28.89,  28.89,  28.89, 
     39.38,  38.38,  37.38,  36.38,  35.38,  34.38,  33.38,  32.38,  31.38,  30.38,  29.38,  28.39,  27.51,  27.03,  26.89,  26.86,  26.86,  26.86,  26.86,  26.86,  26.86, 
     37.96,  36.96,  35.96,  34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.98,  26.13,  25.64,  25.39,  25.33,  25.32,  25.32,  25.32,  25.32,  25.32, 
     36.94,  35.94,  34.94,  33.94,  32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.97,  25.16,  24.59,  24.25,  24.17,  24.16,  24.16,  24.16,  24.16,  24.16, 
     36.16,  35.16,  34.16,  33.16,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.17,  25.22,  24.43,  23.78,  23.44,  23.32,  23.29,  23.29,  23.28,  23.28,  23.28, 
     35.53,  34.53,  33.53,  32.53,  31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.54,  24.61,  23.80,  23.14,  22.76,  22.54,  22.48,  22.48,  22.48,  22.48,  22.48, 
     34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.95,  25.95,  24.97,  24.04,  23.15,  22.43,  21.92,  21.71,  21.64,  21.63,  21.63,  21.63,  21.63, 
     34.38,  33.38,  32.36,  31.36,  30.36,  29.36,  28.36,  27.36,  26.36,  25.36,  24.38,  23.29,  22.19,  21.40,  20.99,  20.73,  20.62,  20.61,  20.60,  20.60,  20.60, 
     33.62,  32.70,  31.70,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.70,  23.65,  22.17,  21.01,  20.44,  20.12,  19.81,  19.72,  19.71,  19.71,  19.71,  19.71, 
     33.03,  32.00,  31.00,  30.00,  29.00,  28.00,  27.00,  26.00,  25.00,  24.00,  22.76,  20.89,  19.98,  19.68,  19.38,  19.13,  19.08,  19.07,  19.07,  19.07,  19.07, 
     32.37,  31.37,  30.37,  29.37,  28.37,  27.37,  26.37,  25.37,  24.37,  23.33,  21.68,  19.76,  19.24,  19.09,  18.83,  18.66,  18.64,  18.64,  18.64,  18.64,  18.64, 
     41.48,  40.48,  39.48,  38.48,  37.48,  36.48,  35.48,  34.48,  33.48,  32.48,  31.48,  30.49,  29.58,  29.04,  28.92,  28.90,  28.90,  28.90,  28.90,  28.90,  28.90, 
     39.39,  38.39,  37.39,  36.39,  35.39,  34.39,  33.39,  32.39,  31.39,  30.39,  29.39,  28.40,  27.52,  27.04,  26.90,  26.87,  26.87,  26.87,  26.87,  26.87,  26.87, 
     37.97,  36.97,  35.97,  34.97,  33.97,  32.97,  31.97,  30.97,  29.97,  28.97,  27.97,  26.99,  26.14,  25.66,  25.42,  25.35,  25.35,  25.35,  25.35,  25.35,  25.35, 
     36.95,  35.95,  34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.96,  25.98,  25.17,  24.62,  24.28,  24.20,  24.19,  24.19,  24.19,  24.19,  24.19, 
     36.18,  35.18,  34.18,  33.18,  32.18,  31.18,  30.18,  29.18,  28.18,  27.18,  26.19,  25.23,  24.45,  23.81,  23.47,  23.35,  23.31,  23.31,  23.31,  23.31,  23.31, 
     35.55,  34.55,  33.55,  32.55,  31.55,  30.55,  29.55,  28.55,  27.55,  26.55,  25.56,  24.64,  23.83,  23.17,  22.78,  22.56,  22.50,  22.50,  22.50,  22.50,  22.50, 
     34.98,  33.98,  32.98,  31.97,  30.98,  29.98,  28.98,  27.98,  26.98,  25.98,  25.00,  24.07,  23.17,  22.44,  21.93,  21.72,  21.66,  21.65,  21.65,  21.65,  21.65, 
     34.41,  33.41,  32.38,  31.38,  30.39,  29.39,  28.39,  27.39,  26.39,  25.39,  24.40,  23.31,  22.20,  21.41,  21.00,  20.75,  20.63,  20.62,  20.62,  20.62,  20.62, 
     33.64,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.67,  22.18,  21.01,  20.45,  20.13,  19.82,  19.73,  19.72,  19.72,  19.72,  19.72, 
     33.05,  32.02,  31.03,  30.03,  29.03,  28.03,  27.03,  26.03,  25.03,  24.02,  22.77,  20.89,  19.99,  19.69,  19.39,  19.13,  19.08,  19.08,  19.08,  19.08,  19.08, 
     32.39,  31.39,  30.39,  29.39,  28.39,  27.39,  26.39,  25.39,  24.39,  23.35,  21.68,  19.77,  19.25,  19.10,  18.83,  18.67,  18.64,  18.64,  18.64,  18.64,  18.64, 
     41.55,  40.55,  39.55,  38.55,  37.55,  36.55,  35.55,  34.55,  33.55,  32.55,  31.56,  30.57,  29.66,  29.12,  28.99,  28.98,  28.97,  28.97,  28.97,  28.97,  28.97, 
     39.46,  38.46,  37.46,  36.46,  35.46,  34.46,  33.46,  32.46,  31.46,  30.46,  29.46,  28.48,  27.60,  27.12,  27.00,  26.98,  26.98,  26.98,  26.98,  26.98,  26.98, 
     38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.06,  27.08,  26.23,  25.79,  25.60,  25.55,  25.54,  25.54,  25.54,  25.54,  25.54, 
     37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.06,  27.07,  26.10,  25.30,  24.82,  24.53,  24.45,  24.43,  24.43,  24.43,  24.43,  24.43, 
     36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.38,  24.63,  24.05,  23.70,  23.55,  23.51,  23.50,  23.50,  23.50,  23.50, 
     35.71,  34.71,  33.71,  32.71,  31.71,  30.71,  29.71,  28.71,  27.71,  26.72,  25.73,  24.81,  24.05,  23.37,  22.94,  22.71,  22.66,  22.65,  22.65,  22.65,  22.65, 
     35.17,  34.17,  33.17,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.17,  25.18,  24.26,  23.34,  22.55,  22.06,  21.87,  21.80,  21.79,  21.79,  21.79,  21.79, 
     34.61,  33.61,  32.58,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.60,  23.44,  22.26,  21.50,  21.12,  20.85,  20.74,  20.72,  20.72,  20.72,  20.72, 
     33.81,  32.91,  31.91,  30.91,  29.91,  28.91,  27.91,  26.91,  25.91,  24.91,  23.82,  22.22,  21.05,  20.54,  20.22,  19.89,  19.80,  19.79,  19.79,  19.79,  19.79, 
     33.22,  32.19,  31.19,  30.19,  29.19,  28.19,  27.19,  26.19,  25.19,  24.17,  22.85,  20.91,  20.05,  19.79,  19.46,  19.18,  19.13,  19.12,  19.12,  19.12,  19.12, 
     32.54,  31.55,  30.54,  29.54,  28.54,  27.54,  26.54,  25.54,  24.54,  23.48,  21.71,  19.80,  19.35,  19.20,  18.89,  18.70,  18.68,  18.67,  18.67,  18.67,  18.67, 
     42.07,  41.07,  40.07,  39.07,  38.07,  37.07,  36.07,  35.07,  34.07,  33.07,  32.07,  31.08,  30.17,  29.63,  29.51,  29.49,  29.49,  29.49,  29.49,  29.49,  29.49, 
     39.98,  38.98,  37.98,  36.98,  35.98,  34.98,  33.98,  32.98,  31.98,  30.98,  29.98,  28.99,  28.11,  27.64,  27.53,  27.52,  27.52,  27.52,  27.52,  27.52,  27.52, 
     38.58,  37.58,  36.58,  35.58,  34.58,  33.58,  32.58,  31.58,  30.58,  29.58,  28.58,  27.60,  26.76,  26.34,  26.22,  26.19,  26.18,  26.18,  26.18,  26.18,  26.18, 
     37.61,  36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.64,  25.84,  25.44,  25.23,  25.14,  25.12,  25.12,  25.12,  25.12,  25.12, 
     36.90,  35.90,  34.90,  33.90,  32.90,  31.90,  30.90,  29.90,  28.90,  27.90,  26.90,  25.94,  25.18,  24.74,  24.42,  24.23,  24.17,  24.16,  24.16,  24.16,  24.16, 
     36.33,  35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.33,  27.33,  26.33,  25.37,  24.63,  24.02,  23.55,  23.34,  23.29,  23.28,  23.28,  23.28,  23.28, 
     35.81,  34.81,  33.81,  32.80,  31.81,  30.81,  29.81,  28.81,  27.81,  26.81,  25.80,  24.76,  23.74,  22.93,  22.56,  22.43,  22.36,  22.34,  22.34,  22.34,  22.34, 
     35.29,  34.29,  33.19,  32.19,  31.19,  30.19,  29.19,  28.19,  27.19,  26.18,  25.12,  23.67,  22.44,  21.86,  21.62,  21.33,  21.19,  21.18,  21.17,  21.17,  21.17, 
     34.17,  33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.34,  24.14,  22.31,  21.29,  20.98,  20.64,  20.26,  20.16,  20.15,  20.15,  20.15,  20.15, 
     33.57,  32.51,  31.51,  30.51,  29.51,  28.51,  27.51,  26.51,  25.51,  24.50,  23.01,  21.04,  20.43,  20.22,  19.77,  19.45,  19.38,  19.38,  19.38,  19.38,  19.38, 
     32.83,  31.84,  30.84,  29.84,  28.84,  27.84,  26.84,  25.83,  24.83,  23.77,  21.80,  20.08,  19.84,  19.61,  19.12,  18.89,  18.85,  18.85,  18.85,  18.85,  18.85, 
     43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.07,  31.16,  30.62,  30.50,  30.48,  30.48,  30.48,  30.48,  30.48,  30.48, 
     40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.97,  33.97,  32.97,  31.97,  30.97,  29.98,  29.10,  28.63,  28.53,  28.51,  28.51,  28.51,  28.51,  28.51,  28.51, 
     39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.57,  32.57,  31.57,  30.57,  29.57,  28.59,  27.75,  27.33,  27.21,  27.17,  27.16,  27.16,  27.16,  27.16,  27.16, 
     38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.60,  30.60,  29.60,  28.60,  27.63,  26.83,  26.42,  26.21,  26.11,  26.08,  26.08,  26.08,  26.08,  26.08, 
     37.89,  36.89,  35.89,  34.89,  33.89,  32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.92,  26.15,  25.68,  25.34,  25.16,  25.11,  25.10,  25.10,  25.10,  25.10, 
     37.32,  36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.31,  27.31,  26.31,  25.36,  24.68,  24.34,  24.21,  24.18,  24.18,  24.18,  24.18,  24.18, 
     36.80,  35.80,  34.80,  33.70,  32.77,  31.77,  30.77,  29.77,  28.76,  27.75,  26.72,  25.27,  24.07,  23.55,  23.37,  23.26,  23.19,  23.17,  23.17,  23.17,  23.17, 
     36.26,  35.26,  33.73,  32.75,  31.75,  30.75,  29.75,  28.75,  27.76,  26.78,  25.68,  23.86,  22.87,  22.60,  22.41,  22.11,  21.97,  21.95,  21.95,  21.95,  21.95, 
     34.31,  33.57,  32.57,  31.57,  30.57,  29.57,  28.57,  27.57,  26.58,  25.69,  24.41,  22.58,  21.93,  21.71,  21.30,  20.93,  20.83,  20.82,  20.82,  20.82,  20.82, 
     33.70,  32.63,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.69,  24.81,  23.21,  21.51,  21.16,  20.85,  20.31,  20.02,  19.96,  19.95,  19.95,  19.95,  19.95, 
     32.93,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.96,  25.06,  24.06,  22.06,  20.84,  20.67,  20.18,  19.61,  19.38,  19.34,  19.34,  19.34,  19.34,  19.34, 
     44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.07,  32.16,  31.62,  31.50,  31.48,  31.48,  31.48,  31.48,  31.48,  31.48, 
     41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.97,  33.97,  32.97,  31.97,  30.98,  30.10,  29.63,  29.53,  29.51,  29.51,  29.51,  29.51,  29.51,  29.51, 
     40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.57,  32.57,  31.57,  30.57,  29.59,  28.75,  28.33,  28.20,  28.17,  28.16,  28.16,  28.16,  28.16,  28.16, 
     39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.60,  30.60,  29.60,  28.63,  27.83,  27.42,  27.21,  27.11,  27.08,  27.07,  27.07,  27.07,  27.07, 
     38.89,  37.89,  36.89,  35.89,  34.89,  33.89,  32.89,  31.89,  30.89,  29.89,  28.89,  27.92,  27.11,  26.57,  26.28,  26.12,  26.07,  26.07,  26.07,  26.07,  26.07, 
     38.32,  37.32,  36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.31,  29.30,  28.30,  27.06,  25.88,  25.45,  25.25,  25.15,  25.12,  25.11,  25.11,  25.11,  25.11, 
     37.80,  36.80,  35.80,  34.22,  33.54,  32.53,  31.52,  30.52,  29.52,  28.57,  27.45,  25.54,  24.65,  24.38,  24.27,  24.16,  24.07,  24.05,  24.05,  24.05,  24.05, 
     37.11,  36.11,  33.85,  32.88,  31.89,  30.89,  29.89,  28.90,  27.95,  27.18,  26.01,  24.28,  23.63,  23.44,  23.20,  22.89,  22.77,  22.75,  22.75,  22.75,  22.75, 
     34.33,  33.60,  32.60,  31.60,  30.60,  29.60,  28.60,  27.62,  26.76,  26.09,  24.72,  23.15,  22.71,  22.38,  21.95,  21.71,  21.64,  21.63,  21.63,  21.63,  21.63, 
     33.72,  32.65,  31.65,  30.65,  29.65,  28.65,  27.66,  26.72,  25.92,  25.10,  23.58,  22.23,  21.93,  21.44,  21.06,  20.82,  20.77,  20.77,  20.77,  20.77,  20.77, 
     32.95,  31.96,  30.95,  29.95,  28.96,  27.96,  26.98,  26.09,  25.30,  24.30,  22.46,  21.75,  21.32,  20.81,  20.38,  20.15,  20.11,  20.11,  20.11,  20.11,  20.11, 
     45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.06,  36.06,  35.06,  34.07,  33.16,  32.62,  32.50,  32.48,  32.48,  32.48,  32.48,  32.48,  32.48, 
     42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.97,  33.97,  32.97,  31.98,  31.10,  30.63,  30.53,  30.51,  30.51,  30.51,  30.51,  30.51,  30.51, 
     41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.57,  32.57,  31.57,  30.59,  29.75,  29.33,  29.20,  29.17,  29.16,  29.16,  29.16,  29.16,  29.16, 
     40.60,  39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.60,  30.60,  29.63,  28.83,  28.42,  28.20,  28.10,  28.07,  28.07,  28.07,  28.07,  28.07, 
     39.89,  38.89,  37.89,  36.89,  35.89,  34.89,  33.89,  32.89,  31.88,  30.88,  29.89,  28.88,  27.90,  27.51,  27.25,  27.11,  27.07,  27.06,  27.06,  27.06,  27.06, 
     39.32,  38.32,  37.32,  36.32,  35.32,  34.32,  33.32,  32.30,  31.29,  30.29,  29.21,  27.43,  26.67,  26.35,  26.20,  26.12,  26.09,  26.08,  26.08,  26.08,  26.08, 
     38.80,  37.80,  36.80,  34.34,  33.83,  32.81,  31.79,  30.81,  29.93,  29.21,  27.81,  26.11,  25.46,  25.30,  25.20,  25.07,  24.96,  24.94,  24.94,  24.94,  24.94, 
     37.52,  36.52,  33.87,  32.90,  31.91,  30.91,  29.91,  28.98,  28.23,  27.83,  26.39,  24.92,  24.47,  24.22,  23.89,  23.69,  23.61,  23.60,  23.60,  23.60,  23.60, 
     34.33,  33.61,  32.61,  31.60,  30.60,  29.60,  28.62,  27.76,  27.15,  26.46,  25.15,  23.86,  23.40,  22.97,  22.77,  22.60,  22.54,  22.53,  22.53,  22.53,  22.53, 
     33.72,  32.65,  31.65,  30.65,  29.65,  28.66,  27.72,  26.93,  26.22,  25.28,  23.98,  22.99,  22.52,  22.21,  21.96,  21.74,  21.69,  21.68,  21.68,  21.68,  21.68, 
     32.95,  31.96,  30.96,  29.96,  28.96,  27.98,  27.10,  26.32,  25.44,  24.45,  22.98,  22.42,  21.93,  21.68,  21.30,  21.06,  21.03,  21.02,  21.02,  21.02,  21.02, 
     46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.06,  36.06,  35.07,  34.16,  33.62,  33.50,  33.48,  33.48,  33.48,  33.48,  33.48,  33.48, 
     43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.97,  33.97,  32.98,  32.10,  31.63,  31.53,  31.51,  31.51,  31.51,  31.51,  31.51,  31.51, 
     42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.57,  32.57,  31.59,  30.75,  30.33,  30.20,  30.17,  30.16,  30.16,  30.16,  30.16,  30.16, 
     41.60,  40.60,  39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.60,  30.63,  29.82,  29.40,  29.20,  29.10,  29.07,  29.07,  29.07,  29.07,  29.07, 
     40.89,  39.89,  38.89,  37.89,  36.89,  35.89,  34.89,  33.88,  32.88,  31.88,  30.88,  29.66,  28.78,  28.44,  28.23,  28.10,  28.06,  28.05,  28.05,  28.05,  28.05, 
     40.32,  39.32,  38.32,  37.32,  36.32,  35.32,  34.31,  33.22,  32.23,  31.27,  29.80,  28.10,  27.45,  27.29,  27.17,  27.09,  27.06,  27.05,  27.05,  27.05,  27.05, 
     39.76,  38.76,  37.76,  34.35,  33.87,  32.85,  31.86,  31.00,  30.35,  30.05,  28.16,  26.82,  26.36,  26.22,  26.09,  25.90,  25.80,  25.78,  25.78,  25.78,  25.78, 
     37.59,  36.59,  33.87,  32.90,  31.91,  30.92,  29.98,  29.23,  28.90,  28.35,  26.92,  25.68,  25.25,  24.88,  24.69,  24.59,  24.52,  24.51,  24.51,  24.51,  24.51, 
     34.33,  33.61,  32.61,  31.60,  30.60,  29.62,  28.76,  28.15,  27.55,  26.61,  25.49,  24.48,  24.01,  23.79,  23.68,  23.53,  23.48,  23.47,  23.47,  23.47,  23.47, 
     33.72,  32.65,  31.65,  30.65,  29.66,  28.72,  27.93,  27.23,  26.33,  25.35,  24.29,  23.61,  23.23,  23.12,  22.89,  22.67,  22.62,  22.62,  22.62,  22.62,  22.62, 
     32.95,  31.96,  30.96,  29.96,  28.98,  28.10,  27.32,  26.46,  25.49,  24.52,  23.52,  22.96,  22.79,  22.62,  22.28,  22.04,  22.00,  22.00,  22.00,  22.00,  22.00, 
     47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.06,  36.07,  35.16,  34.62,  34.50,  34.48,  34.48,  34.48,  34.48,  34.48,  34.48, 
     44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.97,  33.98,  33.10,  32.63,  32.53,  32.51,  32.51,  32.51,  32.51,  32.51,  32.51, 
     43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.57,  32.59,  31.75,  31.33,  31.20,  31.17,  31.16,  31.16,  31.16,  31.16,  31.16, 
     42.60,  41.60,  40.60,  39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.63,  30.79,  30.40,  30.20,  30.10,  30.07,  30.07,  30.07,  30.07,  30.07, 
     41.89,  40.89,  39.89,  38.89,  37.89,  36.89,  35.88,  34.88,  33.88,  32.88,  31.82,  30.25,  29.73,  29.41,  29.23,  29.09,  29.05,  29.04,  29.04,  29.04,  29.04, 
     41.32,  40.32,  39.32,  38.32,  37.32,  36.31,  35.30,  33.92,  33.13,  32.18,  30.10,  28.89,  28.41,  28.24,  28.15,  28.06,  28.02,  28.00,  28.00,  28.00,  28.00, 
     40.51,  39.51,  38.51,  34.36,  33.88,  32.88,  32.00,  31.37,  31.26,  30.60,  28.89,  27.63,  27.26,  27.12,  26.90,  26.77,  26.70,  26.69,  26.69,  26.69,  26.69, 
     37.60,  36.60,  33.87,  32.90,  31.92,  30.98,  30.23,  29.90,  29.51,  28.58,  27.42,  26.38,  25.90,  25.61,  25.58,  25.51,  25.45,  25.44,  25.44,  25.44,  25.44, 
     34.33,  33.61,  32.61,  31.60,  30.62,  29.76,  29.15,  28.56,  27.66,  26.68,  25.73,  25.02,  24.59,  24.70,  24.62,  24.48,  24.42,  24.42,  24.42,  24.42,  24.42, 
     33.72,  32.65,  31.65,  30.66,  29.72,  28.93,  28.23,  27.34,  26.37,  25.44,  24.59,  24.09,  24.14,  24.05,  23.86,  23.65,  23.61,  23.60,  23.60,  23.60,  23.60, 
     32.95,  31.96,  30.96,  29.98,  29.10,  28.33,  27.46,  26.50,  25.53,  24.68,  23.87,  23.45,  23.72,  23.61,  23.27,  23.04,  23.00,  22.99,  22.99,  22.99,  22.99, 
     48.06,  47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.06,  37.07,  36.16,  35.62,  35.50,  35.48,  35.48,  35.48,  35.48,  35.48,  35.48, 
     45.97,  44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.97,  34.98,  34.10,  33.63,  33.53,  33.51,  33.51,  33.51,  33.51,  33.51,  33.51, 
     44.57,  43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.57,  33.59,  32.75,  32.33,  32.20,  32.17,  32.16,  32.16,  32.16,  32.16,  32.16, 
     43.60,  42.60,  41.60,  40.60,  39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.60,  32.60,  31.78,  31.40,  31.20,  31.10,  31.07,  31.07,  31.07,  31.07,  31.07, 
     42.89,  41.89,  40.89,  39.89,  38.89,  37.88,  36.88,  35.88,  34.88,  33.87,  32.48,  31.21,  30.61,  30.41,  30.22,  30.09,  30.05,  30.04,  30.04,  30.04,  30.04, 
     42.32,  41.32,  40.32,  39.32,  38.31,  37.30,  36.30,  34.59,  34.12,  32.75,  30.88,  29.71,  29.32,  29.23,  29.13,  29.02,  28.92,  28.89,  28.89,  28.89,  28.89, 
     40.78,  39.78,  38.78,  34.36,  33.90,  33.02,  32.37,  32.28,  32.17,  30.98,  29.53,  28.52,  28.17,  27.93,  27.73,  27.68,  27.62,  27.61,  27.61,  27.61,  27.61, 
     37.60,  36.60,  33.87,  32.91,  31.98,  31.23,  30.90,  30.52,  29.69,  28.70,  27.71,  26.89,  26.52,  26.50,  26.50,  26.44,  26.39,  26.38,  26.38,  26.38,  26.38, 
     34.33,  33.61,  32.61,  31.62,  30.76,  30.15,  29.56,  28.67,  27.69,  26.78,  25.95,  25.56,  25.32,  25.64,  25.58,  25.44,  25.38,  25.38,  25.37,  25.37,  25.37, 
     33.72,  32.65,  31.66,  30.72,  29.93,  29.23,  28.34,  27.37,  26.44,  25.58,  24.76,  24.28,  25.06,  25.03,  24.86,  24.65,  24.60,  24.59,  24.59,  24.59,  24.59, 
     32.95,  31.96,  30.98,  30.10,  29.33,  28.46,  27.50,  26.54,  25.67,  24.81,  23.95,  23.63,  24.71,  24.61,  24.27,  24.03,  23.99,  23.99,  23.99,  23.99,  23.99, 
     49.06,  48.06,  47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.06,  38.07,  37.16,  36.62,  36.50,  36.48,  36.48,  36.48,  36.48,  36.48,  36.48, 
     46.97,  45.97,  44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.97,  35.98,  35.10,  34.63,  34.53,  34.51,  34.51,  34.51,  34.51,  34.51,  34.51, 
     45.57,  44.57,  43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.57,  34.59,  33.75,  33.33,  33.20,  33.17,  33.16,  33.16,  33.16,  33.16,  33.16, 
     44.60,  43.60,  42.60,  41.60,  40.60,  39.60,  38.60,  37.60,  36.60,  35.60,  34.60,  33.51,  32.78,  32.40,  32.19,  32.10,  32.07,  32.07,  32.07,  32.07,  32.07, 
     43.89,  42.89,  41.89,  40.89,  39.88,  38.88,  37.88,  36.88,  35.88,  34.81,  33.04,  32.13,  31.60,  31.39,  31.21,  31.09,  31.05,  31.04,  31.04,  31.04,  31.04, 
     43.32,  42.32,  41.32,  40.31,  39.30,  38.30,  37.30,  35.58,  35.05,  33.06,  31.71,  30.66,  30.31,  30.20,  30.08,  29.90,  29.81,  29.80,  29.79,  29.79,  29.79, 
     40.81,  39.81,  38.81,  34.38,  34.05,  33.38,  33.28,  33.20,  32.70,  31.56,  30.32,  29.33,  28.94,  28.68,  28.64,  28.62,  28.56,  28.55,  28.55,  28.55,  28.55, 
     37.60,  36.60,  33.88,  32.97,  32.23,  31.90,  31.52,  30.70,  29.72,  28.78,  27.93,  27.42,  26.97,  27.43,  27.45,  27.40,  27.35,  27.34,  27.34,  27.34,  27.34, 
     34.33,  33.61,  32.63,  31.76,  31.15,  30.56,  29.67,  28.70,  27.78,  26.90,  26.08,  25.72,  26.28,  26.59,  26.55,  26.42,  26.37,  26.36,  26.36,  26.36,  26.36, 
     33.72,  32.66,  31.72,  30.93,  30.23,  29.34,  28.37,  27.44,  26.57,  25.62,  24.79,  24.33,  26.03,  26.02,  25.85,  25.64,  25.59,  25.59,  25.59,  25.59,  25.59, 
     32.95,  31.98,  31.10,  30.33,  29.46,  28.50,  27.54,  26.66,  25.73,  24.90,  23.97,  23.66,  25.71,  25.61,  25.26,  25.03,  24.99,  24.99,  24.99,  24.99,  24.99, 
     50.06,  49.06,  48.06,  47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.06,  39.07,  38.16,  37.62,  37.50,  37.48,  37.48,  37.48,  37.48,  37.48,  37.48, 
     47.97,  46.97,  45.97,  44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.97,  36.98,  36.10,  35.63,  35.53,  35.51,  35.51,  35.51,  35.51,  35.51,  35.51, 
     46.57,  45.57,  44.57,  43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.57,  35.59,  34.75,  34.33,  34.20,  34.17,  34.16,  34.16,  34.16,  34.16,  34.16, 
     45.60,  44.60,  43.60,  42.60,  41.60,  40.60,  39.60,  38.60,  37.60,  36.60,  35.55,  34.50,  33.78,  33.40,  33.19,  33.10,  33.07,  33.07,  33.07,  33.07,  33.07, 
     44.89,  43.89,  42.89,  41.88,  40.88,  39.88,  38.88,  37.88,  36.87,  35.45,  34.01,  32.97,  32.60,  32.39,  32.21,  32.09,  32.04,  32.03,  32.03,  32.03,  32.03, 
     44.31,  43.31,  42.30,  41.29,  40.29,  39.29,  38.29,  36.58,  35.68,  33.84,  32.44,  31.55,  31.28,  31.17,  30.99,  30.81,  30.77,  30.75,  30.75,  30.75,  30.75, 
     40.82,  39.82,  38.82,  34.52,  34.41,  34.29,  34.20,  33.78,  32.85,  31.85,  30.81,  30.00,  29.68,  29.45,  29.56,  29.54,  29.48,  29.47,  29.47,  29.47,  29.47, 
     37.60,  36.60,  33.94,  33.22,  32.90,  32.52,  31.70,  30.73,  29.79,  28.92,  28.09,  27.87,  27.56,  28.37,  28.40,  28.36,  28.31,  28.30,  28.30,  28.30,  28.30, 
     34.33,  33.63,  32.77,  32.15,  31.56,  30.67,  29.70,  28.78,  27.89,  26.92,  26.11,  25.75,  27.26,  27.56,  27.54,  27.41,  27.36,  27.35,  27.35,  27.35,  27.35, 
     33.73,  32.71,  31.93,  31.23,  30.34,  29.37,  28.44,  27.57,  26.59,  25.63,  24.79,  24.33,  27.03,  27.02,  26.85,  26.64,  26.59,  26.59,  26.59,  26.59,  26.59, 
     32.97,  32.10,  31.33,  30.46,  29.50,  28.54,  27.66,  26.72,  25.74,  24.92,  23.97,  23.66,  26.71,  26.60,  26.26,  26.03,  25.99,  25.99,  25.99,  25.99,  25.99, 
     51.06,  50.06,  49.06,  48.06,  47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.06,  40.07,  39.16,  38.62,  38.50,  38.48,  38.48,  38.48,  38.48,  38.48,  38.48, 
     48.97,  47.97,  46.97,  45.97,  44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.97,  37.98,  37.10,  36.63,  36.53,  36.51,  36.51,  36.51,  36.51,  36.51,  36.51, 
     47.57,  46.57,  45.57,  44.57,  43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.57,  36.59,  35.75,  35.33,  35.20,  35.17,  35.16,  35.16,  35.16,  35.16,  35.16, 
     46.60,  45.60,  44.60,  43.60,  42.60,  41.60,  40.60,  39.60,  38.60,  37.59,  36.43,  35.50,  34.77,  34.39,  34.19,  34.10,  34.07,  34.07,  34.07,  34.07,  34.07, 
     45.89,  44.89,  43.88,  42.88,  41.88,  40.88,  39.88,  38.88,  37.81,  36.01,  34.94,  33.96,  33.58,  33.38,  33.21,  33.08,  33.03,  33.00,  32.99,  32.99,  32.99, 
     45.23,  44.23,  43.22,  42.22,  41.22,  40.23,  39.23,  37.54,  36.01,  34.72,  33.42,  32.54,  32.24,  32.12,  31.90,  31.78,  31.73,  31.71,  31.71,  31.71,  31.71, 
     40.82,  39.82,  38.82,  34.89,  35.32,  35.21,  34.78,  33.92,  32.93,  32.00,  31.12,  30.46,  30.21,  30.39,  30.48,  30.47,  30.43,  30.42,  30.42,  30.42,  30.42, 
     37.60,  36.60,  34.19,  33.90,  33.52,  32.70,  31.73,  30.79,  29.92,  28.95,  28.16,  27.96,  28.54,  29.33,  29.37,  29.33,  29.28,  29.27,  29.27,  29.27,  29.27, 
     34.36,  33.77,  33.16,  32.56,  31.67,  30.70,  29.78,  28.89,  27.91,  26.92,  26.11,  25.76,  28.24,  28.55,  28.53,  28.41,  28.35,  28.35,  28.35,  28.35,  28.35, 
     33.78,  32.93,  32.23,  31.34,  30.37,  29.44,  28.57,  27.59,  26.60,  25.63,  24.79,  24.33,  28.02,  28.01,  27.85,  27.64,  27.59,  27.59,  27.59,  27.59,  27.59, 
     33.09,  32.33,  31.46,  30.50,  29.54,  28.66,  27.72,  26.73,  25.75,  24.93,  23.97,  23.66,  27.70,  27.60,  27.26,  27.03,  26.99,  26.99,  26.99,  26.99,  26.99, 
     52.06,  51.06,  50.06,  49.06,  48.06,  47.06,  46.06,  45.06,  44.06,  43.06,  42.06,  41.07,  40.16,  39.62,  39.50,  39.48,  39.48,  39.48,  39.48,  39.48,  39.48, 
     49.97,  48.97,  47.97,  46.97,  45.97,  44.97,  43.97,  42.97,  41.97,  40.97,  39.97,  38.98,  38.10,  37.63,  37.53,  37.51,  37.51,  37.51,  37.51,  37.51,  37.51, 
     48.57,  47.57,  46.57,  45.57,  44.57,  43.57,  42.57,  41.57,  40.57,  39.57,  38.57,  37.59,  36.75,  36.33,  36.20,  36.17,  36.16,  36.16,  36.16,  36.16,  36.16, 
     47.60,  46.60,  45.60,  44.60,  43.60,  42.60,  41.60,  40.60,  39.60,  38.54,  37.42,  36.48,  35.77,  35.39,  35.19,  35.10,  35.07,  35.07,  35.07,  35.07,  35.07, 
     46.89,  45.88,  44.88,  43.88,  42.88,  41.88,  40.88,  39.87,  38.45,  36.98,  35.78,  34.96,  34.56,  34.38,  34.20,  34.07,  33.95,  33.87,  33.86,  33.86,  33.86, 
     45.80,  44.80,  43.80,  42.83,  41.85,  40.88,  39.88,  38.28,  36.74,  35.49,  34.35,  33.50,  33.21,  33.03,  32.82,  32.74,  32.70,  32.69,  32.69,  32.69,  32.69, 
     40.82,  39.82,  38.84,  35.81,  36.24,  35.79,  34.93,  33.95,  33.02,  32.14,  31.28,  31.10,  30.74,  31.34,  31.42,  31.42,  31.39,  31.37,  31.37,  31.37,  31.37, 
     37.60,  36.61,  34.87,  34.52,  33.70,  32.73,  31.79,  30.92,  29.95,  28.96,  28.17,  27.98,  29.53,  30.30,  30.34,  30.31,  30.26,  30.25,  30.25,  30.25,  30.25, 
     34.51,  34.16,  33.56,  32.67,  31.70,  30.78,  29.89,  28.91,  27.91,  26.92,  26.11,  25.76,  29.23,  29.54,  29.53,  29.40,  29.35,  29.34,  29.34,  29.34,  29.34, 
     33.98,  33.23,  32.34,  31.37,  30.44,  29.57,  28.59,  27.59,  26.60,  25.63,  24.79,  24.33,  29.02,  29.01,  28.85,  28.64,  28.59,  28.59,  28.59,  28.59,  28.59, 
     33.32,  32.46,  31.50,  30.54,  29.66,  28.72,  27.73,  26.73,  25.75,  24.93,  23.98,  23.66,  28.70,  28.60,  28.26,  28.03,  27.99,  27.99,  27.99,  27.99,  27.99}; 

  setup_cool_interp_grid_(&my_rates->LH2, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_HD(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {16,  16.0, 1.0}, // log10(number-density like)
    {11,   1.6, 0.2}, // log10(temperature)
    {21, -12.0, 1.0} // log10(H2 number density)
  };

  double L[] =
     {37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.76,  23.76,  22.83,  22.20,  22.02,  22.00,  21.99,  21.99, 
      37.14,  36.14,  35.14,  34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.15,  22.24,  21.68,  21.47,  21.41,  21.40,  21.40, 
      36.69,  35.69,  34.69,  33.69,  32.69,  31.69,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.69,  23.69,  22.70,  21.83,  21.26,  20.93,  20.84,  20.83,  20.83, 
      36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.32,  22.34,  21.48,  20.84,  20.45,  20.33,  20.31,  20.31, 
      36.00,  35.00,  34.00,  33.00,  32.00,  31.00,  30.00,  29.00,  28.00,  27.00,  26.00,  25.00,  24.00,  23.00,  22.03,  21.16,  20.47,  20.06,  19.94,  19.92,  19.92, 
      35.71,  34.71,  33.71,  32.71,  31.71,  30.71,  29.71,  28.71,  27.71,  26.71,  25.71,  24.71,  23.71,  22.71,  21.74,  20.86,  20.15,  19.78,  19.69,  19.68,  19.68, 
      35.44,  34.44,  33.44,  32.44,  31.44,  30.44,  29.44,  28.44,  27.44,  26.44,  25.44,  24.44,  23.44,  22.44,  21.47,  20.57,  19.89,  19.60,  19.55,  19.54,  19.54, 
      35.19,  34.19,  33.19,  32.19,  31.19,  30.19,  29.19,  28.19,  27.19,  26.19,  25.19,  24.19,  23.19,  22.19,  21.20,  20.31,  19.70,  19.50,  19.47,  19.46,  19.46, 
      34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.95,  25.95,  24.95,  23.95,  22.95,  21.95,  20.95,  20.07,  19.57,  19.44,  19.42,  19.42,  19.42, 
      34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.73,  24.73,  23.73,  22.72,  21.71,  20.69,  19.87,  19.48,  19.40,  19.39,  19.39,  19.39, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.50,  24.50,  23.50,  22.50,  21.47,  20.43,  19.71,  19.43,  19.38,  19.37,  19.37,  19.37, 
      37.76,  36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.76,  29.76,  28.76,  27.76,  26.76,  25.76,  24.76,  23.76,  22.83,  22.20,  22.02,  22.00,  21.99,  21.99, 
      37.14,  36.14,  35.14,  34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.15,  22.24,  21.68,  21.47,  21.41,  21.40,  21.40, 
      36.69,  35.69,  34.69,  33.69,  32.69,  31.69,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.69,  23.69,  22.71,  21.83,  21.26,  20.93,  20.84,  20.83,  20.83, 
      36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.32,  22.34,  21.48,  20.84,  20.45,  20.33,  20.31,  20.31, 
      36.00,  35.00,  34.00,  33.00,  32.00,  31.00,  30.00,  29.00,  28.00,  27.00,  26.00,  25.00,  24.00,  23.00,  22.03,  21.17,  20.47,  20.06,  19.94,  19.92,  19.92, 
      35.71,  34.71,  33.71,  32.71,  31.71,  30.71,  29.71,  28.71,  27.71,  26.71,  25.71,  24.71,  23.71,  22.71,  21.74,  20.86,  20.15,  19.78,  19.69,  19.68,  19.68, 
      35.44,  34.44,  33.44,  32.44,  31.44,  30.44,  29.44,  28.44,  27.44,  26.44,  25.44,  24.44,  23.44,  22.44,  21.47,  20.57,  19.89,  19.60,  19.55,  19.54,  19.54, 
      35.19,  34.19,  33.19,  32.19,  31.19,  30.19,  29.19,  28.19,  27.19,  26.19,  25.19,  24.19,  23.19,  22.19,  21.21,  20.31,  19.70,  19.50,  19.47,  19.46,  19.46, 
      34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.95,  27.95,  26.95,  25.95,  24.95,  23.95,  22.95,  21.95,  20.95,  20.07,  19.57,  19.44,  19.42,  19.42,  19.42, 
      34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.73,  24.73,  23.73,  22.72,  21.71,  20.69,  19.87,  19.48,  19.40,  19.39,  19.39,  19.39, 
      34.51,  33.51,  32.51,  31.51,  30.51,  29.51,  28.51,  27.51,  26.51,  25.51,  24.51,  23.51,  22.50,  21.47,  20.43,  19.71,  19.43,  19.38,  19.37,  19.37,  19.37, 
      37.77,  36.77,  35.77,  34.77,  33.77,  32.77,  31.77,  30.77,  29.77,  28.77,  27.77,  26.77,  25.77,  24.77,  23.78,  22.84,  22.22,  22.04,  22.01,  22.01,  22.01, 
      37.16,  36.16,  35.16,  34.16,  33.16,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.16,  25.16,  24.16,  23.17,  22.26,  21.69,  21.48,  21.42,  21.41,  21.41, 
      36.70,  35.70,  34.70,  33.70,  32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.71,  22.72,  21.84,  21.26,  20.94,  20.85,  20.83,  20.83, 
      36.34,  35.34,  34.34,  33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.34,  22.36,  21.49,  20.85,  20.46,  20.33,  20.32,  20.31, 
      36.01,  35.01,  34.01,  33.01,  32.01,  31.01,  30.01,  29.01,  28.01,  27.01,  26.01,  25.01,  24.01,  23.02,  22.04,  21.17,  20.48,  20.06,  19.94,  19.93,  19.93, 
      35.72,  34.72,  33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.75,  20.87,  20.16,  19.78,  19.69,  19.68,  19.68, 
      35.45,  34.45,  33.45,  32.45,  31.45,  30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.45,  22.45,  21.47,  20.58,  19.90,  19.60,  19.55,  19.54,  19.54, 
      35.20,  34.20,  33.20,  32.20,  31.20,  30.20,  29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.21,  20.31,  19.70,  19.50,  19.47,  19.46,  19.46, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.96,  25.96,  24.96,  23.96,  22.96,  21.96,  20.95,  20.08,  19.57,  19.44,  19.42,  19.42,  19.42, 
      34.74,  33.74,  32.74,  31.74,  30.74,  29.74,  28.74,  27.74,  26.74,  25.74,  24.74,  23.74,  22.73,  21.72,  20.69,  19.88,  19.48,  19.40,  19.39,  19.39,  19.39, 
      34.51,  33.51,  32.51,  31.51,  30.51,  29.51,  28.51,  27.51,  26.51,  25.51,  24.51,  23.51,  22.51,  21.48,  20.44,  19.72,  19.43,  19.38,  19.37,  19.37,  19.37, 
      37.93,  36.93,  35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.94,  23.00,  22.36,  22.17,  22.14,  22.14,  22.14, 
      37.31,  36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.31,  27.31,  26.31,  25.31,  24.31,  23.32,  22.40,  21.80,  21.56,  21.49,  21.49,  21.48, 
      36.85,  35.85,  34.85,  33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.85,  23.85,  22.87,  21.95,  21.33,  20.99,  20.89,  20.88,  20.87, 
      36.48,  35.48,  34.48,  33.48,  32.48,  31.48,  30.48,  29.48,  28.48,  27.48,  26.48,  25.48,  24.48,  23.48,  22.49,  21.57,  20.90,  20.49,  20.37,  20.35,  20.34, 
      36.14,  35.14,  34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.14,  22.15,  21.23,  20.52,  20.09,  19.97,  19.95,  19.95, 
      35.83,  34.83,  33.83,  32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.83,  21.83,  20.91,  20.19,  19.80,  19.71,  19.70,  19.70, 
      35.55,  34.55,  33.55,  32.55,  31.55,  30.55,  29.55,  28.55,  27.55,  26.55,  25.55,  24.55,  23.55,  22.55,  21.54,  20.62,  19.92,  19.62,  19.56,  19.56,  19.56, 
      35.29,  34.29,  33.29,  32.29,  31.29,  30.29,  29.29,  28.29,  27.29,  26.29,  25.29,  24.29,  23.29,  22.28,  21.26,  20.35,  19.73,  19.51,  19.48,  19.47,  19.47, 
      35.05,  34.05,  33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.05,  23.04,  22.03,  20.99,  20.10,  19.59,  19.44,  19.43,  19.42,  19.42, 
      34.81,  33.81,  32.81,  31.81,  30.81,  29.81,  28.81,  27.81,  26.81,  25.81,  24.81,  23.81,  22.81,  21.78,  20.73,  19.90,  19.50,  19.41,  19.40,  19.39,  19.39, 
      34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.59,  22.58,  21.53,  20.47,  19.74,  19.44,  19.38,  19.38,  19.38,  19.38, 
      38.69,  37.69,  36.69,  35.69,  34.69,  33.69,  32.69,  31.69,  30.69,  29.69,  28.69,  27.69,  26.69,  25.69,  24.70,  23.75,  23.08,  22.84,  22.78,  22.78,  22.77, 
      38.05,  37.05,  36.05,  35.05,  34.05,  33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.06,  23.10,  22.37,  22.03,  21.94,  21.92,  21.92, 
      37.53,  36.53,  35.53,  34.53,  33.53,  32.53,  31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.52,  22.53,  21.77,  21.36,  21.23,  21.21,  21.20, 
      37.03,  36.03,  35.03,  34.03,  33.03,  32.03,  31.03,  30.03,  29.03,  28.03,  27.03,  26.03,  25.03,  24.03,  23.02,  22.04,  21.26,  20.78,  20.62,  20.59,  20.59, 
      36.57,  35.57,  34.57,  33.57,  32.57,  31.57,  30.57,  29.57,  28.57,  27.57,  26.57,  25.57,  24.57,  23.57,  22.58,  21.63,  20.82,  20.32,  20.18,  20.16,  20.16, 
      36.18,  35.18,  34.18,  33.18,  32.18,  31.18,  30.18,  29.18,  28.18,  27.18,  26.18,  25.18,  24.18,  23.18,  22.21,  21.26,  20.45,  20.00,  19.89,  19.88,  19.88, 
      35.84,  34.84,  33.84,  32.84,  31.84,  30.84,  29.84,  28.84,  27.84,  26.84,  25.84,  24.84,  23.84,  22.85,  21.88,  20.91,  20.15,  19.78,  19.70,  19.69,  19.69, 
      35.54,  34.54,  33.54,  32.54,  31.54,  30.54,  29.54,  28.54,  27.54,  26.54,  25.54,  24.54,  23.54,  22.55,  21.58,  20.61,  19.92,  19.63,  19.57,  19.57,  19.57, 
      35.27,  34.27,  33.27,  32.27,  31.27,  30.27,  29.27,  28.27,  27.27,  26.27,  25.27,  24.27,  23.27,  22.28,  21.29,  20.35,  19.75,  19.53,  19.49,  19.49,  19.49, 
      35.02,  34.02,  33.02,  32.02,  31.02,  30.02,  29.02,  28.02,  27.02,  26.02,  25.02,  24.02,  23.02,  22.03,  21.00,  20.13,  19.62,  19.46,  19.44,  19.44,  19.44, 
      34.77,  33.77,  32.77,  31.77,  30.77,  29.77,  28.77,  27.77,  26.77,  25.77,  24.77,  23.77,  22.77,  21.77,  20.72,  19.94,  19.53,  19.42,  19.40,  19.40,  19.40, 
      39.67,  38.67,  37.67,  36.67,  35.67,  34.67,  33.67,  32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.68,  24.72,  24.03,  23.76,  23.68,  23.68,  23.67, 
      38.88,  37.88,  36.88,  35.88,  34.88,  33.88,  32.88,  31.88,  30.88,  29.88,  28.88,  27.88,  26.88,  25.89,  24.90,  24.00,  23.32,  22.95,  22.83,  22.81,  22.81, 
      38.04,  37.04,  36.04,  35.04,  34.04,  33.04,  32.04,  31.04,  30.04,  29.04,  28.04,  27.04,  26.05,  25.06,  24.16,  23.43,  22.70,  22.20,  22.01,  21.98,  21.98, 
      37.31,  36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.31,  27.31,  26.31,  25.32,  24.36,  23.62,  22.93,  22.12,  21.62,  21.44,  21.41,  21.41, 
      36.74,  35.74,  34.74,  33.74,  32.74,  31.74,  30.74,  29.74,  28.74,  27.74,  26.74,  25.74,  24.75,  23.83,  23.24,  22.44,  21.67,  21.18,  21.02,  21.00,  20.99, 
      36.29,  35.29,  34.29,  33.29,  32.29,  31.29,  30.29,  29.29,  28.29,  27.29,  26.29,  25.29,  24.31,  23.45,  22.91,  22.02,  21.31,  20.83,  20.69,  20.67,  20.67, 
      35.92,  34.92,  33.92,  32.92,  31.92,  30.92,  29.92,  28.92,  27.92,  26.92,  25.92,  24.93,  23.95,  23.17,  22.56,  21.71,  20.99,  20.53,  20.41,  20.39,  20.39, 
      35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.66,  22.97,  22.21,  21.44,  20.71,  20.27,  20.17,  20.15,  20.15, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.33,  27.33,  26.33,  25.33,  24.33,  23.40,  22.78,  21.91,  21.18,  20.44,  20.05,  19.96,  19.95,  19.94, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.07,  28.07,  27.07,  26.07,  25.07,  24.08,  23.18,  22.57,  21.67,  20.93,  20.18,  19.86,  19.78,  19.77,  19.77, 
      34.82,  33.82,  32.82,  31.82,  30.82,  29.82,  28.82,  27.82,  26.82,  25.82,  24.82,  23.84,  22.99,  22.31,  21.47,  20.70,  19.99,  19.71,  19.65,  19.64,  19.64, 
      40.53,  39.53,  38.53,  37.53,  36.53,  35.53,  34.53,  33.53,  32.53,  31.53,  30.53,  29.53,  28.53,  27.53,  26.58,  25.70,  25.03,  24.76,  24.68,  24.67,  24.67, 
      39.26,  38.26,  37.26,  36.26,  35.26,  34.26,  33.26,  32.26,  31.26,  30.26,  29.26,  28.26,  27.26,  26.33,  25.69,  24.99,  24.31,  23.91,  23.75,  23.72,  23.71, 
      38.16,  37.16,  36.16,  35.16,  34.16,  33.16,  32.16,  31.16,  30.16,  29.16,  28.16,  27.16,  26.18,  25.38,  25.05,  24.39,  23.65,  23.19,  23.00,  22.96,  22.96, 
      37.35,  36.35,  35.35,  34.35,  33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.36,  26.36,  25.41,  24.80,  24.56,  23.80,  23.11,  22.62,  22.44,  22.41,  22.40, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.76,  29.76,  28.76,  27.76,  26.76,  25.77,  24.87,  24.46,  24.10,  23.35,  22.67,  22.18,  22.01,  21.99,  21.99, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.31,  27.31,  26.31,  25.32,  24.49,  24.22,  23.66,  23.00,  22.30,  21.82,  21.68,  21.66,  21.66, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.94,  24.97,  24.22,  24.01,  23.32,  22.70,  21.98,  21.53,  21.40,  21.39,  21.39, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.62,  24.67,  24.05,  23.78,  23.07,  22.44,  21.69,  21.26,  21.16,  21.14,  21.14, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.33,  27.33,  26.33,  25.34,  24.42,  23.93,  23.52,  22.85,  22.18,  21.41,  21.02,  20.93,  20.91,  20.91, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.07,  28.07,  27.07,  26.07,  25.09,  24.20,  23.83,  23.27,  22.65,  21.90,  21.11,  20.79,  20.71,  20.70,  20.70, 
      34.82,  33.82,  32.82,  31.82,  30.82,  29.82,  28.82,  27.82,  26.82,  25.83,  24.85,  24.02,  23.73,  23.05,  22.47,  21.63,  20.74,  20.57,  20.50,  20.49,  20.49, 
      40.95,  39.95,  38.95,  37.95,  36.95,  35.95,  34.95,  33.95,  32.95,  31.95,  30.95,  29.95,  28.97,  28.11,  27.53,  26.70,  26.03,  25.75,  25.66,  25.65,  25.65, 
      39.32,  38.32,  37.32,  36.32,  35.32,  34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.33,  27.41,  26.95,  26.69,  25.98,  25.29,  24.90,  24.74,  24.71,  24.71, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.20,  26.41,  26.25,  26.03,  25.32,  24.64,  24.17,  23.98,  23.95,  23.95, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.36,  31.36,  30.36,  29.36,  28.36,  27.37,  26.42,  25.83,  25.77,  25.43,  24.78,  24.10,  23.60,  23.41,  23.38,  23.38, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.76,  29.76,  28.76,  27.77,  26.78,  25.87,  25.49,  25.42,  24.97,  24.34,  23.65,  23.15,  22.98,  22.96,  22.96, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.31,  27.31,  26.33,  25.49,  25.27,  25.13,  24.61,  23.99,  23.27,  22.79,  22.65,  22.63,  22.63, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.94,  25.97,  25.23,  25.11,  24.88,  24.32,  23.69,  22.95,  22.50,  22.38,  22.36,  22.36, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.67,  25.06,  24.99,  24.67,  24.06,  23.41,  22.66,  22.24,  22.14,  22.12,  22.12, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.33,  27.33,  26.34,  25.42,  24.95,  24.87,  24.46,  23.84,  23.15,  22.39,  22.01,  21.91,  21.90,  21.90, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.07,  28.07,  27.07,  26.09,  25.20,  24.88,  24.75,  24.25,  23.64,  22.87,  22.10,  21.78,  21.70,  21.69,  21.69, 
      34.82,  33.82,  32.82,  31.82,  30.82,  29.82,  28.82,  27.83,  26.83,  25.85,  25.03,  24.82,  24.62,  24.05,  23.45,  22.60,  21.73,  21.56,  21.49,  21.48,  21.48, 
      41.02,  40.02,  39.02,  38.02,  37.02,  36.02,  35.02,  34.02,  33.02,  32.02,  31.03,  30.05,  29.24,  28.98,  28.53,  27.70,  27.03,  26.73,  26.64,  26.63,  26.63, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.34,  28.42,  27.99,  27.95,  27.68,  26.95,  26.26,  25.85,  25.69,  25.66,  25.66, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.17,  31.17,  30.17,  29.17,  28.20,  27.42,  27.27,  27.24,  26.95,  26.30,  25.58,  25.08,  24.89,  24.86,  24.85, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.36,  31.36,  30.36,  29.36,  28.37,  27.42,  26.83,  26.80,  26.74,  26.39,  25.74,  25.00,  24.46,  24.26,  24.23,  24.22, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.76,  29.76,  28.77,  27.78,  26.87,  26.49,  26.48,  26.38,  25.95,  25.29,  24.52,  23.96,  23.79,  23.76,  23.76, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.31,  27.33,  26.49,  26.27,  26.26,  26.11,  25.59,  24.90,  24.10,  23.58,  23.45,  23.43,  23.43, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.93,  28.93,  27.94,  26.97,  26.23,  26.13,  26.09,  25.87,  25.29,  24.57,  23.76,  23.31,  23.20,  23.19,  23.19, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.67,  26.06,  26.02,  25.97,  25.66,  25.02,  24.26,  23.47,  23.09,  23.00,  22.99,  22.99, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.33,  27.34,  26.42,  25.95,  25.94,  25.85,  25.45,  24.77,  23.97,  23.22,  22.89,  22.81,  22.80,  22.80, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.07,  28.07,  27.09,  26.20,  25.88,  25.86,  25.74,  25.24,  24.54,  23.69,  22.98,  22.70,  22.63,  22.62,  22.62, 
      34.82,  33.82,  32.82,  31.82,  30.82,  29.82,  28.83,  27.83,  26.85,  26.03,  25.83,  25.80,  25.61,  25.02,  24.32,  23.44,  22.67,  22.51,  22.44,  22.43,  22.43, 
      41.03,  40.03,  39.03,  38.03,  37.03,  36.03,  35.03,  34.03,  33.03,  32.04,  31.06,  30.25,  30.07,  29.98,  29.53,  28.70,  28.02,  27.72,  27.63,  27.61,  27.61, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.33,  33.33,  32.33,  31.33,  30.34,  29.42,  28.99,  28.99,  28.95,  28.65,  27.94,  27.23,  26.78,  26.62,  26.59,  26.59, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.17,  31.17,  30.17,  29.20,  28.42,  28.27,  28.27,  28.23,  27.94,  27.27,  26.48,  25.94,  25.74,  25.71,  25.71, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.36,  31.36,  30.36,  29.37,  28.42,  27.83,  27.80,  27.79,  27.73,  27.38,  26.69,  25.83,  25.25,  25.05,  25.02,  25.02, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.76,  29.77,  28.78,  27.87,  27.49,  27.49,  27.47,  27.37,  26.93,  26.19,  25.29,  24.70,  24.53,  24.50,  24.50, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.31,  28.33,  27.49,  27.28,  27.27,  27.25,  27.10,  26.56,  25.76,  24.84,  24.28,  24.14,  24.12,  24.12, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.93,  28.94,  27.97,  27.23,  27.13,  27.13,  27.09,  26.86,  26.23,  25.38,  24.47,  23.96,  23.84,  23.83,  23.82, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.67,  27.06,  27.02,  27.02,  26.96,  26.64,  25.94,  25.04,  24.15,  23.70,  23.61,  23.60,  23.60, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.33,  28.34,  27.42,  26.95,  26.94,  26.93,  26.85,  26.41,  25.66,  24.72,  23.87,  23.51,  23.44,  23.43,  23.43, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.07,  28.09,  27.20,  26.88,  26.88,  26.86,  26.73,  26.18,  25.40,  24.42,  23.62,  23.35,  23.30,  23.29,  23.29, 
      34.82,  33.82,  32.82,  31.82,  30.82,  29.83,  28.83,  27.85,  27.03,  26.83,  26.83,  26.80,  26.60,  25.95,  25.14,  24.14,  23.38,  23.22,  23.17,  23.17,  23.17, 
      41.03,  40.03,  39.03,  38.03,  37.03,  36.03,  35.03,  34.03,  33.04,  32.06,  31.25,  31.08,  31.07,  30.98,  30.52,  29.70,  29.02,  28.72,  28.63,  28.61,  28.61, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.33,  33.33,  32.33,  31.34,  30.42,  30.00,  29.99,  29.99,  29.94,  29.65,  28.94,  28.23,  27.78,  27.62,  27.59,  27.59, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.17,  31.17,  30.20,  29.42,  29.27,  29.27,  29.27,  29.22,  28.94,  28.27,  27.48,  26.93,  26.74,  26.71,  26.70, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.36,  31.36,  30.37,  29.42,  28.83,  28.80,  28.80,  28.79,  28.73,  28.38,  27.69,  26.82,  26.24,  26.04,  26.01,  26.01, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.76,  30.77,  29.78,  28.87,  28.49,  28.49,  28.49,  28.47,  28.37,  27.93,  27.19,  26.27,  25.68,  25.50,  25.48,  25.47, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.31,  29.33,  28.49,  28.28,  28.27,  28.27,  28.25,  28.09,  27.55,  26.76,  25.80,  25.23,  25.09,  25.07,  25.07, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.93,  29.94,  28.97,  28.23,  28.13,  28.13,  28.12,  28.09,  27.85,  27.22,  26.37,  25.39,  24.88,  24.76,  24.75,  24.75, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.67,  28.06,  28.02,  28.02,  28.02,  27.96,  27.62,  26.93,  26.01,  25.04,  24.58,  24.49,  24.48,  24.48, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.33,  29.34,  28.42,  27.95,  27.94,  27.94,  27.93,  27.84,  27.39,  26.66,  25.67,  24.73,  24.33,  24.25,  24.24,  24.24, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.07,  29.09,  28.20,  27.88,  27.88,  27.88,  27.86,  27.72,  27.17,  26.40,  25.33,  24.45,  24.11,  24.03,  24.03,  24.03, 
      34.82,  33.82,  32.82,  31.82,  30.83,  29.83,  28.85,  28.03,  27.83,  27.83,  27.83,  27.80,  27.58,  26.94,  26.13,  25.00,  24.19,  23.89,  23.83,  23.83,  23.83, 
      41.03,  40.03,  39.03,  38.03,  37.03,  36.03,  35.03,  34.04,  33.06,  32.25,  32.08,  32.08,  32.07,  31.98,  31.52,  30.70,  30.02,  29.72,  29.63,  29.61,  29.61, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.33,  33.33,  32.34,  31.42,  31.00,  30.99,  30.99,  30.99,  30.94,  30.65,  29.94,  29.23,  28.78,  28.62,  28.59,  28.59, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.17,  31.20,  30.42,  30.27,  30.27,  30.27,  30.27,  30.22,  29.94,  29.27,  28.48,  27.93,  27.74,  27.71,  27.70, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.36,  31.37,  30.42,  29.83,  29.80,  29.80,  29.80,  29.79,  29.73,  29.38,  28.69,  27.82,  27.24,  27.04,  27.01,  27.01, 
      36.76,  35.76,  34.76,  33.76,  32.76,  31.77,  30.78,  29.87,  29.49,  29.49,  29.49,  29.49,  29.47,  29.37,  28.93,  28.19,  27.27,  26.68,  26.50,  26.48,  26.47, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.31,  30.33,  29.49,  29.28,  29.28,  29.27,  29.27,  29.25,  29.09,  28.55,  27.76,  26.80,  26.23,  26.09,  26.07,  26.07, 
      35.93,  34.93,  33.93,  32.93,  31.93,  30.94,  29.97,  29.23,  29.13,  29.13,  29.13,  29.12,  29.09,  28.85,  28.22,  27.37,  26.39,  25.88,  25.76,  25.75,  25.74, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.67,  29.06,  29.02,  29.02,  29.02,  29.02,  28.96,  28.62,  27.93,  27.01,  26.04,  25.58,  25.48,  25.47,  25.47, 
      35.33,  34.33,  33.33,  32.33,  31.33,  30.34,  29.42,  28.95,  28.94,  28.94,  28.94,  28.93,  28.84,  28.39,  27.66,  26.67,  25.71,  25.31,  25.23,  25.22,  25.22, 
      35.07,  34.07,  33.07,  32.07,  31.07,  30.09,  29.20,  28.88,  28.88,  28.88,  28.88,  28.86,  28.71,  28.17,  27.40,  26.33,  25.41,  25.07,  25.00,  24.99,  24.99, 
      34.82,  33.82,  32.82,  31.83,  30.83,  29.85,  29.03,  28.83,  28.83,  28.83,  28.83,  28.80,  28.57,  27.94,  27.13,  25.99,  25.10,  24.84,  24.78,  24.78,  24.78, 
      41.03,  40.03,  39.03,  38.03,  37.03,  36.03,  35.04,  34.06,  33.25,  33.08,  33.08,  33.08,  33.07,  32.98,  32.52,  31.70,  31.02,  30.72,  30.63,  30.61,  30.61, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.33,  33.34,  32.42,  32.00,  31.99,  31.99,  31.99,  31.99,  31.94,  31.65,  30.94,  30.23,  29.78,  29.62,  29.59,  29.59, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.17,  32.20,  31.42,  31.27,  31.27,  31.27,  31.27,  31.27,  31.22,  30.94,  30.27,  29.48,  28.93,  28.74,  28.71,  28.70, 
      37.36,  36.36,  35.36,  34.36,  33.36,  32.37,  31.42,  30.83,  30.80,  30.80,  30.80,  30.80,  30.79,  30.73,  30.38,  29.69,  28.82,  28.24,  28.04,  28.01,  28.01, 
      36.76,  35.76,  34.76,  33.76,  32.77,  31.78,  30.87,  30.49,  30.49,  30.49,  30.49,  30.49,  30.47,  30.37,  29.93,  29.19,  28.27,  27.68,  27.50,  27.48,  27.47, 
      36.31,  35.31,  34.31,  33.31,  32.31,  31.33,  30.49,  30.28,  30.28,  30.28,  30.27,  30.27,  30.25,  30.09,  29.55,  28.76,  27.80,  27.23,  27.09,  27.07,  27.07, 
      35.93,  34.93,  33.93,  32.93,  31.94,  30.97,  30.23,  30.13,  30.13,  30.13,  30.13,  30.12,  30.09,  29.85,  29.22,  28.37,  27.39,  26.88,  26.76,  26.75,  26.74, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.67,  30.06,  30.02,  30.02,  30.02,  30.02,  30.02,  29.96,  29.62,  28.93,  28.01,  27.04,  26.58,  26.48,  26.47,  26.47, 
      35.33,  34.33,  33.33,  32.33,  31.34,  30.42,  29.95,  29.94,  29.94,  29.94,  29.94,  29.93,  29.84,  29.39,  28.66,  27.67,  26.71,  26.31,  26.23,  26.22,  26.22, 
      35.07,  34.07,  33.07,  32.07,  31.09,  30.20,  29.88,  29.88,  29.88,  29.88,  29.88,  29.86,  29.71,  29.17,  28.40,  27.33,  26.41,  26.07,  26.00,  25.99,  25.99, 
      34.82,  33.82,  32.83,  31.83,  30.85,  30.03,  29.83,  29.83,  29.83,  29.83,  29.82,  29.79,  29.57,  28.94,  28.13,  26.99,  26.10,  25.84,  25.78,  25.78,  25.78, 
      41.03,  40.03,  39.03,  38.03,  37.03,  36.04,  35.06,  34.25,  34.08,  34.08,  34.08,  34.08,  34.07,  33.98,  33.52,  32.70,  32.02,  31.72,  31.63,  31.61,  31.61, 
      39.33,  38.33,  37.33,  36.33,  35.33,  34.34,  33.42,  33.00,  32.99,  32.99,  32.99,  32.99,  32.99,  32.94,  32.65,  31.94,  31.23,  30.78,  30.62,  30.59,  30.59, 
      38.17,  37.17,  36.17,  35.17,  34.17,  33.20,  32.42,  32.27,  32.27,  32.27,  32.27,  32.27,  32.27,  32.22,  31.94,  31.27,  30.48,  29.93,  29.74,  29.71,  29.70, 
      37.36,  36.36,  35.36,  34.36,  33.37,  32.42,  31.83,  31.80,  31.80,  31.80,  31.80,  31.80,  31.79,  31.73,  31.38,  30.69,  29.82,  29.24,  29.04,  29.01,  29.01, 
      36.76,  35.76,  34.76,  33.77,  32.78,  31.87,  31.49,  31.49,  31.49,  31.49,  31.49,  31.49,  31.47,  31.37,  30.93,  30.19,  29.27,  28.68,  28.50,  28.48,  28.47, 
      36.31,  35.31,  34.31,  33.31,  32.33,  31.49,  31.28,  31.28,  31.28,  31.28,  31.27,  31.27,  31.25,  31.09,  30.55,  29.76,  28.80,  28.23,  28.09,  28.07,  28.07, 
      35.93,  34.93,  33.93,  32.94,  31.97,  31.23,  31.13,  31.13,  31.13,  31.13,  31.13,  31.12,  31.09,  30.85,  30.22,  29.37,  28.39,  27.88,  27.76,  27.75,  27.74, 
      35.62,  34.62,  33.62,  32.62,  31.67,  31.06,  31.02,  31.02,  31.02,  31.02,  31.02,  31.02,  30.96,  30.62,  29.93,  29.01,  28.04,  27.58,  27.48,  27.47,  27.47, 
      35.33,  34.33,  33.33,  32.34,  31.42,  30.95,  30.94,  30.94,  30.94,  30.94,  30.94,  30.93,  30.84,  30.39,  29.66,  28.67,  27.71,  27.31,  27.23,  27.22,  27.22, 
      35.07,  34.07,  33.07,  32.09,  31.20,  30.88,  30.88,  30.88,  30.88,  30.88,  30.88,  30.86,  30.71,  30.17,  29.40,  28.33,  27.41,  27.07,  27.00,  26.99,  26.99, 
      34.82,  33.83,  32.83,  31.85,  31.03,  30.83,  30.83,  30.83,  30.83,  30.83,  30.82,  30.79,  30.57,  29.94,  29.13,  27.99,  27.10,  26.84,  26.78,  26.78,  26.78, 
      41.03,  40.03,  39.03,  38.03,  37.04,  36.06,  35.25,  35.08,  35.08,  35.08,  35.08,  35.08,  35.07,  34.98,  34.52,  33.70,  33.02,  32.72,  32.63,  32.61,  32.61, 
      39.33,  38.33,  37.33,  36.33,  35.34,  34.42,  34.00,  33.99,  33.99,  33.99,  33.99,  33.99,  33.99,  33.94,  33.65,  32.94,  32.23,  31.78,  31.62,  31.59,  31.59, 
      38.17,  37.17,  36.17,  35.17,  34.20,  33.42,  33.27,  33.27,  33.27,  33.27,  33.27,  33.27,  33.27,  33.22,  32.94,  32.27,  31.48,  30.93,  30.74,  30.71,  30.70, 
      37.36,  36.36,  35.36,  34.37,  33.42,  32.83,  32.80,  32.80,  32.80,  32.80,  32.80,  32.80,  32.79,  32.73,  32.38,  31.69,  30.82,  30.24,  30.04,  30.01,  30.01, 
      36.76,  35.76,  34.77,  33.78,  32.87,  32.49,  32.49,  32.49,  32.49,  32.49,  32.49,  32.49,  32.47,  32.37,  31.93,  31.19,  30.27,  29.68,  29.50,  29.48,  29.47, 
      36.31,  35.31,  34.31,  33.33,  32.49,  32.28,  32.28,  32.28,  32.28,  32.28,  32.27,  32.27,  32.25,  32.09,  31.55,  30.76,  29.80,  29.23,  29.09,  29.07,  29.07, 
      35.93,  34.93,  33.94,  32.97,  32.23,  32.13,  32.13,  32.13,  32.13,  32.13,  32.13,  32.12,  32.09,  31.85,  31.22,  30.37,  29.39,  28.88,  28.76,  28.75,  28.74, 
      35.62,  34.62,  33.62,  32.67,  32.06,  32.02,  32.02,  32.02,  32.02,  32.02,  32.02,  32.02,  31.96,  31.62,  30.93,  30.01,  29.04,  28.58,  28.48,  28.47,  28.47, 
      35.33,  34.33,  33.34,  32.42,  31.95,  31.94,  31.94,  31.94,  31.94,  31.94,  31.94,  31.93,  31.84,  31.39,  30.66,  29.67,  28.71,  28.31,  28.23,  28.22,  28.22, 
      35.07,  34.07,  33.09,  32.20,  31.88,  31.88,  31.88,  31.88,  31.88,  31.88,  31.88,  31.86,  31.71,  31.17,  30.40,  29.33,  28.41,  28.07,  28.00,  27.99,  27.99, 
      34.83,  33.83,  32.85,  32.03,  31.83,  31.83,  31.83,  31.83,  31.83,  31.83,  31.82,  31.79,  31.57,  30.94,  30.13,  28.99,  28.10,  27.84,  27.78,  27.78,  27.78, 
      41.03,  40.03,  39.03,  38.04,  37.06,  36.25,  36.08,  36.08,  36.08,  36.08,  36.08,  36.08,  36.07,  35.98,  35.52,  34.70,  34.02,  33.72,  33.63,  33.61,  33.61, 
      39.33,  38.33,  37.33,  36.34,  35.42,  35.00,  34.99,  34.99,  34.99,  34.99,  34.99,  34.99,  34.99,  34.94,  34.65,  33.94,  33.23,  32.78,  32.62,  32.59,  32.59, 
      38.17,  37.17,  36.17,  35.20,  34.42,  34.27,  34.27,  34.27,  34.27,  34.27,  34.27,  34.27,  34.27,  34.22,  33.94,  33.27,  32.48,  31.93,  31.74,  31.71,  31.70, 
      37.36,  36.36,  35.37,  34.42,  33.83,  33.80,  33.80,  33.80,  33.80,  33.80,  33.80,  33.80,  33.79,  33.73,  33.38,  32.69,  31.82,  31.24,  31.04,  31.01,  31.01, 
      36.76,  35.77,  34.78,  33.87,  33.49,  33.49,  33.49,  33.49,  33.49,  33.49,  33.49,  33.49,  33.47,  33.37,  32.93,  32.19,  31.27,  30.68,  30.50,  30.48,  30.47, 
      36.31,  35.31,  34.33,  33.49,  33.28,  33.28,  33.28,  33.28,  33.28,  33.28,  33.27,  33.27,  33.25,  33.09,  32.55,  31.76,  30.80,  30.23,  30.09,  30.07,  30.07, 
      35.93,  34.94,  33.97,  33.23,  33.13,  33.13,  33.13,  33.13,  33.13,  33.13,  33.13,  33.12,  33.09,  32.85,  32.22,  31.37,  30.39,  29.88,  29.76,  29.75,  29.74, 
      35.62,  34.62,  33.67,  33.06,  33.02,  33.02,  33.02,  33.02,  33.02,  33.02,  33.02,  33.02,  32.96,  32.62,  31.93,  31.01,  30.04,  29.58,  29.48,  29.47,  29.47, 
      35.33,  34.34,  33.42,  32.95,  32.94,  32.94,  32.94,  32.94,  32.94,  32.94,  32.94,  32.93,  32.84,  32.39,  31.66,  30.67,  29.71,  29.31,  29.23,  29.22,  29.22, 
      35.07,  34.09,  33.20,  32.88,  32.88,  32.88,  32.88,  32.88,  32.88,  32.88,  32.88,  32.86,  32.71,  32.17,  31.40,  30.33,  29.41,  29.07,  29.00,  28.99,  28.99, 
      34.83,  33.85,  33.03,  32.83,  32.83,  32.83,  32.83,  32.83,  32.83,  32.83,  32.82,  32.79,  32.57,  31.94,  31.13,  29.99,  29.10,  28.84,  28.78,  28.78,  28.78};

  setup_cool_interp_grid_(&my_rates->LHD, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_CI(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {13,  15.0, 1.0}, // log10(number-density like)
    {16,   0.6, 0.2}, // log10(temperature)
    {17, -10.0, 1.0}  // log10(H2 number density)
  };

  double L[] =
     {36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.62,  24.67,  23.98,  23.76,  23.73,  23.73, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.62,  24.62,  23.68,  23.01,  22.81,  22.78,  22.78, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.96,  25.96,  24.97,  23.97,  23.03,  22.40,  22.21,  22.19,  22.18, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.50,  24.50,  23.51,  22.58,  21.97,  21.78,  21.76,  21.76, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.15,  22.24,  21.65,  21.48,  21.46,  21.46, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.85,  23.85,  22.87,  21.98,  21.44,  21.30,  21.28,  21.28, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.63,  25.63,  24.63,  23.64,  22.65,  21.80,  21.32,  21.20,  21.18,  21.18, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.47,  25.47,  24.47,  23.47,  22.50,  21.67,  21.24,  21.14,  21.13,  21.13, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.35,  24.35,  23.35,  22.38,  21.58,  21.19,  21.11,  21.10,  21.10, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.26,  22.29,  21.52,  21.16,  21.09,  21.08,  21.08, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.17,  25.17,  24.17,  23.18,  22.22,  21.47,  21.14,  21.07,  21.07,  21.07, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.11,  24.11,  23.11,  22.16,  21.43,  21.12,  21.07,  21.06,  21.06, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.05,  23.05,  22.11,  21.40,  21.11,  21.06,  21.06,  21.06, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  25.99,  24.99,  23.99,  23.00,  22.06,  21.37,  21.10,  21.06,  21.05,  21.05, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.94,  24.94,  23.94,  22.95,  22.01,  21.34,  21.09,  21.06,  21.05,  21.05, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.89,  24.89,  23.89,  22.90,  21.97,  21.32,  21.09,  21.05,  21.05,  21.05, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.62,  24.67,  24.00,  23.80,  23.77,  23.77, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.62,  24.62,  23.68,  23.03,  22.84,  22.82,  22.82, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.96,  25.96,  24.97,  23.97,  23.04,  22.42,  22.23,  22.21,  22.21, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.50,  24.50,  23.51,  22.59,  21.98,  21.80,  21.77,  21.77, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.15,  22.24,  21.66,  21.49,  21.47,  21.47, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.85,  23.85,  22.87,  21.98,  21.45,  21.31,  21.29,  21.29, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.63,  25.63,  24.63,  23.64,  22.66,  21.80,  21.32,  21.20,  21.19,  21.19, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.47,  25.47,  24.47,  23.47,  22.50,  21.67,  21.24,  21.14,  21.13,  21.13, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.35,  24.35,  23.35,  22.38,  21.59,  21.19,  21.11,  21.10,  21.10, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.26,  22.29,  21.52,  21.16,  21.09,  21.08,  21.08, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.17,  25.17,  24.17,  23.18,  22.22,  21.47,  21.14,  21.08,  21.07,  21.07, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.11,  24.11,  23.11,  22.16,  21.43,  21.12,  21.07,  21.06,  21.06, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.05,  23.05,  22.11,  21.40,  21.11,  21.06,  21.06,  21.06, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  25.99,  24.99,  23.99,  23.00,  22.06,  21.37,  21.10,  21.06,  21.05,  21.05, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.94,  24.94,  23.94,  22.95,  22.01,  21.34,  21.10,  21.06,  21.05,  21.05, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.89,  24.89,  23.89,  22.90,  21.97,  21.32,  21.09,  21.05,  21.05,  21.05, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.63,  24.73,  24.22,  24.11,  24.09,  24.09, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.62,  24.63,  23.74,  23.24,  23.13,  23.11,  23.11, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.96,  25.96,  24.97,  23.98,  23.10,  22.58,  22.43,  22.41,  22.41, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.50,  24.50,  23.52,  22.65,  22.08,  21.91,  21.88,  21.88, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.14,  23.16,  22.29,  21.73,  21.56,  21.54,  21.54, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.85,  23.85,  22.88,  22.02,  21.50,  21.35,  21.34,  21.33, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.63,  25.63,  24.63,  23.64,  22.67,  21.84,  21.36,  21.24,  21.22,  21.22, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.47,  25.47,  24.47,  23.48,  22.52,  21.71,  21.27,  21.17,  21.15,  21.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.35,  24.35,  23.35,  22.40,  21.62,  21.22,  21.12,  21.11,  21.11, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.26,  22.32,  21.55,  21.18,  21.10,  21.09,  21.09, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.17,  25.17,  24.17,  23.18,  22.24,  21.50,  21.16,  21.08,  21.07,  21.07, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.11,  24.11,  23.12,  22.18,  21.46,  21.14,  21.07,  21.06,  21.06, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.05,  23.06,  22.13,  21.43,  21.13,  21.07,  21.06,  21.06, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  25.99,  24.99,  23.99,  23.00,  22.08,  21.40,  21.12,  21.06,  21.06,  21.05, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.94,  24.94,  23.94,  22.95,  22.04,  21.37,  21.11,  21.06,  21.05,  21.05, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.89,  24.89,  23.89,  22.90,  22.00,  21.35,  21.10,  21.06,  21.05,  21.05, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.62,  25.71,  25.17,  25.05,  25.03,  25.03,  25.03, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.62,  25.63,  24.73,  24.19,  24.02,  23.96,  23.95,  23.95, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.96,  25.97,  24.98,  24.08,  23.49,  23.18,  23.08,  23.07,  23.07, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.50,  24.52,  23.62,  22.93,  22.56,  22.47,  22.46,  22.46, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.14,  24.16,  23.25,  22.52,  22.13,  22.04,  22.03,  22.03, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.85,  23.87,  22.97,  22.23,  21.84,  21.73,  21.72,  21.72, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.63,  25.63,  24.64,  23.66,  22.77,  22.05,  21.64,  21.51,  21.50,  21.49, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.47,  25.47,  24.48,  23.51,  22.62,  21.91,  21.50,  21.36,  21.34,  21.34, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.35,  24.35,  23.39,  22.51,  21.82,  21.41,  21.26,  21.24,  21.24, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.26,  23.30,  22.43,  21.75,  21.34,  21.20,  21.17,  21.17, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.17,  25.17,  24.18,  23.22,  22.37,  21.70,  21.30,  21.15,  21.13,  21.13, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.11,  24.11,  23.16,  22.31,  21.66,  21.27,  21.12,  21.10,  21.10, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.05,  24.06,  23.11,  22.26,  21.62,  21.24,  21.10,  21.08,  21.08, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  25.99,  24.99,  24.00,  23.05,  22.22,  21.59,  21.22,  21.09,  21.07,  21.07, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.94,  24.94,  23.95,  23.01,  22.18,  21.56,  21.20,  21.08,  21.06,  21.06, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.89,  24.89,  23.90,  22.96,  22.14,  21.53,  21.19,  21.08,  21.06,  21.06, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.61,  27.62,  26.71,  26.17,  26.04,  26.02,  26.00,  26.00,  26.00, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.62,  26.63,  25.73,  25.19,  25.00,  24.88,  24.85,  24.84,  24.84, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.96,  26.97,  25.98,  25.08,  24.49,  24.16,  24.04,  24.03,  24.03,  24.03, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.50,  25.52,  24.62,  23.93,  23.55,  23.46,  23.45,  23.44,  23.44, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.14,  25.16,  24.25,  23.51,  23.12,  23.03,  23.02,  23.02,  23.02, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.85,  24.87,  23.97,  23.23,  22.82,  22.70,  22.68,  22.68,  22.68, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.63,  25.64,  24.66,  23.77,  23.04,  22.60,  22.43,  22.40,  22.40,  22.40, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.47,  25.48,  24.50,  23.62,  22.91,  22.44,  22.21,  22.16,  22.15,  22.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.35,  24.39,  23.51,  22.81,  22.32,  22.02,  21.94,  21.92,  21.92, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.26,  24.30,  23.43,  22.74,  22.23,  21.88,  21.74,  21.72,  21.71, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.17,  25.18,  24.22,  23.36,  22.68,  22.16,  21.76,  21.57,  21.53,  21.53, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.11,  24.16,  23.31,  22.64,  22.10,  21.68,  21.45,  21.39,  21.38, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.05,  25.06,  24.10,  23.26,  22.59,  22.06,  21.61,  21.35,  21.28,  21.27, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  25.99,  25.00,  24.05,  23.22,  22.56,  22.02,  21.56,  21.29,  21.21,  21.20, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.94,  24.95,  24.00,  23.17,  22.52,  21.98,  21.52,  21.24,  21.16,  21.14, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.89,  24.90,  23.96,  23.14,  22.49,  21.95,  21.49,  21.21,  21.12,  21.11, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.61,  28.62,  27.71,  27.17,  27.04,  27.01,  26.98,  26.97,  26.97,  26.97, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.62,  27.63,  26.73,  26.19,  26.00,  25.88,  25.85,  25.84,  25.84,  25.84, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.96,  27.97,  26.98,  26.08,  25.49,  25.16,  25.04,  25.03,  25.03,  25.03,  25.03, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.50,  26.52,  25.62,  24.93,  24.55,  24.46,  24.45,  24.44,  24.44,  24.44, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.14,  26.16,  25.25,  24.51,  24.12,  24.03,  24.02,  24.02,  24.02,  24.02, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.85,  25.87,  24.97,  24.23,  23.82,  23.70,  23.68,  23.68,  23.68,  23.68, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.63,  26.64,  25.66,  24.77,  24.04,  23.60,  23.43,  23.40,  23.40,  23.40,  23.40, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.47,  26.48,  25.50,  24.62,  23.91,  23.44,  23.21,  23.16,  23.15,  23.15,  23.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.35,  25.39,  24.51,  23.81,  23.32,  23.02,  22.93,  22.92,  22.92,  22.92, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.26,  25.30,  24.43,  23.74,  23.23,  22.87,  22.73,  22.70,  22.70,  22.70, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.17,  26.18,  25.22,  24.36,  23.68,  23.16,  22.76,  22.55,  22.50,  22.49,  22.49, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.11,  25.16,  24.31,  23.64,  23.10,  22.67,  22.38,  22.30,  22.28,  22.28, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.05,  26.06,  25.10,  24.26,  23.59,  23.06,  22.60,  22.25,  22.11,  22.08,  22.08, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  26.99,  26.00,  25.05,  24.22,  23.56,  23.02,  22.54,  22.15,  21.93,  21.88,  21.88, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.94,  25.95,  25.00,  24.17,  23.52,  22.98,  22.49,  22.07,  21.79,  21.70,  21.69, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.89,  25.90,  24.96,  24.14,  23.49,  22.95,  22.46,  22.00,  21.67,  21.54,  21.52, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.61,  29.62,  28.71,  28.17,  28.04,  28.01,  27.98,  27.97,  27.97,  27.97,  27.97, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.62,  28.63,  27.73,  27.19,  27.00,  26.88,  26.85,  26.84,  26.84,  26.84,  26.84, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.96,  28.97,  27.98,  27.08,  26.49,  26.16,  26.04,  26.03,  26.03,  26.03,  26.03,  26.03, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.50,  27.52,  26.62,  25.93,  25.55,  25.46,  25.44,  25.44,  25.44,  25.44,  25.44, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.14,  27.16,  26.25,  25.51,  25.12,  25.03,  25.02,  25.02,  25.02,  25.02,  25.02, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.85,  26.87,  25.97,  25.23,  24.82,  24.70,  24.68,  24.68,  24.68,  24.68,  24.68, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.63,  27.64,  26.66,  25.77,  25.04,  24.60,  24.43,  24.40,  24.40,  24.40,  24.40,  24.40, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.47,  27.48,  26.50,  25.62,  24.91,  24.44,  24.21,  24.16,  24.15,  24.15,  24.15,  24.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.35,  26.39,  25.51,  24.81,  24.32,  24.02,  23.93,  23.92,  23.92,  23.92,  23.92, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.26,  26.30,  25.43,  24.74,  24.23,  23.87,  23.73,  23.70,  23.70,  23.70,  23.70, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.17,  27.18,  26.22,  25.36,  24.68,  24.16,  23.76,  23.55,  23.50,  23.49,  23.49,  23.49, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.11,  26.16,  25.31,  24.64,  24.10,  23.67,  23.38,  23.30,  23.28,  23.28,  23.28, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.05,  27.06,  26.10,  25.26,  24.59,  24.06,  23.60,  23.25,  23.11,  23.08,  23.08,  23.08, 
      32.99,  31.99,  30.99,  29.99,  28.99,  27.99,  27.00,  26.05,  25.22,  24.56,  24.02,  23.54,  23.15,  22.93,  22.88,  22.88,  22.88, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.94,  26.95,  26.00,  25.17,  24.52,  23.98,  23.49,  23.07,  22.78,  22.69,  22.68,  22.67, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.89,  26.90,  25.96,  25.14,  24.49,  23.95,  23.46,  23.00,  22.66,  22.50,  22.48,  22.47, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.61,  30.62,  29.71,  29.17,  29.04,  29.01,  28.98,  28.97,  28.97,  28.97,  28.97,  28.97, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.62,  29.63,  28.73,  28.19,  28.00,  27.88,  27.85,  27.84,  27.84,  27.84,  27.84,  27.84, 
      34.96,  33.96,  32.96,  31.96,  30.96,  29.97,  28.98,  28.08,  27.49,  27.16,  27.04,  27.03,  27.03,  27.03,  27.03,  27.03,  27.03, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.50,  28.52,  27.62,  26.93,  26.55,  26.46,  26.44,  26.44,  26.44,  26.44,  26.44,  26.44, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.14,  28.16,  27.25,  26.51,  26.12,  26.03,  26.02,  26.02,  26.01,  26.01,  26.01,  26.01, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.85,  27.87,  26.97,  26.23,  25.82,  25.70,  25.68,  25.68,  25.68,  25.68,  25.68,  25.68, 
      33.63,  32.63,  31.63,  30.63,  29.63,  28.64,  27.66,  26.77,  26.04,  25.60,  25.43,  25.40,  25.40,  25.40,  25.40,  25.40,  25.40, 
      33.47,  32.47,  31.47,  30.47,  29.47,  28.48,  27.50,  26.62,  25.91,  25.44,  25.21,  25.16,  25.15,  25.15,  25.15,  25.15,  25.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.35,  27.39,  26.51,  25.81,  25.32,  25.02,  24.93,  24.92,  24.92,  24.92,  24.92,  24.92, 
      33.25,  32.25,  31.25,  30.25,  29.25,  28.26,  27.30,  26.43,  25.74,  25.23,  24.87,  24.73,  24.70,  24.70,  24.70,  24.70,  24.70, 
      33.17,  32.17,  31.17,  30.17,  29.17,  28.18,  27.22,  26.36,  25.68,  25.16,  24.76,  24.54,  24.50,  24.49,  24.49,  24.49,  24.49, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.11,  27.16,  26.31,  25.64,  25.10,  24.67,  24.38,  24.30,  24.28,  24.28,  24.28,  24.28, 
      33.05,  32.05,  31.05,  30.05,  29.05,  28.06,  27.10,  26.26,  25.59,  25.06,  24.60,  24.25,  24.11,  24.08,  24.08,  24.08,  24.08, 
      32.99,  31.99,  30.99,  29.99,  28.99,  28.00,  27.05,  26.22,  25.56,  25.02,  24.54,  24.15,  23.93,  23.88,  23.88,  23.88,  23.88, 
      32.94,  31.94,  30.94,  29.94,  28.94,  27.95,  27.00,  26.17,  25.52,  24.98,  24.49,  24.07,  23.78,  23.69,  23.67,  23.67,  23.67, 
      32.89,  31.89,  30.89,  29.89,  28.89,  27.90,  26.96,  26.13,  25.49,  24.95,  24.45,  24.00,  23.66,  23.50,  23.48,  23.47,  23.47, 
      36.61,  35.61,  34.61,  33.61,  32.61,  31.62,  30.71,  30.17,  30.04,  30.01,  29.98,  29.97,  29.97,  29.97,  29.97,  29.97,  29.97, 
      35.62,  34.62,  33.62,  32.62,  31.62,  30.63,  29.73,  29.19,  29.00,  28.88,  28.85,  28.84,  28.84,  28.84,  28.84,  28.84,  28.84, 
      34.96,  33.96,  32.96,  31.96,  30.97,  29.98,  29.08,  28.49,  28.15,  28.04,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02, 
      34.50,  33.50,  32.50,  31.50,  30.50,  29.52,  28.62,  27.92,  27.54,  27.45,  27.43,  27.43,  27.43,  27.43,  27.43,  27.43,  27.43, 
      34.14,  33.14,  32.14,  31.14,  30.14,  29.16,  28.25,  27.51,  27.11,  27.02,  27.01,  27.00,  27.00,  27.00,  27.00,  27.00,  27.00, 
      33.85,  32.85,  31.85,  30.85,  29.85,  28.87,  27.97,  27.23,  26.81,  26.69,  26.67,  26.67,  26.67,  26.67,  26.67,  26.67,  26.67, 
      33.63,  32.63,  31.63,  30.63,  29.64,  28.66,  27.76,  27.04,  26.59,  26.43,  26.40,  26.39,  26.39,  26.39,  26.39,  26.39,  26.39, 
      33.47,  32.47,  31.47,  30.47,  29.48,  28.50,  27.62,  26.90,  26.43,  26.20,  26.15,  26.15,  26.15,  26.15,  26.15,  26.15,  26.15, 
      33.35,  32.35,  31.35,  30.35,  29.35,  28.39,  27.51,  26.81,  26.31,  26.02,  25.93,  25.92,  25.92,  25.92,  25.92,  25.92,  25.92, 
      33.25,  32.25,  31.25,  30.25,  29.26,  28.30,  27.43,  26.74,  26.22,  25.87,  25.73,  25.70,  25.70,  25.70,  25.70,  25.70,  25.70, 
      33.17,  32.17,  31.17,  30.17,  29.18,  28.22,  27.36,  26.68,  26.15,  25.75,  25.54,  25.50,  25.49,  25.49,  25.49,  25.49,  25.49, 
      33.11,  32.11,  31.11,  30.11,  29.11,  28.16,  27.31,  26.63,  26.10,  25.66,  25.38,  25.30,  25.28,  25.28,  25.28,  25.28,  25.28, 
      33.05,  32.05,  31.05,  30.05,  29.06,  28.10,  27.26,  26.59,  26.06,  25.59,  25.25,  25.11,  25.08,  25.08,  25.08,  25.08,  25.08, 
      32.99,  31.99,  30.99,  29.99,  29.00,  28.05,  27.21,  26.55,  26.02,  25.54,  25.15,  24.93,  24.88,  24.88,  24.88,  24.87,  24.87, 
      32.94,  31.94,  30.94,  29.94,  28.95,  28.00,  27.17,  26.52,  25.98,  25.49,  25.07,  24.78,  24.69,  24.67,  24.67,  24.67,  24.67, 
      32.89,  31.89,  30.89,  29.89,  28.90,  27.96,  27.13,  26.49,  25.95,  25.45,  25.00,  24.66,  24.50,  24.48,  24.47,  24.47,  24.47, 
      36.61,  35.61,  34.61,  33.61,  32.62,  31.71,  31.17,  31.04,  31.01,  30.98,  30.97,  30.97,  30.97,  30.97,  30.97,  30.97,  30.97, 
      35.62,  34.62,  33.62,  32.62,  31.63,  30.73,  30.19,  30.00,  29.87,  29.84,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83, 
      34.96,  33.96,  32.96,  31.97,  30.98,  30.08,  29.48,  29.13,  29.00,  28.98,  28.98,  28.98,  28.98,  28.98,  28.98,  28.98,  28.98, 
      34.50,  33.50,  32.50,  31.50,  30.52,  29.62,  28.91,  28.48,  28.37,  28.36,  28.36,  28.35,  28.35,  28.35,  28.35,  28.35,  28.35, 
      34.14,  33.14,  32.14,  31.14,  30.16,  29.25,  28.48,  28.04,  27.93,  27.92,  27.91,  27.91,  27.91,  27.91,  27.91,  27.91,  27.91, 
      33.85,  32.85,  31.85,  30.85,  29.87,  28.96,  28.19,  27.75,  27.62,  27.60,  27.59,  27.59,  27.59,  27.59,  27.59,  27.59,  27.59, 
      33.63,  32.63,  31.63,  30.64,  29.66,  28.76,  28.00,  27.54,  27.37,  27.34,  27.34,  27.34,  27.34,  27.34,  27.34,  27.34,  27.34, 
      33.47,  32.47,  31.47,  30.48,  29.50,  28.61,  27.87,  27.39,  27.17,  27.11,  27.11,  27.11,  27.11,  27.11,  27.11,  27.11,  27.11, 
      33.35,  32.35,  31.35,  30.35,  29.39,  28.50,  27.77,  27.28,  26.99,  26.91,  26.89,  26.89,  26.89,  26.89,  26.89,  26.89,  26.89, 
      33.25,  32.25,  31.25,  30.26,  29.30,  28.41,  27.70,  27.20,  26.85,  26.71,  26.69,  26.68,  26.68,  26.68,  26.68,  26.68,  26.68, 
      33.17,  32.17,  31.17,  30.18,  29.22,  28.35,  27.64,  27.13,  26.74,  26.53,  26.48,  26.48,  26.48,  26.48,  26.48,  26.48,  26.48, 
      33.11,  32.11,  31.11,  30.11,  29.16,  28.29,  27.60,  27.08,  26.65,  26.38,  26.29,  26.28,  26.28,  26.28,  26.28,  26.28,  26.28, 
      33.05,  32.05,  31.05,  30.06,  29.10,  28.24,  27.56,  27.04,  26.59,  26.25,  26.10,  26.08,  26.07,  26.07,  26.07,  26.07,  26.07, 
      32.99,  31.99,  30.99,  30.00,  29.05,  28.19,  27.52,  27.00,  26.53,  26.14,  25.93,  25.88,  25.87,  25.87,  25.87,  25.87,  25.87, 
      32.94,  31.94,  30.94,  29.95,  29.00,  28.15,  27.49,  26.97,  26.49,  26.06,  25.78,  25.69,  25.67,  25.67,  25.67,  25.67,  25.67, 
      32.89,  31.89,  30.89,  29.90,  28.95,  28.11,  27.46,  26.94,  26.45,  26.00,  25.66,  25.50,  25.47,  25.47,  25.47,  25.47,  25.47, 
      36.61,  35.61,  34.61,  33.62,  32.71,  32.17,  32.04,  32.01,  31.98,  31.97,  31.97,  31.97,  31.97,  31.97,  31.97,  31.97,  31.97, 
      35.62,  34.62,  33.62,  32.63,  31.73,  31.19,  31.00,  30.86,  30.82,  30.82,  30.82,  30.82,  30.82,  30.82,  30.82,  30.82,  30.82, 
      34.96,  33.96,  32.97,  31.98,  31.08,  30.48,  30.09,  29.92,  29.90,  29.90,  29.90,  29.90,  29.90,  29.90,  29.90,  29.90,  29.90, 
      34.50,  33.50,  32.50,  31.52,  30.62,  29.89,  29.36,  29.19,  29.17,  29.16,  29.16,  29.16,  29.16,  29.16,  29.16,  29.16,  29.16, 
      34.14,  33.14,  32.14,  31.16,  30.24,  29.44,  28.86,  28.67,  28.64,  28.64,  28.64,  28.64,  28.64,  28.64,  28.64,  28.64,  28.64, 
      33.85,  32.85,  31.85,  30.87,  29.96,  29.13,  28.54,  28.33,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29, 
      33.63,  32.63,  31.64,  30.66,  29.75,  28.92,  28.34,  28.10,  28.06,  28.05,  28.05,  28.05,  28.05,  28.05,  28.05,  28.05,  28.05, 
      33.47,  32.47,  31.48,  30.50,  29.59,  28.77,  28.20,  27.94,  27.88,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87, 
      33.35,  32.35,  31.35,  30.38,  29.48,  28.67,  28.11,  27.82,  27.72,  27.71,  27.71,  27.71,  27.71,  27.71,  27.71,  27.71,  27.71, 
      33.25,  32.25,  31.26,  30.29,  29.39,  28.59,  28.04,  27.71,  27.58,  27.55,  27.55,  27.55,  27.55,  27.55,  27.55,  27.55,  27.55, 
      33.17,  32.17,  31.18,  30.22,  29.32,  28.53,  27.98,  27.63,  27.44,  27.39,  27.39,  27.39,  27.39,  27.39,  27.39,  27.39,  27.39, 
      33.11,  32.11,  31.11,  30.15,  29.25,  28.48,  27.94,  27.57,  27.31,  27.23,  27.21,  27.21,  27.21,  27.21,  27.21,  27.21,  27.21, 
      33.05,  32.05,  31.06,  30.10,  29.20,  28.43,  27.91,  27.51,  27.20,  27.06,  27.04,  27.03,  27.03,  27.03,  27.03,  27.03,  27.03, 
      32.99,  31.99,  31.00,  30.05,  29.15,  28.39,  27.87,  27.47,  27.11,  26.90,  26.85,  26.85,  26.85,  26.85,  26.85,  26.85,  26.85, 
      32.94,  31.94,  30.95,  30.00,  29.10,  28.36,  27.85,  27.43,  27.04,  26.76,  26.67,  26.66,  26.65,  26.65,  26.65,  26.65,  26.65, 
      32.89,  31.89,  30.90,  29.95,  29.06,  28.32,  27.82,  27.40,  26.98,  26.64,  26.49,  26.46,  26.46,  26.46,  26.46,  26.46,  26.46, 
      36.61,  35.61,  34.62,  33.71,  33.17,  33.04,  33.01,  32.98,  32.97,  32.97,  32.97,  32.97,  32.97,  32.97,  32.97,  32.97,  32.97, 
      35.62,  34.62,  33.63,  32.73,  32.19,  32.00,  31.86,  31.82,  31.82,  31.82,  31.82,  31.82,  31.82,  31.82,  31.82,  31.82,  31.82, 
      34.96,  33.97,  32.98,  32.08,  31.48,  31.09,  30.92,  30.90,  30.89,  30.89,  30.89,  30.89,  30.89,  30.89,  30.89,  30.89,  30.89, 
      34.50,  33.50,  32.52,  31.62,  30.89,  30.35,  30.18,  30.15,  30.15,  30.15,  30.15,  30.15,  30.15,  30.15,  30.15,  30.15,  30.15, 
      34.14,  33.14,  32.16,  31.24,  30.44,  29.84,  29.63,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59, 
      33.85,  32.85,  31.87,  30.96,  30.12,  29.50,  29.23,  29.17,  29.17,  29.17,  29.17,  29.17,  29.17,  29.17,  29.17,  29.17,  29.17, 
      33.63,  32.64,  31.66,  30.74,  29.91,  29.28,  28.94,  28.85,  28.84,  28.83,  28.83,  28.83,  28.83,  28.83,  28.83,  28.83,  28.83, 
      33.47,  32.48,  31.50,  30.59,  29.76,  29.13,  28.73,  28.59,  28.56,  28.56,  28.56,  28.56,  28.56,  28.56,  28.56,  28.56,  28.56, 
      33.35,  32.35,  31.38,  30.48,  29.65,  29.02,  28.58,  28.39,  28.35,  28.34,  28.34,  28.34,  28.34,  28.34,  28.34,  28.34,  28.34, 
      33.25,  32.26,  31.29,  30.39,  29.57,  28.94,  28.48,  28.23,  28.18,  28.17,  28.17,  28.17,  28.17,  28.17,  28.17,  28.17,  28.17, 
      33.17,  32.18,  31.22,  30.31,  29.51,  28.88,  28.40,  28.12,  28.04,  28.03,  28.03,  28.03,  28.03,  28.03,  28.03,  28.03,  28.03, 
      33.11,  32.11,  31.15,  30.25,  29.45,  28.83,  28.33,  28.03,  27.92,  27.91,  27.90,  27.90,  27.90,  27.90,  27.90,  27.90,  27.90, 
      33.05,  32.06,  31.10,  30.20,  29.41,  28.79,  28.29,  27.96,  27.82,  27.79,  27.79,  27.79,  27.79,  27.79,  27.79,  27.79,  27.79, 
      32.99,  32.00,  31.05,  30.15,  29.37,  28.75,  28.25,  27.90,  27.72,  27.67,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66, 
      32.94,  31.95,  31.00,  30.10,  29.33,  28.71,  28.21,  27.85,  27.62,  27.54,  27.52,  27.52,  27.52,  27.52,  27.52,  27.52,  27.52, 
      32.89,  31.90,  30.95,  30.06,  29.29,  28.68,  28.18,  27.82,  27.54,  27.40,  27.37,  27.37,  27.37,  27.37,  27.37,  27.37,  27.37, 
      36.61,  35.62,  34.71,  34.17,  34.04,  34.01,  33.98,  33.97,  33.97,  33.97,  33.97,  33.97,  33.97,  33.97,  33.97,  33.97,  33.97, 
      35.62,  34.63,  33.73,  33.19,  33.00,  32.86,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82,  32.82, 
      34.97,  33.98,  33.08,  32.48,  32.09,  31.92,  31.90,  31.89,  31.89,  31.89,  31.89,  31.89,  31.89,  31.89,  31.89,  31.89,  31.89, 
      34.50,  33.52,  32.62,  31.89,  31.35,  31.18,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15,  31.15, 
      34.14,  33.16,  32.24,  31.44,  30.84,  30.63,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59, 
      33.85,  32.87,  31.96,  31.12,  30.50,  30.23,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17,  30.17, 
      33.64,  32.66,  31.74,  30.91,  30.28,  29.94,  29.85,  29.84,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83,  29.83, 
      33.48,  32.50,  31.59,  30.76,  30.13,  29.73,  29.58,  29.56,  29.55,  29.55,  29.55,  29.55,  29.55,  29.55,  29.55,  29.55,  29.55, 
      33.35,  32.38,  31.48,  30.65,  30.02,  29.58,  29.36,  29.31,  29.31,  29.31,  29.31,  29.31,  29.31,  29.31,  29.31,  29.31,  29.31, 
      33.26,  32.29,  31.39,  30.57,  29.94,  29.47,  29.18,  29.09,  29.08,  29.08,  29.08,  29.08,  29.08,  29.08,  29.08,  29.08,  29.08, 
      33.18,  32.22,  31.31,  30.51,  29.88,  29.38,  29.04,  28.89,  28.86,  28.86,  28.86,  28.86,  28.86,  28.86,  28.86,  28.86,  28.86, 
      33.11,  32.15,  31.25,  30.45,  29.83,  29.32,  28.92,  28.71,  28.65,  28.65,  28.65,  28.65,  28.65,  28.65,  28.65,  28.65,  28.65, 
      33.06,  32.10,  31.20,  30.41,  29.79,  29.27,  28.83,  28.55,  28.46,  28.45,  28.45,  28.45,  28.45,  28.45,  28.45,  28.45,  28.45, 
      33.00,  32.05,  31.15,  30.37,  29.75,  29.22,  28.77,  28.43,  28.30,  28.28,  28.28,  28.28,  28.28,  28.28,  28.28,  28.28,  28.28, 
      32.95,  32.00,  31.10,  30.33,  29.71,  29.19,  28.71,  28.34,  28.17,  28.13,  28.13,  28.13,  28.13,  28.13,  28.13,  28.13,  28.13, 
      32.90,  31.95,  31.06,  30.29,  29.68,  29.15,  28.67,  28.26,  28.06,  28.01,  28.01,  28.01,  28.01,  28.01,  28.01,  28.01,  28.01}; 

  setup_cool_interp_grid_(&my_rates->LCI, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_CII(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {13,  15.0, 1.0}, // log10(number-density like)
    {16,   0.6, 0.2}, // log10(temperature)
    {17, -10.0, 1.0}  // log10(H2 number density)
  };

  double L[] =
     {42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.75,  30.76,  29.85,  29.31,  29.19,  29.17, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.06,  27.07,  26.17,  25.63,  25.52,  25.50, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.73,  24.74,  23.84,  23.31,  23.20,  23.18, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.26,  22.36,  21.85,  21.74,  21.72, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.32,  22.33,  21.44,  20.94,  20.83,  20.82, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.74,  20.86,  20.40,  20.31,  20.30, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.34,  22.34,  21.36,  20.51,  20.10,  20.02,  20.02, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.09,  23.09,  22.09,  21.12,  20.30,  19.94,  19.88,  19.87, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.93,  22.93,  21.94,  20.96,  20.17,  19.85,  19.80,  19.80, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.83,  21.83,  20.86,  20.09,  19.80,  19.76,  19.76, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.75,  23.75,  22.75,  21.76,  20.79,  20.04,  19.78,  19.74,  19.73, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.70,  22.70,  21.71,  20.75,  20.01,  19.76,  19.72,  19.72, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.67,  22.67,  21.67,  20.71,  19.99,  19.75,  19.71,  19.71, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.64,  22.64,  21.64,  20.68,  19.97,  19.74,  19.71,  19.70, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.61,  22.61,  21.62,  20.66,  19.96,  19.74,  19.70,  19.70, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.59,  22.59,  21.60,  20.65,  19.95,  19.73,  19.70,  19.70, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.75,  30.76,  29.85,  29.32,  29.20,  29.19, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.06,  27.07,  26.17,  25.65,  25.53,  25.52, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.73,  24.74,  23.84,  23.32,  23.21,  23.20, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.26,  22.37,  21.86,  21.75,  21.74, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.32,  22.33,  21.44,  20.95,  20.85,  20.83, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.74,  20.86,  20.41,  20.32,  20.31, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.34,  22.34,  21.36,  20.51,  20.11,  20.03,  20.02, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.09,  23.09,  22.09,  21.12,  20.30,  19.94,  19.88,  19.88, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.93,  22.93,  21.94,  20.96,  20.17,  19.86,  19.81,  19.80, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.83,  21.83,  20.86,  20.10,  19.81,  19.76,  19.76, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.75,  23.75,  22.75,  21.76,  20.79,  20.05,  19.78,  19.74,  19.73, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.70,  22.70,  21.71,  20.75,  20.02,  19.76,  19.72,  19.72, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.67,  22.67,  21.67,  20.71,  19.99,  19.75,  19.71,  19.71, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.64,  22.64,  21.64,  20.69,  19.98,  19.74,  19.71,  19.71, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.61,  22.61,  21.62,  20.67,  19.96,  19.74,  19.71,  19.70, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.59,  22.59,  21.60,  20.65,  19.95,  19.73,  19.70,  19.70, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.75,  30.76,  29.89,  29.43,  29.34,  29.33, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.06,  27.08,  26.20,  25.75,  25.66,  25.65, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.73,  24.75,  23.88,  23.43,  23.34,  23.33, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.25,  23.27,  22.40,  21.97,  21.88,  21.87, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.32,  22.33,  21.48,  21.05,  20.97,  20.96, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.74,  20.90,  20.50,  20.42,  20.41, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.34,  22.34,  21.37,  20.55,  20.18,  20.10,  20.10, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.09,  23.09,  22.10,  21.13,  20.34,  19.99,  19.93,  19.92, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.93,  22.93,  21.94,  20.97,  20.21,  19.89,  19.84,  19.83, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.83,  21.83,  20.87,  20.14,  19.83,  19.78,  19.78, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.75,  23.75,  22.75,  21.76,  20.81,  20.08,  19.80,  19.75,  19.75, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.70,  22.70,  21.71,  20.76,  20.05,  19.78,  19.73,  19.73, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.67,  22.67,  21.67,  20.73,  20.03,  19.76,  19.72,  19.72, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.64,  22.64,  21.64,  20.70,  20.01,  19.75,  19.71,  19.71, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.61,  22.61,  21.62,  20.68,  20.00,  19.75,  19.71,  19.70, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.59,  22.59,  21.60,  20.66,  19.99,  19.74,  19.71,  19.70, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.75,  30.83,  30.22,  30.07,  30.05,  30.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.06,  28.07,  27.14,  26.55,  26.39,  26.38,  26.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.73,  25.74,  24.81,  24.22,  24.08,  24.06,  24.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.25,  24.26,  23.34,  22.76,  22.61,  22.59,  22.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.32,  23.33,  22.41,  21.83,  21.68,  21.66,  21.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.73,  21.82,  21.24,  21.07,  21.05,  21.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.34,  22.35,  21.46,  20.87,  20.65,  20.62,  20.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.09,  23.09,  22.11,  21.23,  20.63,  20.36,  20.31,  20.31, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.93,  22.93,  21.95,  21.08,  20.47,  20.17,  20.10,  20.09, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.83,  21.85,  20.99,  20.37,  20.04,  19.96,  19.95, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.75,  23.75,  22.76,  21.78,  20.93,  20.31,  19.96,  19.87,  19.86, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.70,  22.71,  21.73,  20.89,  20.27,  19.91,  19.81,  19.80, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.67,  22.67,  21.70,  20.86,  20.24,  19.88,  19.77,  19.76, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.64,  22.64,  21.67,  20.83,  20.21,  19.85,  19.75,  19.74, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.61,  22.62,  21.65,  20.81,  20.20,  19.84,  19.74,  19.72, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.59,  22.60,  21.63,  20.80,  20.19,  19.83,  19.73,  19.71, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.75,  31.83,  31.22,  31.07,  31.05,  31.05,  31.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.06,  29.07,  28.14,  27.55,  27.39,  27.38,  27.37,  27.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.73,  26.74,  25.81,  25.22,  25.08,  25.06,  25.06,  25.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.25,  25.26,  24.34,  23.76,  23.61,  23.59,  23.59,  23.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.32,  24.33,  23.41,  22.83,  22.68,  22.66,  22.66,  22.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.73,  22.82,  22.24,  22.07,  22.05,  22.04,  22.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.34,  23.35,  22.46,  21.86,  21.65,  21.61,  21.61,  21.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.09,  23.11,  22.23,  21.62,  21.34,  21.28,  21.27,  21.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.93,  22.95,  22.08,  21.46,  21.11,  21.00,  20.99,  20.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.83,  22.85,  21.99,  21.36,  20.95,  20.77,  20.74,  20.74, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.75,  23.76,  22.78,  21.93,  21.29,  20.84,  20.58,  20.52,  20.51, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.70,  23.71,  22.73,  21.89,  21.25,  20.76,  20.43,  20.32,  20.30, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.67,  22.70,  21.85,  21.21,  20.70,  20.32,  20.15,  20.12, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.64,  22.67,  21.83,  21.19,  20.67,  20.25,  20.03,  19.99, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.61,  23.62,  22.65,  21.81,  21.17,  20.64,  20.19,  19.95,  19.89, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.59,  23.60,  22.63,  21.80,  21.16,  20.62,  20.16,  19.89,  19.83, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.75,  32.83,  32.22,  32.07,  32.05,  32.05,  32.05,  32.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.06,  30.07,  29.14,  28.55,  28.39,  28.38,  28.37,  28.37,  28.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.73,  27.74,  26.81,  26.22,  26.08,  26.06,  26.06,  26.06,  26.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.25,  26.26,  25.34,  24.76,  24.61,  24.59,  24.59,  24.59,  24.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.32,  25.33,  24.41,  23.83,  23.68,  23.66,  23.66,  23.66,  23.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.72,  24.73,  23.82,  23.24,  23.07,  23.05,  23.04,  23.04,  23.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.34,  24.35,  23.46,  22.86,  22.65,  22.61,  22.61,  22.61,  22.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.09,  24.11,  23.23,  22.62,  22.34,  22.28,  22.27,  22.27,  22.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.93,  23.95,  23.08,  22.46,  22.11,  22.00,  21.99,  21.99,  21.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.83,  23.85,  22.99,  22.36,  21.95,  21.77,  21.74,  21.74,  21.74, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.75,  24.76,  23.78,  22.93,  22.29,  21.84,  21.58,  21.51,  21.51,  21.50, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.70,  24.71,  23.73,  22.89,  22.25,  21.76,  21.43,  21.31,  21.29,  21.29, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.67,  23.70,  22.85,  22.21,  21.70,  21.32,  21.12,  21.08,  21.07, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.64,  23.67,  22.83,  22.19,  21.67,  21.24,  20.96,  20.88,  20.87, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.61,  24.62,  23.65,  22.81,  22.17,  21.64,  21.18,  20.84,  20.69,  20.67, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.59,  24.60,  23.63,  22.80,  22.16,  21.62,  21.14,  20.74,  20.52,  20.47, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.75,  33.83,  33.22,  33.07,  33.05,  33.05,  33.05,  33.05,  33.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.06,  31.07,  30.14,  29.55,  29.39,  29.38,  29.37,  29.37,  29.37,  29.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.73,  28.74,  27.81,  27.22,  27.08,  27.06,  27.06,  27.06,  27.06,  27.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.25,  27.26,  26.34,  25.76,  25.61,  25.59,  25.59,  25.59,  25.59,  25.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.32,  26.33,  25.41,  24.83,  24.68,  24.66,  24.66,  24.66,  24.66,  24.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.72,  25.73,  24.82,  24.24,  24.07,  24.05,  24.04,  24.04,  24.04,  24.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.34,  25.35,  24.46,  23.86,  23.65,  23.61,  23.61,  23.61,  23.61,  23.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.09,  25.11,  24.23,  23.62,  23.34,  23.28,  23.27,  23.27,  23.27,  23.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.93,  24.95,  24.08,  23.46,  23.11,  23.00,  22.99,  22.99,  22.99,  22.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.83,  24.85,  23.99,  23.36,  22.95,  22.77,  22.74,  22.74,  22.74,  22.74, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.75,  25.76,  24.78,  23.93,  23.29,  22.84,  22.58,  22.51,  22.51,  22.50,  22.50, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.70,  25.71,  24.73,  23.89,  23.25,  22.76,  22.43,  22.31,  22.29,  22.29,  22.29, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.67,  24.70,  23.85,  23.21,  22.70,  22.32,  22.12,  22.08,  22.07,  22.07, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.64,  24.67,  23.83,  23.19,  22.67,  22.24,  21.96,  21.88,  21.87,  21.87, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.61,  25.62,  24.65,  23.81,  23.17,  22.64,  22.18,  21.84,  21.69,  21.67,  21.66, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.59,  25.60,  24.63,  23.80,  23.16,  22.62,  22.14,  21.74,  21.52,  21.47,  21.46, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.75,  34.83,  34.22,  34.07,  34.05,  34.05,  34.05,  34.05,  34.05,  34.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.06,  32.07,  31.14,  30.55,  30.39,  30.38,  30.37,  30.37,  30.37,  30.37,  30.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.73,  29.74,  28.81,  28.22,  28.08,  28.06,  28.06,  28.06,  28.06,  28.06,  28.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.25,  28.26,  27.34,  26.76,  26.61,  26.59,  26.59,  26.59,  26.59,  26.59,  26.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.32,  27.33,  26.41,  25.83,  25.68,  25.66,  25.66,  25.66,  25.66,  25.66,  25.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.72,  26.73,  25.82,  25.24,  25.07,  25.05,  25.04,  25.04,  25.04,  25.04,  25.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.34,  26.35,  25.46,  24.86,  24.65,  24.61,  24.61,  24.61,  24.61,  24.61,  24.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.09,  26.11,  25.23,  24.62,  24.34,  24.28,  24.27,  24.27,  24.27,  24.27,  24.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.93,  25.95,  25.08,  24.46,  24.11,  24.00,  23.99,  23.99,  23.99,  23.99,  23.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.83,  25.85,  24.99,  24.36,  23.95,  23.77,  23.74,  23.74,  23.74,  23.74,  23.74, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.75,  26.76,  25.78,  24.93,  24.29,  23.84,  23.58,  23.51,  23.51,  23.50,  23.50,  23.50, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.70,  26.71,  25.73,  24.89,  24.25,  23.76,  23.43,  23.31,  23.29,  23.29,  23.29,  23.29, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.67,  25.70,  24.85,  24.21,  23.70,  23.32,  23.12,  23.08,  23.07,  23.07,  23.07, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.64,  25.67,  24.83,  24.19,  23.67,  23.24,  22.96,  22.88,  22.87,  22.87,  22.87, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.61,  26.62,  25.65,  24.81,  24.17,  23.64,  23.18,  22.84,  22.69,  22.67,  22.66,  22.66, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.60,  25.63,  24.80,  24.16,  23.62,  23.14,  22.74,  22.52,  22.47,  22.46,  22.46, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.75,  35.83,  35.22,  35.07,  35.05,  35.05,  35.05,  35.05,  35.05,  35.05,  35.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.06,  33.07,  32.14,  31.55,  31.39,  31.38,  31.37,  31.37,  31.37,  31.37,  31.37,  31.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.73,  30.74,  29.81,  29.22,  29.08,  29.06,  29.06,  29.06,  29.06,  29.06,  29.06,  29.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.25,  29.26,  28.34,  27.76,  27.61,  27.59,  27.59,  27.59,  27.59,  27.59,  27.59,  27.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.32,  28.33,  27.41,  26.83,  26.68,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.72,  27.73,  26.82,  26.24,  26.07,  26.05,  26.04,  26.04,  26.04,  26.04,  26.04,  26.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.34,  27.35,  26.46,  25.86,  25.65,  25.61,  25.61,  25.61,  25.61,  25.61,  25.61,  25.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.09,  27.11,  26.23,  25.62,  25.34,  25.28,  25.27,  25.27,  25.27,  25.27,  25.27,  25.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.93,  26.95,  26.08,  25.46,  25.11,  25.00,  24.99,  24.99,  24.99,  24.99,  24.99,  24.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.83,  26.85,  25.99,  25.36,  24.95,  24.77,  24.74,  24.74,  24.74,  24.74,  24.74,  24.74, 
      32.75,  31.75,  30.75,  29.75,  28.75,  27.76,  26.78,  25.93,  25.29,  24.84,  24.58,  24.51,  24.51,  24.50,  24.50,  24.50,  24.50, 
      32.70,  31.70,  30.70,  29.70,  28.70,  27.71,  26.73,  25.89,  25.25,  24.76,  24.43,  24.31,  24.29,  24.29,  24.29,  24.29,  24.29, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.67,  26.70,  25.85,  25.21,  24.70,  24.32,  24.12,  24.08,  24.07,  24.07,  24.07,  24.07, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.64,  26.67,  25.83,  25.19,  24.67,  24.24,  23.96,  23.88,  23.87,  23.87,  23.87,  23.87, 
      32.61,  31.61,  30.61,  29.61,  28.61,  27.62,  26.65,  25.81,  25.17,  24.64,  24.18,  23.84,  23.69,  23.67,  23.66,  23.66,  23.66, 
      32.59,  31.59,  30.59,  29.59,  28.59,  27.60,  26.63,  25.80,  25.16,  24.62,  24.14,  23.74,  23.52,  23.47,  23.46,  23.46,  23.46, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.75,  36.83,  36.22,  36.07,  36.05,  36.05,  36.05,  36.05,  36.05,  36.05,  36.05,  36.05, 
      39.06,  38.06,  37.06,  36.06,  35.06,  34.07,  33.14,  32.55,  32.39,  32.38,  32.37,  32.37,  32.37,  32.37,  32.37,  32.37,  32.37, 
      36.73,  35.73,  34.73,  33.73,  32.73,  31.74,  30.81,  30.22,  30.08,  30.06,  30.06,  30.06,  30.06,  30.06,  30.06,  30.06,  30.06, 
      35.25,  34.25,  33.25,  32.25,  31.25,  30.26,  29.34,  28.76,  28.61,  28.59,  28.59,  28.59,  28.59,  28.59,  28.59,  28.59,  28.59, 
      34.32,  33.32,  32.32,  31.32,  30.32,  29.33,  28.41,  27.83,  27.68,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66, 
      33.72,  32.72,  31.72,  30.72,  29.72,  28.73,  27.82,  27.24,  27.07,  27.05,  27.04,  27.04,  27.04,  27.04,  27.04,  27.04,  27.04, 
      33.34,  32.34,  31.34,  30.34,  29.34,  28.35,  27.46,  26.86,  26.65,  26.61,  26.61,  26.61,  26.61,  26.61,  26.61,  26.61,  26.61, 
      33.09,  32.09,  31.09,  30.09,  29.09,  28.11,  27.23,  26.62,  26.34,  26.28,  26.27,  26.27,  26.27,  26.27,  26.27,  26.27,  26.27, 
      32.93,  31.93,  30.93,  29.93,  28.93,  27.95,  27.08,  26.46,  26.11,  26.00,  25.99,  25.99,  25.99,  25.99,  25.99,  25.99,  25.99, 
      32.83,  31.83,  30.83,  29.83,  28.83,  27.85,  26.99,  26.36,  25.95,  25.77,  25.74,  25.74,  25.74,  25.74,  25.74,  25.74,  25.74, 
      32.75,  31.75,  30.75,  29.75,  28.76,  27.78,  26.93,  26.29,  25.84,  25.58,  25.51,  25.51,  25.50,  25.50,  25.50,  25.50,  25.50, 
      32.70,  31.70,  30.70,  29.70,  28.71,  27.73,  26.89,  26.25,  25.76,  25.43,  25.31,  25.29,  25.29,  25.29,  25.29,  25.29,  25.29, 
      32.67,  31.67,  30.67,  29.67,  28.67,  27.70,  26.85,  26.21,  25.70,  25.32,  25.12,  25.08,  25.07,  25.07,  25.07,  25.07,  25.07, 
      32.64,  31.64,  30.64,  29.64,  28.64,  27.67,  26.83,  26.19,  25.67,  25.24,  24.96,  24.88,  24.87,  24.87,  24.87,  24.87,  24.87, 
      32.61,  31.61,  30.61,  29.61,  28.62,  27.65,  26.81,  26.17,  25.64,  25.18,  24.84,  24.69,  24.67,  24.66,  24.66,  24.66,  24.66, 
      32.59,  31.59,  30.59,  29.59,  28.60,  27.63,  26.80,  26.16,  25.62,  25.14,  24.74,  24.52,  24.47,  24.46,  24.46,  24.46,  24.46, 
      42.75,  41.75,  40.75,  39.75,  38.75,  37.83,  37.22,  37.07,  37.05,  37.05,  37.05,  37.05,  37.05,  37.05,  37.05,  37.05,  37.05, 
      39.06,  38.06,  37.06,  36.06,  35.07,  34.14,  33.55,  33.39,  33.38,  33.37,  33.37,  33.37,  33.37,  33.37,  33.37,  33.37,  33.37, 
      36.73,  35.73,  34.73,  33.73,  32.74,  31.81,  31.22,  31.08,  31.06,  31.06,  31.06,  31.06,  31.06,  31.06,  31.06,  31.06,  31.06, 
      35.25,  34.25,  33.25,  32.25,  31.26,  30.34,  29.76,  29.61,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59, 
      34.32,  33.32,  32.32,  31.32,  30.33,  29.41,  28.83,  28.68,  28.66,  28.66,  28.66,  28.66,  28.66,  28.66,  28.66,  28.66,  28.66, 
      33.72,  32.72,  31.72,  30.72,  29.73,  28.82,  28.24,  28.07,  28.05,  28.04,  28.04,  28.04,  28.04,  28.04,  28.04,  28.04,  28.04, 
      33.34,  32.34,  31.34,  30.34,  29.35,  28.46,  27.86,  27.65,  27.61,  27.61,  27.61,  27.61,  27.61,  27.61,  27.61,  27.61,  27.61, 
      33.09,  32.09,  31.09,  30.09,  29.11,  28.23,  27.62,  27.34,  27.28,  27.27,  27.27,  27.27,  27.27,  27.27,  27.27,  27.27,  27.27, 
      32.93,  31.93,  30.93,  29.93,  28.95,  28.08,  27.46,  27.11,  27.00,  26.99,  26.99,  26.99,  26.99,  26.99,  26.99,  26.99,  26.99, 
      32.83,  31.83,  30.83,  29.83,  28.85,  27.99,  27.36,  26.95,  26.77,  26.74,  26.74,  26.74,  26.74,  26.74,  26.74,  26.74,  26.74, 
      32.75,  31.75,  30.75,  29.76,  28.78,  27.93,  27.29,  26.84,  26.58,  26.51,  26.51,  26.50,  26.50,  26.50,  26.50,  26.50,  26.50, 
      32.70,  31.70,  30.70,  29.71,  28.73,  27.89,  27.25,  26.76,  26.43,  26.31,  26.29,  26.29,  26.29,  26.29,  26.29,  26.29,  26.29, 
      32.67,  31.67,  30.67,  29.67,  28.70,  27.85,  27.21,  26.70,  26.32,  26.12,  26.08,  26.07,  26.07,  26.07,  26.07,  26.07,  26.07, 
      32.64,  31.64,  30.64,  29.64,  28.67,  27.83,  27.19,  26.67,  26.24,  25.96,  25.88,  25.87,  25.87,  25.87,  25.87,  25.87,  25.87, 
      32.61,  31.61,  30.61,  29.62,  28.65,  27.81,  27.17,  26.64,  26.18,  25.84,  25.69,  25.67,  25.66,  25.66,  25.66,  25.66,  25.66, 
      32.59,  31.59,  30.59,  29.60,  28.63,  27.80,  27.16,  26.62,  26.14,  25.74,  25.52,  25.47,  25.46,  25.46,  25.46,  25.46,  25.46, 
      42.75,  41.75,  40.75,  39.75,  38.83,  38.22,  38.07,  38.05,  38.05,  38.05,  38.05,  38.05,  38.05,  38.05,  38.05,  38.05,  38.05, 
      39.06,  38.06,  37.06,  36.07,  35.14,  34.55,  34.39,  34.38,  34.37,  34.37,  34.37,  34.37,  34.37,  34.37,  34.37,  34.37,  34.37, 
      36.73,  35.73,  34.73,  33.74,  32.81,  32.22,  32.08,  32.06,  32.06,  32.06,  32.06,  32.06,  32.06,  32.06,  32.06,  32.06,  32.06, 
      35.25,  34.25,  33.25,  32.26,  31.34,  30.76,  30.61,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59, 
      34.32,  33.32,  32.32,  31.33,  30.41,  29.83,  29.68,  29.66,  29.66,  29.66,  29.66,  29.66,  29.66,  29.66,  29.66,  29.66,  29.66, 
      33.72,  32.72,  31.72,  30.73,  29.82,  29.24,  29.07,  29.05,  29.04,  29.04,  29.04,  29.04,  29.04,  29.04,  29.04,  29.04,  29.04, 
      33.34,  32.34,  31.34,  30.35,  29.46,  28.86,  28.65,  28.61,  28.61,  28.61,  28.61,  28.61,  28.61,  28.61,  28.61,  28.61,  28.61, 
      33.09,  32.09,  31.09,  30.11,  29.23,  28.62,  28.34,  28.28,  28.27,  28.27,  28.27,  28.27,  28.27,  28.27,  28.27,  28.27,  28.27, 
      32.93,  31.93,  30.93,  29.95,  29.08,  28.46,  28.11,  28.00,  27.99,  27.99,  27.99,  27.99,  27.99,  27.99,  27.99,  27.99,  27.99, 
      32.83,  31.83,  30.83,  29.85,  28.99,  28.36,  27.95,  27.77,  27.74,  27.74,  27.74,  27.74,  27.74,  27.74,  27.74,  27.74,  27.74, 
      32.75,  31.75,  30.76,  29.78,  28.93,  28.29,  27.84,  27.58,  27.51,  27.51,  27.50,  27.50,  27.50,  27.50,  27.50,  27.50,  27.50, 
      32.70,  31.70,  30.71,  29.73,  28.89,  28.25,  27.76,  27.43,  27.31,  27.29,  27.29,  27.29,  27.29,  27.29,  27.29,  27.29,  27.29, 
      32.67,  31.67,  30.67,  29.70,  28.85,  28.21,  27.70,  27.32,  27.12,  27.08,  27.07,  27.07,  27.07,  27.07,  27.07,  27.07,  27.07, 
      32.64,  31.64,  30.64,  29.67,  28.83,  28.19,  27.67,  27.24,  26.96,  26.88,  26.87,  26.87,  26.87,  26.87,  26.87,  26.87,  26.87, 
      32.61,  31.61,  30.62,  29.65,  28.81,  28.17,  27.64,  27.18,  26.84,  26.69,  26.67,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66, 
      32.59,  31.59,  30.60,  29.63,  28.80,  28.16,  27.62,  27.14,  26.74,  26.52,  26.47,  26.46,  26.46,  26.46,  26.46,  26.46,  26.46, 
      42.75,  41.75,  40.75,  39.83,  39.22,  39.07,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05,  39.05, 
      39.06,  38.06,  37.07,  36.14,  35.55,  35.39,  35.38,  35.37,  35.37,  35.37,  35.37,  35.37,  35.37,  35.37,  35.37,  35.37,  35.37, 
      36.73,  35.73,  34.74,  33.81,  33.22,  33.08,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06,  33.06, 
      35.25,  34.25,  33.26,  32.34,  31.76,  31.61,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59,  31.59, 
      34.32,  33.32,  32.33,  31.41,  30.83,  30.68,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66,  30.66, 
      33.72,  32.72,  31.73,  30.82,  30.24,  30.07,  30.05,  30.04,  30.04,  30.04,  30.04,  30.04,  30.04,  30.04,  30.04,  30.04,  30.04, 
      33.34,  32.34,  31.35,  30.46,  29.86,  29.65,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61,  29.61, 
      33.09,  32.09,  31.11,  30.23,  29.62,  29.34,  29.28,  29.27,  29.27,  29.27,  29.27,  29.27,  29.27,  29.27,  29.27,  29.27,  29.27, 
      32.93,  31.93,  30.95,  30.08,  29.46,  29.11,  29.00,  28.99,  28.99,  28.99,  28.99,  28.99,  28.99,  28.99,  28.99,  28.99,  28.99, 
      32.83,  31.83,  30.85,  29.99,  29.36,  28.95,  28.77,  28.74,  28.74,  28.74,  28.74,  28.74,  28.74,  28.74,  28.74,  28.74,  28.74, 
      32.75,  31.76,  30.78,  29.93,  29.29,  28.84,  28.58,  28.51,  28.51,  28.50,  28.50,  28.50,  28.50,  28.50,  28.50,  28.50,  28.50, 
      32.70,  31.71,  30.73,  29.89,  29.25,  28.76,  28.43,  28.31,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29,  28.29, 
      32.67,  31.67,  30.70,  29.85,  29.21,  28.70,  28.32,  28.12,  28.08,  28.07,  28.07,  28.07,  28.07,  28.07,  28.07,  28.07,  28.07, 
      32.64,  31.64,  30.67,  29.83,  29.19,  28.67,  28.24,  27.96,  27.88,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87,  27.87, 
      32.61,  31.62,  30.65,  29.81,  29.17,  28.64,  28.18,  27.84,  27.69,  27.67,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66,  27.66, 
      32.59,  31.60,  30.63,  29.80,  29.16,  28.62,  28.14,  27.74,  27.52,  27.47,  27.46,  27.46,  27.46,  27.46,  27.46,  27.46,  27.46};

  setup_cool_interp_grid_(&my_rates->LCII, rank, params, L, log10(coolunit));

}


extern "C" void initialize_cooling_rate_OI(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {13, 15.0, 1.0}, // log10(number-density like)
    {16,  0.8, 0.2}, // log10(temperature)
    {16, -5.0, 1.0}  // log10(H2 number density)
  };

  double L[] =
     {45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.39,  37.39,  36.39,  35.40,  34.46,  33.81,  33.62,  33.60,  33.60, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.42,  30.42,  29.42,  28.50,  27.92,  27.78,  27.76,  27.76, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.60,  25.61,  24.71,  24.20,  24.08,  24.07,  24.07, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.13,  25.13,  24.14,  23.15,  22.28,  21.84,  21.75,  21.74,  21.74, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.53,  22.53,  21.55,  20.72,  20.35,  20.28,  20.28,  20.28, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.45,  22.46,  21.46,  20.49,  19.71,  19.41,  19.36,  19.36,  19.36, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.72,  20.73,  19.77,  19.06,  18.83,  18.79,  18.79,  18.79, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.20,  20.21,  19.28,  18.65,  18.48,  18.46,  18.46,  18.46, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.82,  20.82,  19.83,  18.94,  18.41,  18.28,  18.27,  18.27,  18.27, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.52,  21.52,  20.53,  19.55,  18.70,  18.26,  18.17,  18.16,  18.16,  18.16, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.28,  22.28,  21.29,  20.29,  19.32,  18.52,  18.17,  18.11,  18.10,  18.10,  18.10, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.08,  21.08,  20.09,  19.13,  18.40,  18.11,  18.07,  18.07,  18.07,  18.07, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.90,  20.90,  19.91,  18.98,  18.30,  18.08,  18.05,  18.04,  18.04,  18.04, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.74,  20.74,  19.75,  18.84,  18.23,  18.06,  18.03,  18.03,  18.03,  18.03, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.60,  18.72,  18.18,  18.04,  18.02,  18.02,  18.02,  18.02, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.43,  20.43,  19.45,  18.61,  18.14,  18.03,  18.02,  18.02,  18.02,  18.02, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.39,  37.39,  36.39,  35.40,  34.46,  33.81,  33.63,  33.61,  33.61, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.42,  30.42,  29.42,  28.50,  27.93,  27.79,  27.77,  27.77, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.60,  25.61,  24.71,  24.20,  24.09,  24.08,  24.08, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.13,  25.13,  24.14,  23.15,  22.29,  21.85,  21.76,  21.75,  21.75, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.53,  22.53,  21.55,  20.72,  20.36,  20.29,  20.29,  20.29, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.45,  22.46,  21.46,  20.49,  19.72,  19.42,  19.37,  19.37,  19.36, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.72,  20.73,  19.77,  19.07,  18.83,  18.80,  18.80,  18.80, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.20,  20.21,  19.28,  18.66,  18.49,  18.47,  18.46,  18.46, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.82,  20.82,  19.83,  18.94,  18.41,  18.29,  18.27,  18.27,  18.27, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.52,  21.52,  20.53,  19.55,  18.70,  18.26,  18.18,  18.17,  18.17,  18.17, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.28,  22.28,  21.29,  20.29,  19.32,  18.53,  18.17,  18.11,  18.10,  18.10,  18.10, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.08,  21.08,  20.09,  19.13,  18.40,  18.12,  18.07,  18.07,  18.07,  18.07, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.90,  20.90,  19.91,  18.98,  18.31,  18.08,  18.05,  18.05,  18.04,  18.04, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.74,  20.74,  19.75,  18.84,  18.23,  18.06,  18.03,  18.03,  18.03,  18.03, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.60,  18.72,  18.18,  18.04,  18.03,  18.02,  18.02,  18.02, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.43,  20.43,  19.45,  18.61,  18.14,  18.03,  18.02,  18.02,  18.02,  18.02, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.39,  37.39,  36.39,  35.40,  34.47,  33.88,  33.73,  33.71,  33.71, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.42,  30.42,  29.43,  28.52,  28.00,  27.88,  27.87,  27.86, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.59,  27.59,  26.60,  25.61,  24.74,  24.28,  24.19,  24.18,  24.18, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.13,  25.13,  24.14,  23.16,  22.32,  21.93,  21.86,  21.85,  21.85, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.53,  22.53,  21.56,  20.76,  20.44,  20.39,  20.38,  20.38, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.45,  22.46,  21.46,  20.50,  19.76,  19.50,  19.46,  19.46,  19.46, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.72,  20.73,  19.78,  19.12,  18.91,  18.88,  18.88,  18.88, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.20,  20.21,  19.30,  18.71,  18.55,  18.53,  18.53,  18.53, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.82,  20.82,  19.84,  18.96,  18.45,  18.33,  18.32,  18.32,  18.32, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.52,  21.52,  20.53,  19.55,  18.72,  18.30,  18.21,  18.20,  18.20,  18.20, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.28,  22.28,  21.29,  20.29,  19.33,  18.55,  18.20,  18.13,  18.12,  18.12,  18.12, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.08,  21.08,  20.09,  19.14,  18.42,  18.13,  18.09,  18.08,  18.08,  18.08, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.90,  20.90,  19.91,  18.98,  18.33,  18.09,  18.06,  18.05,  18.05,  18.05, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.74,  20.74,  19.75,  18.85,  18.26,  18.07,  18.04,  18.04,  18.04,  18.04, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.60,  18.73,  18.20,  18.05,  18.03,  18.03,  18.03,  18.03, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.43,  20.43,  19.46,  18.62,  18.16,  18.04,  18.02,  18.02,  18.02,  18.02, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.39,  37.39,  36.40,  35.43,  34.66,  34.37,  34.33,  34.32,  34.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.42,  30.42,  29.46,  28.75,  28.51,  28.48,  28.48,  28.48, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.59,  27.60,  26.60,  25.66,  25.01,  24.82,  24.79,  24.79,  24.79, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.13,  25.13,  24.14,  23.22,  22.63,  22.48,  22.47,  22.47,  22.47, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.53,  22.54,  21.64,  21.12,  21.01,  20.99,  20.99,  20.99, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.45,  22.46,  21.47,  20.60,  20.15,  20.06,  20.05,  20.05,  20.05, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.72,  20.74,  19.91,  19.51,  19.44,  19.43,  19.43,  19.43, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.20,  20.24,  19.45,  19.08,  19.01,  19.00,  19.00,  19.00, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.82,  20.82,  19.87,  19.13,  18.78,  18.70,  18.69,  18.69,  18.69, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.52,  21.52,  20.53,  19.60,  18.90,  18.56,  18.47,  18.46,  18.46,  18.46, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.28,  22.28,  21.29,  20.30,  19.38,  18.73,  18.40,  18.32,  18.31,  18.30,  18.30, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.08,  21.08,  20.10,  19.21,  18.60,  18.29,  18.21,  18.20,  18.20,  18.20, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.90,  20.90,  19.93,  19.07,  18.49,  18.21,  18.14,  18.13,  18.13,  18.13, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.74,  20.74,  19.77,  18.95,  18.41,  18.15,  18.09,  18.09,  18.09,  18.09, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.63,  18.84,  18.34,  18.11,  18.06,  18.06,  18.06,  18.06, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.43,  20.44,  19.49,  18.75,  18.28,  18.08,  18.04,  18.04,  18.04,  18.04, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.39,  37.40,  36.43,  35.66,  35.37,  35.32,  35.32,  35.32,  35.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.42,  30.46,  29.75,  29.51,  29.48,  29.48,  29.48,  29.48, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.59,  28.60,  27.60,  26.66,  26.00,  25.82,  25.79,  25.79,  25.79,  25.79, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.13,  25.14,  24.22,  23.63,  23.47,  23.45,  23.45,  23.45,  23.45, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.53,  23.54,  22.63,  22.10,  21.97,  21.95,  21.95,  21.95,  21.95, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.45,  23.46,  22.47,  21.59,  21.11,  21.00,  20.98,  20.98,  20.98,  20.98, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.72,  21.74,  20.90,  20.47,  20.39,  20.38,  20.37,  20.37,  20.37, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.20,  21.23,  20.42,  20.04,  19.96,  19.95,  19.95,  19.95,  19.95, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.82,  20.86,  20.10,  19.72,  19.63,  19.62,  19.62,  19.62,  19.62, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.52,  21.53,  20.59,  19.86,  19.48,  19.36,  19.35,  19.34,  19.34,  19.34, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.28,  22.29,  21.30,  20.37,  19.68,  19.28,  19.13,  19.10,  19.10,  19.10,  19.10, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.08,  21.10,  20.20,  19.54,  19.11,  18.92,  18.88,  18.87,  18.87,  18.87, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.90,  20.92,  20.05,  19.41,  18.97,  18.73,  18.67,  18.66,  18.66,  18.66, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.74,  20.77,  19.92,  19.30,  18.84,  18.57,  18.49,  18.48,  18.48,  18.48, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.62,  19.81,  19.20,  18.74,  18.44,  18.35,  18.33,  18.33,  18.33, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.43,  20.48,  19.70,  19.11,  18.65,  18.34,  18.24,  18.23,  18.23,  18.23, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.39,  38.40,  37.43,  36.66,  36.37,  36.32,  36.32,  36.32,  36.32,  36.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.42,  31.46,  30.75,  30.51,  30.48,  30.48,  30.48,  30.48,  30.48, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.59,  29.60,  28.60,  27.66,  27.00,  26.81,  26.78,  26.78,  26.78,  26.78,  26.78, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.13,  26.14,  25.22,  24.62,  24.45,  24.38,  24.36,  24.36,  24.36,  24.36, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.53,  24.54,  23.63,  23.09,  22.88,  22.78,  22.77,  22.77,  22.77,  22.77, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.45,  24.46,  23.47,  22.59,  22.10,  21.96,  21.94,  21.94,  21.94,  21.94,  21.94, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.72,  22.74,  21.89,  21.47,  21.38,  21.37,  21.37,  21.37,  21.37,  21.37, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.20,  22.23,  21.42,  21.04,  20.96,  20.95,  20.95,  20.95,  20.95,  20.95, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.82,  21.86,  21.10,  20.72,  20.63,  20.62,  20.62,  20.62,  20.62,  20.62, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.52,  22.53,  21.59,  20.86,  20.48,  20.36,  20.34,  20.34,  20.34,  20.34,  20.34, 
      28.28,  27.28,  26.28,  25.28,  24.28,  23.29,  22.30,  21.37,  20.68,  20.27,  20.12,  20.10,  20.10,  20.09,  20.09,  20.09, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.08,  22.10,  21.20,  20.54,  20.10,  19.91,  19.87,  19.87,  19.87,  19.87,  19.87, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.90,  21.92,  21.05,  20.41,  19.96,  19.72,  19.66,  19.65,  19.65,  19.65,  19.65, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.74,  21.76,  20.92,  20.30,  19.84,  19.54,  19.45,  19.44,  19.44,  19.44,  19.44, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.58,  21.62,  20.81,  20.20,  19.72,  19.39,  19.25,  19.23,  19.23,  19.23,  19.23, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.43,  21.48,  20.70,  20.11,  19.62,  19.25,  19.07,  19.03,  19.03,  19.03,  19.03, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.39,  39.40,  38.43,  37.66,  37.37,  37.32,  37.32,  37.32,  37.32,  37.32,  37.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.42,  32.46,  31.75,  31.51,  31.48,  31.47,  31.47,  31.47,  31.47,  31.47, 
      35.59,  34.59,  33.59,  32.59,  31.59,  30.60,  29.60,  28.66,  28.00,  27.81,  27.76,  27.70,  27.68,  27.68,  27.68,  27.68, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.13,  27.14,  26.22,  25.62,  25.43,  25.21,  25.00,  24.95,  24.95,  24.94,  24.94, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.53,  25.54,  24.63,  24.09,  23.86,  23.73,  23.70,  23.70,  23.70,  23.70,  23.70, 
      30.45,  29.45,  28.45,  27.45,  26.45,  25.46,  24.47,  23.59,  23.10,  22.96,  22.94,  22.94,  22.94,  22.94,  22.94,  22.94, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.72,  23.74,  22.89,  22.47,  22.38,  22.37,  22.37,  22.37,  22.37,  22.37,  22.37, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.20,  23.23,  22.42,  22.04,  21.96,  21.95,  21.95,  21.95,  21.95,  21.95,  21.95, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.82,  22.86,  22.10,  21.72,  21.63,  21.62,  21.62,  21.62,  21.62,  21.62,  21.62, 
      28.52,  27.52,  26.52,  25.52,  24.52,  23.53,  22.59,  21.86,  21.48,  21.36,  21.34,  21.34,  21.34,  21.34,  21.34,  21.34, 
      28.28,  27.28,  26.28,  25.28,  24.29,  23.30,  22.37,  21.68,  21.27,  21.12,  21.10,  21.10,  21.09,  21.09,  21.09,  21.09, 
      28.08,  27.08,  26.08,  25.08,  24.08,  23.10,  22.20,  21.54,  21.10,  20.91,  20.87,  20.87,  20.87,  20.87,  20.87,  20.87, 
      27.90,  26.90,  25.90,  24.90,  23.90,  22.92,  22.05,  21.41,  20.96,  20.72,  20.66,  20.65,  20.65,  20.65,  20.65,  20.65, 
      27.74,  26.74,  25.74,  24.74,  23.74,  22.76,  21.92,  21.30,  20.84,  20.54,  20.45,  20.44,  20.44,  20.44,  20.44,  20.44, 
      27.58,  26.58,  25.58,  24.58,  23.58,  22.62,  21.81,  21.20,  20.72,  20.39,  20.25,  20.23,  20.23,  20.23,  20.23,  20.23, 
      27.43,  26.43,  25.43,  24.43,  23.43,  22.48,  21.70,  21.11,  20.62,  20.25,  20.07,  20.03,  20.03,  20.02,  20.02,  20.02, 
      45.39,  44.39,  43.39,  42.39,  41.39,  40.40,  39.43,  38.66,  38.37,  38.32,  38.32,  38.32,  38.32,  38.32,  38.32,  38.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.42,  33.46,  32.75,  32.51,  32.48,  32.47,  32.45,  32.44,  32.44,  32.44,  32.44, 
      35.59,  34.59,  33.59,  32.59,  31.60,  30.60,  29.66,  29.00,  28.81,  28.75,  28.55,  28.28,  28.21,  28.20,  28.20,  28.20, 
      33.13,  32.13,  31.13,  30.13,  29.13,  28.14,  27.22,  26.62,  26.43,  26.16,  25.75,  25.62,  25.60,  25.60,  25.60,  25.60, 
      31.53,  30.53,  29.53,  28.53,  27.53,  26.54,  25.63,  25.09,  24.86,  24.73,  24.70,  24.70,  24.70,  24.70,  24.70,  24.70, 
      30.45,  29.45,  28.45,  27.45,  26.46,  25.47,  24.59,  24.10,  23.96,  23.94,  23.93,  23.93,  23.93,  23.93,  23.93,  23.93, 
      29.72,  28.72,  27.72,  26.72,  25.72,  24.74,  23.89,  23.47,  23.38,  23.37,  23.37,  23.37,  23.37,  23.37,  23.37,  23.37, 
      29.20,  28.20,  27.20,  26.20,  25.20,  24.23,  23.42,  23.04,  22.96,  22.95,  22.95,  22.95,  22.95,  22.95,  22.95,  22.95, 
      28.82,  27.82,  26.82,  25.82,  24.82,  23.86,  23.10,  22.72,  22.63,  22.62,  22.62,  22.62,  22.62,  22.62,  22.62,  22.62, 
      28.52,  27.52,  26.52,  25.52,  24.53,  23.59,  22.86,  22.47,  22.36,  22.34,  22.34,  22.34,  22.34,  22.34,  22.34,  22.34, 
      28.28,  27.28,  26.28,  25.29,  24.30,  23.37,  22.68,  22.27,  22.12,  22.10,  22.09,  22.09,  22.09,  22.09,  22.09,  22.09, 
      28.08,  27.08,  26.08,  25.08,  24.10,  23.20,  22.53,  22.10,  21.91,  21.87,  21.86,  21.86,  21.86,  21.86,  21.86,  21.86, 
      27.90,  26.90,  25.90,  24.90,  23.92,  23.05,  22.41,  21.96,  21.72,  21.66,  21.65,  21.65,  21.65,  21.65,  21.65,  21.65, 
      27.74,  26.74,  25.74,  24.74,  23.76,  22.92,  22.30,  21.83,  21.54,  21.45,  21.44,  21.44,  21.44,  21.44,  21.44,  21.44, 
      27.58,  26.58,  25.58,  24.58,  23.62,  22.81,  22.20,  21.72,  21.39,  21.25,  21.23,  21.23,  21.23,  21.23,  21.23,  21.23, 
      27.43,  26.43,  25.43,  24.43,  23.48,  22.70,  22.11,  21.62,  21.25,  21.07,  21.03,  21.02,  21.02,  21.02,  21.02,  21.02, 
      45.39,  44.39,  43.39,  42.39,  41.40,  40.43,  39.66,  39.37,  39.32,  39.32,  39.32,  39.32,  39.32,  39.32,  39.32,  39.32, 
      39.42,  38.42,  37.42,  36.42,  35.42,  34.46,  33.75,  33.51,  33.48,  33.47,  33.41,  33.28,  33.23,  33.23,  33.23,  33.23, 
      35.59,  34.59,  33.59,  32.60,  31.60,  30.66,  30.00,  29.81,  29.75,  29.50,  28.87,  28.45,  28.37,  28.36,  28.35,  28.35, 
      33.13,  32.13,  31.13,  30.13,  29.14,  28.22,  27.62,  27.43,  27.16,  26.75,  26.61,  26.59,  26.59,  26.59,  26.59,  26.59, 
      31.53,  30.53,  29.53,  28.53,  27.54,  26.63,  26.09,  25.86,  25.72,  25.70,  25.69,  25.69,  25.69,  25.69,  25.69,  25.69, 
      30.45,  29.45,  28.45,  27.46,  26.47,  25.59,  25.09,  24.95,  24.92,  24.92,  24.92,  24.92,  24.92,  24.92,  24.92,  24.92, 
      29.72,  28.72,  27.72,  26.72,  25.74,  24.89,  24.45,  24.35,  24.34,  24.34,  24.34,  24.34,  24.34,  24.34,  24.34,  24.34, 
      29.20,  28.20,  27.20,  26.20,  25.23,  24.42,  24.01,  23.93,  23.92,  23.91,  23.91,  23.91,  23.91,  23.91,  23.91,  23.91, 
      28.82,  27.82,  26.82,  25.82,  24.86,  24.09,  23.70,  23.60,  23.59,  23.59,  23.59,  23.59,  23.59,  23.59,  23.59,  23.59, 
      28.52,  27.52,  26.52,  25.53,  24.59,  23.85,  23.45,  23.34,  23.32,  23.32,  23.32,  23.32,  23.32,  23.32,  23.32,  23.32, 
      28.28,  27.28,  26.29,  25.30,  24.37,  23.67,  23.26,  23.10,  23.08,  23.08,  23.08,  23.08,  23.08,  23.08,  23.08,  23.08, 
      28.08,  27.08,  26.08,  25.10,  24.19,  23.52,  23.09,  22.90,  22.86,  22.85,  22.85,  22.85,  22.85,  22.85,  22.85,  22.85, 
      27.90,  26.90,  25.90,  24.92,  24.04,  23.40,  22.95,  22.71,  22.65,  22.64,  22.64,  22.64,  22.64,  22.64,  22.64,  22.64, 
      27.74,  26.74,  25.74,  24.76,  23.91,  23.29,  22.83,  22.54,  22.45,  22.43,  22.43,  22.43,  22.43,  22.43,  22.43,  22.43, 
      27.58,  26.58,  25.58,  24.62,  23.80,  23.19,  22.72,  22.38,  22.25,  22.23,  22.23,  22.23,  22.23,  22.23,  22.23,  22.23, 
      27.43,  26.43,  25.43,  24.48,  23.69,  23.10,  22.62,  22.25,  22.07,  22.03,  22.02,  22.02,  22.02,  22.02,  22.02,  22.02, 
      45.39,  44.39,  43.39,  42.40,  41.43,  40.66,  40.37,  40.32,  40.32,  40.32,  40.32,  40.31,  40.31,  40.31,  40.31,  40.31, 
      39.42,  38.42,  37.42,  36.42,  35.46,  34.75,  34.51,  34.48,  34.47,  34.39,  34.05,  33.65,  33.54,  33.53,  33.53,  33.53, 
      35.59,  34.59,  33.60,  32.60,  31.66,  31.00,  30.81,  30.75,  30.49,  29.78,  29.04,  28.75,  28.71,  28.71,  28.71,  28.71, 
      33.13,  32.13,  31.13,  30.14,  29.22,  28.62,  28.43,  28.15,  27.74,  27.61,  27.59,  27.59,  27.59,  27.59,  27.59,  27.59, 
      31.53,  30.53,  29.53,  28.54,  27.63,  27.09,  26.85,  26.70,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66,  26.66, 
      30.45,  29.45,  28.46,  27.47,  26.59,  26.07,  25.86,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82, 
      29.72,  28.72,  27.72,  26.74,  25.89,  25.37,  25.20,  25.18,  25.18,  25.18,  25.18,  25.18,  25.18,  25.18,  25.18,  25.18, 
      29.20,  28.20,  27.20,  26.23,  25.40,  24.90,  24.75,  24.73,  24.72,  24.72,  24.72,  24.72,  24.72,  24.72,  24.72,  24.72, 
      28.82,  27.82,  26.82,  25.86,  25.05,  24.57,  24.42,  24.40,  24.40,  24.40,  24.40,  24.40,  24.40,  24.40,  24.40,  24.40, 
      28.52,  27.52,  26.53,  25.58,  24.80,  24.33,  24.18,  24.16,  24.16,  24.16,  24.16,  24.16,  24.16,  24.16,  24.16,  24.16, 
      28.28,  27.29,  26.29,  25.36,  24.61,  24.15,  23.98,  23.95,  23.95,  23.95,  23.95,  23.95,  23.95,  23.95,  23.95,  23.95, 
      28.08,  27.08,  26.10,  25.18,  24.45,  24.00,  23.80,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76, 
      27.90,  26.90,  25.92,  25.02,  24.33,  23.88,  23.64,  23.58,  23.57,  23.57,  23.57,  23.57,  23.57,  23.57,  23.57,  23.57, 
      27.74,  26.74,  25.76,  24.89,  24.22,  23.77,  23.49,  23.40,  23.39,  23.39,  23.39,  23.39,  23.39,  23.39,  23.39,  23.39, 
      27.58,  26.58,  25.61,  24.76,  24.12,  23.67,  23.35,  23.22,  23.20,  23.20,  23.20,  23.20,  23.20,  23.20,  23.20,  23.20, 
      27.43,  26.43,  25.47,  24.65,  24.03,  23.58,  23.22,  23.05,  23.01,  23.00,  23.00,  23.00,  23.00,  23.00,  23.00,  23.00, 
      45.39,  44.39,  43.40,  42.43,  41.66,  41.37,  41.32,  41.32,  41.32,  41.32,  41.30,  41.25,  41.22,  41.22,  41.22,  41.22, 
      39.42,  38.42,  37.42,  36.46,  35.75,  35.51,  35.48,  35.47,  35.39,  34.99,  34.22,  33.71,  33.59,  33.58,  33.58,  33.58, 
      35.59,  34.60,  33.60,  32.66,  32.00,  31.81,  31.75,  31.49,  30.77,  30.02,  29.71,  29.66,  29.65,  29.65,  29.65,  29.65, 
      33.13,  32.13,  31.14,  30.22,  29.62,  29.43,  29.15,  28.74,  28.61,  28.59,  28.59,  28.59,  28.59,  28.59,  28.59,  28.59, 
      31.53,  30.53,  29.54,  28.63,  28.09,  27.83,  27.67,  27.63,  27.63,  27.63,  27.63,  27.63,  27.63,  27.63,  27.63,  27.63, 
      30.45,  29.46,  28.47,  27.59,  27.04,  26.79,  26.73,  26.72,  26.72,  26.72,  26.72,  26.72,  26.72,  26.72,  26.72,  26.72, 
      29.72,  28.72,  27.74,  26.88,  26.32,  26.07,  26.03,  26.02,  26.02,  26.02,  26.02,  26.02,  26.02,  26.02,  26.02,  26.02, 
      29.20,  28.20,  27.23,  26.39,  25.81,  25.57,  25.52,  25.52,  25.52,  25.52,  25.52,  25.52,  25.52,  25.52,  25.52,  25.52, 
      28.82,  27.82,  26.86,  26.04,  25.46,  25.20,  25.14,  25.14,  25.14,  25.14,  25.14,  25.14,  25.14,  25.14,  25.14,  25.14, 
      28.52,  27.53,  26.58,  25.78,  25.20,  24.92,  24.85,  24.84,  24.84,  24.84,  24.84,  24.84,  24.84,  24.84,  24.84,  24.84, 
      28.29,  27.29,  26.36,  25.58,  25.00,  24.70,  24.62,  24.61,  24.61,  24.61,  24.61,  24.61,  24.61,  24.61,  24.61,  24.61, 
      28.08,  27.10,  26.17,  25.41,  24.84,  24.52,  24.44,  24.43,  24.43,  24.43,  24.43,  24.43,  24.43,  24.43,  24.43,  24.43, 
      27.90,  26.92,  26.02,  25.27,  24.71,  24.38,  24.29,  24.27,  24.27,  24.27,  24.27,  24.27,  24.27,  24.27,  24.27,  24.27, 
      27.74,  26.76,  25.88,  25.15,  24.59,  24.26,  24.15,  24.14,  24.14,  24.13,  24.13,  24.13,  24.13,  24.13,  24.13,  24.13, 
      27.58,  26.61,  25.75,  25.04,  24.49,  24.16,  24.03,  24.00,  24.00,  24.00,  24.00,  24.00,  24.00,  24.00,  24.00,  24.00, 
      27.43,  26.47,  25.63,  24.94,  24.40,  24.06,  23.90,  23.86,  23.86,  23.86,  23.86,  23.86,  23.86,  23.86,  23.86,  23.86, 
      45.39,  44.40,  43.43,  42.66,  42.37,  42.32,  42.32,  42.32,  42.32,  42.30,  42.16,  41.88,  41.78,  41.77,  41.77,  41.77, 
      39.42,  38.42,  37.46,  36.75,  36.51,  36.48,  36.47,  36.39,  35.98,  35.13,  34.25,  33.72,  33.61,  33.59,  33.59,  33.59, 
      35.60,  34.60,  33.66,  33.00,  32.81,  32.75,  32.49,  31.77,  31.02,  30.71,  30.66,  30.65,  30.65,  30.65,  30.65,  30.65, 
      33.13,  32.14,  31.22,  30.62,  30.43,  30.15,  29.74,  29.61,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59,  29.59, 
      31.53,  30.54,  29.63,  29.09,  28.83,  28.67,  28.63,  28.63,  28.63,  28.63,  28.63,  28.63,  28.63,  28.63,  28.63,  28.63, 
      30.46,  29.47,  28.59,  28.04,  27.79,  27.73,  27.72,  27.72,  27.72,  27.72,  27.72,  27.72,  27.72,  27.72,  27.72,  27.72, 
      29.72,  28.74,  27.88,  27.32,  27.07,  27.03,  27.02,  27.02,  27.02,  27.02,  27.02,  27.02,  27.02,  27.02,  27.02,  27.02, 
      29.20,  28.23,  27.39,  26.81,  26.56,  26.52,  26.51,  26.51,  26.51,  26.51,  26.51,  26.51,  26.51,  26.51,  26.51,  26.51, 
      28.82,  27.86,  27.04,  26.46,  26.19,  26.14,  26.13,  26.13,  26.13,  26.13,  26.13,  26.13,  26.13,  26.13,  26.13,  26.13, 
      28.53,  27.58,  26.78,  26.20,  25.90,  25.83,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82,  25.82, 
      28.29,  27.36,  26.58,  25.99,  25.67,  25.57,  25.56,  25.56,  25.56,  25.56,  25.56,  25.56,  25.56,  25.56,  25.56,  25.56, 
      28.10,  27.17,  26.41,  25.83,  25.47,  25.34,  25.32,  25.32,  25.32,  25.32,  25.32,  25.32,  25.32,  25.32,  25.32,  25.32, 
      27.92,  27.02,  26.27,  26.07,  25.30,  25.13,  25.10,  25.09,  25.09,  25.09,  25.09,  25.09,  25.09,  25.09,  25.09,  25.09, 
      27.76,  26.88,  26.15,  25.57,  25.16,  24.94,  24.89,  24.88,  24.88,  24.88,  24.88,  24.88,  24.88,  24.88,  24.88,  24.88, 
      27.61,  26.75,  26.04,  25.47,  25.03,  24.77,  24.69,  24.68,  24.68,  24.68,  24.68,  24.68,  24.68,  24.68,  24.68,  24.68, 
      27.47,  26.63,  25.93,  25.37,  24.91,  24.62,  24.52,  24.51,  24.51,  24.51,  24.51,  24.51,  24.51,  24.51,  24.51,  24.51, 
      45.40,  44.43,  43.66,  43.37,  43.32,  43.32,  43.32,  43.32,  43.30,  43.13,  42.57,  42.05,  41.92,  41.90,  41.90,  41.90, 
      39.42,  38.46,  37.75,  37.51,  37.48,  37.47,  37.39,  36.98,  36.12,  35.16,  34.27,  33.79,  33.69,  33.68,  33.67,  33.67, 
      35.60,  34.66,  34.00,  33.81,  33.75,  33.49,  32.77,  32.02,  31.71,  31.66,  31.65,  31.65,  31.65,  31.65,  31.65,  31.65, 
      33.14,  32.22,  31.62,  31.43,  31.15,  30.74,  30.61,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59,  30.59, 
      31.54,  30.63,  30.09,  29.83,  29.67,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63,  29.63, 
      30.47,  29.59,  29.04,  28.79,  28.73,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72,  28.72, 
      29.74,  28.88,  28.32,  28.07,  28.03,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02,  28.02, 
      29.23,  28.39,  27.81,  27.56,  27.52,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51,  27.51, 
      28.86,  28.04,  27.46,  27.19,  27.14,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13,  27.13, 
      28.58,  27.78,  27.20,  26.90,  26.83,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82,  26.82, 
      28.36,  27.58,  26.99,  26.67,  26.57,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56,  26.56, 
      28.17,  27.41,  26.83,  26.47,  26.34,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32,  26.32, 
      28.02,  27.27,  26.69,  26.30,  26.13,  26.10,  26.09,  26.09,  26.09,  26.09,  26.09,  26.09,  26.09,  26.09,  26.09,  26.09, 
      27.88,  27.15,  26.57,  26.16,  25.94,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88,  25.88, 
      27.75,  27.04,  26.47,  26.03,  25.76,  25.68,  25.67,  25.67,  25.67,  25.67,  25.67,  25.67,  25.67,  25.67,  25.67,  25.67, 
      27.63,  26.93,  26.37,  25.91,  25.60,  25.48,  25.46,  25.46,  25.46,  25.46,  25.46,  25.46,  25.46,  25.46,  25.46,  25.46}; 

  setup_cool_interp_grid_(&my_rates->LOI, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_CO(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {11, 14.0, 0.5}, // log10(number-density like)
    {11,  1.0, 0.2}, // log10(temperature)
    {14, -3.0, 1.0}  // log10(H2I number density)
  };

  double L[] =
     {27.77,  26.77,  25.77,  24.78,  23.80,  22.84,  21.99,  21.40,  21.16,  21.10,  21.09,  21.08,  21.08,  21.08, 
      27.51,  26.51,  25.52,  24.52,  23.54,  22.58,  21.71,  21.04,  20.72,  20.63,  20.60,  20.60,  20.60,  20.60, 
      27.29,  26.29,  25.29,  24.30,  23.31,  22.36,  21.46,  20.72,  20.31,  20.17,  20.13,  20.12,  20.12,  20.12, 
      27.11,  26.11,  25.12,  24.12,  23.14,  22.17,  21.26,  20.46,  19.94,  19.74,  19.68,  19.67,  19.67,  19.67, 
      26.96,  25.96,  24.96,  23.97,  22.98,  22.01,  21.08,  20.23,  19.61,  19.33,  19.25,  19.24,  19.23,  19.23, 
      26.82,  25.82,  24.82,  23.83,  22.84,  21.86,  20.91,  20.02,  19.31,  18.94,  18.83,  18.81,  18.80,  18.80, 
      26.65,  25.65,  24.65,  23.66,  22.66,  21.68,  20.72,  19.80,  19.03,  18.58,  18.43,  18.39,  18.39,  18.39, 
      26.48,  25.49,  24.49,  23.49,  22.49,  21.50,  20.53,  19.60,  18.77,  18.23,  18.03,  17.98,  17.97,  17.97, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.31,  20.33,  19.38,  18.52,  17.91,  17.66,  17.59,  17.58,  17.57, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.17,  18.29,  17.63,  17.32,  17.22,  17.20,  17.20, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.93,  18.98,  18.09,  17.39,  17.02,  16.90,  16.87,  16.86, 
      27.77,  26.77,  25.77,  24.78,  23.80,  22.85,  21.99,  21.41,  21.17,  21.11,  21.10,  21.09,  21.09,  21.09, 
      27.51,  26.51,  25.52,  24.52,  23.54,  22.59,  21.71,  21.05,  20.73,  20.63,  20.61,  20.60,  20.60,  20.60, 
      27.29,  26.29,  25.29,  24.30,  23.31,  22.36,  21.46,  20.73,  20.31,  20.17,  20.14,  20.13,  20.13,  20.13, 
      27.11,  26.11,  25.12,  24.12,  23.14,  22.17,  21.26,  20.46,  19.94,  19.74,  19.69,  19.68,  19.67,  19.67, 
      26.96,  25.96,  24.97,  23.97,  22.98,  22.01,  21.08,  20.23,  19.61,  19.33,  19.25,  19.24,  19.23,  19.23, 
      26.82,  25.82,  24.82,  23.83,  22.84,  21.86,  20.91,  20.02,  19.31,  18.94,  18.83,  18.81,  18.80,  18.80, 
      26.65,  25.65,  24.65,  23.66,  22.66,  21.68,  20.72,  19.81,  19.03,  18.58,  18.43,  18.39,  18.39,  18.39, 
      26.48,  25.49,  24.49,  23.49,  22.49,  21.50,  20.53,  19.60,  18.77,  18.23,  18.03,  17.98,  17.97,  17.97, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.31,  20.33,  19.38,  18.53,  17.92,  17.66,  17.59,  17.58,  17.57, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.17,  18.29,  17.63,  17.32,  17.22,  17.20,  17.20, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.93,  18.98,  18.09,  17.39,  17.02,  16.90,  16.87,  16.86, 
      27.77,  26.77,  25.77,  24.78,  23.80,  22.85,  22.01,  21.43,  21.19,  21.13,  21.12,  21.11,  21.11,  21.11, 
      27.51,  26.51,  25.52,  24.52,  23.54,  22.59,  21.72,  21.06,  20.74,  20.65,  20.63,  20.62,  20.62,  20.62, 
      27.29,  26.29,  25.29,  24.30,  23.32,  22.36,  21.47,  20.74,  20.32,  20.18,  20.15,  20.14,  20.14,  20.14, 
      27.11,  26.11,  25.12,  24.12,  23.14,  22.18,  21.26,  20.47,  19.95,  19.75,  19.70,  19.69,  19.68,  19.68, 
      26.96,  25.96,  24.97,  23.97,  22.98,  22.01,  21.08,  20.23,  19.61,  19.33,  19.26,  19.24,  19.24,  19.24, 
      26.82,  25.82,  24.82,  23.83,  22.84,  21.86,  20.91,  20.02,  19.31,  18.94,  18.83,  18.81,  18.80,  18.80, 
      26.65,  25.65,  24.65,  23.66,  22.67,  21.68,  20.72,  19.81,  19.03,  18.58,  18.43,  18.39,  18.39,  18.39, 
      26.48,  25.49,  24.49,  23.49,  22.49,  21.50,  20.53,  19.60,  18.77,  18.23,  18.03,  17.98,  17.97,  17.97, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.31,  20.33,  19.38,  18.53,  17.92,  17.66,  17.59,  17.58,  17.57, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.17,  18.29,  17.63,  17.32,  17.22,  17.20,  17.20, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.93,  18.98,  18.09,  17.39,  17.02,  16.90,  16.87,  16.86, 
      27.77,  26.77,  25.78,  24.79,  23.81,  22.88,  22.05,  21.49,  21.26,  21.20,  21.19,  21.18,  21.18,  21.18, 
      27.51,  26.51,  25.52,  24.53,  23.55,  22.61,  21.75,  21.10,  20.79,  20.69,  20.67,  20.66,  20.66,  20.66, 
      27.29,  26.29,  25.29,  24.30,  23.32,  22.37,  21.49,  20.76,  20.35,  20.21,  20.17,  20.17,  20.16,  20.16, 
      27.11,  26.11,  25.12,  24.13,  23.14,  22.19,  21.28,  20.48,  19.96,  19.76,  19.71,  19.70,  19.70,  19.70, 
      26.96,  25.96,  24.97,  23.97,  22.99,  22.02,  21.09,  20.24,  19.62,  19.34,  19.27,  19.25,  19.25,  19.25, 
      26.82,  25.82,  24.82,  23.83,  22.84,  21.87,  20.92,  20.03,  19.32,  18.95,  18.84,  18.82,  18.81,  18.81, 
      26.65,  25.65,  24.65,  23.66,  22.67,  21.68,  20.72,  19.81,  19.04,  18.58,  18.44,  18.40,  18.40,  18.40, 
      26.48,  25.49,  24.49,  23.49,  22.49,  21.51,  20.53,  19.60,  18.78,  18.24,  18.04,  17.99,  17.98,  17.98, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.31,  20.34,  19.39,  18.53,  17.92,  17.67,  17.60,  17.58,  17.58, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.17,  18.29,  17.63,  17.32,  17.23,  17.20,  17.20, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.94,  18.98,  18.09,  17.39,  17.03,  16.91,  16.88,  16.87, 
      27.77,  26.77,  25.78,  24.79,  23.83,  22.93,  22.16,  21.66,  21.45,  21.39,  21.38,  21.37,  21.37,  21.37, 
      27.51,  26.51,  25.52,  24.54,  23.57,  22.65,  21.82,  21.21,  20.92,  20.83,  20.81,  20.80,  20.80,  20.80, 
      27.29,  26.29,  25.29,  24.31,  23.34,  22.40,  21.54,  20.83,  20.43,  20.29,  20.26,  20.25,  20.25,  20.25, 
      27.11,  26.11,  25.12,  24.13,  23.15,  22.20,  21.31,  20.53,  20.02,  19.82,  19.77,  19.76,  19.76,  19.76, 
      26.96,  25.96,  24.97,  23.98,  23.00,  22.03,  21.11,  20.27,  19.66,  19.38,  19.31,  19.29,  19.29,  19.29, 
      26.82,  25.82,  24.83,  23.83,  22.85,  21.87,  20.93,  20.05,  19.34,  18.97,  18.86,  18.84,  18.83,  18.83, 
      26.65,  25.65,  24.66,  23.66,  22.67,  21.69,  20.73,  19.82,  19.05,  18.60,  18.45,  18.42,  18.41,  18.41, 
      26.49,  25.49,  24.49,  23.49,  22.49,  21.51,  20.54,  19.61,  18.79,  18.24,  18.04,  18.00,  17.99,  17.98, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.32,  20.34,  19.39,  18.53,  17.92,  17.67,  17.60,  17.58,  17.58, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.18,  18.30,  17.63,  17.32,  17.23,  17.20,  17.20, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.94,  18.98,  18.09,  17.39,  17.03,  16.91,  16.88,  16.87, 
      27.77,  26.77,  25.78,  24.81,  23.87,  23.01,  22.34,  21.92,  21.75,  21.69,  21.68,  21.67,  21.67,  21.67, 
      27.51,  26.52,  25.52,  24.54,  23.59,  22.71,  21.95,  21.42,  21.16,  21.08,  21.05,  21.05,  21.05,  21.05, 
      27.29,  26.29,  25.30,  24.32,  23.36,  22.45,  21.63,  20.98,  20.62,  20.49,  20.46,  20.45,  20.45,  20.45, 
      27.11,  26.12,  25.12,  24.14,  23.17,  22.24,  21.37,  20.63,  20.15,  19.96,  19.92,  19.91,  19.90,  19.90, 
      26.96,  25.97,  24.97,  23.98,  23.01,  22.06,  21.15,  20.34,  19.74,  19.48,  19.41,  19.39,  19.39,  19.39, 
      26.82,  25.82,  24.83,  23.84,  22.86,  21.89,  20.96,  20.09,  19.40,  19.03,  18.93,  18.91,  18.90,  18.90, 
      26.65,  25.65,  24.66,  23.66,  22.67,  21.70,  20.75,  19.85,  19.09,  18.65,  18.50,  18.47,  18.46,  18.46, 
      26.49,  25.49,  24.49,  23.49,  22.50,  21.51,  20.55,  19.62,  18.81,  18.27,  18.08,  18.03,  18.02,  18.02, 
      26.30,  25.30,  24.30,  23.30,  22.31,  21.32,  20.34,  19.40,  18.55,  17.94,  17.69,  17.62,  17.61,  17.61, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.14,  19.18,  18.31,  17.65,  17.33,  17.24,  17.22,  17.22, 
      25.91,  24.91,  23.91,  22.91,  21.91,  20.92,  19.94,  18.98,  18.10,  17.40,  17.04,  16.92,  16.89,  16.88, 
      27.77,  26.78,  25.79,  24.82,  23.92,  23.15,  22.59,  22.26,  22.12,  22.07,  22.05,  22.04,  22.04,  22.04, 
      27.51,  26.52,  25.53,  24.56,  23.63,  22.81,  22.15,  21.71,  21.49,  21.42,  21.39,  21.39,  21.38,  21.38, 
      27.29,  26.29,  25.30,  24.33,  23.39,  22.52,  21.77,  21.22,  20.91,  20.80,  20.76,  20.76,  20.75,  20.75, 
      27.11,  26.12,  25.13,  24.15,  23.19,  22.29,  21.47,  20.81,  20.39,  20.22,  20.18,  20.16,  20.16,  20.16, 
      26.96,  25.97,  24.97,  23.99,  23.02,  22.09,  21.22,  20.46,  19.93,  19.69,  19.62,  19.60,  19.60,  19.60, 
      26.82,  25.83,  24.83,  23.84,  22.87,  21.91,  21.01,  20.18,  19.53,  19.19,  19.09,  19.07,  19.06,  19.06, 
      26.65,  25.66,  24.66,  23.67,  22.68,  21.72,  20.78,  19.91,  19.19,  18.76,  18.62,  18.59,  18.58,  18.58, 
      26.49,  25.49,  24.49,  23.49,  22.50,  21.53,  20.57,  19.66,  18.87,  18.36,  18.16,  18.12,  18.11,  18.10, 
      26.30,  25.30,  24.30,  23.31,  22.31,  21.33,  20.36,  19.43,  18.59,  18.00,  17.74,  17.68,  17.66,  17.66, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.14,  19.20,  18.33,  17.68,  17.37,  17.27,  17.25,  17.24, 
      25.91,  24.91,  23.91,  22.91,  21.92,  20.92,  19.95,  19.00,  18.12,  17.43,  17.06,  16.94,  16.91,  16.90, 
      27.77,  26.78,  25.80,  24.85,  24.01,  23.35,  22.90,  22.64,  22.52,  22.47,  22.45,  22.44,  22.44,  22.44, 
      27.51,  26.52,  25.53,  24.58,  23.70,  22.96,  22.41,  22.05,  21.87,  21.80,  21.77,  21.77,  21.76,  21.76, 
      27.29,  26.29,  25.31,  24.34,  23.43,  22.63,  21.98,  21.51,  21.26,  21.15,  21.12,  21.11,  21.11,  21.11, 
      27.11,  26.12,  25.13,  24.16,  23.22,  22.36,  21.62,  21.05,  20.70,  20.55,  20.51,  20.49,  20.49,  20.49, 
      26.96,  25.97,  24.97,  23.99,  23.04,  22.14,  21.32,  20.65,  20.19,  19.99,  19.92,  19.90,  19.90,  19.90, 
      26.82,  25.82,  24.83,  23.84,  22.88,  21.94,  21.07,  20.31,  19.75,  19.46,  19.36,  19.34,  19.33,  19.33, 
      26.65,  25.66,  24.66,  23.67,  22.69,  21.73,  20.83,  20.01,  19.35,  18.98,  18.85,  18.82,  18.81,  18.81, 
      26.49,  25.49,  24.49,  23.50,  22.51,  21.54,  20.60,  19.73,  19.00,  18.53,  18.34,  18.30,  18.29,  18.28, 
      26.30,  25.30,  24.30,  23.31,  22.32,  21.33,  20.38,  19.48,  18.69,  18.13,  17.89,  17.82,  17.80,  17.80, 
      26.11,  25.11,  24.11,  23.11,  22.12,  21.13,  20.16,  19.23,  18.40,  17.78,  17.47,  17.37,  17.35,  17.34, 
      25.91,  24.91,  23.91,  22.91,  21.92,  20.93,  19.96,  19.03,  18.17,  17.50,  17.14,  17.01,  16.98,  16.97, 
      27.77,  26.78,  25.81,  24.91,  24.15,  23.61,  23.25,  23.04,  22.94,  22.90,  22.88,  22.87,  22.87,  22.87, 
      27.51,  26.52,  25.55,  24.62,  23.81,  23.17,  22.71,  22.42,  22.27,  22.21,  22.18,  22.18,  22.17,  22.17, 
      27.29,  26.29,  25.31,  24.37,  23.51,  22.79,  22.23,  21.85,  21.63,  21.54,  21.51,  21.50,  21.50,  21.49, 
      27.11,  26.12,  25.13,  24.17,  23.27,  22.48,  21.83,  21.34,  21.05,  20.92,  20.88,  20.86,  20.86,  20.86, 
      26.96,  25.97,  24.98,  24.00,  23.07,  22.21,  21.48,  20.90,  20.52,  20.34,  20.28,  20.26,  20.25,  20.25, 
      26.82,  25.82,  24.83,  23.85,  22.89,  21.99,  21.18,  20.51,  20.03,  19.78,  19.69,  19.67,  19.66,  19.66, 
      26.65,  25.65,  24.66,  23.67,  22.70,  21.76,  20.90,  20.16,  19.60,  19.28,  19.16,  19.12,  19.11,  19.11, 
      26.49,  25.49,  24.49,  23.50,  22.51,  21.55,  20.65,  19.84,  19.20,  18.80,  18.63,  18.57,  18.56,  18.55, 
      26.30,  25.30,  24.30,  23.31,  22.32,  21.34,  20.41,  19.55,  18.85,  18.37,  18.14,  18.06,  18.04,  18.04, 
      26.11,  25.11,  24.11,  23.11,  22.12,  21.14,  20.18,  19.29,  18.53,  17.98,  17.70,  17.59,  17.56,  17.55, 
      25.91,  24.91,  23.91,  22.91,  21.92,  20.94,  19.98,  19.08,  18.28,  17.67,  17.32,  17.20,  17.16,  17.15, 
      27.78,  26.79,  25.84,  25.00,  24.37,  23.93,  23.63,  23.46,  23.37,  23.33,  23.31,  23.30,  23.30,  23.30, 
      27.52,  26.53,  25.57,  24.69,  23.97,  23.44,  23.05,  22.81,  22.68,  22.62,  22.60,  22.59,  22.59,  22.59, 
      27.29,  26.30,  25.33,  24.42,  23.63,  23.01,  22.53,  22.21,  22.02,  21.95,  21.92,  21.91,  21.91,  21.90, 
      27.11,  26.12,  25.14,  24.21,  23.36,  22.64,  22.08,  21.67,  21.42,  21.31,  21.27,  21.26,  21.26,  21.26, 
      26.96,  25.97,  24.98,  24.02,  23.12,  22.33,  21.68,  21.19,  20.86,  20.71,  20.66,  20.64,  20.64,  20.63, 
      26.82,  25.83,  24.83,  23.86,  22.92,  22.07,  21.34,  20.75,  20.35,  20.14,  20.06,  20.04,  20.03,  20.03, 
      26.65,  25.66,  24.66,  23.67,  22.71,  21.81,  21.01,  20.37,  19.90,  19.63,  19.51,  19.48,  19.46,  19.46, 
      26.49,  25.49,  24.49,  23.50,  22.52,  21.58,  20.72,  20.01,  19.47,  19.13,  18.97,  18.92,  18.90,  18.89, 
      26.30,  25.30,  24.30,  23.31,  22.32,  21.36,  20.46,  19.68,  19.08,  18.68,  18.47,  18.39,  18.36,  18.36, 
      26.11,  25.11,  24.11,  23.11,  22.12,  21.15,  20.22,  19.39,  18.72,  18.26,  18.01,  17.90,  17.87,  17.85, 
      25.91,  24.91,  23.91,  22.92,  21.93,  20.95,  20.02,  19.16,  18.44,  17.91,  17.63,  17.52,  17.49,  17.48, 
      27.78,  26.81,  25.90,  25.16,  24.63,  24.27,  24.03,  23.89,  23.81,  23.78,  23.77,  23.76,  23.76,  23.76, 
      27.52,  26.54,  25.61,  24.81,  24.19,  23.74,  23.42,  23.21,  23.10,  23.06,  23.04,  23.03,  23.03,  23.03, 
      27.29,  26.31,  25.36,  24.51,  23.81,  23.27,  22.86,  22.58,  22.43,  22.36,  22.34,  22.33,  22.32,  22.32, 
      27.12,  26.13,  25.16,  24.26,  23.48,  22.86,  22.37,  22.02,  21.81,  21.72,  21.68,  21.67,  21.67,  21.67, 
      26.96,  25.97,  24.99,  24.06,  23.21,  22.50,  21.93,  21.50,  21.23,  21.10,  21.06,  21.04,  21.04,  21.03, 
      26.82,  25.83,  24.84,  23.88,  22.98,  22.19,  21.54,  21.04,  20.70,  20.52,  20.45,  20.43,  20.42,  20.42, 
      26.65,  25.66,  24.66,  23.68,  22.74,  21.89,  21.18,  20.63,  20.23,  20.00,  19.89,  19.86,  19.84,  19.84, 
      26.49,  25.49,  24.49,  23.50,  22.53,  21.63,  20.85,  20.23,  19.78,  19.49,  19.34,  19.28,  19.26,  19.26, 
      26.30,  25.30,  24.30,  23.31,  22.33,  21.39,  20.55,  19.86,  19.36,  19.02,  18.84,  18.76,  18.73,  18.72, 
      26.11,  25.11,  24.11,  23.11,  22.13,  21.17,  20.28,  19.54,  18.97,  18.59,  18.37,  18.28,  18.25,  18.23, 
      25.91,  24.91,  23.91,  22.92,  21.93,  20.97,  20.07,  19.28,  18.65,  18.24,  18.03,  17.96,  17.94,  17.93};

  setup_cool_interp_grid_(&my_rates->LCO, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_OH(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {9, 10.0, 1.0}, // log10(number-density like)
    {6, 1.6, 0.2}, // log10(temperature)
    {14, 1.0, 1.0} // log10(H2I number density)
  };

  double L[] =
     {23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.94,  16.17,  15.85,  15.79,  15.78,  15.78,  15.78, 
      23.28,  22.28,  21.28,  20.28,  19.28,  18.28,  17.29,  16.35,  15.60,  15.28,  15.20,  15.18,  15.18,  15.17, 
      22.82,  21.82,  20.82,  19.82,  18.82,  17.83,  16.84,  15.91,  15.17,  14.82,  14.71,  14.68,  14.67,  14.67, 
      22.54,  21.54,  20.54,  19.54,  18.55,  17.55,  16.58,  15.65,  14.91,  14.51,  14.37,  14.33,  14.32,  14.31, 
      22.27,  21.27,  20.27,  19.27,  18.27,  17.28,  16.31,  15.40,  14.64,  14.19,  14.02,  13.97,  13.96,  13.95, 
      22.05,  21.05,  20.05,  19.06,  18.06,  17.07,  16.10,  15.20,  14.44,  13.97,  13.79,  13.74,  13.72,  13.72, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.93,  16.17,  15.86,  15.80,  15.79,  15.79,  15.78, 
      23.28,  22.28,  21.28,  20.28,  19.28,  18.28,  17.29,  16.35,  15.61,  15.28,  15.20,  15.19,  15.18,  15.18, 
      22.82,  21.82,  20.82,  19.82,  18.82,  17.83,  16.84,  15.91,  15.18,  14.82,  14.71,  14.68,  14.67,  14.67, 
      22.54,  21.54,  20.54,  19.54,  18.55,  17.55,  16.58,  15.65,  14.91,  14.51,  14.37,  14.33,  14.32,  14.31, 
      22.27,  21.27,  20.27,  19.27,  18.27,  17.28,  16.31,  15.40,  14.64,  14.20,  14.02,  13.97,  13.96,  13.95, 
      22.05,  21.05,  20.05,  19.06,  18.06,  17.07,  16.11,  15.20,  14.44,  13.98,  13.79,  13.74,  13.72,  13.72, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.94,  16.18,  15.87,  15.81,  15.80,  15.80,  15.80, 
      23.28,  22.28,  21.28,  20.28,  19.28,  18.28,  17.29,  16.35,  15.61,  15.29,  15.21,  15.20,  15.19,  15.19, 
      22.82,  21.82,  20.82,  19.82,  18.82,  17.83,  16.84,  15.91,  15.18,  14.82,  14.71,  14.68,  14.67,  14.67, 
      22.54,  21.54,  20.54,  19.54,  18.55,  17.55,  16.58,  15.66,  14.91,  14.51,  14.37,  14.33,  14.32,  14.31, 
      22.27,  21.27,  20.27,  19.27,  18.27,  17.28,  16.31,  15.40,  14.64,  14.20,  14.02,  13.97,  13.96,  13.95, 
      22.05,  21.05,  20.05,  19.06,  18.06,  17.07,  16.11,  15.20,  14.44,  13.98,  13.80,  13.74,  13.73,  13.72, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.95,  16.24,  15.98,  15.92,  15.91,  15.91,  15.91, 
      23.28,  22.28,  21.28,  20.28,  19.28,  18.28,  17.30,  16.36,  15.66,  15.36,  15.29,  15.27,  15.26,  15.26, 
      22.82,  21.82,  20.82,  19.82,  18.82,  17.83,  16.85,  15.92,  15.20,  14.86,  14.75,  14.72,  14.71,  14.71, 
      22.54,  21.54,  20.54,  19.54,  18.55,  17.55,  16.58,  15.66,  14.93,  14.54,  14.40,  14.36,  14.35,  14.34, 
      22.27,  21.27,  20.27,  19.27,  18.27,  17.28,  16.31,  15.40,  14.66,  14.21,  14.05,  13.99,  13.98,  13.97, 
      22.05,  21.05,  20.05,  19.06,  18.06,  17.07,  16.11,  15.20,  14.45,  13.99,  13.81,  13.75,  13.74,  13.73, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.89,  17.92,  17.09,  16.64,  16.51,  16.48,  16.46,  16.46,  16.46, 
      23.28,  22.28,  21.28,  20.28,  19.28,  18.29,  17.33,  16.48,  15.93,  15.73,  15.66,  15.64,  15.63,  15.63, 
      22.82,  21.82,  20.82,  19.82,  18.82,  17.84,  16.88,  16.01,  15.39,  15.10,  15.00,  14.97,  14.96,  14.95, 
      22.54,  21.54,  20.54,  19.55,  18.55,  17.56,  16.61,  15.73,  15.07,  14.72,  14.59,  14.55,  14.53,  14.53, 
      22.27,  21.27,  20.27,  19.27,  18.28,  17.29,  16.34,  15.46,  14.75,  14.34,  14.17,  14.12,  14.11,  14.10, 
      22.05,  21.05,  20.06,  19.06,  18.06,  17.08,  16.13,  15.24,  14.52,  14.08,  13.90,  13.85,  13.84,  13.83, 
      23.88,  22.88,  21.88,  20.88,  19.89,  18.93,  18.07,  17.50,  17.25,  17.16,  17.13,  17.12,  17.11,  17.11, 
      23.28,  22.28,  21.28,  20.28,  19.29,  18.32,  17.43,  16.75,  16.40,  16.26,  16.21,  16.20,  16.19,  16.19, 
      22.82,  21.82,  20.82,  19.82,  18.83,  17.86,  16.94,  16.20,  15.75,  15.56,  15.49,  15.46,  15.45,  15.45, 
      22.54,  21.54,  20.54,  19.55,  18.56,  17.58,  16.66,  15.89,  15.39,  15.16,  15.09,  15.06,  15.05,  15.05, 
      22.27,  21.27,  20.27,  19.27,  18.28,  17.31,  16.38,  15.58,  15.03,  14.77,  14.68,  14.66,  14.65,  14.65, 
      22.05,  21.05,  20.06,  19.06,  18.07,  17.10,  16.17,  15.36,  14.77,  14.49,  14.39,  14.36,  14.35,  14.35, 
      23.88,  22.88,  21.88,  20.89,  19.93,  19.07,  18.47,  18.16,  18.02,  17.96,  17.93,  17.92,  17.92,  17.92, 
      23.28,  22.28,  21.28,  20.29,  19.32,  18.43,  17.73,  17.30,  17.09,  17.00,  16.97,  16.95,  16.95,  16.95, 
      22.82,  21.82,  20.82,  19.83,  18.86,  17.94,  17.17,  16.63,  16.35,  16.23,  16.18,  16.17,  16.16,  16.16, 
      22.54,  21.55,  20.55,  19.56,  18.59,  17.66,  16.84,  16.24,  15.90,  15.77,  15.73,  15.72,  15.72,  15.72, 
      22.27,  21.27,  20.27,  19.29,  18.31,  17.38,  16.53,  15.85,  15.47,  15.32,  15.28,  15.27,  15.27,  15.27, 
      22.06,  21.06,  20.06,  19.08,  18.10,  17.17,  16.30,  15.59,  15.17,  15.02,  14.98,  14.97,  14.97,  14.97, 
      23.88,  22.88,  21.89,  20.93,  20.07,  19.45,  19.10,  18.93,  18.85,  18.81,  18.80,  18.79,  18.79,  18.79, 
      23.28,  22.28,  21.29,  20.32,  19.43,  18.73,  18.30,  18.06,  17.94,  17.89,  17.86,  17.85,  17.85,  17.85, 
      22.82,  21.82,  20.83,  19.86,  18.94,  18.17,  17.62,  17.31,  17.15,  17.09,  17.06,  17.05,  17.05,  17.05, 
      22.55,  21.55,  20.56,  19.59,  18.66,  17.84,  17.21,  16.82,  16.64,  16.57,  16.55,  16.55,  16.54,  16.54, 
      22.27,  21.27,  20.29,  19.32,  18.38,  17.53,  16.82,  16.34,  16.13,  16.06,  16.04,  16.04,  16.04,  16.04, 
      22.06,  21.06,  20.08,  19.11,  18.17,  17.30,  16.55,  16.01,  15.76,  15.69,  15.67,  15.66,  15.66,  15.66, 
      23.88,  22.89,  21.92,  21.07,  20.47,  20.14,  19.96,  19.85,  19.79,  19.76,  19.75,  19.74,  19.74,  19.73, 
      23.28,  22.29,  21.32,  20.42,  19.73,  19.31,  19.06,  18.92,  18.85,  18.81,  18.80,  18.79,  18.79,  18.79, 
      22.82,  21.83,  20.86,  19.94,  19.17,  18.63,  18.30,  18.14,  18.06,  18.04,  18.03,  18.02,  18.02,  18.02, 
      22.55,  21.56,  20.58,  19.66,  18.84,  18.21,  17.81,  17.62,  17.55,  17.53,  17.52,  17.52,  17.52,  17.52, 
      22.27,  21.29,  20.32,  19.38,  18.53,  17.82,  17.34,  17.12,  17.05,  17.03,  17.02,  17.02,  17.02,  17.02, 
      22.06,  21.08,  20.11,  19.17,  18.30,  17.55,  17.01,  16.74,  16.66,  16.64,  16.64,  16.63,  16.63,  16.63}; 

  setup_cool_interp_grid_(&my_rates->LOH, rank, params, L, log10(coolunit));
}


extern "C" void initialize_cooling_rate_H2O(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit)
{
  const int rank = 3;
  const struct regular_range_ params[3] = {
    {10, 10.0, 1.0}, // log10(number-density like)
    {12, 1.0, 0.2}, // log10(temperature)
    {16, -1.0, 1.0} // log10(H2I number density)
  };

  double L[] =
     {27.85,  26.85,  25.85,  24.85,  23.85,  22.85,  21.85,  20.85,  19.86,  18.91,  18.21,  17.94,  17.89,  17.88,  17.88,  17.87, 
      27.18,  26.18,  25.18,  24.18,  23.18,  22.18,  21.18,  20.19,  19.20,  18.25,  17.50,  17.18,  17.11,  17.09,  17.08,  17.08, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.58,  18.60,  17.65,  16.88,  16.50,  16.39,  16.36,  16.35,  16.35, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.11,  20.11,  19.12,  18.14,  17.20,  16.42,  15.95,  15.79,  15.74,  15.72,  15.72, 
      25.71,  24.71,  23.71,  22.71,  21.71,  20.71,  19.72,  18.72,  17.75,  16.82,  16.02,  15.47,  15.24,  15.17,  15.15,  15.14, 
      25.35,  24.35,  23.35,  22.35,  21.36,  20.36,  19.36,  18.37,  17.40,  16.48,  15.66,  15.06,  14.75,  14.64,  14.61,  14.60, 
      25.03,  24.03,  23.03,  22.03,  21.03,  20.03,  19.04,  18.05,  17.08,  16.16,  15.33,  14.68,  14.30,  14.16,  14.11,  14.10, 
      24.72,  23.72,  22.72,  21.72,  20.72,  19.73,  18.73,  17.74,  16.77,  15.85,  15.00,  14.32,  13.88,  13.70,  13.64,  13.63, 
      24.42,  23.42,  22.42,  21.42,  20.42,  19.43,  18.43,  17.45,  16.48,  15.54,  14.68,  13.96,  13.47,  13.25,  13.19,  13.17, 
      24.15,  23.15,  22.15,  21.15,  20.15,  19.16,  18.16,  17.17,  16.20,  15.26,  14.39,  13.63,  13.09,  12.84,  12.77,  12.75, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.90,  15.93,  14.98,  14.09,  13.30,  12.72,  12.43,  12.35,  12.33, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.63,  17.64,  16.65,  15.67,  14.72,  13.83,  13.02,  12.42,  12.12,  12.04,  12.02, 
      27.85,  26.85,  25.85,  24.85,  23.85,  22.85,  21.85,  20.85,  19.86,  18.92,  18.22,  17.96,  17.91,  17.90,  17.90,  17.90, 
      27.18,  26.18,  25.18,  24.18,  23.18,  22.18,  21.18,  20.19,  19.20,  18.26,  17.52,  17.20,  17.13,  17.11,  17.10,  17.10, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.58,  18.60,  17.66,  16.89,  16.52,  16.40,  16.37,  16.37,  16.36, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.11,  20.11,  19.12,  18.14,  17.20,  16.42,  15.96,  15.79,  15.75,  15.73,  15.73, 
      25.71,  24.71,  23.71,  22.71,  21.71,  20.71,  19.72,  18.72,  17.75,  16.82,  16.02,  15.48,  15.24,  15.17,  15.15,  15.15, 
      25.35,  24.35,  23.35,  22.35,  21.36,  20.36,  19.36,  18.37,  17.40,  16.48,  15.66,  15.06,  14.75,  14.64,  14.61,  14.60, 
      25.03,  24.03,  23.03,  22.03,  21.03,  20.03,  19.04,  18.05,  17.08,  16.16,  15.33,  14.68,  14.31,  14.16,  14.12,  14.11, 
      24.72,  23.72,  22.72,  21.72,  20.72,  19.73,  18.73,  17.74,  16.77,  15.85,  15.00,  14.32,  13.89,  13.70,  13.65,  13.63, 
      24.42,  23.42,  22.42,  21.42,  20.42,  19.43,  18.43,  17.45,  16.48,  15.54,  14.68,  13.96,  13.47,  13.25,  13.19,  13.17, 
      24.15,  23.15,  22.15,  21.15,  20.15,  19.16,  18.16,  17.17,  16.20,  15.26,  14.39,  13.63,  13.09,  12.84,  12.77,  12.75, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.90,  15.93,  14.98,  14.09,  13.30,  12.72,  12.43,  12.35,  12.33, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.63,  17.64,  16.65,  15.67,  14.72,  13.83,  13.02,  12.42,  12.12,  12.04,  12.02, 
      27.85,  26.85,  25.85,  24.85,  23.85,  22.85,  21.85,  20.85,  19.87,  18.95,  18.34,  18.16,  18.13,  18.12,  18.12,  18.12, 
      27.18,  26.18,  25.18,  24.18,  23.18,  22.18,  21.19,  20.19,  19.21,  18.28,  17.61,  17.36,  17.30,  17.29,  17.28,  17.28, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.59,  18.60,  17.68,  16.97,  16.63,  16.53,  16.50,  16.49,  16.49, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.11,  20.11,  19.12,  18.15,  17.22,  16.47,  16.03,  15.87,  15.82,  15.80,  15.80, 
      25.71,  24.71,  23.71,  22.71,  21.71,  20.71,  19.72,  18.73,  17.76,  16.84,  16.04,  15.51,  15.28,  15.21,  15.19,  15.18, 
      25.35,  24.35,  23.35,  22.35,  21.36,  20.36,  19.36,  18.37,  17.41,  16.49,  15.67,  15.07,  14.77,  14.66,  14.63,  14.62, 
      25.03,  24.03,  23.03,  22.03,  21.03,  20.03,  19.04,  18.05,  17.08,  16.16,  15.33,  14.68,  14.31,  14.16,  14.12,  14.11, 
      24.72,  23.72,  22.72,  21.72,  20.72,  19.73,  18.73,  17.74,  16.78,  15.85,  15.00,  14.32,  13.88,  13.70,  13.65,  13.63, 
      24.42,  23.42,  22.42,  21.42,  20.42,  19.43,  18.43,  17.45,  16.48,  15.54,  14.68,  13.96,  13.47,  13.25,  13.19,  13.17, 
      24.15,  23.15,  22.15,  21.15,  20.15,  19.16,  18.16,  17.17,  16.20,  15.26,  14.39,  13.63,  13.09,  12.84,  12.77,  12.75, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.90,  15.93,  14.98,  14.09,  13.30,  12.72,  12.43,  12.35,  12.33, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.63,  17.64,  16.65,  15.67,  14.73,  13.83,  13.02,  12.42,  12.13,  12.04,  12.02, 
      27.85,  26.85,  25.85,  24.85,  23.85,  22.85,  21.85,  20.86,  19.90,  19.15,  18.85,  18.80,  18.79,  18.79,  18.78,  18.78, 
      27.18,  26.18,  25.18,  24.18,  23.18,  22.18,  21.19,  20.20,  19.24,  18.44,  18.03,  17.90,  17.86,  17.85,  17.85,  17.85, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.58,  19.59,  18.64,  17.81,  17.28,  17.05,  16.97,  16.95,  16.94,  16.93, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.11,  20.12,  19.13,  18.18,  17.33,  16.67,  16.31,  16.17,  16.13,  16.11,  16.11, 
      25.71,  24.71,  23.71,  22.71,  21.71,  20.71,  19.72,  18.74,  17.79,  16.91,  16.17,  15.69,  15.47,  15.40,  15.38,  15.38, 
      25.35,  24.35,  23.35,  22.35,  21.36,  20.36,  19.37,  18.39,  17.44,  16.54,  15.76,  15.17,  14.87,  14.77,  14.74,  14.73, 
      25.03,  24.03,  23.03,  22.03,  21.03,  20.04,  19.04,  18.06,  17.10,  16.19,  15.37,  14.72,  14.33,  14.17,  14.13,  14.12, 
      24.72,  23.72,  22.72,  21.72,  20.72,  19.73,  18.73,  17.75,  16.79,  15.87,  15.03,  14.33,  13.88,  13.69,  13.63,  13.62, 
      24.42,  23.42,  22.42,  21.42,  20.43,  19.43,  18.44,  17.45,  16.49,  15.56,  14.70,  13.97,  13.48,  13.26,  13.19,  13.18, 
      24.15,  23.15,  22.15,  21.15,  20.15,  19.16,  18.16,  17.18,  16.21,  15.27,  14.40,  13.64,  13.09,  12.84,  12.77,  12.75, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.88,  17.89,  16.90,  15.93,  14.99,  14.10,  13.31,  12.72,  12.43,  12.35,  12.33, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.63,  17.64,  16.65,  15.68,  14.73,  13.83,  13.02,  12.42,  12.12,  12.04,  12.02, 
      27.85,  26.85,  25.85,  24.85,  23.85,  22.85,  21.86,  20.90,  20.13,  19.80,  19.73,  19.71,  19.70,  19.70,  19.70,  19.70, 
      27.18,  26.18,  25.18,  24.18,  23.18,  22.19,  21.20,  20.24,  19.42,  18.95,  18.76,  18.69,  18.66,  18.65,  18.65,  18.64, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.58,  20.59,  19.64,  18.79,  18.20,  17.88,  17.73,  17.67,  17.65,  17.64,  17.63, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.12,  20.13,  19.17,  18.30,  17.60,  17.13,  16.88,  16.77,  16.73,  16.72,  16.71, 
      25.71,  24.71,  23.71,  22.71,  21.71,  20.72,  19.73,  18.78,  17.88,  17.11,  16.51,  16.14,  15.97,  15.90,  15.88,  15.88, 
      25.35,  24.35,  23.35,  22.36,  21.36,  20.36,  19.38,  18.42,  17.51,  16.69,  16.00,  15.51,  15.25,  15.15,  15.12,  15.11, 
      25.03,  24.03,  23.03,  22.03,  21.03,  20.04,  19.05,  18.08,  17.15,  16.28,  15.52,  14.92,  14.56,  14.42,  14.38,  14.37, 
      24.72,  23.72,  22.72,  21.72,  20.73,  19.73,  18.74,  17.77,  16.83,  15.94,  15.14,  14.47,  14.04,  13.86,  13.80,  13.79, 
      24.42,  23.42,  22.42,  21.43,  20.43,  19.43,  18.45,  17.47,  16.52,  15.61,  14.78,  14.06,  13.55,  13.33,  13.27,  13.26, 
      24.15,  23.15,  22.15,  21.15,  20.16,  19.16,  18.17,  17.19,  16.23,  15.31,  14.45,  13.69,  13.14,  12.89,  12.82,  12.80, 
      23.88,  22.88,  21.88,  20.88,  19.88,  18.89,  17.90,  16.91,  15.95,  15.01,  14.13,  13.34,  12.74,  12.45,  12.36,  12.35, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.64,  17.65,  16.68,  15.73,  14.82,  13.96,  13.19,  12.56,  12.19,  12.07,  12.04, 
      27.85,  26.85,  25.85,  24.85,  23.85,  22.86,  21.90,  21.13,  20.79,  20.70,  20.66,  20.64,  20.63,  20.63,  20.63,  20.62, 
      27.18,  26.18,  25.18,  24.18,  23.19,  22.19,  21.24,  20.43,  19.94,  19.72,  19.61,  19.56,  19.53,  19.51,  19.51,  19.50, 
      26.58,  25.58,  24.58,  23.58,  22.58,  21.59,  20.63,  19.79,  19.19,  18.84,  18.64,  18.53,  18.48,  18.46,  18.45,  18.44, 
      26.11,  25.11,  24.11,  23.11,  22.11,  21.13,  20.17,  19.29,  18.59,  18.11,  17.79,  17.62,  17.53,  17.50,  17.48,  17.48, 
      25.71,  24.71,  23.71,  22.71,  21.72,  20.73,  19.77,  18.87,  18.10,  17.49,  17.05,  16.78,  16.65,  16.60,  16.58,  16.57, 
      25.35,  24.35,  23.36,  22.36,  21.36,  20.38,  19.41,  18.50,  17.67,  16.98,  16.42,  16.04,  15.84,  15.76,  15.73,  15.72, 
      25.03,  24.03,  23.03,  22.03,  21.04,  20.05,  19.07,  18.13,  17.26,  16.47,  15.81,  15.29,  14.99,  14.87,  14.83,  14.82, 
      24.72,  23.72,  22.72,  21.73,  20.73,  19.74,  18.77,  17.82,  16.92,  16.09,  15.36,  14.76,  14.37,  14.21,  14.17,  14.16, 
      24.42,  23.42,  22.42,  21.43,  20.43,  19.44,  18.46,  17.51,  16.59,  15.73,  14.95,  14.28,  13.82,  13.61,  13.55,  13.54, 
      24.15,  23.15,  22.15,  21.15,  20.16,  19.17,  18.18,  17.22,  16.28,  15.40,  14.58,  13.87,  13.35,  13.10,  13.03,  13.02, 
      23.88,  22.88,  21.88,  20.88,  19.89,  18.89,  17.91,  16.93,  15.98,  15.08,  14.23,  13.47,  12.89,  12.60,  12.51,  12.50, 
      23.63,  22.63,  21.63,  20.63,  19.63,  18.64,  17.65,  16.67,  15.72,  14.80,  13.93,  13.14,  12.54,  12.25,  12.17,  12.15, 
      27.85,  26.85,  25.85,  24.85,  23.86,  22.90,  22.13,  21.78,  21.68,  21.63,  21.59,  21.57,  21.55,  21.54,  21.53,  21.53, 
      27.18,  26.18,  25.18,  24.19,  23.19,  22.23,  21.42,  20.94,  20.72,  20.59,  20.51,  20.46,  20.43,  20.41,  20.40,  20.40, 
      26.58,  25.58,  24.58,  23.58,  22.59,  21.63,  20.78,  20.19,  19.83,  19.61,  19.47,  19.39,  19.35,  19.33,  19.32,  19.31, 
      26.11,  25.11,  24.11,  23.11,  22.12,  21.16,  20.29,  19.60,  19.12,  18.78,  18.55,  18.42,  18.36,  18.33,  18.32,  18.31, 
      25.71,  24.71,  23.71,  22.72,  21.73,  20.76,  19.87,  19.10,  18.50,  18.05,  17.72,  17.51,  17.42,  17.38,  17.36,  17.36, 
      25.35,  24.36,  23.36,  22.36,  21.37,  20.41,  19.49,  18.67,  17.98,  17.42,  16.98,  16.68,  16.53,  16.48,  16.46,  16.45, 
      25.03,  24.03,  23.03,  22.04,  21.05,  20.07,  19.13,  18.25,  17.47,  16.80,  16.23,  15.80,  15.57,  15.48,  15.45,  15.45, 
      24.72,  23.72,  22.72,  21.73,  20.74,  19.76,  18.81,  17.91,  17.09,  16.35,  15.71,  15.19,  14.88,  14.76,  14.72,  14.71, 
      24.42,  23.42,  22.43,  21.43,  20.44,  19.46,  18.50,  17.58,  16.73,  15.94,  15.23,  14.64,  14.25,  14.09,  14.04,  14.03, 
      24.15,  23.15,  22.15,  21.16,  20.17,  19.18,  18.21,  17.28,  16.39,  15.56,  14.81,  14.16,  13.71,  13.52,  13.46,  13.45, 
      23.88,  22.88,  21.88,  20.89,  19.89,  18.90,  17.93,  16.98,  16.06,  15.20,  14.41,  13.71,  13.20,  12.96,  12.89,  12.87, 
      23.63,  22.63,  21.63,  20.63,  19.64,  18.65,  17.67,  16.71,  15.78,  14.90,  14.07,  13.34,  12.82,  12.60,  12.54,  12.53, 
      27.85,  26.85,  25.85,  24.86,  23.90,  23.13,  22.78,  22.66,  22.60,  22.56,  22.53,  22.51,  22.49,  22.48,  22.47,  22.47, 
      27.18,  26.18,  25.19,  24.19,  23.23,  22.42,  21.94,  21.71,  21.57,  21.48,  21.41,  21.37,  21.34,  21.33,  21.32,  21.31, 
      26.58,  25.58,  24.58,  23.59,  22.63,  21.78,  21.19,  20.84,  20.60,  20.43,  20.33,  20.27,  20.24,  20.22,  20.21,  20.21, 
      26.11,  25.11,  24.11,  23.12,  22.16,  21.29,  20.60,  20.12,  19.77,  19.51,  19.35,  19.26,  19.22,  19.20,  19.19,  19.19, 
      25.71,  24.71,  23.71,  22.73,  21.76,  20.86,  20.10,  19.52,  19.06,  18.70,  18.45,  18.31,  18.24,  18.22,  18.21,  18.20, 
      25.35,  24.36,  23.36,  22.37,  21.40,  20.48,  19.67,  18.99,  18.45,  17.99,  17.64,  17.42,  17.31,  17.27,  17.26,  17.25, 
      25.03,  24.03,  23.04,  22.04,  21.06,  20.12,  19.24,  18.47,  17.81,  17.24,  16.75,  16.42,  16.25,  16.19,  16.17,  16.16, 
      24.72,  23.72,  22.73,  21.74,  20.76,  19.80,  18.90,  18.08,  17.35,  16.70,  16.14,  15.71,  15.49,  15.41,  15.38,  15.38, 
      24.42,  23.43,  22.43,  21.44,  20.46,  19.50,  18.58,  17.72,  16.94,  16.22,  15.58,  15.07,  14.78,  14.67,  14.65,  14.64, 
      24.15,  23.15,  22.16,  21.16,  20.18,  19.21,  18.27,  17.38,  16.55,  15.79,  15.10,  14.54,  14.21,  14.09,  14.06,  14.05, 
      23.88,  22.88,  21.89,  20.89,  19.90,  18.93,  17.97,  17.06,  16.19,  15.39,  14.65,  14.04,  13.65,  13.51,  13.47,  13.46, 
      23.63,  22.63,  21.63,  20.64,  19.65,  18.67,  17.71,  16.78,  15.89,  15.05,  14.28,  13.64,  13.28,  13.17,  13.14,  13.14, 
      27.85,  26.85,  25.86,  24.89,  24.13,  23.79,  23.69,  23.64,  23.59,  23.56,  23.53,  23.51,  23.49,  23.48,  23.46,  23.46, 
      27.18,  26.18,  25.19,  24.23,  23.42,  22.95,  22.73,  22.59,  22.49,  22.41,  22.36,  22.32,  22.29,  22.28,  22.27,  22.26, 
      26.58,  25.58,  24.59,  23.62,  22.78,  22.20,  21.85,  21.61,  21.43,  21.31,  21.23,  21.18,  21.16,  21.14,  21.14,  21.13, 
      26.11,  25.11,  24.12,  23.16,  22.28,  21.61,  21.15,  20.80,  20.53,  20.34,  20.21,  20.15,  20.12,  20.10,  20.09,  20.09, 
      25.71,  24.71,  23.72,  22.76,  21.86,  21.10,  20.53,  20.09,  19.72,  19.44,  19.25,  19.15,  19.11,  19.09,  19.08,  19.08, 
      25.36,  24.36,  23.37,  22.40,  21.48,  20.67,  20.00,  19.47,  19.01,  18.63,  18.36,  18.21,  18.14,  18.11,  18.11,  18.10, 
      25.03,  24.03,  23.04,  22.06,  21.11,  20.24,  19.47,  18.83,  18.26,  17.76,  17.36,  17.11,  17.00,  16.96,  16.95,  16.95, 
      24.72,  23.73,  22.73,  21.75,  20.80,  19.89,  19.08,  18.36,  17.72,  17.14,  16.65,  16.32,  16.18,  16.13,  16.12,  16.11, 
      24.43,  23.43,  22.44,  21.45,  20.49,  19.57,  18.71,  17.93,  17.21,  16.56,  15.99,  15.59,  15.41,  15.35,  15.33,  15.33, 
      24.15,  23.16,  22.16,  21.18,  20.21,  19.26,  18.37,  17.54,  16.77,  16.06,  15.45,  15.01,  14.82,  14.76,  14.75,  14.74, 
      23.88,  22.89,  21.89,  20.90,  19.93,  18.97,  18.05,  17.18,  16.36,  15.60,  14.93,  14.45,  14.23,  14.18,  14.16,  14.16, 
      23.63,  22.63,  21.64,  20.65,  19.67,  18.71,  17.77,  16.88,  16.02,  15.23,  14.54,  14.08,  13.90,  13.86,  13.85,  13.85, 
      27.85,  26.86,  25.90,  25.13,  24.78,  24.67,  24.61,  24.56,  24.52,  24.49,  24.47,  24.44,  24.43,  24.41,  24.40,  24.39, 
      27.18,  26.19,  25.23,  24.42,  23.95,  23.73,  23.58,  23.47,  23.38,  23.31,  23.27,  23.24,  23.22,  23.20,  23.20,  23.19, 
      26.58,  25.59,  24.62,  23.78,  23.20,  22.86,  22.62,  22.43,  22.29,  22.19,  22.14,  22.10,  22.08,  22.07,  22.07,  22.06, 
      26.11,  25.12,  24.15,  23.28,  22.61,  22.16,  21.81,  21.53,  21.32,  21.18,  21.09,  21.05,  21.03,  21.02,  21.01,  21.01, 
      25.71,  24.72,  23.75,  22.85,  22.11,  21.55,  21.11,  20.75,  20.44,  20.22,  20.09,  20.03,  20.00,  19.99,  19.98,  19.98, 
      25.36,  24.37,  23.40,  22.48,  21.67,  21.01,  20.49,  20.04,  19.64,  19.34,  19.14,  19.05,  19.01,  19.00,  18.99,  18.99, 
      25.03,  24.04,  23.06,  22.11,  21.24,  20.48,  19.84,  19.28,  18.78,  18.34,  18.03,  17.87,  17.81,  17.79,  17.78,  17.78, 
      24.73,  23.73,  22.75,  21.79,  20.89,  20.08,  19.37,  18.74,  18.16,  17.64,  17.25,  17.02,  16.94,  16.91,  16.91,  16.90, 
      24.43,  23.43,  22.45,  21.48,  20.56,  19.71,  18.94,  18.24,  17.60,  17.01,  16.53,  16.24,  16.12,  16.09,  16.08,  16.08, 
      24.16,  23.16,  22.17,  21.20,  20.26,  19.36,  18.53,  17.77,  17.07,  16.43,  15.92,  15.64,  15.54,  15.52,  15.51,  15.51, 
      23.89,  22.89,  21.90,  20.93,  19.97,  19.04,  18.17,  17.35,  16.58,  15.89,  15.34,  15.05,  14.97,  14.95,  14.94,  14.94, 
      23.63,  22.64,  21.65,  20.67,  19.71,  18.77,  17.86,  17.01,  16.20,  15.48,  14.94,  14.70,  14.64,  14.62,  14.62,  14.62}; 

  setup_cool_interp_grid_(&my_rates->LH2O, rank, params, L, log10(coolunit));
}


extern "C" void initialize_primordial_opacity(chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{
  const int rank = 2;
  const struct regular_range_ params[2] = {
    {15, -16.0, 1.0}, // log10(mass-density)
    {29,   1.8, 0.1}  // log10(temperature)
  };

  double kp[29][15] =
   {{-14.39, -13.39, -12.39, -11.39, -10.39, -9.39, -8.39, -7.39,  -6.39,  -5.39,  -4.39,  -3.39,  -2.39,  -1.39,  -0.39}
   ,{-14.18, -13.18, -12.18, -11.18, -10.18, -9.18, -8.18, -7.18,  -6.18,  -5.18,  -4.18,  -3.18,  -2.18,  -1.18,  -0.18}
   ,{-14.03, -13.03, -12.03, -11.03, -10.03, -9.03, -8.03, -7.03,  -6.03,  -5.03,  -4.03,  -3.03,  -2.03,  -1.03,  -0.03}
   ,{-13.92, -12.92, -11.92, -10.92,  -9.92, -8.92, -7.92, -6.92,  -5.92,  -4.92,  -3.92,  -2.92,  -1.92,  -0.92,   0.08}
   ,{-13.84, -12.84, -11.84, -10.84,  -9.84, -8.84, -7.84, -6.84,  -5.84,  -4.84,  -3.84,  -2.84,  -1.84,  -0.84,   0.16}
   ,{-13.79, -12.79, -11.79, -10.79,  -9.79, -8.79, -7.79, -6.79,  -5.79,  -4.79,  -3.79,  -2.79,  -1.79,  -0.79,   0.21}
   ,{-13.78, -12.78, -11.78, -10.78,  -9.78, -8.78, -7.78, -6.78,  -5.78,  -4.78,  -3.78,  -2.78,  -1.78,  -0.78,   0.22}
   ,{-13.80, -12.80, -11.80, -10.80,  -9.80, -8.80, -7.80, -6.80,  -5.80,  -4.80,  -3.80,  -2.80,  -1.80,  -0.80,   0.20}
   ,{-13.84, -12.84, -11.84, -10.84,  -9.84, -8.84, -7.84, -6.84,  -5.84,  -4.84,  -3.84,  -2.84,  -1.84,  -0.84,   0.16}
   ,{-13.87, -12.87, -11.87, -10.87,  -9.87, -8.87, -7.87, -6.87,  -5.87,  -4.87,  -3.87,  -2.87,  -1.87,  -0.87,   0.13}
   ,{-13.88, -12.88, -11.88, -10.88,  -9.89, -8.89, -7.89, -6.89,  -5.89,  -4.89,  -3.88,  -2.88,  -1.88,  -0.88,   0.12}
   ,{-13.85, -12.85, -11.85, -10.85,  -9.85, -8.85, -7.85, -6.85,  -5.85,  -4.85,  -3.85,  -2.85,  -1.85,  -0.85,   0.15}
   ,{-13.75, -12.75, -11.75, -10.75,  -9.75, -8.75, -7.75, -6.75,  -5.75,  -4.75,  -3.75,  -2.75,  -1.75,  -0.75,   0.25}
   ,{-13.59, -12.61, -11.61, -10.60,  -9.60, -8.60, -7.60, -6.60,  -5.60,  -4.60,  -3.60,  -2.60,  -1.60,  -0.60,   0.40}
   ,{-13.80, -12.76, -11.64, -10.54,  -9.51, -8.51, -7.50, -6.50,  -5.50,  -4.50,  -3.50,  -2.50,  -1.50,  -0.50,   0.50}
   ,{-14.06, -13.06, -12.06, -11.05,  -9.97, -8.72, -7.53, -6.46,  -5.44,  -4.43,  -3.43,  -2.43,  -1.43,  -0.43,   0.57}
   ,{-12.33, -11.82, -11.30, -10.73, -10.05, -9.19, -8.22, -7.18,  -5.82,  -4.56,  -3.47,  -2.44,  -1.43,  -0.43,   0.57}
   ,{ -9.70,  -9.22,  -8.73,  -8.23,  -7.73, -7.23, -6.72, -6.19,  -5.51,  -4.60,  -3.69,  -2.57,  -1.50,  -0.47,   0.54}
   ,{ -6.65,  -6.57,  -6.39,  -6.08,  -5.66, -5.19, -4.70, -4.21,  -3.71,  -3.21,  -2.72,  -2.23,  -1.52,  -0.53,   0.49}
   ,{ -3.71,  -3.52,  -3.46,  -3.42,  -3.36, -3.23, -2.96, -2.58,  -2.12,  -1.63,  -1.14,  -0.69,  -0.33,   0.12,   0.71}
   ,{ -3.49,  -2.52,  -1.66,  -1.18,  -1.00, -0.94, -0.89, -0.81,  -0.63,  -0.31,   0.11,   0.57,   1.01,   1.37,   1.76}
   ,{ -3.93,  -2.94,  -1.94,  -0.95,  -0.02,  0.62,  0.89,  0.99,   1.04,   1.12,   1.30,   1.61,   2.02,   2.40,   2.74}
   ,{ -4.29,  -3.30,  -2.34,  -1.38,  -0.38,  0.60,  1.52,  2.13,   2.38,   2.47,   2.54,   2.66,   2.89,   3.23,   3.54}
   ,{ -4.74,  -3.74,  -2.74,  -1.74,  -0.75,  0.20,  1.17,  2.14,   2.95,   3.39,   3.55,  10.00,  10.00,  10.00,  10.00}
   ,{ -5.18,  -4.18,  -3.18,  -2.18,  -1.18, -0.18,  0.81,  1.77,   2.72,   3.60,   4.14,  10.00,  10.00,  10.00,  10.00}
   ,{ -5.40,  -4.54,  -3.60,  -2.61,  -1.61, -0.61,  0.39,  1.39,   2.37,   3.31,   4.18,  10.00,  10.00,  10.00,  10.00}
   ,{ -5.75,  -4.75,  -3.76,  -2.82,  -1.95, -1.00,  0.00,  1.00,   2.00,   2.99,   3.93,  10.00,  10.00,  10.00,  10.00}
   ,{ -6.13,  -5.13,  -4.13,  -3.13,  -2.13, -1.15, -0.25,  0.68,   1.67,   2.67,   3.66,  10.00,  10.00,  10.00,  10.00}
   ,{ -6.45,  -5.45,  -4.45,  -3.45,  -2.45, -1.45, -0.45,  0.53,   1.46,   2.42,   3.41,  10.00,  10.00,  10.00,  10.00}};

  setup_generic_grid_props_(&my_rates->alphap.props, rank, params);

  my_rates->alphap.data = (double*)malloc(my_rates->alphap.props.data_size * sizeof(double));
  for(int iD=0; iD<params[0].count; iD++) {
    double log_rho = params[0].start + iD*params[0].step;
    for(int iT=0; iT<params[1].count; iT++) {
      int itab = iD * params[1].count + iT;
      my_rates->alphap.data[itab] = kp[iT][iD] + log_rho;
    }
  }

}


static int allocate_rates_metal(chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{
    my_rates->k125 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k129 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k130 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k131 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k132 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k133 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k134 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k135 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k136 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k137 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k148 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k149 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k150 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k151 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k152 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->k153 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->kz15 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz16 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz17 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz18 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz19 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz20 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz21 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz22 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz23 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz24 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz25 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz26 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz27 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz28 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz29 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz30 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz31 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz32 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz33 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz34 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz35 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz36 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz37 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz38 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz39 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz40 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz41 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz42 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz43 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz44 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz45 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz46 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz47 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz48 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz49 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz50 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz51 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz52 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz53 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));
    my_rates->kz54 = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    my_rates->cieY06  = (double*)malloc(my_chemistry->NumberOfTemperatureBins * sizeof(double));

    return SUCCESS;
}
