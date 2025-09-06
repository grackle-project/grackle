//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines the function to initialize extra collisional rates
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

#ifndef INIT_EXTRA_COLLISIONAL_RATES_HPP
#define INIT_EXTRA_COLLISIONAL_RATES_HPP

#include "grackle.h"
#include "grackle_macros.h"    // tiny
#include "internal_types.hpp"  // CollisionalRxnRateCollection
#include "internal_units.h"    // InternalGrUnits
#include "LUT.hpp"             // CollisionalRxnLUT
#include "opaque_storage.hpp"  // gr_opaque_storage

namespace grackle::impl {

/// initialize misc primordial_chemistry == 4 and metal chemistry rates
///
/// @todo
/// we should refactor this logic so it is implemented that all "standard"
/// collisional rate initialization logic is implemented in a consistent way
int init_extra_collisional_rates(
  chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
  code_units *my_units)
{
  if (my_chemistry->primordial_chemistry == 0) {
    return GR_SUCCESS;
  }

  // temporarily construct the InternalGrUnits struct
  // -> the construction logic deduplicates a lot of logic that was
  //    previously copied and pasted across a lot of fortran files
  InternalGrUnits internalu = new_internalu_legacy_C_(my_units);
  const double kunit = internalu_calc_kunit_(internalu);

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

  // when we allocated all of the rate buffers, we also initialized all rate
  // constants to have values of tiny
  grackle::impl::CollisionalRxnRateCollection* kcol_rate_tables =
    my_rates->opaque_storage->kcol_rate_tables;

  // Fill in tables of
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

    kcol_rate_tables->data[CollisionalRxnLUT::k125][i] = 6.4e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::k129][i] = 3.9e-19 * pow(ttt300, 1.8) * exp(20.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::k130][i] = 3.9e-19 * pow(ttt300, 1.8) * exp(20.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::k131][i] = 3.4e-9 * pow(ttt300, -0.4);
    kcol_rate_tables->data[CollisionalRxnLUT::k132][i] = 3.0e-16 * pow(ttt300, 0.95) * exp(-ttt/9320.0);
    kcol_rate_tables->data[CollisionalRxnLUT::k133][i] = 5.7e-8 * pow(ttt300, -0.50);
    kcol_rate_tables->data[CollisionalRxnLUT::k134][i] = 4.6e-8 * pow(ttt300, -0.50);
    kcol_rate_tables->data[CollisionalRxnLUT::k135][i] = 4.6e-8 * pow(ttt300, -0.50);
    kcol_rate_tables->data[CollisionalRxnLUT::k136][i] = 6.4e-9 * pow(ttt300, 0.41);
    kcol_rate_tables->data[CollisionalRxnLUT::k137][i] = 1.5e-9 * pow(ttt300, -0.1);
  
    kcol_rate_tables->data[CollisionalRxnLUT::k148][i] = 5.0e-21;
    if(ttt < 1000.0) {
      kcol_rate_tables->data[CollisionalRxnLUT::k149][i] = 7.60e-18 * pow(ttt, -0.50);
    } else {
      kcol_rate_tables->data[CollisionalRxnLUT::k149][i] = 3.45e-16 * pow(ttt, -1.06);
    }
    kcol_rate_tables->data[CollisionalRxnLUT::k150][i] = 3.0e-10 * exp(-6717.0/ttt);
    if(ttt < 4000.0) {
      kcol_rate_tables->data[CollisionalRxnLUT::k151][i] = 1.6e-14 * pow(ttt, -0.33);
    } else {
      kcol_rate_tables->data[CollisionalRxnLUT::k151][i] = 1.0e-15;
    }
    kcol_rate_tables->data[CollisionalRxnLUT::k152][i] = 9.1e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::k153][i] = 1.7e-7 * pow(ttt, -0.5);
  
    kcol_rate_tables->data[CollisionalRxnLUT::kz15][i] = 4.98e-11;
    kcol_rate_tables->data[CollisionalRxnLUT::kz16][i] = 2.70e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz17][i] = 7.00e-14 * pow(ttt300, 2.80) * exp(-1950.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz18][i] = 6.83e-12 * pow(ttt300, 1.60) * exp(-9720.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz19][i] = 3.30e-10 * exp(-8460.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz20][i] = 6.64e-10 * exp(-11700.0/ttt);
    if(ttt < 1.0e7) {
      kcol_rate_tables->data[CollisionalRxnLUT::kz21][i] = 3.43e-13 * pow(ttt300, 2.67) * exp(-3160.0/ttt);
    } else {
      kcol_rate_tables->data[CollisionalRxnLUT::kz21][i] = 3.43e-13 * pow(1.0e7/300.0, 2.67) * exp(-3160.0/1.0e7);
    }
    // The rate comes from an experiment (297-3532 K).
    // We refrain to extrapolate it to high temperatures.
    kcol_rate_tables->data[CollisionalRxnLUT::kz22][i] = 7.00e-10 * exp(-232.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz23][i] = 2.38e-10 * exp(-1760.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz24][i] = 1.55e-12 * pow(ttt300, 1.60) * exp(-1660.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz25][i] = 1.65e-12 * pow(ttt300, 1.14) * exp(-50.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz26][i] = 1.0e-13;
    kcol_rate_tables->data[CollisionalRxnLUT::kz27][i] = 1.0e-17;
    kcol_rate_tables->data[CollisionalRxnLUT::kz28][i] = 1.1e-10 * pow(ttt300, 0.5);
    kcol_rate_tables->data[CollisionalRxnLUT::kz29][i] = 3.3e-11;
    kcol_rate_tables->data[CollisionalRxnLUT::kz30][i] = 9.9e-19 * pow(ttt300, -0.38);
    kcol_rate_tables->data[CollisionalRxnLUT::kz31][i] = 4.9e-20 * pow(ttt300, 1.58);
    kcol_rate_tables->data[CollisionalRxnLUT::kz32][i] = 6.6e-11;
    kcol_rate_tables->data[CollisionalRxnLUT::kz33][i] = 4.34e-11 * pow(ttt300, -0.5) * exp(-30.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz34][i] = 2.1e-9;
    kcol_rate_tables->data[CollisionalRxnLUT::kz35][i] = 6.9e-9;
    kcol_rate_tables->data[CollisionalRxnLUT::kz36][i] = 2.0e-9;
    kcol_rate_tables->data[CollisionalRxnLUT::kz37][i] = 7.7e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz38][i] = 6.2e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz39][i] = 6.8e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz40][i] = 1.7e-9;
    kcol_rate_tables->data[CollisionalRxnLUT::kz41][i] = 1.01e-9;
    kcol_rate_tables->data[CollisionalRxnLUT::kz42][i] = 8.3e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz43][i] = 7.5e-10;
    kcol_rate_tables->data[CollisionalRxnLUT::kz44][i] = 4.4e-12 * pow(ttt300, -0.61);
    kcol_rate_tables->data[CollisionalRxnLUT::kz45][i] = 3.4e-12 * pow(ttt300, -0.63);
    kcol_rate_tables->data[CollisionalRxnLUT::kz46][i] = 1.6e-7 * pow(ttt300, -0.5);
    kcol_rate_tables->data[CollisionalRxnLUT::kz47][i] = 2.0e-7 * pow(ttt300, -0.5);
    kcol_rate_tables->data[CollisionalRxnLUT::kz48][i] = 3.5e-7 * pow(ttt300, -0.5);
    kcol_rate_tables->data[CollisionalRxnLUT::kz49][i] = 6.5e-7 * pow(ttt300, -0.5);
    kcol_rate_tables->data[CollisionalRxnLUT::kz50][i] = 1.95e-7 * pow(ttt300, -0.7);
    kcol_rate_tables->data[CollisionalRxnLUT::kz51][i] = 1.0e-17;
    kcol_rate_tables->data[CollisionalRxnLUT::kz52][i]  = 3.00e-11;
    kcol_rate_tables->data[CollisionalRxnLUT::kz53][i]  = 1.30e-11 * exp(-111.0/ttt);
    kcol_rate_tables->data[CollisionalRxnLUT::kz54][i]  = 2.00e-13;
  }

  // at this point, all values in kcol_rate_tables initialized by this
  // function have cgs units. We need to convert to code units & enforce a
  // minimum value
  int tmp_rateid_list[] = {
    // primordial chemistry >= 4 rates
    CollisionalRxnLUT::k125, CollisionalRxnLUT::k129, CollisionalRxnLUT::k130,
    CollisionalRxnLUT::k131, CollisionalRxnLUT::k132, CollisionalRxnLUT::k133,
    CollisionalRxnLUT::k134, CollisionalRxnLUT::k135, CollisionalRxnLUT::k136,
    CollisionalRxnLUT::k137, CollisionalRxnLUT::k148, CollisionalRxnLUT::k149,
    CollisionalRxnLUT::k150, CollisionalRxnLUT::k151, CollisionalRxnLUT::k152,
    CollisionalRxnLUT::k153,
    // metal chemistry rates
    CollisionalRxnLUT::kz15, CollisionalRxnLUT::kz16, CollisionalRxnLUT::kz17,
    CollisionalRxnLUT::kz18, CollisionalRxnLUT::kz19, CollisionalRxnLUT::kz20,
    CollisionalRxnLUT::kz21, CollisionalRxnLUT::kz22, CollisionalRxnLUT::kz23,
    CollisionalRxnLUT::kz24, CollisionalRxnLUT::kz25, CollisionalRxnLUT::kz26,
    CollisionalRxnLUT::kz27, CollisionalRxnLUT::kz28, CollisionalRxnLUT::kz29,
    CollisionalRxnLUT::kz30, CollisionalRxnLUT::kz31, CollisionalRxnLUT::kz32,
    CollisionalRxnLUT::kz33, CollisionalRxnLUT::kz34, CollisionalRxnLUT::kz35,
    CollisionalRxnLUT::kz36, CollisionalRxnLUT::kz37, CollisionalRxnLUT::kz38,
    CollisionalRxnLUT::kz39, CollisionalRxnLUT::kz40, CollisionalRxnLUT::kz41,
    CollisionalRxnLUT::kz42, CollisionalRxnLUT::kz43, CollisionalRxnLUT::kz44,
    CollisionalRxnLUT::kz45, CollisionalRxnLUT::kz46, CollisionalRxnLUT::kz47,
    CollisionalRxnLUT::kz48, CollisionalRxnLUT::kz49, CollisionalRxnLUT::kz50,
    CollisionalRxnLUT::kz51, CollisionalRxnLUT::kz52, CollisionalRxnLUT::kz53,
    CollisionalRxnLUT::kz54
  };

  size_t tmp_rateid_list_len = sizeof(tmp_rateid_list) / sizeof(int);

  for (size_t list_idx = 0; list_idx < tmp_rateid_list_len; list_idx++) {
    // get the lut_idx (lookup table index)
    int lut_idx = tmp_rateid_list[list_idx];
    // load the data corresponding to the specified rate
    double* ptr = kcol_rate_tables->data[lut_idx];

    // now do the "heavy lifting" of setting up appropriate units (and
    // enforcing a minimum value)
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++) {
      // we use std::fmax, rather than fmax to avoid overloading issues with
      // the custom fmax introduced in utils-cpp.hpp that operates on 3
      // arguments (that function is used in code transcribed from Fortran)
      ptr[i] = std::fmax(ptr[i], tiny) / kunit;
    }
  }

  return GR_SUCCESS;
}

} // namespace grackle::impl

#endif /* INIT_EXTRA_COLLISIONAL_RATES_HPP */

