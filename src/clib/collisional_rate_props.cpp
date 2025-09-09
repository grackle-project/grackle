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
/// "Classic" Rate coefficients. All rates are labelled as in Abel et al., 1996
/// (astro-ph/9608040)
///
/// | rate-name | Reaction                         |
/// | --------- | -------------------------------- |
/// |  k1       | HI + e --> HII + 2e              |
/// |  k2       | HII + e --> HI + photon          |
/// |  k3       | HeI + e --> HeII + 2e            |
/// |  k4       | HeII + e --> HeI + photon        |
/// |  k5       | HeII + e --> HeIII + 2e          |
/// |  k6       | HeIII + e --> HeII + photon      |
/// |  k7       | HI + e --> HM + photon           |
/// |  k8       | HI + HM --> H2I* + e             |
/// |  k9       | HI + HII --> H2II + photon       |
/// |  k10      | H2II + HI --> H2I* + HII         |
/// |  k11      | H2I + HII --> H2II + HI          |
/// |  k12      | H2I + e --> 2HI + e              |
/// |  k13      | H2I + HI --> 3HI                 |
/// |  k14      | HM + e --> HI + 2e               |
/// |  k15      | HM + HI --> 2HI + e              |
/// |  k16      | HM + HI --> 2HI                  |
/// |  k17      | HM + HI --> H2I + e              |
/// |  k18      | H2I + e --> 2HI                  |
/// |  k19      | H2I + HM --> H2I + HI            |
/// |  k20      | Not Used                         |
/// |  k21      | 2HI + H2I --> H2I + H2I          |
/// |  k22      | 2HI + HI --> H2I + HI            |
/// |  k24      | HI + p --> HII + e               |
/// |  k25      | HeIII + p --> HeII + e           |
/// |  k26      | HeI + p --> HeII + e             |
/// |  k27      | HM + p --> HI + e                |
/// |  k28      | H2II + p --> HI + HII            |
/// |  k29      | H2I + p --> H2II + e             |
/// |  k30      | H2II + p --> 2HII + e            |
/// |  k31      | H2I + p --> 2HI                  |
/// |  k50      | HII + DI --> HI + DII            |
/// |  k51      | HI + DII --> HII + DI            |
/// |  k52      | H2I + DII --> HDI + HII          |
/// |  k53      | HDI + HII --> H2I + DII          |
/// |  k54      | H2I + DI --> HDI + HI            |
/// |  k55      | HDI + HI --> H2I + DI            |
/// |  k56      | DI + HM --> HDI + e              |
/// |           | DM + HI --> HDI + e  (SEE BELOW) |
/// |  k57      | HI + HI --> HII + HI + e         |
/// |  k58      | HI + HeI --> HII + HeI + e       |
///
/// DM + HI --> HDI + e is included implicitly by multiplying k56 by two as
/// they are assumed to have the same rate.
///
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

#include <cmath>

#include "grackle.h"
#include "grackle_rate_functions.h"
#include "collisional_rate_props.hpp"  // forward declaration
#include "grackle_macros.h"            // tiny
#include "internal_types.hpp"          // CollisionalRxnRateCollection
#include "internal_units.h"            // InternalGrUnits
#include "LUT.hpp"                     // CollisionalRxnLUT
#include "opaque_storage.hpp"          // gr_opaque_storage

// define rates for primordial_chemistry == 4 and metal_chemistry == 1
// -------------------------------------------------------------------
//
// At some point, we probably want to treat these all a little more consistently
// with the rates for primordial_chemistry == 1, 2, and 3
// - personally, I'm not sold on the fact that we want to expose all of these
//   as individual public functions (see GitHub issue #414)
//
// If we do decide to expose all of these as public functions, we need to do a
// few things first:
// 1. We should add docstrings to **EVERY** single function that includes the
//    citation to the paper where the reaction info was taken from
// 2. We **MUST** figure out what the comment before kz22_rate_ actually means
//    and update the docstrings of all rate functions that the comment is
//    talking about accordingly
// 3. delete all occurrences of GRIMPL_NOEXPORT in this file


/// Annotating a function with this macro means that "symbol" is hidden when
/// compiling shared libraries.
///
/// If the gnu::visibility attribute is not known to the compiler (e.g. you
/// aren't using g++ or clang++), this doesn't do anything.
///
/// When to use
/// ===========
/// This is currently intended to be used over the `static` keyword when we
/// want a function to have C linkage but not be publicly exposed. To 0th
/// order, you can think of this as being equivalent to the `static` keyword.
/// (In reality it means something somewhat different).
///
/// More Context
/// ============
/// If you aren't familiar with the concept of symbol visibility, this paper
/// may be useful to read
/// https://cs.dartmouth.edu/~sergey/cs258/ABI/UlrichDrepper-How-To-Write-Shared-Libraries.pdf
///
/// But, let's give a quick overview. When you create a shared C library,
/// all functions and global variables are associated with "symbols". By
/// default, in unix environments, all of the symbols (except functions/global
/// variables declared with `static`) are fully exposed (i.e. visible) to
/// external programs in a table (stored in the shared library).
/// - An external C program typically accesses these functions by forward
///   declaring the signatures. The convention is to only access the symbols
///   with declarations provided in public headers and expressly described as
///   part of the public API (but in principle, a user could access any visible
///   symbol).
/// - An external C program might dynamically load symbols at runtime (via
///   dlopen), but that's a topic for another time.
///
/// The macro described here explicitly indicates that all annotated functions
/// should not be visible (i.e. they are hidden from external programs).
///
/// If we want to apply this more broadly
/// =====================================
/// In an ideal world, every single Grackle function, other than public API,
/// would be hidden. The benefits include:
/// - reducing the size of the Grackle shared library
/// - improving startup time of programs linked agains a Grackle shared library
/// - enforcing that private functions shouldn't be used everywhere
///
/// If we choose to go this path, then it would be better to make all functions
/// have hidden visibility, by default, and explicitly annotate the public
/// functions. This exact strategy is used on Windows. CMake has support to
/// help do this across compiler vendors
///   https://cmake.org/cmake/help/latest/module/GenerateExportHeader.html
/// Here is a guide that applies to gcc (and probably clang)
///   https://gcc.gnu.org/wiki/Visibility
#define GRIMPL_NOEXPORT [[gnu::visibility("hidden")]]

extern "C" {

[[gnu::visibility("hidden")]] double k125_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k125 = 6.4e-10;
  return std::fmax(k125, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k129_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k129 = 3.9e-19 * pow((T/300.0), 1.8) * exp(20.0/T);
  return std::fmax(k129, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k130_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k130 = 3.9e-19 * pow((T/300.0), 1.8) * exp(20.0/T);
  return std::fmax(k130, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k131_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k131 = 3.4e-9 * pow((T/300.0), -0.4);
  return std::fmax(k131, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k132_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k132 = 3.0e-16 * pow((T/300.0), 0.95) * exp(-T/9320.0);
  return std::fmax(k132, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k133_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k133 = 5.7e-8 * pow((T/300.0), -0.50);
  return std::fmax(k133, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k134_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k134 = 4.6e-8 * pow((T/300.0), -0.50);
  return std::fmax(k134, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k135_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k135 = 4.6e-8 * pow((T/300.0), -0.50);
  return std::fmax(k135, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k136_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k136 = 6.4e-9 * pow((T/300.0), 0.41);
  return std::fmax(k136, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k137_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k137 = 1.5e-9 * pow((T/300.0), -0.1);
  return std::fmax(k137, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k148_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k148 = 5.0e-21;
  return std::fmax(k148, tiny) / kunit;
}

double k149_rate_(double T, double kunit, chemistry_data *my_chemistry) {
  double k149;
  if(T < 1000.0) {
    k149 = 7.60e-18 * pow(T, -0.50);
  } else {
    k149 = 3.45e-16 * pow(T, -1.06);
  }
  return std::fmax(k149, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k150_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k150 = 3.0e-10 * exp(-6717.0/T);
  return std::fmax(k150, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k151_rate_(double T, double kunit, chemistry_data *my_chemistry) {
  double k151;
  if(T < 4000.0) {
    k151 = 1.6e-14 * pow(T, -0.33);
  } else {
    k151 = 1.0e-15;
  }
  return std::fmax(k151, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k152_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k152 = 9.1e-10;
  return std::fmax(k152, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double k153_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double k153 = 1.7e-7 * pow(T, -0.5);
  return std::fmax(k153, tiny) / kunit;
}

  
[[gnu::visibility("hidden")]] double kz15_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz15 = 4.98e-11;
  return std::fmax(kz15, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz16_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz16 = 2.70e-10;
  return std::fmax(kz16, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz17_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz17 = 7.00e-14 * pow((T/300.0), 2.80) * exp(-1950.0/T);
  return std::fmax(kz17, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz18_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz18 = 6.83e-12 * pow((T/300.0), 1.60) * exp(-9720.0/T);
  return std::fmax(kz18, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz19_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz19 = 3.30e-10 * exp(-8460.0/T);
  return std::fmax(kz19, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz20_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz20 = 6.64e-10 * exp(-11700.0/T);
  return std::fmax(kz20, tiny) / kunit;
}


[[gnu::visibility("hidden")]] double kz21_rate_(double T, double kunit, chemistry_data *my_chemistry) {
  double kz21;
  if(T < 1.0e7) {
    kz21 = 3.43e-13 * pow((T / 300.0), 2.67) * exp(-3160.0/T);
  } else {
    kz21 = 3.43e-13 * pow(1.0e7/300.0, 2.67) * exp(-3160.0/1.0e7);
  }
  return std::fmax(kz21, tiny) / kunit;
}

// The following comment has been adapted from earlier versions of Grackle
//    The rate comes from an experiment (297-3532 K).
//    We refrain to extrapolate it to high temperatures.
//
// There are 2 issues here (that were not apparent in older Grackle versions):
// 1. Does this comment just apply to kz22? Or does it apply to kz22 and some
//    of the subsequent rates? (which ones?) Or all of the subsequent rates?
// 2. The statement "We refrain to extrapolate it to high temperatures" needs
//    to be clarified...
//    - I assume that it should read "we refrain from extrapolating it to high
//      temperatures"
//    - but, we don't obviously do that. It seems like we use the values at ALL
//      temperatures.

[[gnu::visibility("hidden")]] double kz22_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz22 = 7.00e-10 * exp(-232.0/T);
  return std::fmax(kz22, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz23_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz23 = 2.38e-10 * exp(-1760.0/T);
  return std::fmax(kz23, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz24_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz24 = 1.55e-12 * pow((T/300.0), 1.60) * exp(-1660.0/T);
  return std::fmax(kz24, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz25_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz25 = 1.65e-12 * pow((T/300.0), 1.14) * exp(-50.0/T);
  return std::fmax(kz25, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz26_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz26 = 1.0e-13;
  return std::fmax(kz26, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz27_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz27 = 1.0e-17;
  return std::fmax(kz27, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz28_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz28 = 1.1e-10 * pow((T/300.0), 0.5);
  return std::fmax(kz28, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz29_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz29 = 3.3e-11;
  return std::fmax(kz29, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz30_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz30 = 9.9e-19 * pow((T/300.0), -0.38);
  return std::fmax(kz30, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz31_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz31 = 4.9e-20 * pow((T/300.0), 1.58);
  return std::fmax(kz31, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz32_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz32 = 6.6e-11;
  return std::fmax(kz32, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz33_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz33 = 4.34e-11 * pow((T/300.0), -0.5) * exp(-30.0/T);
  return std::fmax(kz33, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz34_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz34 = 2.1e-9;
  return std::fmax(kz34, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz35_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz35 = 6.9e-9;
  return std::fmax(kz35, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz36_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz36 = 2.0e-9;
  return std::fmax(kz36, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz37_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz37 = 7.7e-10;
  return std::fmax(kz37, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz38_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz38 = 6.2e-10;
  return std::fmax(kz38, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz39_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz39 = 6.8e-10;
  return std::fmax(kz39, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz40_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz40 = 1.7e-9;
  return std::fmax(kz40, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz41_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz41 = 1.01e-9;
  return std::fmax(kz41, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz42_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz42 = 8.3e-10;
  return std::fmax(kz42, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz43_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz43 = 7.5e-10;
  return std::fmax(kz43, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz44_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz44 = 4.4e-12 * pow((T/300.0), -0.61);
  return std::fmax(kz44, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz45_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz45 = 3.4e-12 * pow((T/300.0), -0.63);
  return std::fmax(kz45, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz46_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz46 = 1.6e-7 * pow((T/300.0), -0.5);
  return std::fmax(kz46, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz47_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz47 = 2.0e-7 * pow((T/300.0), -0.5);
  return std::fmax(kz47, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz48_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz48 = 3.5e-7 * pow((T/300.0), -0.5);
  return std::fmax(kz48, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz49_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz49 = 6.5e-7 * pow((T/300.0), -0.5);
  return std::fmax(kz49, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz50_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz50 = 1.95e-7 * pow((T/300.0), -0.7);
  return std::fmax(kz50, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz51_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz51 = 1.0e-17;
  return std::fmax(kz51, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz52_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz52 = 3.00e-11;
  return std::fmax(kz52, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz53_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz53 = 1.30e-11 * exp(-111.0/T);
  return std::fmax(kz53, tiny) / kunit;
}

[[gnu::visibility("hidden")]] double kz54_rate_(double T, double kunit, chemistry_data* my_chemistry) {
  double kz54 = 2.00e-13;
  return std::fmax(kz54, tiny) / kunit;
}

} // extern "C"

// define functions directly exposed to other parts of Grackle
// -----------------------------------------------------------

int grackle::impl::init_extra_collisional_rates(
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
    (std::log(my_chemistry->TemperatureEnd) -
     std::log(my_chemistry->TemperatureStart)) /
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
    double logT = std::log(my_chemistry->TemperatureStart) + (double)(i)*dlogtem;
    double T = std::exp(logT);

    kcol_rate_tables->data[CollisionalRxnLUT::k125][i] = k125_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k129][i] = k129_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k130][i] = k130_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k131][i] = k131_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k132][i] = k132_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k133][i] = k133_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k134][i] = k134_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k135][i] = k135_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k136][i] = k136_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k137][i] = k137_rate_(T, kunit, my_chemistry);
  
    kcol_rate_tables->data[CollisionalRxnLUT::k148][i] = k148_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k149][i] = k149_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k150][i] = k150_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k151][i] = k151_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k152][i] = k152_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::k153][i] = k153_rate_(T, kunit, my_chemistry);
  
    kcol_rate_tables->data[CollisionalRxnLUT::kz15][i] = kz15_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz16][i] = kz16_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz17][i] = kz17_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz18][i] = kz18_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz19][i] = kz19_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz20][i] = kz20_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz21][i] = kz21_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz22][i] = kz22_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz23][i] = kz23_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz24][i] = kz24_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz25][i] = kz25_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz26][i] = kz26_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz27][i] = kz27_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz28][i] = kz28_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz29][i] = kz29_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz30][i] = kz30_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz31][i] = kz31_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz32][i] = kz32_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz33][i] = kz33_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz34][i] = kz34_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz35][i] = kz35_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz36][i] = kz36_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz37][i] = kz37_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz38][i] = kz38_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz39][i] = kz39_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz40][i] = kz40_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz41][i] = kz41_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz42][i] = kz42_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz43][i] = kz43_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz44][i] = kz44_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz45][i] = kz45_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz46][i] = kz46_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz47][i] = kz47_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz48][i] = kz48_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz49][i] = kz49_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz50][i] = kz50_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz51][i] = kz51_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz52][i] = kz52_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz53][i] = kz53_rate_(T, kunit, my_chemistry);
    kcol_rate_tables->data[CollisionalRxnLUT::kz54][i] = kz54_rate_(T, kunit, my_chemistry);
  }

  return GR_SUCCESS;
}


int grackle::impl::visit_rate_props(const chemistry_data* my_chemistry,
                                    grackle::impl::visit_kcol_prop_callback* cb,
                                    void* ctx)
{
  // shorten `grackle::impl` to `grimpl` within this function
  namespace grimpl = ::grackle::impl;

  if (my_chemistry->primordial_chemistry > 0) {
    cb(grimpl::KColProp{CollisionalRxnLUT::k1, &k1_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k3, &k3_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k4, &k4_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k2, &k2_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k5, &k5_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k6, &k6_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k7, &k7_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k8, &k8_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k9, &k9_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k10, &k10_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k11, &k11_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k12, &k12_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k14, &k14_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k15, &k15_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k16, &k16_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k17, &k17_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k18, &k18_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k19, &k19_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k20, &k20_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k23, &k23_rate, true}, ctx);

    //--------Calculate 3-body H2 rate--------

    // Calculated by the same method as done in the original code. First is the
    // fit to A.E. Orel 1987, J.Chem.Phys., 87, 314, which is matched to the 1/T
    // of Palla e tal (1983) -- which is four times smaller than the Palla rate.

    // Varying threebody and corresponding collisional dissociation rates from
    // Simon.
    cb(grimpl::KColProp{CollisionalRxnLUT::k13, &k13_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k21, &k21_rate, false}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k22, &k22_rate, false}, ctx);

    //--------Deuterium Rates--------
    cb(grimpl::KColProp{CollisionalRxnLUT::k50, &k50_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k51, &k51_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k52, &k52_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k53, &k53_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k54, &k54_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k55, &k55_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k56, &k56_rate, true}, ctx);

    //--------New H Ionization Rates--------
    cb(grimpl::KColProp{CollisionalRxnLUT::k57, &k57_rate, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k58, &k58_rate, true}, ctx);

    // primordial_chem > 3:
    cb(grimpl::KColProp{CollisionalRxnLUT::k125, &k125_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k129, &k129_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k130, &k130_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k131, &k131_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k132, &k132_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k133, &k133_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k134, &k134_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k135, &k135_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k136, &k136_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k137, &k137_rate_, true}, ctx);
  
    cb(grimpl::KColProp{CollisionalRxnLUT::k148, &k148_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k149, &k149_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k150, &k150_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k151, &k151_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k152, &k152_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::k153, &k153_rate_, true}, ctx);
    
    // primordial_chem > 3 && metal_chemistry == 1
    cb(grimpl::KColProp{CollisionalRxnLUT::kz15, &kz15_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz16, &kz16_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz17, &kz17_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz18, &kz18_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz19, &kz19_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz20, &kz20_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz21, &kz21_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz22, &kz22_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz23, &kz23_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz24, &kz24_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz25, &kz25_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz26, &kz26_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz27, &kz27_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz28, &kz28_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz29, &kz29_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz30, &kz30_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz31, &kz31_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz32, &kz32_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz33, &kz33_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz34, &kz34_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz35, &kz35_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz36, &kz36_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz37, &kz37_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz38, &kz38_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz39, &kz39_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz40, &kz40_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz41, &kz41_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz42, &kz42_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz43, &kz43_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz44, &kz44_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz45, &kz45_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz46, &kz46_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz47, &kz47_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz48, &kz48_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz49, &kz49_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz50, &kz50_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz51, &kz51_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz52, &kz52_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz53, &kz53_rate_, true}, ctx);
    cb(grimpl::KColProp{CollisionalRxnLUT::kz54, &kz54_rate_, true}, ctx);
  }

  return GR_SUCCESS;
}
