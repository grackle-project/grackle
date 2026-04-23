//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares rate functions
///
//===----------------------------------------------------------------------===//

#include "grackle.h"

#ifndef GRACKLE_RATE_FUNCTIONS_H
#define GRACKLE_RATE_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

double k1_rate(double T, double units,  chemistry_data *my_chemistry);
double k2_rate(double T, double units,  chemistry_data *my_chemistry);
double k3_rate(double T, double units,  chemistry_data *my_chemistry);
double k4_rate(double T, double units,  chemistry_data *my_chemistry);
double k5_rate(double T, double units,  chemistry_data *my_chemistry);
double k6_rate(double T, double units,  chemistry_data *my_chemistry);
double k7_rate(double T, double units,  chemistry_data *my_chemistry);
double k8_rate(double T, double units,  chemistry_data *my_chemistry);
double k9_rate(double T, double units,  chemistry_data *my_chemistry);
double k10_rate(double T, double units,  chemistry_data *my_chemistry);
double k11_rate(double T, double units,  chemistry_data *my_chemistry);
double k12_rate(double T, double units,  chemistry_data *my_chemistry);
double k13_rate(double T, double units,  chemistry_data *my_chemistry);

void k13dd_rate(double T, double units, double *k13dd_results, chemistry_data *my_chemistry);

double k14_rate(double T, double units,  chemistry_data *my_chemistry);
double k15_rate(double T, double units,  chemistry_data *my_chemistry);
double k16_rate(double T, double units,  chemistry_data *my_chemistry);
double k17_rate(double T, double units,  chemistry_data *my_chemistry);
double k18_rate(double T, double units,  chemistry_data *my_chemistry);
double k19_rate(double T, double units,  chemistry_data *my_chemistry);
double k20_rate(double T, double units,  chemistry_data *my_chemistry);
double k21_rate(double T, double units,  chemistry_data *my_chemistry);
double k22_rate(double T, double units,  chemistry_data *my_chemistry);
double k23_rate(double T, double units,  chemistry_data *my_chemistry);
double k50_rate(double T, double units,  chemistry_data *my_chemistry);
double k51_rate(double T, double units,  chemistry_data *my_chemistry);
double k52_rate(double T, double units,  chemistry_data *my_chemistry);
double k53_rate(double T, double units,  chemistry_data *my_chemistry);
double k54_rate(double T, double units,  chemistry_data *my_chemistry);
double k55_rate(double T, double units,  chemistry_data *my_chemistry);
double k56_rate(double T, double units,  chemistry_data *my_chemistry);
double k57_rate(double T, double units,  chemistry_data *my_chemistry);
double k58_rate(double T, double units,  chemistry_data *my_chemistry);

double h2dust_rate(double T, double T_dust, double units, chemistry_data *my_chemistry);
double h2dust_C_rate(double T, double T_dust, double units, chemistry_data *my_chemistry);
double h2dust_S_rate(double T, double T_dust, double units, chemistry_data *my_chemistry);

double n_cr_n_rate(double T, double units,  chemistry_data *my_chemistry);
double n_cr_d1_rate(double T, double units,  chemistry_data *my_chemistry);
double n_cr_d2_rate(double T, double units,  chemistry_data *my_chemistry);

double ceHI_rate(double T, double units, chemistry_data *my_chemistry);
double ceHeI_rate(double T, double units, chemistry_data *my_chemistry);
double ceHeII_rate(double T, double units, chemistry_data *my_chemistry);

double ciHeIS_rate(double T, double units, chemistry_data *my_chemistry);
double ciHI_rate(double T, double units, chemistry_data *my_chemistry);
double ciHeI_rate(double T, double units, chemistry_data *my_chemistry);
double ciHeII_rate(double T, double units, chemistry_data *my_chemistry);

double reHII_rate(double T, double units, chemistry_data *my_chemistry);
double reHeII1_rate(double T, double units, chemistry_data *my_chemistry);
double reHeII2_rate(double T, double units, chemistry_data *my_chemistry);
double reHeIII_rate(double T, double units, chemistry_data *my_chemistry);

double brem_rate(double T, double units, chemistry_data *my_chemistry);

double vibh_rate(double T, double units, chemistry_data *my_chemistry);
double hyd01k_rate(double T, double units, chemistry_data *my_chemistry);
double h2k01_rate(double T, double units, chemistry_data *my_chemistry);
double rotl_rate(double T, double units, chemistry_data *my_chemistry);
double roth_rate(double T, double units, chemistry_data *my_chemistry);

double GP99LowDensityLimit_rate(double T, double units, chemistry_data *my_chemistry);
double GP99HighDensityLimit_rate(double T, double units, chemistry_data *my_chemistry);
double GAHI_rate(double T, double units, chemistry_data *my_chemistry);
double GAH2_rate(double T, double units, chemistry_data *my_chemistry);
double GAHe_rate(double T, double units, chemistry_data *my_chemistry);
double GAHp_rate(double T, double units, chemistry_data *my_chemistry);
double GAel_rate(double T, double units, chemistry_data *my_chemistry);
double H2LTE_rate(double T, double units, chemistry_data *my_chemistry);
double HDlte_rate(double T, double units, chemistry_data *my_chemistry);

double HDlow_rate(double T, double units, chemistry_data *my_chemistry);
double cie_thin_cooling_rate(double T);
double cieco_rate(double T, double units, chemistry_data *my_chemistry);

/// Calculate gas_grain, the Gas/grain energy transfer rate
double gasGrain_rate(double T, double units, chemistry_data *my_chemistry);

/// Calculate gas_grain2
///
/// The resulting value is similar to the value returned by @ref gasGrain_rate.
/// It is used to compute the Gas/grain energy transfer rate for arbitrary size
/// distributions.
double gasGrain2_rate(double T, double units, chemistry_data *my_chemistry);

/// Calculate regr (as in GRain REcombination cooling)
///
/// This is the value of a cluster of variables taken from equation 9 from
/// [Wolfire+95](https://ui.adsabs.harvard.edu/abs/1995ApJ...443..152W/abstract).
double regr_rate(double T, double units, chemistry_data *my_chemistry);

double grain_growth_rate(double T, double units, chemistry_data *my_chemistry);

double comp_rate(double units, chemistry_data *my_chemistry);

/// Calculate quantity closely related to Γ for photo-electric heating from
/// photo-electric heating by dust-grains
///
/// For added context, we currently use equation 1 (and possibly eqn 2) of
/// [Wolfire+95](https://ui.adsabs.harvard.edu/abs/1995ApJ...443..152W/abstract)
/// to implement photo-electric heating
double gammah_rate(double units, chemistry_data *my_chemistry);
double gamma_isrf_rate(double units, chemistry_data *my_chemistry);
double gamma_isrf2_rate(double units, chemistry_data *my_chemistry);

// given that all of the following are interpolation tables, I don't think we
// should be exposing these functions as part of the public api. These are
// fundamentally different from the other functions: they are initializing
// entries tracked by chemistry_data_storage

void initialize_cooling_rate_H2(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_HD(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_primordial_opacity(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

void initialize_cooling_rate_CI (chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_CII(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_OI (chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_CO (chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_OH (chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);
void initialize_cooling_rate_H2O(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double coolunit);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* GRACKLE_RATE_FUNCTIONS_H */
