// Header file containing all rate function declarations.
#include "grackle_chemistry_data.h"

#ifndef RATE_FUNCTIONS_H
#define RATE_FUNCTIONS_H

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

double gasGrain_rate(double T, double units, chemistry_data *my_chemistry);
double regr_rate(double T, double units, chemistry_data *my_chemistry);

double comp_rate(double units, chemistry_data *my_chemistry);
double gammah_rate(double units, chemistry_data *my_chemistry);
double gamma_isrf_rate(double units, chemistry_data *my_chemistry);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RATE_FUNCTIONS_H */
