#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_chemistry_data.h"

int calc_rates_dust_loc(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_C13(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_C20(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_C25(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_C30(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_F13(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_F15(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_F50(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_F80(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_P170(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_P200(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);
int calc_rates_dust_Y19(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

typedef int calc_yield_rate_fn(int, chemistry_data*, chemistry_data_storage*);

int initialize_dust_yields(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           code_units *my_units)
{

  if (my_chemistry->metal_chemistry == 0)
    return SUCCESS;

//-------kdMgSiO3  :    Mg   +    SiO + 2H2O -> MgSiO3  + 2H2I
//-------kdAC      :    C                    -> AC      
                                 
//-------kdSiM     :    Si                   -> SiM     
//-------kdFeM     :    Fe                   -> FeM     
//-------kdMg2SiO4 :   2Mg   +    SiO + 3H2O -> Mg2SiO4 + 3H2I
//-------kdFe3O4   :   3Fe   +   4H2O        -> Fe3O4   + 4H2I
//-------kdSiO2D   :    SiO2                 -> SiO2D   
//-------kdMgO     :    Mg   +    H2O        -> MgO     +  H2I
//-------kdFeS     :    Fe   +    S          -> FeS     
//-------kdAl2O3   :   2Al   +   3H2O        -> Al2O3   + 3H2I
                     
//-------kdreforg  : 0.5CO   + 0.5CH2 + 1.2N -> reforg (C:H:O:N = 1:1:0.5:1.2)
//-------kdvolorg  :    CO   +   2H2I        -> volorg (CH3OH)
//-------kdH2Oice  :    H2O                  -> H2Oice

      int NSN, NTd, Nmom;
      double Td0, dTd;
      int iSN, iTd, imom, itab;

      NSN = 12;
      my_rates->SN0_N = NSN;

      my_rates->SN0_XC  = malloc(NSN * sizeof(double));
      my_rates->SN0_XO  = malloc(NSN * sizeof(double));
      my_rates->SN0_XMg = malloc(NSN * sizeof(double));
      my_rates->SN0_XAl = malloc(NSN * sizeof(double));
      my_rates->SN0_XSi = malloc(NSN * sizeof(double));
      my_rates->SN0_XS  = malloc(NSN * sizeof(double));
      my_rates->SN0_XFe = malloc(NSN * sizeof(double));

      my_rates->SN0_fC  = malloc(NSN * sizeof(double));
      my_rates->SN0_fO  = malloc(NSN * sizeof(double));
      my_rates->SN0_fMg = malloc(NSN * sizeof(double));
      my_rates->SN0_fAl = malloc(NSN * sizeof(double));
      my_rates->SN0_fSi = malloc(NSN * sizeof(double));
      my_rates->SN0_fS  = malloc(NSN * sizeof(double));
      my_rates->SN0_fFe = malloc(NSN * sizeof(double));

      my_rates->SN0_fSiM      = malloc(NSN * sizeof(double));
      my_rates->SN0_fFeM      = malloc(NSN * sizeof(double));
      my_rates->SN0_fMg2SiO4  = malloc(NSN * sizeof(double));
      my_rates->SN0_fMgSiO3   = malloc(NSN * sizeof(double));
      my_rates->SN0_fFe3O4    = malloc(NSN * sizeof(double));
      my_rates->SN0_fAC       = malloc(NSN * sizeof(double));
      my_rates->SN0_fSiO2D    = malloc(NSN * sizeof(double));
      my_rates->SN0_fMgO      = malloc(NSN * sizeof(double));
      my_rates->SN0_fFeS      = malloc(NSN * sizeof(double));
      my_rates->SN0_fAl2O3    = malloc(NSN * sizeof(double));
      my_rates->SN0_freforg   = malloc(NSN * sizeof(double));
      my_rates->SN0_fvolorg   = malloc(NSN * sizeof(double));
      my_rates->SN0_fH2Oice   = malloc(NSN * sizeof(double));

      for(iSN = 0; iSN < NSN; iSN++) {
        my_rates->SN0_XC [iSN] = 0.0;
        my_rates->SN0_XO [iSN] = 0.0;
        my_rates->SN0_XMg[iSN] = 0.0;
        my_rates->SN0_XAl[iSN] = 0.0;
        my_rates->SN0_XSi[iSN] = 0.0;
        my_rates->SN0_XS [iSN] = 0.0;
        my_rates->SN0_XFe[iSN] = 0.0;

        my_rates->SN0_fC [iSN] = 0.0;
        my_rates->SN0_fO [iSN] = 0.0;
        my_rates->SN0_fMg[iSN] = 0.0;
        my_rates->SN0_fAl[iSN] = 0.0;
        my_rates->SN0_fSi[iSN] = 0.0;
        my_rates->SN0_fS [iSN] = 0.0;
        my_rates->SN0_fFe[iSN] = 0.0;

        my_rates->SN0_fSiM     [iSN] = 0.0;
        my_rates->SN0_fFeM     [iSN] = 0.0;
        my_rates->SN0_fMg2SiO4 [iSN] = 0.0;
        my_rates->SN0_fMgSiO3  [iSN] = 0.0;
        my_rates->SN0_fFe3O4   [iSN] = 0.0;
        my_rates->SN0_fAC      [iSN] = 0.0;
        my_rates->SN0_fSiO2D   [iSN] = 0.0;
        my_rates->SN0_fMgO     [iSN] = 0.0;
        my_rates->SN0_fFeS     [iSN] = 0.0;
        my_rates->SN0_fAl2O3   [iSN] = 0.0;
        my_rates->SN0_freforg  [iSN] = 0.0;
        my_rates->SN0_fvolorg  [iSN] = 0.0;
        my_rates->SN0_fH2Oice  [iSN] = 0.0;
      }

      my_rates->SN0_r0SiM      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0FeM      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Mg2SiO4  = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0MgSiO3   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Fe3O4    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0AC       = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0SiO2D    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0MgO      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0FeS      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Al2O3    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0reforg   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0volorg   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0H2Oice   = malloc(NSN * 3 * sizeof(double));

      itab = 0;
      for(iSN = 0; iSN < NSN; iSN++) {
        for(imom = 0; imom < 3; imom++) {
          my_rates->SN0_r0SiM     [itab] = 0.0;
          my_rates->SN0_r0FeM     [itab] = 0.0;
          my_rates->SN0_r0Mg2SiO4 [itab] = 0.0;
          my_rates->SN0_r0MgSiO3  [itab] = 0.0;
          my_rates->SN0_r0Fe3O4   [itab] = 0.0;
          my_rates->SN0_r0AC      [itab] = 0.0;
          my_rates->SN0_r0SiO2D   [itab] = 0.0;
          my_rates->SN0_r0MgO     [itab] = 0.0;
          my_rates->SN0_r0FeS     [itab] = 0.0;
          my_rates->SN0_r0Al2O3   [itab] = 0.0;
          my_rates->SN0_r0reforg  [itab] = 0.0;
          my_rates->SN0_r0volorg  [itab] = 0.0;
          my_rates->SN0_r0H2Oice  [itab] = 0.0;
          itab++;
        }
      }

      NTd =            35;
      Td0 =     0.0000000;
      dTd =     0.1000000;
     Nmom =             4;

      my_rates->gr_Td = malloc(NTd * Nmom * sizeof(double));
      my_rates->gr_Size = NTd * Nmom;
      my_rates->gr_N[0] = Nmom;
      my_rates->gr_N[1] = NTd;
      my_rates->gr_dT   = dTd;
      for(iTd = 0; iTd < NTd; iTd++)
        my_rates->gr_Td[iTd] = Td0 + (double)iTd * dTd;

      my_rates->SN0_kpSiM      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeM      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMg2SiO4  = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgSiO3   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFe3O4    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAC       = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpSiO2D    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgO      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeS      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAl2O3    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpreforg   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpvolorg   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpH2Oice   = malloc(NSN * Nmom * NTd * sizeof(double));

      itab = 0;
      for(iSN = 0; iSN < NSN; iSN++) {
        for(imom = 0; imom < Nmom; imom++) {
          for(iTd = 0; iTd < NTd; iTd++) {
            my_rates->SN0_kpSiM     [itab] = 0.0;
            my_rates->SN0_kpFeM     [itab] = 0.0;
            my_rates->SN0_kpMg2SiO4 [itab] = 0.0;
            my_rates->SN0_kpMgSiO3  [itab] = 0.0;
            my_rates->SN0_kpFe3O4   [itab] = 0.0;
            my_rates->SN0_kpAC      [itab] = 0.0;
            my_rates->SN0_kpSiO2D   [itab] = 0.0;
            my_rates->SN0_kpMgO     [itab] = 0.0;
            my_rates->SN0_kpFeS     [itab] = 0.0;
            my_rates->SN0_kpAl2O3   [itab] = 0.0;
            my_rates->SN0_kpreforg  [itab] = 0.0;
            my_rates->SN0_kpvolorg  [itab] = 0.0;
            my_rates->SN0_kpH2Oice  [itab] = 0.0;
            itab++;
          }
        }
      }

      calc_yield_rate_fn* fn_list[] = {
        &calc_rates_dust_loc, &calc_rates_dust_C13, &calc_rates_dust_C20,
        &calc_rates_dust_C25, &calc_rates_dust_C30, &calc_rates_dust_F13,
        &calc_rates_dust_F15, &calc_rates_dust_F50, &calc_rates_dust_F80,
        &calc_rates_dust_P170, &calc_rates_dust_P200, &calc_rates_dust_Y19
      };

      int n_funcs = (int)(sizeof(fn_list) / sizeof(calc_yield_rate_fn*));

      for (int i = 0; i < n_funcs; i++) {
        int rv = fn_list[i](i, my_chemistry, my_rates);
        if (rv != SUCCESS) { return rv; }
      }

  return SUCCESS;
}

int local_free_dust_yields(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates)
{

  if (my_chemistry->metal_chemistry == 0)
    return SUCCESS;

  GRACKLE_FREE(my_rates->SN0_XC);
  GRACKLE_FREE(my_rates->SN0_XO);
  GRACKLE_FREE(my_rates->SN0_XMg);
  GRACKLE_FREE(my_rates->SN0_XAl);
  GRACKLE_FREE(my_rates->SN0_XSi);
  GRACKLE_FREE(my_rates->SN0_XS);
  GRACKLE_FREE(my_rates->SN0_XFe);

  GRACKLE_FREE(my_rates->SN0_fC);
  GRACKLE_FREE(my_rates->SN0_fO);
  GRACKLE_FREE(my_rates->SN0_fMg);
  GRACKLE_FREE(my_rates->SN0_fAl);
  GRACKLE_FREE(my_rates->SN0_fSi);
  GRACKLE_FREE(my_rates->SN0_fS);
  GRACKLE_FREE(my_rates->SN0_fFe);

  GRACKLE_FREE(my_rates->SN0_fSiM);
  GRACKLE_FREE(my_rates->SN0_fFeM);
  GRACKLE_FREE(my_rates->SN0_fMg2SiO4);
  GRACKLE_FREE(my_rates->SN0_fMgSiO3);
  GRACKLE_FREE(my_rates->SN0_fFe3O4);
  GRACKLE_FREE(my_rates->SN0_fAC);
  GRACKLE_FREE(my_rates->SN0_fSiO2D);
  GRACKLE_FREE(my_rates->SN0_fMgO);
  GRACKLE_FREE(my_rates->SN0_fFeS);
  GRACKLE_FREE(my_rates->SN0_fAl2O3);
  GRACKLE_FREE(my_rates->SN0_freforg);
  GRACKLE_FREE(my_rates->SN0_fvolorg);
  GRACKLE_FREE(my_rates->SN0_fH2Oice);

  GRACKLE_FREE(my_rates->SN0_r0SiM);
  GRACKLE_FREE(my_rates->SN0_r0FeM);
  GRACKLE_FREE(my_rates->SN0_r0Mg2SiO4);
  GRACKLE_FREE(my_rates->SN0_r0MgSiO3);
  GRACKLE_FREE(my_rates->SN0_r0Fe3O4);
  GRACKLE_FREE(my_rates->SN0_r0AC);
  GRACKLE_FREE(my_rates->SN0_r0SiO2D);
  GRACKLE_FREE(my_rates->SN0_r0MgO);
  GRACKLE_FREE(my_rates->SN0_r0FeS);
  GRACKLE_FREE(my_rates->SN0_r0Al2O3);
  GRACKLE_FREE(my_rates->SN0_r0reforg);
  GRACKLE_FREE(my_rates->SN0_r0volorg);
  GRACKLE_FREE(my_rates->SN0_r0H2Oice);

  GRACKLE_FREE(my_rates->gr_Td);

  GRACKLE_FREE(my_rates->SN0_kpSiM);
  GRACKLE_FREE(my_rates->SN0_kpFeM);
  GRACKLE_FREE(my_rates->SN0_kpMg2SiO4);
  GRACKLE_FREE(my_rates->SN0_kpMgSiO3);
  GRACKLE_FREE(my_rates->SN0_kpFe3O4);
  GRACKLE_FREE(my_rates->SN0_kpAC);
  GRACKLE_FREE(my_rates->SN0_kpSiO2D);
  GRACKLE_FREE(my_rates->SN0_kpMgO);
  GRACKLE_FREE(my_rates->SN0_kpFeS);
  GRACKLE_FREE(my_rates->SN0_kpAl2O3);
  GRACKLE_FREE(my_rates->SN0_kpreforg);
  GRACKLE_FREE(my_rates->SN0_kpvolorg);
  GRACKLE_FREE(my_rates->SN0_kpH2Oice);

  return SUCCESS;
}

int calc_rates_dust_loc(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   1.79042e-01;
  my_rates->SN0_XO [iSN] =   5.11524e-01;
  my_rates->SN0_XMg[iSN] =   3.46246e-02;
  my_rates->SN0_XAl[iSN] =   3.07922e-03;
  my_rates->SN0_XSi[iSN] =   3.76121e-02;
  my_rates->SN0_XS [iSN] =   2.21374e-02;
  my_rates->SN0_XFe[iSN] =   6.77017e-02;

  my_rates->SN0_fC [iSN] =   5.01317e-02;
  my_rates->SN0_fO [iSN] =   2.78491e-01;
  my_rates->SN0_fMg[iSN] =   0.00000e+00;
  my_rates->SN0_fAl[iSN] =   3.07922e-03;
  my_rates->SN0_fSi[iSN] =   3.50813e-03;
  my_rates->SN0_fS [iSN] =   0.00000e+00;
  my_rates->SN0_fFe[iSN] =   1.66568e-04;

  my_rates->SN0_fFeM     [iSN] =   1.35403e-02;
  my_rates->SN0_fMg2SiO4 [iSN] =   1.36165e-01;
  my_rates->SN0_fMgSiO3  [iSN] =   3.84003e-02;
  my_rates->SN0_fFeS     [iSN] =   3.04389e-02;
  my_rates->SN0_freforg  [iSN] =   1.86114e-01;
  my_rates->SN0_fvolorg  [iSN] =   3.81956e-02;
  my_rates->SN0_fH2Oice  [iSN] =   6.33011e-02;

  itab0 = 3 * iSN;
  my_rates->SN0_r0FeM     [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0FeS     [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0reforg  [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0volorg  [itab0 + 0] =   8.33039e-07;
  my_rates->SN0_r0H2Oice  [itab0 + 0] =   8.33039e-07;

  my_rates->SN0_r0FeM     [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0FeS     [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0reforg  [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0volorg  [itab0 + 1] =   1.16161e-12;
  my_rates->SN0_r0H2Oice  [itab0 + 1] =   1.16161e-12;

  my_rates->SN0_r0FeM     [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0FeS     [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0reforg  [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0volorg  [itab0 + 2] =   8.21384e-18;
  my_rates->SN0_r0H2Oice  [itab0 + 2] =   8.21384e-18;

  NTd =            35;
 Nmom =             4;

  double loc_kpFeM[] = 
  {  3.03937e-04,   1.23816e-09,   4.62094e-14,   3.72497e-18,
     5.33931e-04,   1.94666e-09,   6.42306e-14,   4.92267e-18,
     8.23085e-04,   2.82715e-09,   8.64655e-14,   6.40508e-18,
     1.18689e-03,   3.92926e-09,   1.14201e-13,   8.25547e-18,
     1.93090e-03,   5.86828e-09,   1.50433e-13,   1.02764e-17,
     2.95251e-03,   8.43657e-09,   1.95467e-13,   1.27089e-17,
     4.55033e-03,   1.22013e-08,   2.52591e-13,   1.54880e-17,
     7.04895e-03,   1.77533e-08,   3.26595e-13,   1.87358e-17,
     1.09337e-02,   2.59554e-08,   4.24054e-13,   2.26060e-17,
     1.67787e-02,   3.77909e-08,   5.51941e-13,   2.72576e-17,
     2.56790e-02,   5.50228e-08,   7.20688e-13,   3.28090e-17,
     3.89780e-02,   7.95600e-08,   9.39516e-13,   3.93285e-17,
     5.82741e-02,   1.13372e-07,   1.21633e-12,   4.68721e-17,
     8.58228e-02,   1.58873e-07,   1.55917e-12,   5.54857e-17,
     1.24511e-01,   2.18724e-07,   1.97549e-12,   6.52470e-17,
     1.78164e-01,   2.96038e-07,   2.47229e-12,   7.62244e-17,
     2.51965e-01,   3.94501e-07,   3.05495e-12,   8.84292e-17,
     3.53706e-01,   5.19228e-07,   3.72904e-12,   1.01872e-16,
     4.96690e-01,   6.79268e-07,   4.51008e-12,   1.16807e-16,
     7.06391e-01,   8.93298e-07,   5.44240e-12,   1.34078e-16,
     1.03345e+00,   1.19958e-06,   6.62440e-12,   1.55520e-16,
     1.58058e+00,   1.67481e-06,   8.24715e-12,   1.84604e-16,
     2.55750e+00,   2.47088e-06,   1.06606e-11,   2.27517e-16,
     4.39505e+00,   3.88891e-06,   1.44836e-11,   2.94492e-16,
     7.97338e+00,   6.52370e-06,   2.07697e-11,   4.00081e-16,
     1.50676e+01,   1.15419e-05,   3.12837e-11,   5.61771e-16,
     2.91750e+01,   2.11959e-05,   4.90431e-11,   7.99872e-16,
     5.69546e+01,   3.97177e-05,   7.93362e-11,   1.14279e-15,
     1.10473e+02,   7.47151e-05,   1.31340e-10,   1.63629e-15,
     2.10337e+02,   1.39130e-04,   2.20345e-10,   2.35251e-15,
     3.89661e+02,   2.53753e-04,   3.70677e-10,   3.39830e-15,
     6.99152e+02,   4.50456e-04,   6.19407e-10,   4.92457e-15,
     1.21457e+03,   7.76739e-04,   1.02057e-09,   7.13090e-15,
     2.05022e+03,   1.30302e-03,   1.64745e-09,   1.02322e-14,
     3.38793e+03,   2.13731e-03,   2.59306e-09,   1.43690e-14  };

  double loc_kpMg2SiO4[] = 
  {  2.45237e-01,   1.46287e-07,   1.48927e-13,   1.14540e-18,
     3.08964e-01,   1.84301e-07,   1.87627e-13,   1.44311e-18,
     3.89192e-01,   2.32158e-07,   2.36348e-13,   1.81790e-18,
     4.90193e-01,   2.92406e-07,   2.97684e-13,   2.28974e-18,
     6.33568e-01,   3.77932e-07,   3.84757e-13,   2.95995e-18,
     8.28859e-01,   4.94426e-07,   5.03359e-13,   3.87306e-18,
     1.13002e+00,   6.74071e-07,   6.86261e-13,   5.28188e-18,
     1.63064e+00,   9.72699e-07,   9.90315e-13,   7.62608e-18,
     2.46680e+00,   1.47148e-06,   1.49820e-12,   1.15475e-17,
     3.79604e+00,   2.26439e-06,   2.30565e-12,   1.77932e-17,
     5.92492e+00,   3.53430e-06,   3.59906e-12,   2.78297e-17,
     9.23920e+00,   5.51134e-06,   5.61323e-12,   4.35401e-17,
     1.42293e+01,   8.48802e-06,   8.64701e-12,   6.73899e-17,
     2.16427e+01,   1.29103e-05,   1.31570e-11,   1.03293e-16,
     3.24523e+01,   1.93586e-05,   1.97391e-11,   1.56579e-16,
     4.78595e+01,   2.85497e-05,   2.91302e-11,   2.33879e-16,
     7.00593e+01,   4.17932e-05,   4.26757e-11,   3.46486e-16,
     1.06051e+02,   6.32653e-05,   6.46560e-11,   5.28640e-16,
     1.74267e+02,   1.03964e-04,   1.06329e-10,   8.68729e-16,
     3.02053e+02,   1.80203e-04,   1.84359e-10,   1.49313e-15,
     5.00593e+02,   2.98655e-04,   3.05537e-10,   2.44653e-15,
     7.45698e+02,   4.44894e-04,   4.55158e-10,   3.61190e-15,
     1.00149e+03,   5.97517e-04,   6.11446e-10,   4.82507e-15,
     1.23701e+03,   7.38051e-04,   7.55475e-10,   5.93849e-15,
     1.39749e+03,   8.33819e-04,   8.53626e-10,   6.68589e-15,
     1.41344e+03,   8.43342e-04,   8.63374e-10,   6.74095e-15,
     1.26599e+03,   7.55370e-04,   7.73283e-10,   6.02512e-15,
     1.01032e+03,   6.02819e-04,   6.17127e-10,   4.80541e-15,
     7.30148e+02,   4.35656e-04,   4.46054e-10,   3.47626e-15,
     4.87070e+02,   2.90623e-04,   2.97653e-10,   2.32554e-15,
     3.05625e+02,   1.82365e-04,   1.86919e-10,   1.46881e-15,
     1.84266e+02,   1.09966e-04,   1.12993e-10,   9.02608e-16,
     1.10340e+02,   6.58826e-05,   6.82482e-11,   5.72358e-16,
     6.91207e+01,   4.13300e-05,   4.36913e-11,   4.08010e-16,
     4.77899e+01,   2.86629e-05,   3.13943e-11,   3.40524e-16  };

  double loc_kpMgSiO3[] = 
  {  5.12401e-02,   3.05654e-08,   3.11170e-14,   2.39327e-19,
     9.10229e-02,   5.42964e-08,   5.52765e-14,   4.25193e-19,
     1.41106e-01,   8.41719e-08,   8.56916e-14,   6.59185e-19,
     2.04158e-01,   1.21783e-07,   1.23982e-13,   9.53762e-19,
     3.33897e-01,   1.99174e-07,   2.02773e-13,   1.56024e-18,
     5.10947e-01,   3.04787e-07,   3.10298e-13,   2.38807e-18,
     7.83563e-01,   4.67406e-07,   4.75863e-13,   3.66337e-18,
     1.19854e+00,   7.14943e-07,   7.27902e-13,   5.60698e-18,
     1.85772e+00,   1.10816e-06,   1.12831e-12,   8.70099e-18,
     2.92247e+00,   1.74330e-06,   1.77515e-12,   1.37124e-17,
     4.74091e+00,   2.82803e-06,   2.88010e-12,   2.23096e-17,
     7.79816e+00,   4.65174e-06,   4.73847e-12,   3.68720e-17,
     1.27207e+01,   7.58818e-06,   7.73243e-12,   6.06033e-17,
     2.05556e+01,   1.22619e-05,   1.25020e-11,   9.90816e-17,
     3.30490e+01,   1.97147e-05,   2.01153e-11,   1.61650e-16,
     5.32277e+01,   3.17524e-05,   3.24225e-11,   2.63830e-16,
     8.64944e+01,   5.15982e-05,   5.27236e-11,   4.32162e-16,
     1.43144e+02,   8.53944e-05,   8.73076e-11,   7.16006e-16,
     2.41992e+02,   1.44367e-04,   1.47663e-10,   1.20379e-15,
     4.08980e+02,   2.43996e-04,   2.49602e-10,   2.01287e-15,
     6.57282e+02,   3.92138e-04,   4.01100e-10,   3.19447e-15,
     9.65951e+02,   5.76298e-04,   5.89371e-10,   4.64306e-15,
     1.30489e+03,   7.78527e-04,   7.96032e-10,   6.21292e-15,
     1.65692e+03,   9.88559e-04,   1.01028e-09,   7.80402e-15,
     1.95936e+03,   1.16899e-03,   1.19360e-09,   9.11030e-15,
     2.08639e+03,   1.24475e-03,   1.26965e-09,   9.58003e-15,
     1.95845e+03,   1.16840e-03,   1.19068e-09,   8.90255e-15,
     1.62353e+03,   9.68572e-04,   9.86379e-10,   7.32855e-15,
     1.20759e+03,   7.20421e-04,   7.33336e-10,   5.42763e-15,
     8.22295e+02,   4.90559e-04,   4.99228e-10,   3.68795e-15,
     5.22496e+02,   3.11708e-04,   3.17198e-10,   2.34267e-15,
     3.14915e+02,   1.87872e-04,   1.91211e-10,   1.41447e-15,
     1.82496e+02,   1.08876e-04,   1.10875e-10,   8.23864e-16,
     1.02901e+02,   6.13959e-05,   6.26109e-11,   4.69856e-16,
     5.73186e+01,   3.42174e-05,   3.50508e-11,   2.68356e-16  };

  double loc_kpFeS[] = 
  {  1.20726e-01,   7.20182e-08,   7.35280e-14,   5.95267e-19,
     2.32766e-01,   1.38854e-07,   1.41735e-13,   1.14369e-18,
     3.73816e-01,   2.22995e-07,   2.27602e-13,   1.83408e-18,
     5.51388e-01,   3.28923e-07,   3.35702e-13,   2.70320e-18,
     8.55868e-01,   5.10557e-07,   5.21179e-13,   4.21168e-18,
     1.24954e+00,   7.45400e-07,   7.61051e-13,   6.17189e-18,
     1.78078e+00,   1.06231e-06,   1.08497e-12,   8.85200e-18,
     2.44607e+00,   1.45920e-06,   1.49119e-12,   1.23029e-17,
     3.21763e+00,   1.91949e-06,   1.96357e-12,   1.65254e-17,
     4.06334e+00,   2.42404e-06,   2.48333e-12,   2.14914e-17,
     4.90080e+00,   2.92371e-06,   3.00112e-12,   2.69365e-17,
     5.64261e+00,   3.36638e-06,   3.46359e-12,   3.23217e-17,
     6.21740e+00,   3.70948e-06,   3.82669e-12,   3.70749e-17,
     6.56246e+00,   3.91565e-06,   4.05361e-12,   4.08425e-17,
     6.75396e+00,   4.03063e-06,   4.20072e-12,   4.47151e-17,
     7.17709e+00,   4.28498e-06,   4.53461e-12,   5.31551e-17,
     8.60795e+00,   5.14313e-06,   5.58433e-12,   7.44632e-17,
     1.17764e+01,   7.04212e-06,   7.86423e-12,   1.17214e-16,
     1.64787e+01,   9.86444e-06,   1.14471e-11,   1.95038e-16,
     2.14143e+01,   1.28476e-05,   1.62848e-11,   3.61550e-16,
     2.51760e+01,   1.51720e-05,   2.28477e-11,   7.25431e-16,
     2.71833e+01,   1.64950e-05,   3.14267e-11,   1.35273e-15,
     2.77234e+01,   1.69657e-05,   4.09753e-11,   2.14916e-15,
     2.74221e+01,   1.69259e-05,   4.94400e-11,   2.89997e-15,
     2.68227e+01,   1.66816e-05,   5.53227e-11,   3.43309e-15,
     2.62247e+01,   1.64266e-05,   5.86985e-11,   3.72279e-15,
     2.58704e+01,   1.64543e-05,   6.22578e-11,   3.90038e-15,
     2.71342e+01,   1.84130e-05,   7.74409e-11,   4.33089e-15,
     3.69065e+01,   2.84578e-05,   1.28908e-10,   5.42457e-15,
     7.78838e+01,   6.43091e-05,   2.58063e-10,   7.39401e-15,
     1.96086e+02,   1.59699e-04,   5.23332e-10,   1.03910e-14,
     4.47749e+02,   3.52292e-04,   9.69286e-10,   1.44002e-14,
     8.75250e+02,   6.66973e-04,   1.60576e-09,   1.92445e-14,
     1.50099e+03,   1.11385e-03,   2.41814e-09,   2.47015e-14,
     2.33693e+03,   1.69566e-03,   3.38124e-09,   3.05222e-14  };

  double loc_kpreforg[] = 
  {  4.68555e-02,   2.79499e-08,   2.84543e-14,   2.18837e-19,
     6.66016e-02,   3.97287e-08,   4.04457e-14,   3.11065e-19,
     9.14604e-02,   5.45574e-08,   5.55420e-14,   4.27173e-19,
     1.22756e-01,   7.32255e-08,   7.45471e-14,   5.73344e-19,
     2.39395e-01,   1.42802e-07,   1.45382e-13,   1.11857e-18,
     4.33638e-01,   2.58671e-07,   2.63349e-13,   2.02683e-18,
     8.54269e-01,   5.09583e-07,   5.18807e-13,   3.99441e-18,
     1.71749e+00,   1.02450e-06,   1.04308e-12,   8.03545e-18,
     3.29447e+00,   1.96520e-06,   2.00091e-12,   1.54275e-17,
     5.86107e+00,   3.49621e-06,   3.55997e-12,   2.74796e-17,
     9.98591e+00,   5.95675e-06,   6.06595e-12,   4.69112e-17,
     1.64428e+01,   9.80839e-06,   9.98976e-12,   7.74960e-17,
     2.63091e+01,   1.56939e-05,   1.59880e-11,   1.24634e-16,
     4.13343e+01,   2.46568e-05,   2.51287e-11,   1.97408e-16,
     6.41346e+01,   3.82579e-05,   3.90120e-11,   3.09791e-16,
     9.87636e+01,   5.89156e-05,   6.01167e-11,   4.83037e-16,
     1.51958e+02,   9.06488e-05,   9.25564e-11,   7.50644e-16,
     2.34417e+02,   1.39841e-04,   1.42858e-10,   1.16365e-15,
     3.57697e+02,   2.13387e-04,   2.18076e-10,   1.77561e-15,
     5.28819e+02,   3.15479e-04,   3.22528e-10,   2.61825e-15,
     7.52235e+02,   4.48774e-04,   4.59023e-10,   3.71332e-15,
     1.02632e+03,   6.12308e-04,   6.26695e-10,   5.05587e-15,
     1.32111e+03,   7.88218e-04,   8.07423e-10,   6.50989e-15,
     1.57574e+03,   9.40185e-04,   9.64291e-10,   7.80218e-15,
     1.73868e+03,   1.03748e-03,   1.06621e-09,   8.71453e-15,
     1.80162e+03,   1.07516e-03,   1.10855e-09,   9.23861e-15,
     1.79360e+03,   1.07058e-03,   1.10965e-09,   9.55041e-15,
     1.76388e+03,   1.05316e-03,   1.10040e-09,   9.92949e-15,
     1.76147e+03,   1.05218e-03,   1.11104e-09,   1.06225e-14,
     1.80062e+03,   1.07614e-03,   1.14939e-09,   1.16303e-14,
     1.85114e+03,   1.10702e-03,   1.19518e-09,   1.26854e-14,
     1.87614e+03,   1.12283e-03,   1.22470e-09,   1.35171e-14,
     1.86911e+03,   1.11985e-03,   1.23456e-09,   1.40796e-14,
     1.85783e+03,   1.11511e-03,   1.24487e-09,   1.45610e-14,
     1.90737e+03,   1.14893e-03,   1.30468e-09,   1.53873e-14  };

  double loc_kpvolorg[] = 
  {  7.02832e-02,   4.19249e-08,   4.26814e-14,   3.28255e-19,
     9.99024e-02,   5.95931e-08,   6.06685e-14,   4.66597e-19,
     1.37191e-01,   8.18361e-08,   8.33130e-14,   6.40759e-19,
     1.84134e-01,   1.09838e-07,   1.11821e-13,   8.60016e-19,
     3.59093e-01,   2.14204e-07,   2.18074e-13,   1.67786e-18,
     6.50459e-01,   3.88007e-07,   3.95023e-13,   3.04025e-18,
     1.28140e+00,   7.64376e-07,   7.78212e-13,   5.99162e-18,
     2.57622e+00,   1.53675e-06,   1.56462e-12,   1.20532e-17,
     4.94170e+00,   2.94779e-06,   3.00137e-12,   2.31413e-17,
     8.79162e+00,   5.24433e-06,   5.33996e-12,   4.12194e-17,
     1.49789e+01,   8.93514e-06,   9.09893e-12,   7.03669e-17,
     2.46641e+01,   1.47126e-05,   1.49846e-11,   1.16244e-16,
     3.94637e+01,   2.35409e-05,   2.39820e-11,   1.86951e-16,
     6.20013e+01,   3.69851e-05,   3.76930e-11,   2.96113e-16,
     9.62020e+01,   5.73870e-05,   5.85180e-11,   4.64688e-16,
     1.48145e+02,   8.83735e-05,   9.01751e-11,   7.24557e-16,
     2.27937e+02,   1.35973e-04,   1.38835e-10,   1.12597e-15,
     3.51626e+02,   2.09762e-04,   2.14288e-10,   1.74548e-15,
     5.36545e+02,   3.20081e-04,   3.27114e-10,   2.66341e-15,
     7.93228e+02,   4.73218e-04,   4.83793e-10,   3.92738e-15,
     1.12835e+03,   6.73161e-04,   6.88534e-10,   5.56999e-15,
     1.53947e+03,   9.18463e-04,   9.40042e-10,   7.58380e-15,
     1.98167e+03,   1.18233e-03,   1.21114e-09,   9.76484e-15,
     2.36362e+03,   1.41028e-03,   1.44644e-09,   1.17033e-14,
     2.60802e+03,   1.55622e-03,   1.59931e-09,   1.30718e-14,
     2.70242e+03,   1.61274e-03,   1.66283e-09,   1.38579e-14,
     2.69039e+03,   1.60586e-03,   1.66448e-09,   1.43256e-14,
     2.64583e+03,   1.57974e-03,   1.65060e-09,   1.48942e-14,
     2.64222e+03,   1.57827e-03,   1.66657e-09,   1.59338e-14,
     2.70093e+03,   1.61421e-03,   1.72408e-09,   1.74454e-14,
     2.77670e+03,   1.66053e-03,   1.79277e-09,   1.90281e-14,
     2.81420e+03,   1.68424e-03,   1.83704e-09,   2.02756e-14,
     2.80366e+03,   1.67977e-03,   1.85184e-09,   2.11194e-14,
     2.78675e+03,   1.67267e-03,   1.86731e-09,   2.18416e-14,
     2.86105e+03,   1.72339e-03,   1.95702e-09,   2.30809e-14  };

  double loc_kpH2Oice[] = 
  {  6.30862e-02,   3.76318e-08,   3.83108e-14,   2.94648e-19,
     1.09691e-01,   6.54321e-08,   6.66131e-14,   5.12363e-19,
     1.68363e-01,   1.00431e-07,   1.02243e-13,   7.86449e-19,
     2.42226e-01,   1.44491e-07,   1.47099e-13,   1.13150e-18,
     3.92991e-01,   2.34424e-07,   2.38657e-13,   1.83595e-18,
     6.03026e-01,   3.59713e-07,   3.66210e-13,   2.81745e-18,
     9.41509e-01,   5.61623e-07,   5.71771e-13,   4.39959e-18,
     1.54206e+00,   9.19862e-07,   9.36500e-13,   7.20849e-18,
     2.75198e+00,   1.64160e-06,   1.67135e-12,   1.28735e-17,
     5.12965e+00,   3.05991e-06,   3.11550e-12,   2.40170e-17,
     9.67703e+00,   5.77249e-06,   5.87767e-12,   4.53565e-17,
     1.70842e+01,   1.01910e-05,   1.03775e-11,   8.02010e-17,
     2.87858e+01,   1.71712e-05,   1.74882e-11,   1.35552e-16,
     5.22985e+01,   3.11971e-05,   3.17818e-11,   2.47506e-16,
     1.08421e+02,   6.46756e-05,   6.59042e-11,   5.15083e-16,
     2.29694e+02,   1.37018e-04,   1.39634e-10,   1.09171e-15,
     4.35888e+02,   2.60019e-04,   2.64976e-10,   2.06803e-15,
     6.94501e+02,   4.14288e-04,   4.22159e-10,   3.28761e-15,
     9.21763e+02,   5.49857e-04,   5.60297e-10,   4.35726e-15,
     1.06208e+03,   6.33565e-04,   6.45702e-10,   5.02330e-15,
     1.17460e+03,   7.00699e-04,   7.14492e-10,   5.57174e-15,
     1.41297e+03,   8.42922e-04,   8.60130e-10,   6.71698e-15,
     1.86068e+03,   1.11004e-03,   1.13311e-09,   8.82123e-15,
     2.38031e+03,   1.42005e-03,   1.44954e-09,   1.12234e-14,
     2.71131e+03,   1.61753e-03,   1.65107e-09,   1.27373e-14,
     2.72120e+03,   1.62346e-03,   1.65757e-09,   1.27900e-14,
     2.53189e+03,   1.51058e-03,   1.54284e-09,   1.19191e-14,
     2.48249e+03,   1.48118e-03,   1.51103e-09,   1.15227e-14,
     2.89111e+03,   1.72503e-03,   1.75260e-09,   1.28501e-14,
     3.66301e+03,   2.18553e-03,   2.21015e-09,   1.55022e-14,
     4.27842e+03,   2.55260e-03,   2.57302e-09,   1.74959e-14,
     4.29745e+03,   2.56387e-03,   2.57944e-09,   1.72213e-14,
     3.71668e+03,   2.21732e-03,   2.22839e-09,   1.47266e-14,
     2.83563e+03,   1.69168e-03,   1.69914e-09,   1.11685e-14,
     1.96048e+03,   1.16958e-03,   1.17438e-09,   7.69912e-15  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpFeM     [itab0] = loc_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = loc_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = loc_kpMgSiO3  [itab];
      my_rates->SN0_kpFeS     [itab0] = loc_kpFeS     [itab];
      my_rates->SN0_kpreforg  [itab0] = loc_kpreforg  [itab];
      my_rates->SN0_kpvolorg  [itab0] = loc_kpvolorg  [itab];
      my_rates->SN0_kpH2Oice  [itab0] = loc_kpH2Oice  [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_C13(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   2.65314e-01;
  my_rates->SN0_XO [iSN] =   3.00982e-01;
  my_rates->SN0_XMg[iSN] =   3.06651e-02;
  my_rates->SN0_XAl[iSN] =   2.47296e-04;
  my_rates->SN0_XSi[iSN] =   6.38319e-02;
  my_rates->SN0_XS [iSN] =   3.40910e-02;
  my_rates->SN0_XFe[iSN] =   9.62448e-02;

  my_rates->SN0_fC [iSN] =   2.16731e-01;
  my_rates->SN0_fO [iSN] =   2.99231e-01;
  my_rates->SN0_fMg[iSN] =   3.03586e-02;
  my_rates->SN0_fAl[iSN] =   2.47296e-04;
  my_rates->SN0_fSi[iSN] =   4.59041e-02;
  my_rates->SN0_fS [iSN] =   3.40903e-02;
  my_rates->SN0_fFe[iSN] =   7.22586e-02;

  my_rates->SN0_fSiM     [iSN] =   1.65746e-02;
  my_rates->SN0_fFeM     [iSN] =   2.39849e-02;
  my_rates->SN0_fMg2SiO4 [iSN] =   8.69522e-04;
  my_rates->SN0_fMgSiO3  [iSN] =   2.87802e-06;
  my_rates->SN0_fAC      [iSN] =   4.85826e-02;
  my_rates->SN0_fSiO2D   [iSN] =   2.52534e-03;
  my_rates->SN0_fMgO     [iSN] =   1.28672e-05;
  my_rates->SN0_fFeS     [iSN] =   2.09730e-06;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   1.68557e-06;
  my_rates->SN0_r0FeM     [itab0 + 0] =   4.62542e-06;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   1.82163e-06;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   7.26303e-07;
  my_rates->SN0_r0AC      [itab0 + 0] =   4.82296e-06;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   1.33530e-06;
  my_rates->SN0_r0MgO     [itab0 + 0] =   1.59029e-06;
  my_rates->SN0_r0FeS     [itab0 + 0] =   6.16010e-07;

  my_rates->SN0_r0SiM     [itab0 + 1] =   9.75226e-12;
  my_rates->SN0_r0FeM     [itab0 + 1] =   3.82292e-11;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   5.83823e-12;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   7.49856e-13;
  my_rates->SN0_r0AC      [itab0 + 1] =   3.91353e-11;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   5.91862e-12;
  my_rates->SN0_r0MgO     [itab0 + 1] =   7.21459e-12;
  my_rates->SN0_r0FeS     [itab0 + 1] =   4.56500e-13;

  my_rates->SN0_r0SiM     [itab0 + 2] =   1.74046e-16;
  my_rates->SN0_r0FeM     [itab0 + 2] =   4.68445e-16;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   3.61356e-17;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   1.57511e-18;
  my_rates->SN0_r0AC      [itab0 + 2] =   5.15140e-16;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   5.31739e-17;
  my_rates->SN0_r0MgO     [itab0 + 2] =   4.84120e-17;
  my_rates->SN0_r0FeS     [itab0 + 2] =   4.16699e-19;

  NTd =            35;
 Nmom =             4;

  double C13_kpSiM[] = 
  {  1.54619e-01,   2.60128e-07,   1.49475e-12,   2.65148e-17,
     1.94656e-01,   3.27554e-07,   1.88361e-12,   3.34351e-17,
     2.45059e-01,   4.12439e-07,   2.37316e-12,   4.21472e-17,
     3.08513e-01,   5.19301e-07,   2.98947e-12,   5.31150e-17,
     3.88404e-01,   6.53894e-07,   3.76675e-12,   6.69643e-17,
     4.88982e-01,   8.23350e-07,   4.74559e-12,   8.44087e-17,
     6.15605e-01,   1.03673e-06,   5.97899e-12,   1.06403e-16,
     7.75006e-01,   1.30539e-06,   7.53301e-12,   1.34131e-16,
     9.75366e-01,   1.64315e-06,   9.48791e-12,   1.69031e-16,
     1.22493e+00,   2.06392e-06,   1.19246e-11,   2.12553e-16,
     1.52116e+00,   2.56344e-06,   1.48192e-11,   2.64285e-16,
     1.83681e+00,   3.09586e-06,   1.79072e-11,   3.19511e-16,
     2.15662e+00,   3.63549e-06,   2.10413e-11,   3.75635e-16,
     2.55502e+00,   4.30820e-06,   2.49580e-11,   4.45921e-16,
     3.22790e+00,   5.44497e-06,   3.15890e-11,   5.65115e-16,
     4.33126e+00,   7.30995e-06,   4.24882e-11,   7.61352e-16,
     5.81498e+00,   9.82030e-06,   5.72103e-11,   1.02723e-15,
     7.48294e+00,   1.26484e-05,   7.39233e-11,   1.33106e-15,
     9.21324e+00,   1.55973e-05,   9.16644e-11,   1.65851e-15,
     1.11943e+01,   1.90119e-05,   1.13021e-10,   2.06536e-15,
     1.39990e+01,   2.39150e-05,   1.45136e-10,   2.69938e-15,
     1.78867e+01,   3.08010e-05,   1.92146e-10,   3.65687e-15,
     2.17798e+01,   3.79266e-05,   2.45744e-10,   4.82338e-15,
     2.38104e+01,   4.23268e-05,   2.93339e-10,   6.06677e-15,
     2.31200e+01,   4.26489e-05,   3.30423e-10,   7.40538e-15,
     2.05061e+01,   3.98990e-05,   3.56330e-10,   8.76054e-15,
     1.74011e+01,   3.59849e-05,   3.69684e-10,   9.86673e-15,
     1.47227e+01,   3.22080e-05,   3.69774e-10,   1.04731e-14,
     1.27456e+01,   2.91368e-05,   3.59951e-10,   1.05598e-14,
     1.18944e+01,   2.82426e-05,   3.62521e-10,   1.07315e-14,
     1.43878e+01,   3.74123e-05,   4.87177e-10,   1.37994e-14,
     2.85377e+01,   9.23385e-05,   1.20077e-09,   3.00478e-14,
     7.93256e+01,   3.11569e-04,   3.91124e-09,   8.56317e-14,
     2.33021e+02,   9.78251e-04,   1.15601e-08,   2.28453e-13,
     6.42052e+02,   2.55294e-03,   2.77058e-08,   5.06000e-13  };

  double C13_kpFeM[] = 
  {  5.88921e-03,   6.87770e-08,   1.06227e-12,   2.02349e-17,
     1.01153e-02,   1.15444e-07,   1.74157e-12,   3.24501e-17,
     1.53919e-02,   1.73359e-07,   2.57952e-12,   4.74426e-17,
     2.20096e-02,   2.45791e-07,   3.62468e-12,   6.61001e-17,
     3.50047e-02,   3.82752e-07,   5.51942e-12,   9.85221e-17,
     5.25286e-02,   5.65069e-07,   8.00632e-12,   1.40508e-16,
     7.91090e-02,   8.36026e-07,   1.16179e-11,   2.00078e-16,
     1.19191e-01,   1.23677e-06,   1.68430e-11,   2.84365e-16,
     1.79251e-01,   1.82633e-06,   2.43704e-11,   4.03226e-16,
     2.66584e-01,   2.66955e-06,   3.49349e-11,   5.66851e-16,
     3.94058e-01,   3.87941e-06,   4.97989e-11,   7.92464e-16,
     5.74681e-01,   5.56462e-06,   7.01108e-11,   1.09474e-15,
     8.20409e-01,   7.81996e-06,   9.68130e-11,   1.48491e-15,
     1.14378e+00,   1.07400e-05,   1.30804e-10,   1.97313e-15,
     1.55566e+00,   1.44006e-05,   1.72756e-10,   2.56649e-15,
     2.06561e+00,   1.88620e-05,   2.23158e-10,   3.26978e-15,
     2.68069e+00,   2.41572e-05,   2.82187e-10,   4.08370e-15,
     3.40683e+00,   3.03006e-05,   3.49792e-10,   5.00618e-15,
     4.25797e+00,   3.73606e-05,   4.26499e-10,   6.04358e-15,
     5.27551e+00,   4.56092e-05,   5.14976e-10,   7.23179e-15,
     6.55438e+00,   5.57066e-05,   6.21908e-10,   8.66078e-15,
     8.27812e+00,   6.89310e-05,   7.60263e-10,   1.05052e-14,
     1.07753e+01,   8.75306e-05,   9.52770e-10,   1.30714e-14,
     1.46156e+01,   1.15279e-04,   1.23703e-09,   1.68635e-14,
     2.07765e+01,   1.58349e-04,   1.67242e-09,   2.26556e-14,
     3.09403e+01,   2.26776e-04,   2.35074e-09,   3.15878e-14,
     4.80249e+01,   3.37152e-04,   3.41710e-09,   4.53799e-14,
     7.70965e+01,   5.17359e-04,   5.10858e-09,   6.67630e-14,
     1.26831e+02,   8.14201e-04,   7.81618e-09,   1.00160e-13,
     2.11638e+02,   1.30438e-03,   1.21704e-08,   1.52564e-13,
     3.54644e+02,   2.11006e-03,   1.91598e-08,   2.34688e-13,
     5.91556e+02,   3.41851e-03,   3.02724e-08,   3.62228e-13,
     9.74370e+02,   5.49677e-03,   4.75739e-08,   5.56235e-13,
     1.56930e+03,   8.65788e-03,   7.33192e-08,   8.38100e-13,
     2.44406e+03,   1.31389e-02,   1.08763e-07,   1.21599e-12  };

  double C13_kpMg2SiO4[] = 
  {  1.05240e-01,   1.91709e-07,   6.14415e-13,   3.80291e-18,
     1.32588e-01,   2.41526e-07,   7.74078e-13,   4.79114e-18,
     1.67016e-01,   3.04243e-07,   9.75080e-13,   6.03524e-18,
     2.10360e-01,   3.83198e-07,   1.22813e-12,   7.60148e-18,
     2.71887e-01,   4.95279e-07,   1.58734e-12,   9.82485e-18,
     3.55694e-01,   6.47944e-07,   2.07663e-12,   1.28533e-17,
     4.84932e-01,   8.83369e-07,   2.83116e-12,   1.75235e-17,
     6.99767e-01,   1.27472e-06,   4.08543e-12,   2.52870e-17,
     1.05860e+00,   1.92838e-06,   6.18042e-12,   3.82543e-17,
     1.62902e+00,   2.96748e-06,   9.51075e-12,   5.88683e-17,
     2.54260e+00,   4.63171e-06,   1.48447e-11,   9.18851e-17,
     3.96490e+00,   7.22268e-06,   2.31492e-11,   1.43293e-16,
     6.10635e+00,   1.11237e-05,   3.56530e-11,   2.20699e-16,
     9.28776e+00,   1.69193e-05,   5.42303e-11,   3.35716e-16,
     1.39267e+01,   2.53704e-05,   8.13220e-11,   5.03475e-16,
     2.05388e+01,   3.74163e-05,   1.19943e-10,   7.42685e-16,
     3.00662e+01,   5.47749e-05,   1.75612e-10,   1.08767e-15,
     4.55134e+01,   8.29222e-05,   2.65919e-10,   1.64778e-15,
     7.47928e+01,   1.36281e-04,   4.37190e-10,   2.71088e-15,
     1.29641e+02,   2.36241e-04,   7.58090e-10,   4.70339e-15,
     2.14861e+02,   3.91558e-04,   1.25676e-09,   7.80036e-15,
     3.20074e+02,   5.83333e-04,   1.87266e-09,   1.16272e-14,
     4.29885e+02,   7.83522e-04,   2.51590e-09,   1.56273e-14,
     5.31008e+02,   9.67906e-04,   3.10866e-09,   1.93168e-14,
     5.99926e+02,   1.09359e-03,   3.51289e-09,   2.18342e-14,
     6.06790e+02,   1.10614e-03,   3.55349e-09,   2.20891e-14,
     5.43501e+02,   9.90793e-04,   3.18304e-09,   1.97873e-14,
     4.33743e+02,   7.90717e-04,   2.54035e-09,   1.57927e-14,
     3.13466e+02,   5.71463e-04,   1.83604e-09,   1.14151e-14,
     2.09115e+02,   3.81238e-04,   1.22500e-09,   7.61738e-15,
     1.31224e+02,   2.39260e-04,   7.69004e-10,   4.78403e-15,
     7.91470e+01,   1.44373e-04,   4.64502e-10,   2.89404e-15,
     4.74663e+01,   8.67304e-05,   2.80053e-10,   1.75354e-15,
     2.98668e+01,   5.48299e-05,   1.78713e-10,   1.13252e-15,
     2.08636e+01,   3.86952e-05,   1.28439e-10,   8.30592e-16  };

  double C13_kpMgSiO3[] = 
  {  2.19890e-02,   1.59707e-08,   1.64886e-14,   3.46350e-20,
     3.90612e-02,   2.83703e-08,   2.92903e-14,   6.15261e-20,
     6.05539e-02,   4.39805e-08,   4.54068e-14,   9.53799e-20,
     8.76116e-02,   6.36326e-08,   6.56961e-14,   1.37999e-19,
     1.43288e-01,   1.04070e-07,   1.07445e-13,   2.25699e-19,
     2.19266e-01,   1.59254e-07,   1.64418e-13,   3.45380e-19,
     3.36256e-01,   2.44223e-07,   2.52144e-13,   5.29662e-19,
     5.14336e-01,   3.73564e-07,   3.85679e-13,   8.10191e-19,
     7.97216e-01,   5.79021e-07,   5.97800e-13,   1.25585e-18,
     1.25414e+00,   9.10886e-07,   9.40430e-13,   1.97579e-18,
     2.03450e+00,   1.47766e-06,   1.52560e-12,   3.20555e-18,
     3.34648e+00,   2.43056e-06,   2.50942e-12,   5.27367e-18,
     5.45894e+00,   3.96485e-06,   4.09354e-12,   8.60504e-18,
     8.82118e+00,   6.40687e-06,   6.61493e-12,   1.39111e-17,
     1.41825e+01,   1.03009e-05,   1.06356e-11,   2.23816e-17,
     2.28420e+01,   1.65903e-05,   1.71302e-11,   3.60870e-17,
     3.71180e+01,   2.69593e-05,   2.78383e-11,   5.87421e-17,
     6.14285e+01,   4.46166e-05,   4.60756e-11,   9.74574e-17,
     1.03848e+02,   7.54274e-05,   7.79023e-11,   1.65248e-16,
     1.75510e+02,   1.27478e-04,   1.31675e-10,   2.80026e-16,
     2.82066e+02,   2.04875e-04,   2.11635e-10,   4.50867e-16,
     4.14529e+02,   3.01090e-04,   3.11041e-10,   6.63478e-16,
     5.59986e+02,   4.06746e-04,   4.20208e-10,   8.97036e-16,
     7.11059e+02,   5.16484e-04,   5.33582e-10,   1.13854e-15,
     8.40851e+02,   6.10764e-04,   6.30953e-10,   1.34372e-15,
     8.95368e+02,   6.50365e-04,   6.71811e-10,   1.42700e-15,
     8.40461e+02,   6.10482e-04,   6.30563e-10,   1.33611e-15,
     6.96732e+02,   5.06082e-04,   5.22696e-10,   1.10542e-15,
     5.18234e+02,   3.76427e-04,   3.88766e-10,   8.21088e-16,
     3.52885e+02,   2.56323e-04,   2.64719e-10,   5.58648e-16,
     2.24228e+02,   1.62872e-04,   1.68207e-10,   3.54873e-16,
     1.35145e+02,   9.81654e-05,   1.01382e-10,   2.13970e-16,
     7.83182e+01,   5.68892e-05,   5.87587e-11,   1.24192e-16,
     4.41610e+01,   3.20799e-05,   3.31421e-11,   7.03026e-17,
     2.46026e+01,   1.78814e-05,   1.84973e-11,   3.95672e-17  };

  double C13_kpAC[] = 
  {  3.27960e-01,   1.58173e-06,   1.28346e-11,   1.68940e-16,
     4.38754e-01,   2.11613e-06,   1.71717e-11,   2.26045e-16,
     5.78236e-01,   2.78890e-06,   2.26319e-11,   2.97935e-16,
     7.53833e-01,   3.63586e-06,   2.95057e-11,   3.88440e-16,
     1.04018e+00,   5.01714e-06,   4.07179e-11,   5.36095e-16,
     1.41744e+00,   6.83702e-06,   5.54910e-11,   7.30661e-16,
     1.95305e+00,   9.42082e-06,   7.64677e-11,   1.00698e-15,
     2.71551e+00,   1.30993e-05,   1.06337e-10,   1.40054e-15,
     3.79716e+00,   1.83183e-05,   1.48729e-10,   1.95931e-15,
     5.29823e+00,   2.55621e-05,   2.07584e-10,   2.73540e-15,
     7.37977e+00,   3.56090e-05,   2.89250e-10,   3.81296e-15,
     1.02196e+01,   4.93202e-05,   4.00773e-10,   5.28580e-15,
     1.40471e+01,   6.78071e-05,   5.51275e-10,   7.27585e-15,
     1.92118e+01,   9.27674e-05,   7.54736e-10,   9.97084e-15,
     2.61798e+01,   1.26469e-04,   1.02994e-09,   1.36251e-14,
     3.55647e+01,   1.71912e-04,   1.40197e-09,   1.85822e-14,
     4.82256e+01,   2.33320e-04,   1.90655e-09,   2.53388e-14,
     6.54391e+01,   3.17003e-04,   2.59771e-09,   3.46585e-14,
     8.90003e+01,   4.31918e-04,   3.55359e-09,   4.76706e-14,
     1.21150e+02,   5.89447e-04,   4.87725e-09,   6.59296e-14,
     1.64482e+02,   8.03274e-04,   6.70152e-09,   9.15912e-14,
     2.22238e+02,   1.09138e-03,   9.21653e-09,   1.27994e-13,
     2.99304e+02,   1.48143e-03,   1.27247e-08,   1.80612e-13,
     4.02975e+02,   2.01378e-03,   1.76517e-08,   2.56967e-13,
     5.41945e+02,   2.73383e-03,   2.44284e-08,   3.63965e-13,
     7.24087e+02,   3.67841e-03,   3.33149e-08,   5.04207e-13,
     9.57733e+02,   4.88178e-03,   4.44353e-08,   6.76055e-13,
     1.25789e+03,   6.40887e-03,   5.80974e-08,   8.78814e-13,
     1.65268e+03,   8.39077e-03,   7.51368e-08,   1.11835e-12,
     2.18798e+03,   1.10499e-02,   9.71231e-08,   1.40975e-12,
     2.93400e+03,   1.47364e-02,   1.26617e-07,   1.77953e-12,
     3.99924e+03,   2.00041e-02,   1.67709e-07,   2.27029e-12,
     5.55224e+03,   2.77223e-02,   2.26777e-07,   2.94699e-12,
     7.83754e+03,   3.91330e-02,   3.12657e-07,   3.89444e-12,
     1.11477e+04,   5.56050e-02,   4.34205e-07,   5.18557e-12  };

  double C13_kpSiO2D[] = 
  {  7.60358e-02,   1.01529e-07,   4.49989e-13,   4.04236e-18,
     9.07205e-02,   1.21137e-07,   5.36902e-13,   4.82320e-18,
     1.09207e-01,   1.45823e-07,   6.46320e-13,   5.80621e-18,
     1.32481e-01,   1.76900e-07,   7.84068e-13,   7.04376e-18,
     1.58907e-01,   2.12188e-07,   9.40480e-13,   8.44905e-18,
     1.91565e-01,   2.55795e-07,   1.13377e-12,   1.01857e-17,
     2.30490e-01,   3.07773e-07,   1.36416e-12,   1.22557e-17,
     2.76795e-01,   3.69605e-07,   1.63824e-12,   1.47181e-17,
     3.33074e-01,   4.44757e-07,   1.97136e-12,   1.77112e-17,
     4.05326e-01,   5.41239e-07,   2.39905e-12,   2.15541e-17,
     5.08162e-01,   6.78559e-07,   3.00776e-12,   2.70236e-17,
     6.72474e-01,   8.97974e-07,   3.98043e-12,   3.57640e-17,
     9.48554e-01,   1.26665e-06,   5.61486e-12,   5.04526e-17,
     1.41789e+00,   1.89344e-06,   8.39396e-12,   7.54321e-17,
     2.19504e+00,   2.93130e-06,   1.29962e-11,   1.16808e-16,
     3.46727e+00,   4.63056e-06,   2.05341e-11,   1.84616e-16,
     5.76879e+00,   7.70562e-06,   3.41883e-11,   3.07633e-16,
     1.17206e+01,   1.56620e-05,   6.95719e-11,   6.27184e-16,
     3.16503e+01,   4.23092e-05,   1.88132e-10,   1.69866e-15,
     8.68466e+01,   1.16104e-04,   5.16371e-10,   4.66366e-15,
     1.92342e+02,   2.57132e-04,   1.14345e-09,   1.03252e-14,
     3.36302e+02,   4.49562e-04,   1.99868e-09,   1.80397e-14,
     5.05933e+02,   6.76162e-04,   3.00294e-09,   2.70549e-14,
     7.20757e+02,   9.62603e-04,   4.26374e-09,   3.82441e-14,
     9.77487e+02,   1.30417e-03,   5.75586e-09,   5.13190e-14,
     1.18647e+03,   1.58150e-03,   6.95698e-09,   6.16931e-14,
     1.23844e+03,   1.64960e-03,   7.23901e-09,   6.39380e-14,
     1.11176e+03,   1.48010e-03,   6.48454e-09,   5.71201e-14,
     8.76276e+02,   1.16622e-03,   5.10394e-09,   4.48799e-14,
     6.22103e+02,   8.27770e-04,   3.62020e-09,   3.17970e-14,
     4.07211e+02,   5.41759e-04,   2.36830e-09,   2.07862e-14,
     2.50522e+02,   3.33268e-04,   1.45648e-09,   1.27774e-14,
     1.47042e+02,   1.95598e-04,   8.54671e-10,   7.49580e-15,
     8.32958e+01,   1.10803e-04,   4.84133e-10,   4.24560e-15,
     4.59515e+01,   6.11325e-05,   2.67142e-10,   2.34298e-15  };

  double C13_kpMgO[] = 
  {  2.25390e-04,   3.58434e-10,   1.62608e-15,   1.09114e-20,
     4.04968e-04,   6.44015e-10,   2.92166e-15,   1.96051e-20,
     6.31043e-04,   1.00354e-09,   4.55270e-15,   3.05498e-20,
     9.15654e-04,   1.45615e-09,   6.60605e-15,   4.43284e-20,
     1.52197e-03,   2.42038e-09,   1.09804e-14,   7.36816e-20,
     2.37407e-03,   3.77546e-09,   1.71280e-14,   1.14934e-19,
     3.77208e-03,   5.99871e-09,   2.72141e-14,   1.82615e-19,
     6.14347e-03,   9.76990e-09,   4.43228e-14,   2.97421e-19,
     1.01907e-02,   1.62062e-08,   7.35223e-14,   4.93362e-19,
     1.68896e-02,   2.68594e-08,   1.21853e-13,   8.17686e-19,
     2.96122e-02,   4.70930e-08,   2.13650e-13,   1.43371e-18,
     6.10599e-02,   9.71067e-08,   4.40563e-13,   2.95653e-18,
     1.43395e-01,   2.28052e-07,   1.03468e-12,   6.94381e-18,
     3.27370e-01,   5.20656e-07,   2.36231e-12,   1.58544e-17,
     6.39417e-01,   1.01697e-06,   4.61437e-12,   3.09706e-17,
     1.05137e+00,   1.67234e-06,   7.58927e-12,   5.09499e-17,
     1.55802e+00,   2.48015e-06,   1.12675e-11,   7.57604e-17,
     2.94085e+00,   4.69147e-06,   2.13703e-11,   1.44197e-16,
     1.32769e+01,   2.11739e-05,   9.63657e-11,   6.49212e-16,
     7.10101e+01,   1.13066e-04,   5.13371e-10,   3.44675e-15,
     2.54818e+02,   4.05363e-04,   1.83807e-09,   1.23170e-14,
     6.00377e+02,   9.54596e-04,   4.32539e-09,   2.89547e-14,
     9.94365e+02,   1.58056e-03,   7.15868e-09,   4.78926e-14,
     1.24259e+03,   1.97471e-03,   8.94139e-09,   5.97968e-14,
     1.24573e+03,   1.97945e-03,   8.96134e-09,   5.99160e-14,
     1.05442e+03,   1.67534e-03,   7.58381e-09,   5.06979e-14,
     7.85160e+02,   1.24742e-03,   5.64619e-09,   3.77407e-14,
     5.31062e+02,   8.43677e-04,   3.81852e-09,   2.55222e-14,
     3.34296e+02,   5.31069e-04,   2.40355e-09,   1.60640e-14,
     1.99427e+02,   3.16805e-04,   1.43378e-09,   9.58224e-15,
     1.14271e+02,   1.81526e-04,   8.21530e-10,   5.49033e-15,
     6.35097e+01,   1.00887e-04,   4.56573e-10,   3.05124e-15,
     3.44902e+01,   5.47882e-05,   2.47947e-10,   1.65699e-15,
     1.84022e+01,   2.92323e-05,   1.32292e-10,   8.84081e-16,
     9.68612e+00,   1.53865e-05,   6.96321e-11,   4.65335e-16  };

  double C13_kpFeS[] = 
  {  5.18081e-02,   3.19144e-08,   2.36505e-14,   2.15886e-20,
     9.98885e-02,   6.15325e-08,   4.55993e-14,   4.16239e-20,
     1.60418e-01,   9.88194e-08,   7.32312e-14,   6.68467e-20,
     2.36621e-01,   1.45761e-07,   1.08018e-13,   9.86004e-20,
     3.67284e-01,   2.26251e-07,   1.67666e-13,   1.53048e-19,
     5.36223e-01,   3.30320e-07,   2.44787e-13,   2.23446e-19,
     7.64198e-01,   4.70755e-07,   3.48859e-13,   3.18445e-19,
     1.04970e+00,   6.46627e-07,   4.79190e-13,   4.37414e-19,
     1.38080e+00,   8.50591e-07,   6.30340e-13,   5.75387e-19,
     1.74373e+00,   1.07416e-06,   7.96018e-13,   7.26623e-19,
     2.10311e+00,   1.29554e-06,   9.60079e-13,   8.76383e-19,
     2.42144e+00,   1.49164e-06,   1.10541e-12,   1.00905e-18,
     2.66810e+00,   1.64360e-06,   1.21802e-12,   1.11185e-18,
     2.81617e+00,   1.73481e-06,   1.28563e-12,   1.17358e-18,
     2.89834e+00,   1.78545e-06,   1.32317e-12,   1.20789e-18,
     3.07987e+00,   1.89732e-06,   1.40614e-12,   1.28370e-18,
     3.69380e+00,   2.27562e-06,   1.68664e-12,   1.53995e-18,
     5.05331e+00,   3.11331e-06,   2.30770e-12,   2.10726e-18,
     7.07089e+00,   4.35648e-06,   3.22940e-12,   2.94918e-18,
     9.18826e+00,   5.66121e-06,   4.19684e-12,   3.83303e-18,
     1.08014e+01,   6.65544e-06,   4.93433e-12,   4.50719e-18,
     1.16611e+01,   7.18584e-06,   5.32843e-12,   4.86833e-18,
     1.18910e+01,   7.32863e-06,   5.43584e-12,   4.96854e-18,
     1.17598e+01,   7.24978e-06,   5.38006e-12,   4.92122e-18,
     1.15004e+01,   7.09317e-06,   5.26837e-12,   4.82521e-18,
     1.12405e+01,   6.93914e-06,   5.16238e-12,   4.73954e-18,
     1.10769e+01,   6.85905e-06,   5.13138e-12,   4.75017e-18,
     1.15509e+01,   7.24181e-06,   5.54113e-12,   5.30090e-18,
     1.54818e+01,   9.90289e-06,   7.85008e-12,   7.89232e-18,
     3.23248e+01,   2.09189e-05,   1.69099e-11,   1.74499e-17,
     8.13876e+01,   5.27406e-05,   4.27128e-11,   4.41701e-17,
     1.86528e+02,   1.20594e-04,   9.72543e-11,   9.99819e-17,
     3.66028e+02,   2.35977e-04,   1.89337e-10,   1.93263e-16,
     6.29868e+02,   4.05043e-04,   3.23450e-10,   3.27929e-16,
     9.83739e+02,   6.31245e-04,   5.01948e-10,   5.05703e-16  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = C13_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = C13_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = C13_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = C13_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = C13_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = C13_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = C13_kpMgO     [itab];
      my_rates->SN0_kpFeS     [itab0] = C13_kpFeS     [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_C20(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   1.00183e-01;
  my_rates->SN0_XO [iSN] =   6.06515e-01;
  my_rates->SN0_XMg[iSN] =   2.75968e-02;
  my_rates->SN0_XAl[iSN] =   1.87118e-04;
  my_rates->SN0_XSi[iSN] =   1.00051e-01;
  my_rates->SN0_XS [iSN] =   6.02208e-02;
  my_rates->SN0_XFe[iSN] =   3.07560e-02;

  my_rates->SN0_fC [iSN] =   8.74563e-02;
  my_rates->SN0_fO [iSN] =   6.04383e-01;
  my_rates->SN0_fMg[iSN] =   2.63753e-02;
  my_rates->SN0_fAl[iSN] =   1.87118e-04;
  my_rates->SN0_fSi[iSN] =   6.44592e-02;
  my_rates->SN0_fS [iSN] =   6.02018e-02;
  my_rates->SN0_fFe[iSN] =   2.69505e-02;

  my_rates->SN0_fSiM     [iSN] =   3.44388e-02;
  my_rates->SN0_fFeM     [iSN] =   3.77223e-03;
  my_rates->SN0_fMg2SiO4 [iSN] =   1.90086e-03;
  my_rates->SN0_fMgSiO3  [iSN] =   2.57266e-06;
  my_rates->SN0_fAC      [iSN] =   1.27270e-02;
  my_rates->SN0_fSiO2D   [iSN] =   1.65484e-03;
  my_rates->SN0_fMgO     [iSN] =   9.48713e-04;
  my_rates->SN0_fFeS     [iSN] =   5.23050e-05;
  my_rates->SN0_fAl2O3   [iSN] =   1.31693e-29;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   1.24861e-05;
  my_rates->SN0_r0FeM     [itab0 + 0] =   6.67024e-06;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   1.41253e-06;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   1.01138e-06;
  my_rates->SN0_r0AC      [itab0 + 0] =   7.95099e-07;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   1.40285e-06;
  my_rates->SN0_r0MgO     [itab0 + 0] =   1.29303e-06;
  my_rates->SN0_r0FeS     [itab0 + 0] =   1.68897e-06;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   9.21063e-08;

  my_rates->SN0_r0SiM     [itab0 + 1] =   2.86508e-10;
  my_rates->SN0_r0FeM     [itab0 + 1] =   7.50596e-11;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   4.77566e-12;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.31688e-12;
  my_rates->SN0_r0AC      [itab0 + 1] =   2.51133e-12;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   3.98828e-12;
  my_rates->SN0_r0MgO     [itab0 + 1] =   1.06240e-11;
  my_rates->SN0_r0FeS     [itab0 + 1] =   3.16618e-12;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   9.03508e-15;

  my_rates->SN0_r0SiM     [itab0 + 2] =   1.01028e-14;
  my_rates->SN0_r0FeM     [itab0 + 2] =   1.22752e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   3.08016e-17;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   2.89696e-18;
  my_rates->SN0_r0AC      [itab0 + 2] =   4.21640e-17;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   1.93974e-17;
  my_rates->SN0_r0MgO     [itab0 + 2] =   1.57687e-16;
  my_rates->SN0_r0FeS     [itab0 + 2] =   6.72598e-18;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   9.36936e-22;

  NTd =            35;
 Nmom =             4;

  double C20_kpSiM[] = 
  {  1.53894e-01,   1.90916e-06,   4.34900e-11,   1.52207e-15,
     1.93844e-01,   2.40648e-06,   5.48634e-11,   1.92178e-15,
     2.44138e-01,   3.03256e-06,   6.91797e-11,   2.42474e-15,
     3.07454e-01,   3.82073e-06,   8.72020e-11,   3.05783e-15,
     3.87243e-01,   4.81526e-06,   1.09978e-10,   3.85936e-15,
     4.87709e-01,   6.06778e-06,   1.38669e-10,   4.86916e-15,
     6.14251e-01,   7.64642e-06,   1.74856e-10,   6.14383e-15,
     7.73625e-01,   9.63590e-06,   2.20495e-10,   7.75268e-15,
     9.74036e-01,   1.21392e-05,   2.77959e-10,   9.78002e-15,
     1.22376e+00,   1.52600e-05,   3.49642e-10,   1.23105e-14,
     1.52029e+00,   1.89682e-05,   4.34881e-10,   1.53223e-14,
     1.83650e+00,   2.29254e-05,   5.25941e-10,   1.85447e-14,
     2.15714e+00,   2.69443e-05,   6.18616e-10,   2.18354e-14,
     2.55729e+00,   3.19712e-05,   7.34846e-10,   2.59753e-14,
     3.23398e+00,   4.04866e-05,   9.32032e-10,   3.30022e-14,
     4.34499e+00,   5.44911e-05,   1.25683e-09,   4.45846e-14,
     5.84259e+00,   7.34292e-05,   1.69748e-09,   6.03345e-14,
     7.53509e+00,   9.49825e-05,   2.20255e-09,   7.84868e-14,
     9.31292e+00,   1.17996e-04,   2.75076e-09,   9.84428e-14,
     1.14053e+01,   1.46050e-04,   3.44207e-09,   1.24265e-13,
     1.44699e+01,   1.88865e-04,   4.53802e-09,   1.66375e-13,
     1.88525e+01,   2.52409e-04,   6.22006e-09,   2.32633e-13,
     2.35897e+01,   3.27106e-04,   8.33967e-09,   3.20293e-13,
     2.71065e+01,   3.99935e-04,   1.08043e-08,   4.33797e-13,
     2.87408e+01,   4.69399e-04,   1.38602e-08,   5.94648e-13,
     2.88428e+01,   5.34587e-04,   1.74785e-08,   8.07129e-13,
     2.80514e+01,   5.86377e-04,   2.09655e-08,   1.03092e-12,
     2.67893e+01,   6.13530e-04,   2.33810e-08,   1.20137e-12,
     2.53851e+01,   6.14601e-04,   2.43011e-08,   1.28034e-12,
     2.53848e+01,   6.23633e-04,   2.48207e-08,   1.31287e-12,
     3.52560e+01,   8.08661e-04,   3.01974e-08,   1.52398e-12,
     9.64731e+01,   1.84326e-03,   5.73265e-08,   2.51143e-12,
     3.56496e+02,   5.70937e-03,   1.44867e-07,   5.33930e-12,
     1.16716e+03,   1.65654e-02,   3.60906e-07,   1.15456e-11,
     3.02792e+03,   3.93985e-02,   7.71226e-07,   2.22345e-11  };

  double C20_kpFeM[] = 
  {  1.10506e-02,   1.69983e-07,   3.50089e-12,   9.16404e-17,
     1.85837e-02,   2.76971e-07,   5.53036e-12,   1.40642e-16,
     2.79449e-02,   4.08949e-07,   8.01809e-12,   2.00458e-16,
     3.96605e-02,   5.73574e-07,   1.11125e-11,   2.74730e-16,
     6.19392e-02,   8.69904e-07,   1.63567e-11,   3.93140e-16,
     9.16817e-02,   1.25867e-06,   2.31128e-11,   5.43193e-16,
     1.36080e-01,   1.82251e-06,   3.26005e-11,   7.47226e-16,
     2.02056e-01,   2.63811e-06,   4.59129e-11,   1.02492e-15,
     2.99593e-01,   3.81390e-06,   6.45602e-11,   1.40281e-15,
     4.39758e-01,   5.46636e-06,   9.01024e-11,   1.90718e-15,
     6.41903e-01,   7.79581e-06,   1.25163e-10,   2.58082e-15,
     9.24985e-01,   1.09870e-05,   1.71974e-10,   3.45658e-15,
     1.30585e+00,   1.51941e-05,   2.32250e-10,   4.55695e-15,
     1.80164e+00,   2.05666e-05,   3.07564e-10,   5.90101e-15,
     2.42648e+00,   2.72192e-05,   3.99035e-10,   7.50126e-15,
     3.19207e+00,   3.52384e-05,   5.07449e-10,   9.36601e-15,
     4.10565e+00,   4.46618e-05,   6.32978e-10,   1.14947e-14,
     5.17172e+00,   5.54917e-05,   7.75385e-10,   1.38820e-14,
     6.40478e+00,   6.78244e-05,   9.35760e-10,   1.65481e-14,
     7.85614e+00,   8.21045e-05,   1.11980e-09,   1.95943e-14,
     9.64785e+00,   9.94326e-05,   1.34166e-09,   2.32660e-14,
     1.20159e+01,   1.21941e-04,   1.62875e-09,   2.80374e-14,
     1.53778e+01,   1.53371e-04,   2.02922e-09,   3.47445e-14,
     2.04423e+01,   1.99938e-04,   2.62238e-09,   4.47607e-14,
     2.83887e+01,   2.71546e-04,   3.53030e-09,   6.01372e-14,
     4.11773e+01,   3.83718e-04,   4.93303e-09,   8.37180e-14,
     6.21107e+01,   5.61284e-04,   7.10254e-09,   1.19520e-13,
     9.68132e+01,   8.45164e-04,   1.04719e-08,   1.73689e-13,
     1.54804e+02,   1.30324e-03,   1.57446e-08,   2.55980e-13,
     2.51780e+02,   2.04562e-03,   2.40381e-08,   3.81486e-13,
     4.12838e+02,   3.24595e-03,   3.70714e-08,   5.72650e-13,
     6.76592e+02,   5.16747e-03,   5.73792e-08,   8.61246e-13,
     1.09863e+03,   8.17937e-03,   8.83910e-08,   1.28821e-12,
     1.74646e+03,   1.26954e-02,   1.33666e-07,   1.89228e-12,
     2.67882e+03,   1.89754e-02,   1.94753e-07,   2.68226e-12  };

  double C20_kpMg2SiO4[] = 
  {  1.05240e-01,   1.48654e-07,   5.02591e-13,   3.24156e-18,
     1.32588e-01,   1.87284e-07,   6.33194e-13,   4.08391e-18,
     1.67016e-01,   2.35915e-07,   7.97614e-13,   5.14437e-18,
     2.10360e-01,   2.97139e-07,   1.00461e-12,   6.47941e-18,
     2.71887e-01,   3.84048e-07,   1.29844e-12,   8.37459e-18,
     3.55694e-01,   5.02428e-07,   1.69868e-12,   1.09560e-17,
     4.84932e-01,   6.84981e-07,   2.31589e-12,   1.49369e-17,
     6.99767e-01,   9.88442e-07,   3.34188e-12,   2.15544e-17,
     1.05860e+00,   1.49530e-06,   5.05559e-12,   3.26077e-17,
     1.62902e+00,   2.30104e-06,   7.77981e-12,   5.01791e-17,
     2.54260e+00,   3.59152e-06,   1.21430e-11,   7.83229e-17,
     3.96490e+00,   5.60062e-06,   1.89363e-11,   1.22144e-16,
     6.10634e+00,   8.62558e-06,   2.91646e-11,   1.88128e-16,
     9.28774e+00,   1.31196e-05,   4.43615e-11,   2.86175e-16,
     1.39267e+01,   1.96729e-05,   6.65240e-11,   4.29191e-16,
     2.05387e+01,   2.90138e-05,   9.81191e-11,   6.33137e-16,
     3.00660e+01,   4.24748e-05,   1.43665e-10,   9.27312e-16,
     4.55128e+01,   6.43028e-05,   2.17560e-10,   1.40505e-15,
     7.47912e+01,   1.05684e-04,   3.57721e-10,   2.31202e-15,
     1.29637e+02,   1.83206e-04,   6.20346e-10,   4.01209e-15,
     2.14853e+02,   3.03661e-04,   1.02847e-09,   6.65466e-15,
     3.20060e+02,   4.52395e-04,   1.53258e-09,   9.92043e-15,
     4.29862e+02,   6.07661e-04,   2.05913e-09,   1.33350e-14,
     5.30972e+02,   7.50675e-04,   2.54443e-09,   1.64851e-14,
     5.99878e+02,   8.48166e-04,   2.87540e-09,   1.86348e-14,
     6.06737e+02,   8.57905e-04,   2.90867e-09,   1.88528e-14,
     5.43451e+02,   7.68444e-04,   2.60547e-09,   1.68884e-14,
     4.33701e+02,   6.13270e-04,   2.07941e-09,   1.34791e-14,
     3.13436e+02,   4.43222e-04,   1.50292e-09,   9.74306e-15,
     2.09092e+02,   2.95687e-04,   1.00276e-09,   6.50192e-15,
     1.31208e+02,   1.85574e-04,   6.29536e-10,   4.08396e-15,
     7.91293e+01,   1.11988e-04,   3.80343e-10,   2.47142e-15,
     4.74361e+01,   6.72978e-05,   2.29477e-10,   1.49906e-15,
     2.98121e+01,   4.25812e-05,   1.46690e-10,   9.70392e-16,
     2.07674e+01,   3.01013e-05,   1.05718e-10,   7.13724e-16  };

  double C20_kpMgSiO3[] = 
  {  2.19890e-02,   2.22393e-08,   2.89570e-14,   6.37012e-20,
     3.90612e-02,   3.95059e-08,   5.14391e-14,   1.13159e-19,
     6.05539e-02,   6.12433e-08,   7.97425e-14,   1.75423e-19,
     8.76116e-02,   8.86090e-08,   1.15374e-13,   2.53808e-19,
     1.43288e-01,   1.44919e-07,   1.88693e-13,   4.15103e-19,
     2.19266e-01,   2.21762e-07,   2.88749e-13,   6.35214e-19,
     3.36256e-01,   3.40084e-07,   4.42810e-13,   9.74136e-19,
     5.14336e-01,   5.20191e-07,   6.77322e-13,   1.49005e-18,
     7.97217e-01,   8.06292e-07,   1.04985e-12,   2.30962e-18,
     1.25414e+00,   1.26842e-06,   1.65157e-12,   3.63350e-18,
     2.03450e+00,   2.05766e-06,   2.67923e-12,   5.89467e-18,
     3.34648e+00,   3.38458e-06,   4.40703e-12,   9.69677e-18,
     5.45894e+00,   5.52111e-06,   7.18906e-12,   1.58199e-17,
     8.82120e+00,   8.92167e-06,   1.16172e-11,   2.55686e-17,
     1.41826e+01,   1.43441e-05,   1.86785e-11,   4.11217e-17,
     2.28421e+01,   2.31025e-05,   3.00847e-11,   6.62628e-17,
     3.71183e+01,   3.75417e-05,   4.88915e-11,   1.07761e-16,
     6.14292e+01,   6.21307e-05,   8.09232e-11,   1.78540e-16,
     1.03850e+02,   1.05037e-04,   1.36825e-10,   3.02237e-16,
     1.75513e+02,   1.77524e-04,   2.31276e-10,   5.11420e-16,
     2.82073e+02,   2.85307e-04,   3.71726e-10,   8.22603e-16,
     4.14541e+02,   4.19298e-04,   5.46336e-10,   1.20964e-15,
     5.60007e+02,   5.66440e-04,   7.38095e-10,   1.63478e-15,
     7.11090e+02,   7.19267e-04,   9.37240e-10,   2.07555e-15,
     8.40892e+02,   8.50566e-04,   1.10826e-09,   2.45241e-15,
     8.95414e+02,   9.05715e-04,   1.18001e-09,   2.60842e-15,
     8.40504e+02,   8.50170e-04,   1.10753e-09,   2.44576e-15,
     6.96768e+02,   7.04778e-04,   9.18057e-10,   2.02574e-15,
     5.18260e+02,   5.24216e-04,   6.82816e-10,   1.50584e-15,
     3.52903e+02,   3.56958e-04,   4.64940e-10,   1.02501e-15,
     2.24241e+02,   2.26818e-04,   2.95430e-10,   6.51237e-16,
     1.35153e+02,   1.36707e-04,   1.78065e-10,   3.92585e-16,
     7.83236e+01,   7.92263e-05,   1.03204e-10,   2.27687e-16,
     4.41657e+01,   4.46783e-05,   5.82153e-11,   1.28647e-16,
     2.46133e+01,   2.49140e-05,   3.25047e-11,   7.21799e-17  };

  double C20_kpAC[] = 
  {  3.27960e-01,   2.60760e-07,   8.23594e-13,   1.38274e-17,
     4.38752e-01,   3.48855e-07,   1.10197e-12,   1.85031e-17,
     5.78230e-01,   4.59761e-07,   1.45242e-12,   2.43895e-17,
     7.53824e-01,   5.99382e-07,   1.89361e-12,   3.18000e-17,
     1.04013e+00,   8.27053e-07,   2.61333e-12,   4.38935e-17,
     1.41735e+00,   1.12702e-06,   3.56171e-12,   5.98308e-17,
     1.95293e+00,   1.55292e-06,   4.90855e-12,   8.24705e-17,
     2.71532e+00,   2.15922e-06,   6.82677e-12,   1.14729e-16,
     3.79677e+00,   3.01935e-06,   9.54991e-12,   1.60553e-16,
     5.29747e+00,   4.21303e-06,   1.33318e-11,   2.24238e-16,
     7.37841e+00,   5.86846e-06,   1.85820e-11,   3.12737e-16,
     1.02169e+01,   8.12703e-06,   2.57565e-11,   4.33854e-16,
     1.40423e+01,   1.11717e-05,   3.54480e-11,   5.97791e-16,
     1.92026e+01,   1.52804e-05,   4.85668e-11,   8.20344e-16,
     2.61626e+01,   2.08251e-05,   6.63455e-11,   1.12316e-15,
     3.55324e+01,   2.82955e-05,   9.04435e-11,   1.53595e-15,
     4.81644e+01,   3.83784e-05,   1.23253e-10,   2.10252e-15,
     6.53233e+01,   5.20969e-05,   1.68441e-10,   2.89162e-15,
     8.87789e+01,   7.08925e-05,   2.31406e-10,   4.00785e-15,
     1.20729e+02,   9.65780e-05,   3.19528e-10,   5.60280e-15,
     1.63666e+02,   1.31274e-04,   4.42922e-10,   7.90380e-15,
     2.20664e+02,   1.77696e-04,   6.17150e-10,   1.12922e-14,
     2.96276e+02,   2.39957e-04,   8.67756e-10,   1.64168e-14,
     3.97347e+02,   3.24130e-04,   1.23022e-09,   2.41672e-14,
     5.32072e+02,   4.37186e-04,   1.73755e-09,   3.52998e-14,
     7.07827e+02,   5.84876e-04,   2.40229e-09,   4.99094e-14,
     9.32126e+02,   7.72477e-04,   3.21611e-09,   6.73490e-14,
     1.21808e+03,   1.00957e-03,   4.17439e-09,   8.68091e-14,
     1.58941e+03,   1.31461e-03,   5.30260e-09,   1.07984e-13,
     2.08259e+03,   1.71725e-03,   6.66833e-09,   1.31338e-13,
     2.74876e+03,   2.26047e-03,   8.38813e-09,   1.58137e-13,
     3.65886e+03,   3.00609e-03,   1.06436e-08,   1.90525e-13,
     4.91536e+03,   4.04510e-03,   1.37047e-08,   2.31589e-13,
     6.67260e+03,   5.51150e-03,   1.79301e-08,   2.84849e-13,
     9.16963e+03,   7.59478e-03,   2.36873e-08,   3.52355e-13  };

  double C20_kpSiO2D[] = 
  {  7.60358e-02,   1.06666e-07,   3.03247e-13,   1.47482e-18,
     9.07206e-02,   1.27267e-07,   3.61815e-13,   1.75967e-18,
     1.09208e-01,   1.53201e-07,   4.35546e-13,   2.11827e-18,
     1.32481e-01,   1.85851e-07,   5.28369e-13,   2.56972e-18,
     1.58907e-01,   2.22922e-07,   6.33765e-13,   3.08233e-18,
     1.91565e-01,   2.68735e-07,   7.64012e-13,   3.71581e-18,
     2.30490e-01,   3.23342e-07,   9.19258e-13,   4.47087e-18,
     2.76795e-01,   3.88301e-07,   1.10394e-12,   5.36908e-18,
     3.33074e-01,   4.67253e-07,   1.32840e-12,   6.46082e-18,
     4.05326e-01,   5.68613e-07,   1.61658e-12,   7.86245e-18,
     5.08161e-01,   7.12876e-07,   2.02672e-12,   9.85731e-18,
     6.72474e-01,   9.43381e-07,   2.68206e-12,   1.30448e-17,
     9.48552e-01,   1.33068e-06,   3.78319e-12,   1.84008e-17,
     1.41789e+00,   1.98912e-06,   5.65528e-12,   2.75074e-17,
     2.19503e+00,   3.07934e-06,   8.75500e-12,   4.25864e-17,
     3.46724e+00,   4.86415e-06,   1.38300e-11,   6.72794e-17,
     5.76869e+00,   8.09313e-06,   2.30131e-11,   1.11983e-16,
     1.17202e+01,   1.64442e-05,   4.67706e-11,   2.27727e-16,
     3.16486e+01,   4.44090e-05,   1.26336e-10,   6.15454e-16,
     8.68419e+01,   1.21859e-04,   3.46686e-10,   1.68908e-15,
     1.92333e+02,   2.69887e-04,   7.67805e-10,   3.74058e-15,
     3.36287e+02,   4.71894e-04,   1.34247e-09,   6.53939e-15,
     5.05924e+02,   7.09951e-04,   2.01945e-09,   9.83185e-15,
     7.20787e+02,   1.01144e-03,   2.87594e-09,   1.39829e-14,
     9.77602e+02,   1.37169e-03,   3.89813e-09,   1.89179e-14,
     1.18669e+03,   1.66488e-03,   4.72879e-09,   2.29108e-14,
     1.23874e+03,   1.73774e-03,   4.93372e-09,   2.38742e-14,
     1.11204e+03,   1.55988e-03,   4.42750e-09,   2.14067e-14,
     8.76525e+02,   1.22945e-03,   3.48896e-09,   1.68597e-14,
     6.22287e+02,   8.72816e-04,   2.47658e-09,   1.19634e-14,
     4.07334e+02,   5.71310e-04,   1.62094e-09,   7.82836e-15,
     2.50598e+02,   3.51474e-04,   9.97166e-10,   4.81513e-15,
     1.47088e+02,   2.06293e-04,   5.85254e-10,   2.82583e-15,
     8.33219e+01,   1.16862e-04,   3.31541e-10,   1.60077e-15,
     4.59660e+01,   6.44731e-05,   1.82924e-10,   8.83269e-16  };

  double C20_kpMgO[] = 
  {  2.25389e-04,   2.91423e-10,   2.39426e-15,   3.55346e-20,
     4.04967e-04,   5.23622e-10,   4.30206e-15,   6.38511e-20,
     6.31042e-04,   8.15942e-10,   6.70384e-15,   9.94996e-20,
     9.15653e-04,   1.18395e-09,   9.72751e-15,   1.44378e-19,
     1.52197e-03,   1.96795e-09,   1.61693e-14,   2.39993e-19,
     2.37407e-03,   3.06977e-09,   2.52225e-14,   3.74371e-19,
     3.77209e-03,   4.87751e-09,   4.00763e-14,   5.94853e-19,
     6.14348e-03,   7.94395e-09,   6.52733e-14,   9.68874e-19,
     1.01907e-02,   1.31776e-08,   1.08281e-13,   1.60732e-18,
     1.68897e-02,   2.18408e-08,   1.79477e-13,   2.66430e-18,
     2.96125e-02,   3.82967e-08,   3.14746e-13,   4.67293e-18,
     6.10610e-02,   7.89818e-08,   6.49304e-13,   9.64262e-18,
     1.43400e-01,   1.85524e-07,   1.52569e-12,   2.26652e-17,
     3.27383e-01,   4.23645e-07,   3.48510e-12,   5.17901e-17,
     6.39452e-01,   8.27698e-07,   6.81196e-12,   1.01271e-16,
     1.05149e+00,   1.36261e-06,   1.12352e-11,   1.67335e-16,
     1.55882e+00,   2.03505e-06,   1.69775e-11,   2.55740e-16,
     2.94474e+00,   3.91135e-06,   3.34978e-11,   5.17225e-16,
     1.32877e+01,   1.75173e-05,   1.48167e-10,   2.26035e-15,
     7.10035e+01,   9.20831e-05,   7.58876e-10,   1.12902e-14,
     2.54672e+02,   3.27285e-04,   2.65806e-09,   3.89823e-14,
     5.99878e+02,   7.67225e-04,   6.18279e-09,   8.99815e-14,
     9.93398e+02,   1.26701e-03,   1.01648e-08,   1.47279e-13,
     1.24125e+03,   1.58042e-03,   1.26439e-08,   1.82695e-13,
     1.24435e+03,   1.58257e-03,   1.26385e-08,   1.82292e-13,
     1.05321e+03,   1.33853e-03,   1.06770e-08,   1.53820e-13,
     7.84243e+02,   9.96171e-04,   7.93973e-09,   1.14294e-13,
     5.30434e+02,   6.73540e-04,   5.36536e-09,   7.71941e-14,
     3.33897e+02,   4.23882e-04,   3.37536e-09,   4.85451e-14,
     1.99187e+02,   2.52825e-04,   2.01271e-09,   2.89397e-14,
     1.14133e+02,   1.44851e-04,   1.15293e-09,   1.65742e-14,
     6.34329e+01,   8.04978e-05,   6.40626e-10,   9.20833e-15,
     3.44484e+01,   4.37133e-05,   3.47851e-10,   4.99954e-15,
     1.83799e+01,   2.33224e-05,   1.85579e-10,   2.66709e-15,
     9.67434e+00,   1.22756e-05,   9.76752e-11,   1.40372e-15  };

  double C20_kpFeS[] = 
  {  5.18099e-02,   8.75057e-08,   1.64042e-13,   3.48481e-19,
     9.98914e-02,   1.68714e-07,   3.16278e-13,   6.71882e-19,
     1.60422e-01,   2.70949e-07,   5.07932e-13,   1.07902e-18,
     2.36626e-01,   3.99656e-07,   7.49210e-13,   1.59157e-18,
     3.67294e-01,   6.20351e-07,   1.16293e-12,   2.47046e-18,
     5.36239e-01,   9.05695e-07,   1.69785e-12,   3.60681e-18,
     7.64225e-01,   1.29076e-06,   2.41970e-12,   5.14028e-18,
     1.04974e+00,   1.77298e-06,   3.32371e-12,   7.06072e-18,
     1.38085e+00,   2.33223e-06,   4.37210e-12,   9.28786e-18,
     1.74382e+00,   2.94529e-06,   5.52140e-12,   1.17294e-17,
     2.10325e+00,   3.55236e-06,   6.65945e-12,   1.41471e-17,
     2.42167e+00,   4.09020e-06,   7.66778e-12,   1.62893e-17,
     2.66848e+00,   4.50708e-06,   8.44938e-12,   1.79499e-17,
     2.81672e+00,   4.75754e-06,   8.91905e-12,   1.89480e-17,
     2.89933e+00,   4.89720e-06,   9.18118e-12,   1.95057e-17,
     3.08200e+00,   5.20610e-06,   9.76116e-12,   2.07399e-17,
     3.69875e+00,   6.24871e-06,   1.17179e-11,   2.49020e-17,
     5.06356e+00,   8.55561e-06,   1.60466e-11,   3.41079e-17,
     7.08905e+00,   1.19792e-05,   2.24709e-11,   4.77704e-17,
     9.21663e+00,   1.55761e-05,   2.92218e-11,   6.21311e-17,
     1.08429e+01,   1.83272e-05,   3.43894e-11,   7.31341e-17,
     1.17217e+01,   1.98178e-05,   3.71989e-11,   7.91392e-17,
     1.19809e+01,   2.02656e-05,   3.80618e-11,   8.10292e-17,
     1.18982e+01,   2.01424e-05,   3.78698e-11,   8.07160e-17,
     1.17194e+01,   1.98680e-05,   3.74208e-11,   7.99219e-17,
     1.16101e+01,   1.97358e-05,   3.72971e-11,   7.99631e-17,
     1.19784e+01,   2.05514e-05,   3.92888e-11,   8.53362e-17,
     1.49136e+01,   2.65380e-05,   5.30405e-11,   1.20988e-16,
     2.57980e+01,   4.84122e-05,   1.02966e-10,   2.50769e-16,
     6.12919e+01,   1.18083e-04,   2.58706e-10,   6.49359e-16,
     1.55902e+02,   3.00837e-04,   6.60279e-10,   1.66030e-15,
     3.47036e+02,   6.65409e-04,   1.45009e-09,   3.62025e-15,
     6.57260e+02,   1.25056e-03,   2.70176e-09,   6.68619e-15,
     1.09327e+03,   2.06455e-03,   4.42254e-09,   1.08501e-14,
     1.65394e+03,   3.10062e-03,   6.58692e-09,   1.60229e-14  };

  double C20_kpAl2O3[] = 
  {  9.93250e-04,   9.14846e-11,   8.97410e-18,   9.30612e-25,
     1.81240e-03,   1.66933e-10,   1.63752e-17,   1.69810e-24,
     2.84365e-03,   2.61918e-10,   2.56926e-17,   2.66432e-24,
     4.14191e-03,   3.81496e-10,   3.74225e-17,   3.88071e-24,
     7.18271e-03,   6.61573e-10,   6.48964e-17,   6.72974e-24,
     1.13364e-02,   1.04415e-09,   1.02425e-16,   1.06215e-23,
     1.77361e-02,   1.63360e-09,   1.60247e-16,   1.66176e-23,
     2.59477e-02,   2.38995e-09,   2.34440e-16,   2.43114e-23,
     3.45425e-02,   3.18159e-09,   3.12095e-16,   3.23642e-23,
     4.22006e-02,   3.88695e-09,   3.81286e-16,   3.95393e-23,
     4.71420e-02,   4.34208e-09,   4.25932e-16,   4.41691e-23,
     4.91934e-02,   4.53102e-09,   4.44466e-16,   4.60911e-23,
     5.05162e-02,   4.65286e-09,   4.56418e-16,   4.73304e-23,
     5.78201e-02,   5.32560e-09,   5.22410e-16,   5.41738e-23,
     8.84237e-02,   8.14438e-09,   7.98916e-16,   8.28474e-23,
     1.78786e-01,   1.64673e-08,   1.61535e-15,   1.67511e-22,
     4.36404e-01,   4.01956e-08,   3.94295e-15,   4.08884e-22,
     1.63796e+00,   1.50867e-07,   1.47992e-14,   1.53467e-21,
     8.50819e+00,   7.83659e-07,   7.68723e-14,   7.97165e-21,
     3.92751e+01,   3.61749e-06,   3.54854e-13,   3.67984e-20,
     1.41439e+02,   1.30275e-05,   1.27792e-12,   1.32520e-19,
     3.83709e+02,   3.53420e-05,   3.46684e-12,   3.59511e-19,
     7.70411e+02,   7.09598e-05,   6.96073e-12,   7.21827e-19,
     1.16399e+03,   1.07211e-04,   1.05167e-11,   1.09058e-18,
     1.37566e+03,   1.26707e-04,   1.24292e-11,   1.28891e-18,
     1.33070e+03,   1.22566e-04,   1.20230e-11,   1.24678e-18,
     1.09978e+03,   1.01297e-04,   9.93663e-12,   1.03043e-18,
     8.05638e+02,   7.42044e-05,   7.27901e-12,   7.54832e-19,
     5.38690e+02,   4.96167e-05,   4.86711e-12,   5.04718e-19,
     3.36338e+02,   3.09789e-05,   3.03884e-12,   3.15127e-19,
     1.99460e+02,   1.83715e-05,   1.80214e-12,   1.86881e-19,
     1.13787e+02,   1.04805e-05,   1.02808e-12,   1.06611e-19,
     6.30411e+01,   5.80648e-06,   5.69582e-13,   5.90655e-20,
     3.41529e+01,   3.14570e-06,   3.08575e-13,   3.19991e-20,
     1.81893e+01,   1.67535e-06,   1.64342e-13,   1.70422e-20  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = C20_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = C20_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = C20_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = C20_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = C20_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = C20_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = C20_kpMgO     [itab];
      my_rates->SN0_kpFeS     [itab0] = C20_kpFeS     [itab];
      my_rates->SN0_kpAl2O3   [itab0] = C20_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_C25(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   1.75488e-01;
  my_rates->SN0_XO [iSN] =   5.69674e-01;
  my_rates->SN0_XMg[iSN] =   3.12340e-02;
  my_rates->SN0_XAl[iSN] =   2.98415e-04;
  my_rates->SN0_XSi[iSN] =   8.33205e-02;
  my_rates->SN0_XS [iSN] =   4.73930e-02;
  my_rates->SN0_XFe[iSN] =   1.98197e-02;

  my_rates->SN0_fC [iSN] =   1.34092e-01;
  my_rates->SN0_fO [iSN] =   5.53726e-01;
  my_rates->SN0_fMg[iSN] =   2.48100e-02;
  my_rates->SN0_fAl[iSN] =   2.98415e-04;
  my_rates->SN0_fSi[iSN] =   3.47760e-02;
  my_rates->SN0_fS [iSN] =   4.72556e-02;
  my_rates->SN0_fFe[iSN] =   1.46955e-02;

  my_rates->SN0_fSiM     [iSN] =   3.83373e-02;
  my_rates->SN0_fFeM     [iSN] =   4.88366e-03;
  my_rates->SN0_fMg2SiO4 [iSN] =   1.68068e-02;
  my_rates->SN0_fMgSiO3  [iSN] =   2.49736e-05;
  my_rates->SN0_fAC      [iSN] =   4.13961e-02;
  my_rates->SN0_fSiO2D   [iSN] =   1.46546e-02;
  my_rates->SN0_fMgO     [iSN] =   1.09289e-03;
  my_rates->SN0_fFeS     [iSN] =   3.77935e-04;
  my_rates->SN0_fAl2O3   [iSN] =   1.65550e-31;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   1.72153e-05;
  my_rates->SN0_r0FeM     [itab0 + 0] =   1.96666e-05;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   2.33213e-06;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   1.55439e-06;
  my_rates->SN0_r0AC      [itab0 + 0] =   7.93494e-07;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   2.56804e-06;
  my_rates->SN0_r0MgO     [itab0 + 0] =   3.58420e-06;
  my_rates->SN0_r0FeS     [itab0 + 0] =   9.61035e-07;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   1.99526e-08;

  my_rates->SN0_r0SiM     [itab0 + 1] =   6.33208e-10;
  my_rates->SN0_r0FeM     [itab0 + 1] =   5.88305e-10;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   2.48648e-11;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   4.30058e-12;
  my_rates->SN0_r0AC      [itab0 + 1] =   3.53402e-12;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   4.82971e-11;
  my_rates->SN0_r0MgO     [itab0 + 1] =   3.09713e-11;
  my_rates->SN0_r0FeS     [itab0 + 1] =   2.46507e-12;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   3.98107e-16;

  my_rates->SN0_r0SiM     [itab0 + 2] =   4.04318e-14;
  my_rates->SN0_r0FeM     [itab0 + 2] =   2.42323e-14;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   4.29427e-16;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   1.92568e-17;
  my_rates->SN0_r0AC      [itab0 + 2] =   1.04050e-16;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   2.53766e-15;
  my_rates->SN0_r0MgO     [itab0 + 2] =   4.03929e-16;
  my_rates->SN0_r0FeS     [itab0 + 2] =   1.42549e-17;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   7.94328e-24;

  NTd =            35;
 Nmom =             4;

  double C25_kpSiM[] = 
  {  1.53307e-01,   2.58151e-06,   8.97185e-11,   5.13410e-15,
     1.93187e-01,   3.26103e-06,   1.14053e-10,   6.60852e-15,
     2.43381e-01,   4.11484e-06,   1.44443e-10,   8.42753e-15,
     3.06566e-01,   5.18903e-06,   1.82603e-10,   1.07022e-14,
     3.86268e-01,   6.55217e-06,   2.31853e-10,   1.37374e-14,
     4.86630e-01,   8.26891e-06,   2.93869e-10,   1.75576e-14,
     6.13093e-01,   1.04376e-05,   3.72751e-10,   2.24815e-14,
     7.72441e-01,   1.31779e-05,   4.73205e-10,   2.88483e-14,
     9.72908e-01,   1.66348e-05,   6.00927e-10,   3.70670e-14,
     1.22279e+00,   2.09532e-05,   7.61455e-10,   4.75186e-14,
     1.51967e+00,   2.61010e-05,   9.54694e-10,   6.03376e-14,
     1.83660e+00,   3.16418e-05,   1.16836e-09,   7.52509e-14,
     2.15883e+00,   3.73832e-05,   1.40366e-09,   9.34936e-14,
     2.56188e+00,   4.46730e-05,   1.71564e-09,   1.19386e-13,
     3.24331e+00,   5.69556e-05,   2.23297e-09,   1.61351e-13,
     4.36192e+00,   7.70449e-05,   3.06366e-09,   2.26659e-13,
     5.87089e+00,   1.04185e-04,   4.17746e-09,   3.12654e-13,
     7.58000e+00,   1.35161e-04,   5.44737e-09,   4.09330e-13,
     9.38530e+00,   1.68565e-04,   6.82760e-09,   5.12932e-13,
     1.15371e+01,   2.10285e-04,   8.59830e-09,   6.44799e-13,
     1.47388e+01,   2.75860e-04,   1.14755e-08,   8.58985e-13,
     1.93877e+01,   3.75887e-04,   1.59953e-08,   1.19600e-12,
     2.46008e+01,   5.00971e-04,   2.20360e-08,   1.65571e-12,
     2.90408e+01,   6.45456e-04,   3.01610e-08,   2.31085e-12,
     3.23701e+01,   8.27330e-04,   4.23216e-08,   3.35508e-12,
     3.49305e+01,   1.05238e-03,   5.92405e-08,   4.87268e-12,
     3.68020e+01,   1.28200e-03,   7.79075e-08,   6.59988e-12,
     3.76194e+01,   1.45174e-03,   9.27975e-08,   8.01634e-12,
     3.72617e+01,   1.52446e-03,   1.00177e-07,   8.74632e-12,
     3.76576e+01,   1.55784e-03,   1.02687e-07,   8.96074e-12,
     4.91245e+01,   1.84545e-03,   1.14505e-07,   9.64022e-12,
     1.17100e+02,   3.32870e-03,   1.70057e-07,   1.27411e-11,
     3.93519e+02,   8.31574e-03,   3.27938e-07,   2.07522e-11,
     1.23314e+03,   2.12665e-02,   6.72506e-07,   3.61691e-11,
     3.12736e+03,   4.71942e-02,   1.26825e-06,   5.97225e-11  };

  double C25_kpFeM[] = 
  {  7.05387e-02,   2.70513e-06,   1.33741e-10,   7.96280e-15,
     1.06564e-01,   3.93018e-06,   1.88633e-10,   1.09924e-14,
     1.50619e-01,   5.42584e-06,   2.55756e-10,   1.47090e-14,
     2.05371e-01,   7.28371e-06,   3.39206e-10,   1.93364e-14,
     2.90528e-01,   9.90213e-06,   4.47050e-10,   2.49139e-14,
     3.98319e-01,   1.31460e-05,   5.78557e-10,   3.16454e-14,
     5.44029e-01,   1.73143e-05,   7.40062e-10,   3.96165e-14,
     7.41852e-01,   2.27123e-05,   9.40496e-10,   4.91736e-14,
     1.01111e+00,   2.97427e-05,   1.19128e-09,   6.07485e-14,
     1.37135e+00,   3.87934e-05,   1.50307e-09,   7.47433e-14,
     1.85427e+00,   5.04362e-05,   1.88890e-09,   9.15116e-14,
     2.48538e+00,   6.50566e-05,   2.35529e-09,   1.11141e-13,
     3.28338e+00,   8.28882e-05,   2.90473e-09,   1.33599e-13,
     4.26495e+00,   1.04109e-04,   3.53805e-09,   1.58791e-13,
     5.44209e+00,   1.28847e-04,   4.25649e-09,   1.86722e-13,
     6.82357e+00,   1.57201e-04,   5.06195e-09,   2.17471e-13,
     8.41130e+00,   1.89169e-04,   5.95453e-09,   2.51081e-13,
     1.02033e+01,   2.24714e-04,   6.93482e-09,   2.87663e-13,
     1.22164e+01,   2.64256e-04,   8.01879e-09,   3.27999e-13,
     1.45282e+01,   3.09507e-04,   9.26157e-09,   3.74467e-13,
     1.73262e+01,   3.64447e-04,   1.07851e-08,   4.32120e-13,
     2.09720e+01,   4.36715e-04,   1.28220e-08,   5.10543e-13,
     2.61035e+01,   5.39852e-04,   1.57872e-08,   6.26923e-13,
     3.37725e+01,   6.96147e-04,   2.03646e-08,   8.09666e-13,
     4.55764e+01,   9.38083e-04,   2.75201e-08,   1.09797e-12,
     6.37992e+01,   1.30793e-03,   3.84014e-08,   1.53493e-12,
     9.17796e+01,   1.86084e-03,   5.43135e-08,   2.16276e-12,
     1.34749e+02,   2.67751e-03,   7.69998e-08,   3.03160e-12,
     2.01156e+02,   3.88432e-03,   1.09109e-07,   4.21570e-12,
     3.04328e+02,   5.67395e-03,   1.54550e-07,   5.82179e-12,
     4.64551e+02,   8.32593e-03,   2.18657e-07,   7.98599e-12,
     7.11390e+02,   1.22244e-02,   3.08194e-07,   1.08633e-11,
     1.08442e+03,   1.78487e-02,   4.30797e-07,   1.46043e-11,
     1.62442e+03,   2.56300e-02,   5.92006e-07,   1.92777e-11,
     2.34943e+03,   3.56242e-02,   7.89629e-07,   2.47479e-11  };

  double C25_kpMg2SiO4[] = 
  {  1.05240e-01,   2.45433e-07,   2.61677e-12,   4.51929e-17,
     1.32588e-01,   3.09211e-07,   3.29676e-12,   5.69367e-17,
     1.67016e-01,   3.89504e-07,   4.15283e-12,   7.17213e-17,
     2.10360e-01,   4.90585e-07,   5.23055e-12,   9.03341e-17,
     2.71887e-01,   6.34079e-07,   6.76050e-12,   1.16758e-16,
     3.55694e-01,   8.29533e-07,   8.84446e-12,   1.52750e-16,
     4.84933e-01,   1.13094e-06,   1.20582e-11,   2.08253e-16,
     6.99770e-01,   1.63200e-06,   1.74006e-11,   3.00524e-16,
     1.05860e+00,   2.46891e-06,   2.63246e-11,   4.54655e-16,
     1.62903e+00,   3.79938e-06,   4.05116e-11,   6.99697e-16,
     2.54264e+00,   5.93041e-06,   6.32372e-11,   1.09224e-15,
     3.96499e+00,   9.24857e-06,   9.86269e-11,   1.70359e-15,
     6.10655e+00,   1.42452e-05,   1.51927e-10,   2.62446e-15,
     9.28824e+00,   2.16706e-05,   2.31158e-10,   3.99359e-15,
     1.39278e+01,   3.25027e-05,   3.46790e-10,   5.99240e-15,
     2.05413e+01,   4.79530e-05,   5.11839e-10,   8.84697e-15,
     3.00722e+01,   7.02477e-05,   7.50339e-10,   1.29761e-14,
     4.55290e+01,   1.06478e-04,   1.13881e-09,   1.97131e-14,
     7.48333e+01,   1.75301e-04,   1.87829e-09,   3.25565e-14,
     1.29734e+02,   3.04341e-04,   3.26603e-09,   5.66753e-14,
     2.15039e+02,   5.04949e-04,   5.42460e-09,   9.42056e-14,
     3.20373e+02,   7.52945e-04,   8.09632e-09,   1.40698e-13,
     4.30336e+02,   1.01240e-03,   1.08979e-08,   1.89530e-13,
     5.31620e+02,   1.25190e-03,   1.34895e-08,   2.34771e-13,
     6.00659e+02,   1.41538e-03,   1.52608e-08,   2.65718e-13,
     6.07548e+02,   1.43203e-03,   1.54444e-08,   2.68965e-13,
     5.44184e+02,   1.28282e-03,   1.38365e-08,   2.40975e-13,
     4.34292e+02,   1.02387e-03,   1.10444e-08,   1.92359e-13,
     3.13872e+02,   7.40120e-04,   7.98523e-09,   1.39097e-13,
     2.09392e+02,   4.93970e-04,   5.33184e-09,   9.29066e-14,
     1.31415e+02,   3.10355e-04,   3.35361e-09,   5.84810e-14,
     7.92901e+01,   1.87957e-04,   2.03808e-09,   3.56245e-14,
     4.76038e+01,   1.14252e-04,   1.25235e-09,   2.20472e-14,
     3.00283e+01,   7.42825e-05,   8.34278e-10,   1.49154e-14,
     2.10539e+01,   5.48540e-05,   6.37875e-10,   1.16338e-14  };

  double C25_kpMgSiO3[] = 
  {  2.19890e-02,   3.41795e-08,   9.45655e-14,   4.23439e-19,
     3.90612e-02,   6.07164e-08,   1.67986e-13,   7.52197e-19,
     6.05539e-02,   9.41245e-08,   2.60417e-13,   1.16608e-18,
     8.76116e-02,   1.36183e-07,   3.76781e-13,   1.68713e-18,
     1.43288e-01,   2.22725e-07,   6.16221e-13,   2.75928e-18,
     2.19266e-01,   3.40825e-07,   9.42974e-13,   4.22240e-18,
     3.36256e-01,   5.22673e-07,   1.44610e-12,   6.47526e-18,
     5.14336e-01,   7.99479e-07,   2.21195e-12,   9.90458e-18,
     7.97217e-01,   1.23919e-06,   3.42851e-12,   1.53521e-17,
     1.25414e+00,   1.94943e-06,   5.39358e-12,   2.41515e-17,
     2.03450e+00,   3.16241e-06,   8.74964e-12,   3.91798e-17,
     3.34649e+00,   5.20178e-06,   1.43922e-11,   6.44481e-17,
     5.45897e+00,   8.48547e-06,   2.34778e-11,   1.05137e-16,
     8.82126e+00,   1.37119e-05,   3.79391e-11,   1.69905e-16,
     1.41827e+01,   2.20461e-05,   6.10001e-11,   2.73202e-16,
     2.28425e+01,   3.55077e-05,   9.82519e-11,   4.40095e-16,
     3.71193e+01,   5.77019e-05,   1.59676e-10,   7.15351e-16,
     6.14319e+01,   9.54994e-05,   2.64297e-10,   1.18434e-15,
     1.03856e+02,   1.61458e-04,   4.46893e-10,   2.00314e-15,
     1.75529e+02,   2.72897e-04,   7.55426e-10,   3.38696e-15,
     2.82103e+02,   4.38604e-04,   1.21423e-09,   5.44494e-15,
     4.14591e+02,   6.44612e-04,   1.78465e-09,   8.00387e-15,
     5.60087e+02,   8.70868e-04,   2.41120e-09,   1.08148e-14,
     7.11212e+02,   1.10590e-03,   3.06203e-09,   1.37338e-14,
     8.41053e+02,   1.30782e-03,   3.62100e-09,   1.62382e-14,
     8.95593e+02,   1.39263e-03,   3.85557e-09,   1.72860e-14,
     8.40674e+02,   1.30722e-03,   3.61884e-09,   1.62210e-14,
     6.96909e+02,   1.08365e-03,   2.99976e-09,   1.34436e-14,
     5.18364e+02,   8.06018e-04,   2.23111e-09,   9.99751e-15,
     3.52974e+02,   5.48846e-04,   1.51920e-09,   6.80696e-15,
     2.24287e+02,   3.48752e-04,   9.65345e-10,   4.32526e-15,
     1.35182e+02,   2.10204e-04,   5.81862e-10,   2.60719e-15,
     7.83444e+01,   1.21833e-04,   3.37284e-10,   1.51159e-15,
     4.41840e+01,   6.87277e-05,   1.90332e-10,   8.53442e-16,
     2.46537e+01,   3.84248e-05,   1.06643e-10,   4.79266e-16  };

  double C25_kpAC[] = 
  {  3.27960e-01,   2.60233e-07,   1.15896e-12,   3.41207e-17,
     4.38752e-01,   3.48153e-07,   1.55085e-12,   4.56711e-17,
     5.78230e-01,   4.58837e-07,   2.04421e-12,   6.02121e-17,
     7.53824e-01,   5.98180e-07,   2.66530e-12,   7.85179e-17,
     1.04013e+00,   8.25404e-07,   3.67884e-12,   1.08416e-16,
     1.41735e+00,   1.12479e-06,   5.01451e-12,   1.47828e-16,
     1.95293e+00,   1.54986e-06,   6.91189e-12,   2.03855e-16,
     2.71532e+00,   2.15499e-06,   9.61533e-12,   2.83773e-16,
     3.79678e+00,   3.01350e-06,   1.34554e-11,   3.97470e-16,
     5.29747e+00,   4.20498e-06,   1.87919e-11,   5.55745e-16,
     7.37842e+00,   5.85746e-06,   2.62073e-11,   7.76231e-16,
     1.02169e+01,   8.11222e-06,   3.63532e-11,   1.07888e-15,
     1.40424e+01,   1.11520e-05,   5.00793e-11,   1.48986e-15,
     1.92027e+01,   1.52549e-05,   6.86959e-11,   2.04983e-15,
     2.61627e+01,   2.07929e-05,   9.39920e-11,   2.81528e-15,
     3.55327e+01,   2.82565e-05,   1.28407e-10,   3.86530e-15,
     4.81650e+01,   3.83349e-05,   1.75518e-10,   5.31971e-15,
     6.53244e+01,   5.20561e-05,   2.40904e-10,   7.37173e-15,
     8.87808e+01,   7.08733e-05,   3.33003e-10,   1.03268e-14,
     1.20733e+02,   9.66248e-05,   4.63924e-10,   1.46567e-14,
     1.63674e+02,   1.31491e-04,   6.51645e-10,   2.11366e-14,
     2.20677e+02,   1.78321e-04,   9.26445e-10,   3.11929e-14,
     2.96305e+02,   2.41503e-04,   1.34108e-09,   4.74137e-14,
     3.97405e+02,   3.27522e-04,   1.97072e-09,   7.35006e-14,
     5.32189e+02,   4.43695e-04,   2.88200e-09,   1.12542e-13,
     7.08046e+02,   5.95675e-04,   4.08283e-09,   1.64175e-13,
     9.32526e+02,   7.88039e-04,   5.50980e-09,   2.23600e-13,
     1.21879e+03,   1.02941e-03,   7.08050e-09,   2.84124e-13,
     1.59074e+03,   1.33766e-03,   8.76094e-09,   3.40697e-13,
     2.08510e+03,   1.74270e-03,   1.05986e-08,   3.91613e-13,
     2.75362e+03,   2.28917e-03,   1.27282e-08,   4.38283e-13,
     3.66839e+03,   3.04259e-03,   1.53805e-08,   4.84412e-13,
     4.93379e+03,   4.10065e-03,   1.89004e-08,   5.35144e-13,
     6.70658e+03,   5.60645e-03,   2.37440e-08,   5.95776e-13,
     9.22668e+03,   7.75494e-03,   3.03854e-08,   6.69355e-13  };

  double C25_kpSiO2D[] = 
  {  7.60344e-02,   1.95196e-07,   3.66716e-12,   1.92423e-16,
     9.07191e-02,   2.32906e-07,   4.37632e-12,   2.29682e-16,
     1.09206e-01,   2.80380e-07,   5.26909e-12,   2.76586e-16,
     1.32480e-01,   3.40146e-07,   6.39301e-12,   3.35635e-16,
     1.58906e-01,   4.08019e-07,   7.66999e-12,   4.02759e-16,
     1.91564e-01,   4.91897e-07,   9.24810e-12,   4.85715e-16,
     2.30489e-01,   5.91872e-07,   1.11292e-11,   5.84611e-16,
     2.76795e-01,   7.10808e-07,   1.33674e-11,   7.02310e-16,
     3.33075e-01,   8.55378e-07,   1.60886e-11,   8.45443e-16,
     4.05328e-01,   1.04100e-06,   1.95833e-11,   1.02932e-15,
     5.08167e-01,   1.30521e-06,   2.45600e-11,   1.29133e-15,
     6.72485e-01,   1.72749e-06,   3.25210e-11,   1.71095e-15,
     9.48580e-01,   2.43730e-06,   4.59219e-11,   2.41877e-15,
     1.41796e+00,   3.64482e-06,   6.87751e-11,   3.63027e-15,
     2.19521e+00,   5.64613e-06,   1.06783e-10,   5.65494e-15,
     3.46773e+00,   8.92808e-06,   1.69393e-10,   9.00627e-15,
     5.77034e+00,   1.48887e-05,   2.83967e-10,   1.51653e-14,
     1.17273e+01,   3.03917e-05,   5.84644e-10,   3.13500e-14,
     3.16762e+01,   8.23694e-05,   1.59258e-09,   8.53491e-14,
     8.69213e+01,   2.26080e-04,   4.36458e-09,   2.32841e-13,
     1.92500e+02,   5.00291e-04,   9.62913e-09,   5.11437e-13,
     3.36556e+02,   8.73556e-04,   1.67531e-08,   8.86507e-13,
     5.06180e+02,   1.30807e-03,   2.48422e-08,   1.30609e-12,
     7.20631e+02,   1.84289e-03,   3.42245e-08,   1.77586e-12,
     9.76452e+02,   2.46285e-03,   4.43878e-08,   2.26306e-12,
     1.18428e+03,   2.95079e-03,   5.17620e-08,   2.59653e-12,
     1.23545e+03,   3.05120e-03,   5.24657e-08,   2.59989e-12,
     1.10864e+03,   2.72192e-03,   4.61791e-08,   2.26949e-12,
     8.73605e+02,   2.13674e-03,   3.59368e-08,   1.75670e-12,
     6.20107e+02,   1.51303e-03,   2.53058e-08,   1.23288e-12,
     4.05864e+02,   9.88756e-04,   1.64791e-08,   8.01183e-13,
     2.49677e+02,   6.07668e-04,   1.01055e-08,   4.90701e-13,
     1.46540e+02,   3.56444e-04,   5.92003e-09,   2.87265e-13,
     8.30108e+01,   2.01877e-04,   3.35130e-09,   1.62587e-13,
     4.57957e+01,   1.11407e-04,   1.85045e-09,   8.98118e-14  };

  double C25_kpMgO[] = 
  {  2.25388e-04,   8.07807e-10,   6.97998e-15,   9.10286e-20,
     4.04965e-04,   1.45145e-09,   1.25417e-14,   1.63564e-19,
     6.31040e-04,   2.26174e-09,   1.95435e-14,   2.54881e-19,
     9.15651e-04,   3.28184e-09,   2.83582e-14,   3.69842e-19,
     1.52197e-03,   5.45504e-09,   4.71373e-14,   6.14765e-19,
     2.37408e-03,   8.50921e-09,   7.35292e-14,   9.58978e-19,
     3.77210e-03,   1.35201e-08,   1.16831e-13,   1.52374e-18,
     6.14351e-03,   2.20201e-08,   1.90284e-13,   2.48178e-18,
     1.01908e-02,   3.65275e-08,   3.15657e-13,   4.11707e-18,
     1.68899e-02,   6.05411e-08,   5.23195e-13,   6.82425e-18,
     2.96134e-02,   1.06156e-07,   9.17483e-13,   1.19682e-17,
     6.10648e-02,   2.18932e-07,   1.89256e-12,   2.46925e-17,
     1.43413e-01,   5.14257e-07,   4.44656e-12,   5.80289e-17,
     3.27427e-01,   1.17431e-06,   1.01561e-11,   1.32572e-16,
     6.39567e-01,   2.29431e-06,   1.98486e-11,   2.59169e-16,
     1.05188e+00,   3.77698e-06,   3.27186e-11,   4.27780e-16,
     1.56137e+00,   5.64036e-06,   4.92651e-11,   6.49403e-16,
     2.95878e+00,   1.08379e-05,   9.64135e-11,   1.29373e-15,
     1.33369e+01,   4.85546e-05,   4.28194e-10,   5.69626e-15,
     7.10715e+01,   2.55311e-04,   2.21076e-09,   2.88829e-14,
     2.54519e+02,   9.07495e-04,   7.77715e-09,   1.00557e-13,
     5.99024e+02,   2.12736e-03,   1.81309e-08,   2.33128e-13,
     9.91502e+02,   3.51313e-03,   2.98462e-08,   3.82528e-13,
     1.23853e+03,   4.38206e-03,   3.71546e-08,   4.75243e-13,
     1.24134e+03,   4.38799e-03,   3.71575e-08,   4.74666e-13,
     1.05053e+03,   3.71132e-03,   3.14012e-08,   4.00789e-13,
     7.82169e+02,   2.76204e-03,   2.33559e-08,   2.97932e-13,
     5.28995e+02,   1.86748e-03,   1.57854e-08,   2.01282e-13,
     3.32978e+02,   1.17527e-03,   9.93166e-09,   1.26606e-13,
     1.98632e+02,   7.00988e-04,   5.92262e-09,   7.54858e-14,
     1.13813e+02,   4.01617e-04,   3.39280e-09,   4.32365e-14,
     6.32536e+01,   2.23189e-04,   1.88528e-09,   2.40230e-14,
     3.43507e+01,   1.21200e-04,   1.02371e-09,   1.30436e-14,
     1.83277e+01,   6.46641e-05,   5.46160e-10,   6.95860e-15,
     9.64683e+00,   3.40355e-05,   2.87461e-10,   3.66245e-15  };

  double C25_kpFeS[] = 
  {  5.18089e-02,   4.97944e-08,   1.27767e-13,   7.39409e-19,
     9.98898e-02,   9.60047e-08,   2.46329e-13,   1.42543e-18,
     1.60420e-01,   1.54180e-07,   3.95589e-13,   2.28909e-18,
     2.36623e-01,   2.27418e-07,   5.83496e-13,   3.37637e-18,
     3.67289e-01,   3.53003e-07,   9.05730e-13,   5.24118e-18,
     5.36230e-01,   5.15376e-07,   1.32237e-12,   7.65247e-18,
     7.64209e-01,   7.34494e-07,   1.88465e-12,   1.09071e-17,
     1.04972e+00,   1.00891e-06,   2.58893e-12,   1.49849e-17,
     1.38083e+00,   1.32717e-06,   3.40583e-12,   1.97167e-17,
     1.74377e+00,   1.67607e-06,   4.30174e-12,   2.49099e-17,
     2.10317e+00,   2.02159e-06,   5.18938e-12,   3.00615e-17,
     2.42155e+00,   2.32777e-06,   5.97696e-12,   3.46445e-17,
     2.66826e+00,   2.56519e-06,   6.58920e-12,   3.82277e-17,
     2.81642e+00,   2.70807e-06,   6.96112e-12,   4.04503e-17,
     2.89877e+00,   2.78835e-06,   7.17920e-12,   4.18739e-17,
     3.08081e+00,   2.96627e-06,   7.66777e-12,   4.51307e-17,
     3.69596e+00,   3.56434e-06,   9.27486e-12,   5.54010e-17,
     5.05776e+00,   4.88601e-06,   1.28025e-11,   7.76522e-17,
     7.07898e+00,   6.85096e-06,   1.80991e-11,   1.11855e-16,
     9.20178e+00,   8.93168e-06,   2.39573e-11,   1.53397e-16,
     1.08226e+01,   1.05533e-05,   2.89938e-11,   1.96013e-16,
     1.16924e+01,   1.14675e-05,   3.24020e-11,   2.32497e-16,
     1.19355e+01,   1.17825e-05,   3.42098e-11,   2.58640e-16,
     1.18230e+01,   1.17613e-05,   3.49907e-11,   2.75275e-16,
     1.15929e+01,   1.16487e-05,   3.54453e-11,   2.86649e-16,
     1.13868e+01,   1.16310e-05,   3.63433e-11,   2.99851e-16,
     1.14199e+01,   1.23167e-05,   4.12834e-11,   3.50696e-16,
     1.29371e+01,   1.75436e-05,   7.43145e-11,   6.75730e-16,
     2.05440e+01,   3.98903e-05,   2.16132e-10,   2.04334e-15,
     4.78755e+01,   1.10160e-04,   6.39249e-10,   5.90057e-15,
     1.21252e+02,   2.78758e-04,   1.57840e-09,   1.39121e-14,
     2.69644e+02,   5.88705e-04,   3.18490e-09,   2.68674e-14,
     5.11453e+02,   1.05406e-03,   5.44914e-09,   4.43131e-14,
     8.53209e+02,   1.66585e-03,   8.25940e-09,   6.51454e-14,
     1.29539e+03,   2.40455e-03,   1.14623e-08,   8.80401e-14  };

  double C25_kpAl2O3[] = 
  {  9.93250e-04,   1.98179e-11,   3.95420e-19,   7.88967e-27,
     1.81240e-03,   3.61621e-11,   7.21529e-19,   1.43964e-26,
     2.84365e-03,   5.67382e-11,   1.13208e-18,   2.25879e-26,
     4.14191e-03,   8.26420e-11,   1.64892e-18,   3.29004e-26,
     7.18271e-03,   1.43314e-10,   2.85949e-18,   5.70543e-26,
     1.13364e-02,   2.26190e-10,   4.51309e-18,   9.00479e-26,
     1.77361e-02,   3.53881e-10,   7.06085e-18,   1.40883e-25,
     2.59477e-02,   5.17725e-10,   1.03300e-17,   2.06110e-25,
     3.45425e-02,   6.89214e-10,   1.37516e-17,   2.74381e-25,
     4.22006e-02,   8.42014e-10,   1.68004e-17,   3.35212e-25,
     4.71420e-02,   9.40607e-10,   1.87676e-17,   3.74462e-25,
     4.91934e-02,   9.81537e-10,   1.95842e-17,   3.90757e-25,
     5.05162e-02,   1.00793e-09,   2.01109e-17,   4.01264e-25,
     5.78201e-02,   1.15366e-09,   2.30186e-17,   4.59282e-25,
     8.84237e-02,   1.76428e-09,   3.52021e-17,   7.02374e-25,
     1.78786e-01,   3.56725e-09,   7.11761e-17,   1.42015e-24,
     4.36404e-01,   8.70740e-09,   1.73736e-16,   3.46648e-24,
     1.63796e+00,   3.26816e-08,   6.52083e-16,   1.30108e-23,
     8.50817e+00,   1.69760e-07,   3.38716e-15,   6.75828e-23,
     3.92751e+01,   7.83641e-07,   1.56357e-14,   3.11973e-22,
     1.41433e+02,   2.82196e-06,   5.63055e-14,   1.12344e-21,
     3.83709e+02,   7.65599e-06,   1.52757e-13,   3.04791e-21,
     7.70411e+02,   1.53717e-05,   3.06706e-13,   6.11959e-21,
     1.16399e+03,   2.32246e-05,   4.63392e-13,   9.24589e-21,
     1.37566e+03,   2.74481e-05,   5.47662e-13,   1.09273e-20,
     1.33070e+03,   2.65509e-05,   5.29761e-13,   1.05701e-20,
     1.09978e+03,   2.19435e-05,   4.37830e-13,   8.73585e-21,
     8.05638e+02,   1.60746e-05,   3.20730e-13,   6.39941e-21,
     5.38690e+02,   1.07483e-05,   2.14456e-13,   4.27897e-21,
     3.36338e+02,   6.71083e-06,   1.33899e-13,   2.67163e-21,
     1.99460e+02,   3.97975e-06,   7.94065e-14,   1.58437e-21,
     1.13787e+02,   2.27035e-06,   4.52995e-14,   9.03844e-22,
     6.30411e+01,   1.25784e-06,   2.50971e-14,   5.00753e-22,
     3.41529e+01,   6.81441e-07,   1.35965e-14,   2.71286e-22,
     1.81893e+01,   3.62924e-07,   7.24128e-15,   1.44483e-22  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = C25_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = C25_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = C25_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = C25_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = C25_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = C25_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = C25_kpMgO     [itab];
      my_rates->SN0_kpFeS     [itab0] = C25_kpFeS     [itab];
      my_rates->SN0_kpAl2O3   [itab0] = C25_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_C30(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   4.99965e-02;
  my_rates->SN0_XO [iSN] =   7.32832e-01;
  my_rates->SN0_XMg[iSN] =   3.87430e-02;
  my_rates->SN0_XAl[iSN] =   8.61678e-04;
  my_rates->SN0_XSi[iSN] =   7.18810e-02;
  my_rates->SN0_XS [iSN] =   3.70455e-02;
  my_rates->SN0_XFe[iSN] =   1.45822e-02;

  my_rates->SN0_fC [iSN] =   4.93773e-02;
  my_rates->SN0_fO [iSN] =   7.29130e-01;
  my_rates->SN0_fMg[iSN] =   3.76731e-02;
  my_rates->SN0_fAl[iSN] =   8.61678e-04;
  my_rates->SN0_fSi[iSN] =   4.01269e-02;
  my_rates->SN0_fS [iSN] =   3.68812e-02;
  my_rates->SN0_fFe[iSN] =   1.23641e-02;

  my_rates->SN0_fSiM     [iSN] =   2.91389e-02;
  my_rates->SN0_fFeM     [iSN] =   1.93065e-03;
  my_rates->SN0_fMg2SiO4 [iSN] =   7.73041e-04;
  my_rates->SN0_fMgSiO3  [iSN] =   4.17376e-06;
  my_rates->SN0_fAC      [iSN] =   6.19235e-04;
  my_rates->SN0_fSiO2D   [iSN] =   5.27016e-03;
  my_rates->SN0_fMgO     [iSN] =   1.33978e-03;
  my_rates->SN0_fFeS     [iSN] =   4.51744e-04;
  my_rates->SN0_fAl2O3   [iSN] =   5.79251e-12;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   2.56305e-05;
  my_rates->SN0_r0FeM     [itab0 + 0] =   2.05800e-05;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   4.70227e-07;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   1.02156e-06;
  my_rates->SN0_r0AC      [itab0 + 0] =   1.17005e-06;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   1.62875e-06;
  my_rates->SN0_r0MgO     [itab0 + 0] =   2.32229e-06;
  my_rates->SN0_r0FeS     [itab0 + 0] =   1.69769e-06;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   7.63588e-08;

  my_rates->SN0_r0SiM     [itab0 + 1] =   1.02092e-09;
  my_rates->SN0_r0FeM     [itab0 + 1] =   5.92424e-10;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   6.29420e-13;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.31765e-12;
  my_rates->SN0_r0AC      [itab0 + 1] =   2.37154e-12;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   1.12314e-11;
  my_rates->SN0_r0MgO     [itab0 + 1] =   1.39783e-11;
  my_rates->SN0_r0FeS     [itab0 + 1] =   6.40794e-12;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   8.88224e-15;

  my_rates->SN0_r0SiM     [itab0 + 2] =   5.78476e-14;
  my_rates->SN0_r0FeM     [itab0 + 2] =   2.26690e-14;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   1.71079e-18;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   2.63083e-18;
  my_rates->SN0_r0AC      [itab0 + 2] =   7.59875e-18;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   1.91031e-16;
  my_rates->SN0_r0MgO     [itab0 + 2] =   1.49800e-16;
  my_rates->SN0_r0FeS     [itab0 + 2] =   4.40126e-17;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   1.42247e-21;

  NTd =            35;
 Nmom =             4;

  double C30_kpSiM[] = 
  {  1.52613e-01,   3.87036e-06,   1.51475e-10,   8.34686e-15,
     1.92410e-01,   4.88554e-06,   1.91604e-10,   1.05937e-14,
     2.42503e-01,   6.16264e-06,   2.42001e-10,   1.34057e-14,
     3.05564e-01,   7.77004e-06,   3.05397e-10,   1.69392e-14,
     3.85174e-01,   9.80462e-06,   3.86046e-10,   2.14732e-14,
     4.85438e-01,   1.23675e-05,   4.87634e-10,   2.71828e-14,
     6.11833e-01,   1.56019e-05,   6.16102e-10,   3.44270e-14,
     7.71157e-01,   1.96838e-05,   7.78585e-10,   4.36231e-14,
     9.71669e-01,   2.48269e-05,   9.83743e-10,   5.52763e-14,
     1.22169e+00,   3.12457e-05,   1.24021e-09,   6.98829e-14,
     1.51883e+00,   3.88844e-05,   1.54619e-09,   8.73862e-14,
     1.83610e+00,   4.70620e-05,   1.87584e-09,   1.06471e-13,
     2.15867e+00,   5.54249e-05,   2.21779e-09,   1.26800e-13,
     2.56252e+00,   6.59468e-05,   2.65228e-09,   1.53070e-13,
     3.24657e+00,   8.37697e-05,   3.38561e-09,   1.97020e-13,
     4.37158e+00,   1.13090e-04,   4.58835e-09,   2.68525e-13,
     5.89356e+00,   1.52861e-04,   6.22099e-09,   3.65325e-13,
     7.62770e+00,   1.98484e-04,   8.10272e-09,   4.76950e-13,
     9.48448e+00,   2.48132e-04,   1.01763e-08,   6.00571e-13,
     1.17614e+01,   3.11117e-04,   1.28785e-08,   7.63878e-13,
     1.52604e+01,   4.11692e-04,   1.73236e-08,   1.03683e-12,
     2.04857e+01,   5.67120e-04,   2.43755e-08,   1.47603e-12,
     2.67121e+01,   7.65926e-04,   3.38756e-08,   2.08532e-12,
     3.30214e+01,   1.00600e-03,   4.67023e-08,   2.95970e-12,
     3.94796e+01,   1.32148e-03,   6.59032e-08,   4.35947e-12,
     4.60630e+01,   1.71870e-03,   9.26057e-08,   6.41018e-12,
     5.17888e+01,   2.12529e-03,   1.22060e-07,   8.76593e-12,
     5.52863e+01,   2.42692e-03,   1.45553e-07,   1.07154e-11,
     5.60930e+01,   2.55895e-03,   1.57221e-07,   1.17347e-11,
     5.70683e+01,   2.61927e-03,   1.61476e-07,   1.20713e-11,
     7.25749e+01,   3.09971e-03,   1.82166e-07,   1.32162e-11,
     1.59194e+02,   5.49929e-03,   2.77545e-07,   1.82852e-11,
     4.83739e+02,   1.31658e-02,   5.43977e-07,   3.13196e-11,
     1.39550e+03,   3.19882e-02,   1.11499e-06,   5.65773e-11,
     3.29912e+03,   6.76125e-02,   2.08154e-06,   9.54501e-11  };

  double C30_kpFeM[] = 
  {  7.21495e-02,   2.58300e-06,   1.19134e-10,   6.77090e-15,
     1.09797e-01,   3.79107e-06,   1.69411e-10,   9.39020e-15,
     1.55805e-01,   5.26522e-06,   2.30859e-10,   1.26025e-14,
     2.12966e-01,   7.09577e-06,   3.07229e-10,   1.66015e-14,
     3.02975e-01,   9.74053e-06,   4.08459e-10,   2.15040e-14,
     4.17015e-01,   1.30295e-05,   5.32431e-10,   2.74392e-14,
     5.71703e-01,   1.73021e-05,   6.86687e-10,   3.45374e-14,
     7.82114e-01,   2.28860e-05,   8.80449e-10,   4.31318e-14,
     1.06869e+00,   3.02137e-05,   1.12556e-09,   5.36402e-14,
     1.45194e+00,   3.97008e-05,   1.43310e-09,   6.64513e-14,
     1.96519e+00,   5.19736e-05,   1.81745e-09,   8.19495e-14,
     2.63479e+00,   6.74573e-05,   2.28640e-09,   1.00267e-13,
     3.47952e+00,   8.64064e-05,   2.84327e-09,   1.21406e-13,
     4.51582e+00,   1.09014e-04,   3.48958e-09,   1.45307e-13,
     5.75506e+00,   1.35409e-04,   4.22673e-09,   1.71980e-13,
     7.20522e+00,   1.65686e-04,   5.05643e-09,   2.01492e-13,
     8.86719e+00,   1.99826e-04,   5.97834e-09,   2.33868e-13,
     1.07378e+01,   2.37775e-04,   6.99239e-09,   2.69183e-13,
     1.28337e+01,   2.79956e-04,   8.11382e-09,   3.08137e-13,
     1.52347e+01,   3.28163e-04,   9.39756e-09,   3.52929e-13,
     1.81342e+01,   3.86588e-04,   1.09665e-08,   4.08289e-13,
     2.19061e+01,   4.63288e-04,   1.30553e-08,   4.83202e-13,
     2.72090e+01,   5.72531e-04,   1.60822e-08,   5.93756e-13,
     3.51270e+01,   7.37785e-04,   2.07358e-08,   7.66520e-13,
     4.72945e+01,   9.93265e-04,   2.79935e-08,   1.03837e-12,
     6.60183e+01,   1.38363e-03,   3.90367e-08,   1.45073e-12,
     9.46224e+01,   1.96740e-03,   5.52478e-08,   2.04604e-12,
     1.38263e+02,   2.83025e-03,   7.85071e-08,   2.87658e-12,
     2.05210e+02,   4.10585e-03,   1.11674e-07,   4.02015e-12,
     3.08399e+02,   5.99700e-03,   1.58965e-07,   5.58896e-12,
     4.67338e+02,   8.79616e-03,   2.26155e-07,   7.72827e-12,
     7.10129e+02,   1.29022e-02,   3.20597e-07,   1.06077e-11,
     1.07381e+03,   1.88070e-02,   4.50609e-07,   1.43975e-11,
     1.59524e+03,   2.69407e-02,   6.22216e-07,   1.91854e-11,
     2.28747e+03,   3.73289e-02,   8.32996e-07,   2.48404e-11  };

  double C30_kpMg2SiO4[] = 
  {  1.05240e-01,   4.94867e-08,   6.62401e-14,   1.80043e-19,
     1.32588e-01,   6.23464e-08,   8.34533e-14,   2.26830e-19,
     1.67016e-01,   7.85357e-08,   1.05123e-13,   2.85730e-19,
     2.10360e-01,   9.89168e-08,   1.32405e-13,   3.59881e-19,
     2.71887e-01,   1.27849e-07,   1.71131e-13,   4.65142e-19,
     3.55694e-01,   1.67257e-07,   2.23881e-13,   6.08518e-19,
     4.84932e-01,   2.28028e-07,   3.05226e-13,   8.29619e-19,
     6.99767e-01,   3.29050e-07,   4.40448e-13,   1.19716e-18,
     1.05860e+00,   4.97781e-07,   6.66304e-13,   1.81105e-18,
     1.62902e+00,   7.66009e-07,   1.02534e-12,   2.78694e-18,
     2.54260e+00,   1.19560e-06,   1.60037e-12,   4.34992e-18,
     3.96488e+00,   1.86440e-06,   2.49561e-12,   6.78334e-18,
     6.10630e+00,   2.87136e-06,   3.84349e-12,   1.04472e-17,
     9.28767e+00,   4.36734e-06,   5.84598e-12,   1.58906e-17,
     1.39265e+01,   6.54868e-06,   8.76595e-12,   2.38284e-17,
     2.05383e+01,   9.65780e-06,   1.29279e-11,   3.51435e-17,
     3.00651e+01,   1.41377e-05,   1.89253e-11,   5.14515e-17,
     4.55105e+01,   2.14011e-05,   2.86497e-11,   7.79011e-17,
     7.47848e+01,   3.51681e-05,   4.70833e-11,   1.28054e-16,
     1.29623e+02,   6.09573e-05,   8.16150e-11,   2.22014e-16,
     2.14824e+02,   1.01026e-04,   1.35269e-10,   3.68020e-16,
     3.20010e+02,   1.50495e-04,   2.01516e-10,   5.48326e-16,
     4.29781e+02,   2.02124e-04,   2.70665e-10,   7.36593e-16,
     5.30850e+02,   2.49665e-04,   3.34350e-10,   9.10050e-16,
     5.99723e+02,   2.82064e-04,   3.77759e-10,   1.02832e-15,
     6.06569e+02,   2.85288e-04,   3.82090e-10,   1.04017e-15,
     5.43292e+02,   2.55531e-04,   3.42243e-10,   9.31723e-16,
     4.33571e+02,   2.03926e-04,   2.73131e-10,   7.43589e-16,
     3.13340e+02,   1.47377e-04,   1.97394e-10,   5.37419e-16,
     2.09022e+02,   9.83146e-05,   1.31685e-10,   3.58547e-16,
     1.31159e+02,   6.16930e-05,   8.26409e-11,   2.25056e-16,
     7.90789e+01,   3.72049e-05,   4.98596e-11,   1.35886e-16,
     4.73566e+01,   2.23007e-05,   2.99370e-11,   8.18168e-17,
     2.96727e+01,   1.40105e-05,   1.88993e-11,   5.20353e-17,
     2.05268e+01,   9.75253e-06,   1.32995e-11,   3.71721e-17  };

  double C30_kpMgSiO3[] = 
  {  2.19890e-02,   2.24631e-08,   2.89738e-14,   5.78493e-20,
     3.90612e-02,   3.99034e-08,   5.14691e-14,   1.02764e-19,
     6.05539e-02,   6.18594e-08,   7.97889e-14,   1.59308e-19,
     8.76116e-02,   8.95004e-08,   1.15442e-13,   2.30492e-19,
     1.43288e-01,   1.46377e-07,   1.88803e-13,   3.76970e-19,
     2.19266e-01,   2.23993e-07,   2.88916e-13,   5.76861e-19,
     3.36256e-01,   3.43505e-07,   4.43068e-13,   8.84648e-19,
     5.14336e-01,   5.25424e-07,   6.77716e-13,   1.35317e-18,
     7.97217e-01,   8.14404e-07,   1.05046e-12,   2.09745e-18,
     1.25414e+00,   1.28118e-06,   1.65253e-12,   3.29971e-18,
     2.03450e+00,   2.07836e-06,   2.68078e-12,   5.35315e-18,
     3.34648e+00,   3.41863e-06,   4.40956e-12,   8.80595e-18,
     5.45894e+00,   5.57665e-06,   7.19317e-12,   1.43665e-17,
     8.82120e+00,   9.01141e-06,   1.16237e-11,   2.32195e-17,
     1.41826e+01,   1.44884e-05,   1.86888e-11,   3.73434e-17,
     2.28421e+01,   2.33348e-05,   3.01009e-11,   6.01739e-17,
     3.71183e+01,   3.79191e-05,   4.89167e-11,   9.78570e-17,
     6.14292e+01,   6.27552e-05,   8.09620e-11,   1.62128e-16,
     1.03850e+02,   1.06093e-04,   1.36885e-10,   2.74448e-16,
     1.75514e+02,   1.79307e-04,   2.31368e-10,   4.64387e-16,
     2.82073e+02,   2.88172e-04,   3.71865e-10,   7.46941e-16,
     4.14542e+02,   4.23507e-04,   5.46530e-10,   1.09837e-15,
     5.60007e+02,   5.72126e-04,   7.38347e-10,   1.48437e-15,
     7.11091e+02,   7.26489e-04,   9.37565e-10,   1.88456e-15,
     8.40894e+02,   8.59108e-04,   1.10867e-09,   2.22674e-15,
     8.95416e+02,   9.14814e-04,   1.18049e-09,   2.36840e-15,
     8.40506e+02,   8.58713e-04,   1.10802e-09,   2.22074e-15,
     6.96770e+02,   7.11863e-04,   9.18488e-10,   1.83938e-15,
     5.18262e+02,   5.29487e-04,   6.83149e-10,   1.36732e-15,
     3.52905e+02,   3.60548e-04,   4.65172e-10,   9.30726e-16,
     2.24242e+02,   2.29099e-04,   2.95579e-10,   5.91332e-16,
     1.35153e+02,   1.38081e-04,   1.78153e-10,   3.56469e-16,
     7.83239e+01,   8.00228e-05,   1.03253e-10,   2.06735e-16,
     4.41660e+01,   4.51272e-05,   5.82391e-11,   1.16797e-16,
     2.46140e+01,   2.51643e-05,   3.25109e-11,   6.54871e-17  };

  double C30_kpAC[] = 
  {  3.27960e-01,   3.83729e-07,   7.77768e-13,   2.49208e-18,
     4.38752e-01,   5.13360e-07,   1.04052e-12,   3.33400e-18,
     5.78230e-01,   6.76557e-07,   1.37130e-12,   4.39392e-18,
     7.53823e-01,   8.82009e-07,   1.78773e-12,   5.72828e-18,
     1.04013e+00,   1.21701e-06,   2.46677e-12,   7.90427e-18,
     1.41736e+00,   1.65839e-06,   3.36142e-12,   1.07712e-17,
     1.95293e+00,   2.28504e-06,   4.63158e-12,   1.48415e-17,
     2.71531e+00,   3.17707e-06,   6.43967e-12,   2.06359e-17,
     3.79677e+00,   4.44245e-06,   9.00459e-12,   2.88564e-17,
     5.29746e+00,   6.19839e-06,   1.25639e-11,   4.02650e-17,
     7.37839e+00,   8.63325e-06,   1.74996e-11,   5.60866e-17,
     1.02169e+01,   1.19546e-05,   2.42326e-11,   7.76741e-17,
     1.40423e+01,   1.64308e-05,   3.33069e-11,   1.06774e-16,
     1.92025e+01,   2.24692e-05,   4.55495e-11,   1.46048e-16,
     2.61625e+01,   3.06137e-05,   6.20638e-11,   1.99051e-16,
     3.55322e+01,   4.15788e-05,   8.43006e-11,   2.70467e-16,
     4.81640e+01,   5.63626e-05,   1.14288e-10,   3.66871e-16,
     6.53224e+01,   7.64461e-05,   1.55039e-10,   4.98051e-16,
     8.87770e+01,   1.03904e-04,   2.10775e-10,   6.77818e-16,
     1.20724e+02,   1.41310e-04,   2.86753e-10,   9.23530e-16,
     1.63658e+02,   1.91596e-04,   3.88991e-10,   1.25555e-15,
     2.20645e+02,   2.58372e-04,   5.24942e-10,   1.69987e-15,
     2.96236e+02,   3.47007e-04,   7.05751e-10,   2.29590e-15,
     3.97277e+02,   4.65581e-04,   9.48153e-10,   3.10191e-15,
     5.32000e+02,   6.23846e-04,   1.27229e-09,   4.18567e-15,
     7.07873e+02,   8.30738e-04,   1.69662e-09,   5.60591e-15,
     9.32623e+02,   1.09569e-03,   2.24083e-09,   7.42175e-15,
     1.21977e+03,   1.43532e-03,   2.93992e-09,   9.74128e-15,
     1.59377e+03,   1.87988e-03,   3.85840e-09,   1.27719e-14,
     2.09261e+03,   2.47712e-03,   5.09977e-09,   1.68549e-14,
     2.77023e+03,   3.29669e-03,   6.81899e-09,   2.25141e-14,
     3.70299e+03,   4.44048e-03,   9.24971e-09,   3.05579e-14,
     5.00261e+03,   6.06097e-03,   1.27486e-08,   4.22365e-14,
     6.83605e+03,   8.38300e-03,   1.78358e-08,   5.93501e-14,
     9.45120e+03,   1.17154e-02,   2.51744e-08,   8.40301e-14  };

  double C30_kpSiO2D[] = 
  {  7.60354e-02,   1.23833e-07,   8.53754e-13,   1.45191e-17,
     9.07201e-02,   1.47751e-07,   1.01868e-12,   1.73242e-17,
     1.09207e-01,   1.77861e-07,   1.22631e-12,   2.08557e-17,
     1.32481e-01,   2.15768e-07,   1.48770e-12,   2.53015e-17,
     1.58907e-01,   2.58811e-07,   1.78454e-12,   3.03507e-17,
     1.91564e-01,   3.12003e-07,   2.15137e-12,   3.65904e-17,
     2.30490e-01,   3.75405e-07,   2.58860e-12,   4.40274e-17,
     2.76795e-01,   4.50827e-07,   3.10875e-12,   5.28751e-17,
     3.33074e-01,   5.42499e-07,   3.74099e-12,   6.36301e-17,
     4.05326e-01,   6.60191e-07,   4.55275e-12,   7.74389e-17,
     5.08163e-01,   8.27703e-07,   5.70817e-12,   9.70950e-17,
     6.72477e-01,   1.09537e-06,   7.55468e-12,   1.28511e-16,
     9.48561e-01,   1.54515e-06,   1.06581e-11,   1.81319e-16,
     1.41791e+00,   2.30989e-06,   1.59363e-11,   2.71154e-16,
     2.19508e+00,   3.57636e-06,   2.46813e-11,   4.20045e-16,
     3.46738e+00,   5.65063e-06,   3.90198e-11,   6.64369e-16,
     5.76921e+00,   9.40776e-06,   6.50678e-11,   1.10919e-15,
     1.17225e+01,   1.91430e-05,   1.32874e-10,   2.27105e-15,
     3.16576e+01,   5.17610e-05,   3.60367e-10,   6.17301e-15,
     8.68678e+01,   1.42065e-04,   9.89621e-10,   1.69586e-14,
     1.92388e+02,   3.14591e-04,   2.19062e-09,   3.75294e-14,
     3.36374e+02,   5.49866e-04,   3.82570e-09,   6.54980e-14,
     5.05999e+02,   8.26100e-04,   5.72796e-09,   9.78085e-14,
     7.20701e+02,   1.17289e-03,   8.06418e-09,   1.36816e-13,
     9.77142e+02,   1.58336e-03,   1.07627e-08,   1.81003e-13,
     1.18576e+03,   1.91391e-03,   1.28756e-08,   2.14810e-13,
     1.23749e+03,   1.99168e-03,   1.32966e-08,   2.20518e-13,
     1.11076e+03,   1.78423e-03,   1.18503e-08,   1.95740e-13,
     8.75424e+02,   1.40444e-03,   9.29642e-09,   1.53150e-13,
     6.21466e+02,   9.96201e-04,   6.57984e-09,   1.08212e-13,
     4.06780e+02,   6.51722e-04,   4.29860e-09,   7.06171e-14,
     2.50252e+02,   4.00808e-04,   2.64130e-09,   4.33611e-14,
     1.46882e+02,   2.35201e-04,   1.54913e-09,   2.54207e-14,
     8.32048e+01,   1.33228e-04,   8.77332e-10,   1.43944e-14,
     4.59015e+01,   7.35091e-05,   4.84195e-10,   7.94538e-15  };

  double C30_kpMgO[] = 
  {  2.25389e-04,   5.23412e-10,   3.15038e-15,   3.37599e-20,
     4.04967e-04,   9.40443e-10,   5.66055e-15,   6.06603e-20,
     6.31042e-04,   1.46545e-09,   8.82067e-15,   9.45260e-20,
     9.15653e-04,   2.12641e-09,   1.27990e-14,   1.37160e-19,
     1.52197e-03,   3.53447e-09,   2.12745e-14,   2.27991e-19,
     2.37407e-03,   5.51331e-09,   3.31857e-14,   3.55642e-19,
     3.77209e-03,   8.75995e-09,   5.27284e-14,   5.65083e-19,
     6.14348e-03,   1.42671e-08,   8.58786e-14,   9.20363e-19,
     1.01907e-02,   2.36663e-08,   1.42459e-13,   1.52677e-18,
     1.68897e-02,   3.92242e-08,   2.36116e-13,   2.53063e-18,
     2.96125e-02,   6.87744e-08,   4.14030e-13,   4.43787e-18,
     6.10613e-02,   1.41824e-07,   8.53932e-13,   9.15481e-18,
     1.43400e-01,   3.33097e-07,   2.00598e-12,   2.15106e-17,
     3.27386e-01,   7.60540e-07,   4.58100e-12,   4.91345e-17,
     6.39460e-01,   1.48568e-06,   8.95096e-12,   9.60338e-17,
     1.05151e+00,   2.44421e-06,   1.47414e-11,   1.58361e-16,
     1.55897e+00,   3.63523e-06,   2.20702e-11,   2.39002e-16,
     2.94595e+00,   6.92078e-06,   4.26541e-11,   4.70201e-16,
     1.32939e+01,   3.11374e-05,   1.90589e-10,   2.08296e-15,
     7.10280e+01,   1.65213e-04,   9.96517e-10,   1.06988e-14,
     2.54734e+02,   5.90224e-04,   3.53091e-09,   3.75292e-14,
     5.99997e+02,   1.38734e-03,   8.26333e-09,   8.73603e-14,
     9.93558e+02,   2.29460e-03,   1.36329e-08,   1.43684e-13,
     1.24146e+03,   2.86493e-03,   1.69946e-08,   1.78772e-13,
     1.24450e+03,   2.87058e-03,   1.70110e-08,   1.78725e-13,
     1.05332e+03,   2.42888e-03,   1.43842e-08,   1.51004e-13,
     7.84315e+02,   1.80814e-03,   1.07031e-08,   1.12298e-13,
     5.30478e+02,   1.22276e-03,   7.23576e-09,   7.58903e-14,
     3.33924e+02,   7.69623e-04,   4.55334e-09,   4.77443e-14,
     1.99202e+02,   4.59084e-04,   2.71568e-09,   2.84703e-14,
     1.14141e+02,   2.63040e-04,   1.55583e-09,   1.63087e-14,
     6.34373e+01,   1.46185e-04,   8.64588e-10,   9.06207e-15,
     3.44508e+01,   7.93863e-05,   4.69493e-10,   4.92062e-15,
     1.83812e+01,   4.23560e-05,   2.50488e-10,   2.62517e-15,
     9.67501e+00,   2.22941e-05,   1.31841e-10,   1.38170e-15  };

  double C30_kpFeS[] = 
  {  5.18102e-02,   8.79700e-08,   3.32172e-13,   2.28308e-18,
     9.98920e-02,   1.69607e-07,   6.40403e-13,   4.40131e-18,
     1.60423e-01,   2.72381e-07,   1.02844e-12,   7.06800e-18,
     2.36628e-01,   4.01766e-07,   1.51696e-12,   1.04252e-17,
     3.67296e-01,   6.23631e-07,   2.35471e-12,   1.61832e-17,
     5.36243e-01,   9.10491e-07,   3.43791e-12,   2.36285e-17,
     7.64229e-01,   1.29761e-06,   4.89979e-12,   3.36782e-17,
     1.04975e+00,   1.78244e-06,   6.73094e-12,   4.62695e-17,
     1.38087e+00,   2.34473e-06,   8.85504e-12,   6.08809e-17,
     1.74385e+00,   2.96123e-06,   1.11849e-11,   7.69178e-17,
     2.10328e+00,   3.57180e-06,   1.34937e-11,   9.28279e-17,
     2.42172e+00,   4.11304e-06,   1.55431e-11,   1.06984e-16,
     2.66854e+00,   4.53297e-06,   1.71379e-11,   1.18057e-16,
     2.81685e+00,   4.78622e-06,   1.81100e-11,   1.24934e-16,
     2.89955e+00,   4.92990e-06,   1.86888e-11,   1.29362e-16,
     3.08249e+00,   5.24916e-06,   1.99900e-11,   1.39496e-16,
     3.69977e+00,   6.31698e-06,   2.42367e-11,   1.71353e-16,
     5.06558e+00,   8.67290e-06,   3.35347e-11,   2.40280e-16,
     7.09332e+00,   1.21816e-05,   4.75377e-11,   3.46171e-16,
     9.22649e+00,   1.59270e-05,   6.32195e-11,   4.74390e-16,
     1.08624e+01,   1.89007e-05,   7.70051e-11,   6.04198e-16,
     1.17514e+01,   2.06434e-05,   8.66021e-11,   7.12327e-16,
     1.20172e+01,   2.13237e-05,   9.19108e-11,   7.86989e-16,
     1.19350e+01,   2.14078e-05,   9.44200e-11,   8.32977e-16,
     1.17499e+01,   2.13514e-05,   9.60658e-11,   8.64823e-16,
     1.16255e+01,   2.15518e-05,   9.90976e-11,   9.04821e-16,
     1.19573e+01,   2.36139e-05,   1.14559e-10,   1.06630e-15,
     1.51193e+01,   3.80077e-05,   2.17630e-10,   2.11441e-15,
     2.88990e+01,   9.93223e-05,   6.64795e-10,   6.60023e-15,
     7.41785e+01,   2.87865e-04,   1.99265e-09,   1.93332e-14,
     1.88656e+02,   7.26466e-04,   4.90072e-09,   4.57065e-14,
     4.09136e+02,   1.50900e-03,   9.80616e-09,   8.81118e-14,
     7.53739e+02,   2.65184e-03,   1.66314e-08,   1.44869e-13,
     1.22335e+03,   4.11549e-03,   2.50007e-08,   2.12244e-13,
     1.81062e+03,   5.83587e-03,   3.44214e-08,   2.85851e-13  };

  double C30_kpAl2O3[] = 
  {  9.93250e-04,   7.58434e-11,   8.82228e-18,   1.41287e-24,
     1.81240e-03,   1.38393e-10,   1.60982e-17,   2.57809e-24,
     2.84365e-03,   2.17138e-10,   2.52580e-17,   4.04502e-24,
     4.14191e-03,   3.16271e-10,   3.67894e-17,   5.89176e-24,
     7.18271e-03,   5.48463e-10,   6.37986e-17,   1.02172e-23,
     1.13364e-02,   8.65631e-10,   1.00692e-16,   1.61257e-23,
     1.77361e-02,   1.35430e-09,   1.57536e-16,   2.52291e-23,
     2.59477e-02,   1.98134e-09,   2.30474e-16,   3.69100e-23,
     3.45425e-02,   2.63763e-09,   3.06815e-16,   4.91359e-23,
     4.22006e-02,   3.22239e-09,   3.74836e-16,   6.00294e-23,
     4.71420e-02,   3.59971e-09,   4.18726e-16,   6.70583e-23,
     4.91934e-02,   3.75635e-09,   4.36947e-16,   6.99764e-23,
     5.05162e-02,   3.85735e-09,   4.48697e-16,   7.18580e-23,
     5.78201e-02,   4.41508e-09,   5.13572e-16,   8.22477e-23,
     8.84237e-02,   6.75193e-09,   7.85401e-16,   1.25781e-22,
     1.78786e-01,   1.36519e-08,   1.58802e-15,   2.54319e-22,
     4.36405e-01,   3.33234e-08,   3.87627e-15,   6.20781e-22,
     1.63797e+00,   1.25074e-07,   1.45490e-14,   2.33003e-21,
     8.50820e+00,   6.49679e-07,   7.55727e-14,   1.21029e-20,
     3.92752e+01,   2.99901e-06,   3.48854e-13,   5.58686e-20,
     1.41437e+02,   1.08001e-05,   1.25630e-12,   2.01195e-19,
     3.83709e+02,   2.92996e-05,   3.40821e-12,   5.45820e-19,
     7.70412e+02,   5.88278e-05,   6.84300e-12,   1.09590e-18,
     1.16399e+03,   8.88808e-05,   1.03388e-11,   1.65574e-18,
     1.37566e+03,   1.05044e-04,   1.22190e-11,   1.95685e-18,
     1.33070e+03,   1.01611e-04,   1.18196e-11,   1.89289e-18,
     1.09978e+03,   8.39785e-05,   9.76864e-12,   1.56444e-18,
     8.05639e+02,   6.15177e-05,   7.15589e-12,   1.14600e-18,
     5.38690e+02,   4.11337e-05,   4.78477e-12,   7.66273e-19,
     3.36338e+02,   2.56824e-05,   2.98743e-12,   4.78432e-19,
     1.99460e+02,   1.52305e-05,   1.77165e-12,   2.83727e-19,
     1.13787e+02,   8.68865e-06,   1.01068e-12,   1.61859e-19,
     6.30411e+01,   4.81374e-06,   5.59946e-13,   8.96744e-20,
     3.41529e+01,   2.60788e-06,   3.03354e-13,   4.85817e-20,
     1.81893e+01,   1.38891e-06,   1.61561e-13,   2.58738e-20  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = C30_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = C30_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = C30_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = C30_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = C30_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = C30_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = C30_kpMgO     [itab];
      my_rates->SN0_kpFeS     [itab0] = C30_kpFeS     [itab];
      my_rates->SN0_kpAl2O3   [itab0] = C30_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_F13(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   6.69235e-01;
  my_rates->SN0_XO [iSN] =   3.30556e-01;
  my_rates->SN0_XMg[iSN] =   1.86824e-04;
  my_rates->SN0_XAl[iSN] =   1.97017e-07;
  my_rates->SN0_XSi[iSN] =   1.30184e-05;
  my_rates->SN0_XS [iSN] =   0.00000e+00;
  my_rates->SN0_XFe[iSN] =   8.90341e-06;

  my_rates->SN0_fC [iSN] =   4.93693e-01;
  my_rates->SN0_fO [iSN] =   3.30556e-01;
  my_rates->SN0_fMg[iSN] =   1.86824e-04;
  my_rates->SN0_fAl[iSN] =   1.97017e-07;
  my_rates->SN0_fSi[iSN] =   1.30184e-05;
  my_rates->SN0_fS [iSN] =   0.00000e+00;
  my_rates->SN0_fFe[iSN] =   8.90341e-06;

  my_rates->SN0_fFeM     [iSN] =   6.31648e-26;
  my_rates->SN0_fMg2SiO4 [iSN] =   2.06081e-16;
  my_rates->SN0_fMgSiO3  [iSN] =   3.19262e-15;
  my_rates->SN0_fFe3O4   [iSN] =   4.37192e-15;
  my_rates->SN0_fAC      [iSN] =   1.75542e-01;
  my_rates->SN0_fSiO2D   [iSN] =   1.92019e-16;
  my_rates->SN0_fAl2O3   [iSN] =   6.23283e-17;

  itab0 = 3 * iSN;
  my_rates->SN0_r0FeM     [itab0 + 0] =   4.02937e-08;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   4.03307e-08;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   4.03157e-08;
  my_rates->SN0_r0Fe3O4   [itab0 + 0] =   4.03312e-08;
  my_rates->SN0_r0AC      [itab0 + 0] =   6.60867e-06;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   4.03146e-08;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   4.03146e-08;

  my_rates->SN0_r0FeM     [itab0 + 1] =   1.67044e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   1.67330e-15;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.67182e-15;
  my_rates->SN0_r0Fe3O4   [itab0 + 1] =   1.67336e-15;
  my_rates->SN0_r0AC      [itab0 + 1] =   5.49310e-11;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   1.67171e-15;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   1.67171e-15;

  my_rates->SN0_r0FeM     [itab0 + 2] =   7.11477e-23;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   7.13316e-23;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   7.12190e-23;
  my_rates->SN0_r0Fe3O4   [itab0 + 2] =   7.13357e-23;
  my_rates->SN0_r0AC      [itab0 + 2] =   5.25955e-16;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   7.12105e-23;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   7.12106e-23;

  NTd =            35;
 Nmom =             4;

  double F13_kpFeM[] = 
  {  1.23621e-05,   4.98941e-13,   2.07173e-20,   8.83710e-28,
     2.19539e-05,   8.86065e-13,   3.67916e-20,   1.56937e-27,
     3.40291e-05,   1.37342e-12,   5.70280e-20,   2.43256e-27,
     4.92310e-05,   1.98698e-12,   8.25041e-20,   3.51925e-27,
     8.08514e-05,   3.26317e-12,   1.35494e-19,   5.77953e-27,
     1.25020e-04,   5.04574e-12,   2.09508e-19,   8.93651e-27,
     1.96586e-04,   7.93387e-12,   3.29418e-19,   1.40509e-26,
     3.14491e-04,   1.26917e-11,   5.26937e-19,   2.24747e-26,
     5.06850e-04,   2.04532e-11,   8.49130e-19,   3.62146e-26,
     8.07286e-04,   3.25749e-11,   1.35229e-18,   5.76710e-26,
     1.28668e-03,   5.19155e-11,   2.15506e-18,   9.19009e-26,
     2.05241e-03,   8.28056e-11,   3.43709e-18,   1.46562e-25,
     3.27026e-03,   1.31928e-10,   5.47555e-18,   2.33466e-25,
     5.23898e-03,   2.11325e-10,   8.76988e-18,   3.73889e-25,
     8.45023e-03,   3.40811e-10,   1.41417e-17,   6.02834e-25,
     1.37158e-02,   5.53101e-10,   2.29473e-17,   9.78076e-25,
     2.24100e-02,   9.03572e-10,   3.74827e-17,   1.59741e-24,
     3.70042e-02,   1.49181e-09,   6.18765e-17,   2.63669e-24,
     6.21585e-02,   2.50559e-09,   1.03913e-16,   4.42746e-24,
     1.07033e-01,   4.31400e-09,   1.78894e-16,   7.62148e-24,
     1.90089e-01,   7.66095e-09,   3.17659e-16,   1.35323e-23,
     3.49470e-01,   1.40834e-08,   5.83928e-16,   2.48739e-23,
     6.64947e-01,   2.67957e-08,   1.11096e-15,   4.73220e-23,
     1.30413e+00,   5.25515e-08,   2.17873e-15,   9.28021e-23,
     2.61640e+00,   1.05429e-07,   4.37088e-15,   1.86172e-22,
     5.31791e+00,   2.14284e-07,   8.88372e-15,   3.78386e-22,
     1.08366e+01,   4.36654e-07,   1.81025e-14,   7.71039e-22,
     2.19132e+01,   8.82975e-07,   3.66056e-14,   1.55913e-21,
     4.35354e+01,   1.75422e-06,   7.27247e-14,   3.09753e-21,
     8.42362e+01,   3.39422e-06,   1.40714e-13,   5.99335e-21,
     1.57704e+02,   6.35452e-06,   2.63439e-13,   1.12205e-20,
     2.84822e+02,   1.14766e-05,   4.75784e-13,   2.02648e-20,
     4.96653e+02,   2.00121e-05,   8.29639e-13,   3.53363e-20,
     8.39966e+02,   3.38455e-05,   1.40313e-12,   5.97626e-20,
     1.38932e+03,   5.59813e-05,   2.32081e-12,   9.88487e-20  };

  double F13_kpMg2SiO4[] = 
  {  1.05240e-01,   4.24440e-09,   1.76098e-16,   7.50693e-24,
     1.32588e-01,   5.34735e-09,   2.21859e-16,   9.45769e-24,
     1.67016e-01,   6.73589e-09,   2.79469e-16,   1.19135e-23,
     2.10360e-01,   8.48395e-09,   3.51995e-16,   1.50053e-23,
     2.71887e-01,   1.09654e-08,   4.54949e-16,   1.93941e-23,
     3.55694e-01,   1.43454e-08,   5.95184e-16,   2.53722e-23,
     4.84932e-01,   1.95577e-08,   8.11439e-16,   3.45910e-23,
     6.99767e-01,   2.82221e-08,   1.17092e-15,   4.99155e-23,
     1.05860e+00,   4.26939e-08,   1.77135e-15,   7.55113e-23,
     1.62902e+00,   6.56994e-08,   2.72584e-15,   1.16200e-22,
     2.54260e+00,   1.02545e-07,   4.25454e-15,   1.81368e-22,
     3.96488e+00,   1.59906e-07,   6.63444e-15,   2.82821e-22,
     6.10630e+00,   2.46271e-07,   1.02177e-14,   4.35572e-22,
     9.28766e+00,   3.74578e-07,   1.55411e-14,   6.62503e-22,
     1.39265e+01,   5.61664e-07,   2.33032e-14,   9.93397e-22,
     2.05382e+01,   8.28321e-07,   3.43667e-14,   1.46502e-21,
     3.00649e+01,   1.21254e-06,   5.03077e-14,   2.14458e-21,
     4.55102e+01,   1.83546e-06,   7.61523e-14,   3.24631e-21,
     7.47839e+01,   3.01609e-06,   1.25136e-13,   5.33445e-21,
     1.29621e+02,   5.22769e-06,   2.16895e-13,   9.24605e-21,
     2.14820e+02,   8.66384e-06,   3.59459e-13,   1.53234e-20,
     3.20002e+02,   1.29059e-05,   5.35460e-13,   2.28262e-20,
     4.29768e+02,   1.73329e-05,   7.19133e-13,   3.06560e-20,
     5.30827e+02,   2.14086e-05,   8.88234e-13,   3.78647e-20,
     5.99694e+02,   2.41861e-05,   1.00347e-12,   4.27771e-20,
     6.06537e+02,   2.44620e-05,   1.01492e-12,   4.32652e-20,
     5.43262e+02,   2.19101e-05,   9.09042e-13,   3.87517e-20,
     4.33545e+02,   1.74852e-05,   7.25453e-13,   3.09255e-20,
     3.13324e+02,   1.26366e-05,   5.24285e-13,   2.23499e-20,
     2.09006e+02,   8.42935e-06,   3.49730e-13,   1.49087e-20,
     1.31150e+02,   5.28937e-06,   2.19454e-13,   9.35513e-21,
     7.90681e+01,   3.18887e-06,   1.32305e-13,   5.64006e-21,
     4.73389e+01,   1.90921e-06,   7.92124e-14,   3.37676e-21,
     2.96409e+01,   1.19544e-06,   4.95982e-14,   2.11433e-21,
     2.04708e+01,   8.25601e-07,   3.42539e-14,   1.46021e-21  };

  double F13_kpMgSiO3[] = 
  {  2.19890e-02,   8.86503e-10,   3.67618e-17,   1.56604e-24,
     3.90612e-02,   1.57478e-09,   6.53036e-17,   2.78190e-24,
     6.05539e-02,   2.44128e-09,   1.01236e-16,   4.31259e-24,
     8.76116e-02,   3.53213e-09,   1.46471e-16,   6.23961e-24,
     1.43288e-01,   5.77674e-09,   2.39552e-16,   1.02048e-23,
     2.19266e-01,   8.83988e-09,   3.66575e-16,   1.56159e-23,
     3.36256e-01,   1.35564e-08,   5.62160e-16,   2.39478e-23,
     5.14336e-01,   2.07358e-08,   8.59879e-16,   3.66305e-23,
     7.97216e-01,   3.21404e-08,   1.33281e-15,   5.67770e-23,
     1.25414e+00,   5.05616e-08,   2.09670e-15,   8.93186e-23,
     2.03450e+00,   8.20224e-08,   3.40133e-15,   1.44895e-22,
     3.34648e+00,   1.34916e-07,   5.59472e-15,   2.38333e-22,
     5.45893e+00,   2.20081e-07,   9.12638e-15,   3.88780e-22,
     8.82117e+00,   3.55632e-07,   1.47474e-14,   6.28235e-22,
     1.41825e+01,   5.71778e-07,   2.37106e-14,   1.01006e-21,
     2.28419e+01,   9.20889e-07,   3.81877e-14,   1.62678e-21,
     3.71178e+01,   1.49643e-06,   6.20544e-14,   2.64349e-21,
     6.14272e+01,   2.47648e-06,   1.02696e-13,   4.37479e-21,
     1.03847e+02,   4.18665e-06,   1.73613e-13,   7.39586e-21,
     1.75507e+02,   7.07567e-06,   2.93416e-13,   1.24994e-20,
     2.82060e+02,   1.13715e-05,   4.71555e-13,   2.00880e-20,
     4.14519e+02,   1.67116e-05,   6.93003e-13,   2.95216e-20,
     5.59961e+02,   2.25752e-05,   9.36157e-13,   3.98799e-20,
     7.11024e+02,   2.86655e-05,   1.18871e-12,   5.06385e-20,
     8.40805e+02,   3.38977e-05,   1.40568e-12,   5.98813e-20,
     8.95312e+02,   3.60952e-05,   1.49681e-12,   6.37633e-20,
     8.40415e+02,   3.38819e-05,   1.40503e-12,   5.98535e-20,
     6.96693e+02,   2.80877e-05,   1.16475e-12,   4.96178e-20,
     5.18202e+02,   2.08917e-05,   8.66344e-13,   3.69059e-20,
     3.52864e+02,   1.42260e-05,   5.89927e-13,   2.51306e-20,
     2.24210e+02,   9.03919e-06,   3.74840e-13,   1.59680e-20,
     1.35138e+02,   5.44818e-06,   2.25927e-13,   9.62439e-21,
     7.83119e+01,   3.15720e-06,   1.30924e-13,   5.57730e-21,
     4.41553e+01,   1.78015e-06,   7.38199e-14,   3.14470e-21,
     2.45888e+01,   9.91317e-07,   4.11082e-14,   1.75119e-21  };

  double F13_kpFe3O4[] = 
  {  1.47700e-02,   5.95693e-10,   2.47155e-17,   1.05363e-24,
     2.47694e-02,   9.98982e-10,   4.14481e-17,   1.76695e-24,
     3.73580e-02,   1.50669e-09,   6.25133e-17,   2.66496e-24,
     5.32060e-02,   2.14587e-09,   8.90327e-17,   3.79549e-24,
     8.50036e-02,   3.42830e-09,   1.42241e-16,   6.06380e-24,
     1.29213e-01,   5.21132e-09,   2.16220e-16,   9.21750e-24,
     2.00170e-01,   8.07309e-09,   3.34956e-16,   1.42793e-23,
     3.15560e-01,   1.27269e-08,   5.28045e-16,   2.25107e-23,
     5.01384e-01,   2.02214e-08,   8.38995e-16,   3.57666e-23,
     7.88907e-01,   3.18176e-08,   1.32012e-15,   5.62773e-23,
     1.24250e+00,   5.01116e-08,   2.07915e-15,   8.86347e-23,
     1.95225e+00,   7.87365e-08,   3.26681e-15,   1.39265e-22,
     3.04002e+00,   1.22608e-07,   5.08704e-15,   2.16862e-22,
     4.68918e+00,   1.89120e-07,   7.84668e-15,   3.34506e-22,
     7.12599e+00,   2.87400e-07,   1.19243e-14,   5.08338e-22,
     1.05834e+01,   4.26842e-07,   1.77098e-14,   7.54974e-22,
     1.52356e+01,   6.14471e-07,   2.54946e-14,   1.08684e-21,
     2.13345e+01,   8.60449e-07,   3.57003e-14,   1.52192e-21,
     2.98061e+01,   1.20212e-06,   4.98762e-14,   2.12624e-21,
     4.27642e+01,   1.72473e-06,   7.15598e-14,   3.05062e-21,
     6.30370e+01,   2.54236e-06,   1.05483e-13,   4.49679e-21,
     9.29361e+01,   3.74823e-06,   1.55515e-13,   6.62967e-21,
     1.32987e+02,   5.36353e-06,   2.22535e-13,   9.48673e-21,
     1.82150e+02,   7.34635e-06,   3.04803e-13,   1.29938e-20,
     2.40388e+02,   9.69513e-06,   4.02254e-13,   1.71482e-20,
     3.12065e+02,   1.25860e-05,   5.22197e-13,   2.22614e-20,
     4.08414e+02,   1.64718e-05,   6.83423e-13,   2.91345e-20,
     5.49591e+02,   2.21657e-05,   9.19662e-13,   3.92055e-20,
     7.67451e+02,   3.09523e-05,   1.28422e-12,   5.47467e-20,
     1.10725e+03,   4.46570e-05,   1.85283e-12,   7.89869e-20,
     1.62060e+03,   6.53608e-05,   2.71184e-12,   1.15607e-19,
     2.33999e+03,   9.43747e-05,   3.91564e-12,   1.66925e-19,
     3.24367e+03,   1.30821e-04,   5.42783e-12,   2.31390e-19,
     4.25716e+03,   1.71697e-04,   7.12376e-12,   3.03688e-19,
     5.34010e+03,   2.15373e-04,   8.93591e-12,   3.80941e-19  };

  double F13_kpAC[] = 
  {  3.27960e-01,   2.16737e-06,   1.80151e-11,   1.72491e-16,
     4.38754e-01,   2.89959e-06,   2.41015e-11,   2.30770e-16,
     5.78235e-01,   3.82140e-06,   3.17638e-11,   3.04139e-16,
     7.53832e-01,   4.98189e-06,   4.14101e-11,   3.96504e-16,
     1.04018e+00,   6.87442e-06,   5.71418e-11,   5.47146e-16,
     1.41746e+00,   9.36786e-06,   7.78688e-11,   7.45623e-16,
     1.95306e+00,   1.29077e-05,   1.07295e-10,   1.02741e-15,
     2.71551e+00,   1.79470e-05,   1.49187e-10,   1.42858e-15,
     3.79717e+00,   2.50964e-05,   2.08623e-10,   1.99781e-15,
     5.29825e+00,   3.50184e-05,   2.91115e-10,   2.78791e-15,
     7.37979e+00,   4.87780e-05,   4.05522e-10,   3.88380e-15,
     1.02197e+01,   6.75525e-05,   5.61646e-10,   5.37953e-15,
     1.40472e+01,   9.28592e-05,   7.72124e-10,   7.39643e-15,
     1.92121e+01,   1.27015e-04,   1.05627e-09,   1.01201e-14,
     2.61803e+01,   1.73108e-04,   1.43985e-09,   1.37985e-14,
     3.55654e+01,   2.35211e-04,   1.95692e-09,   1.87600e-14,
     4.82268e+01,   3.19039e-04,   2.65535e-09,   2.54677e-14,
     6.54410e+01,   4.33093e-04,   3.60654e-09,   3.46142e-14,
     8.90033e+01,   5.89369e-04,   4.91159e-09,   4.71850e-14,
     1.21154e+02,   8.02913e-04,   6.69826e-09,   6.44364e-14,
     1.64484e+02,   1.09135e-03,   9.11847e-09,   8.78899e-14,
     2.22227e+02,   1.47698e-03,   1.23680e-08,   1.19551e-13,
     2.99250e+02,   1.99360e-03,   1.67456e-08,   1.62501e-13,
     4.02841e+02,   2.69134e-03,   2.26888e-08,   2.21192e-13,
     5.41754e+02,   3.62947e-03,   3.07031e-08,   3.00611e-13,
     7.24111e+02,   4.86183e-03,   4.12300e-08,   4.04886e-13,
     9.58822e+02,   6.44664e-03,   5.47294e-08,   5.38055e-13,
     1.26198e+03,   8.49034e-03,   7.20550e-08,   7.07781e-13,
     1.66376e+03,   1.11958e-02,   9.48672e-08,   9.29402e-13,
     2.21388e+03,   1.49013e-02,   1.25959e-07,   1.22898e-12,
     2.99023e+03,   2.01435e-02,   1.69784e-07,   1.64820e-12,
     4.11644e+03,   2.77858e-02,   2.33536e-07,   2.25436e-12,
     5.78855e+03,   3.92074e-02,   3.28708e-07,   3.15462e-12,
     8.29050e+03,   5.63946e-02,   4.71738e-07,   4.50079e-12,
     1.19408e+04,   8.14776e-02,   6.79738e-07,   6.44637e-12  };

  double F13_kpSiO2D[] = 
  {  7.60360e-02,   3.06536e-09,   1.27110e-16,   5.41456e-24,
     9.07207e-02,   3.65737e-09,   1.51659e-16,   6.46027e-24,
     1.09208e-01,   4.40266e-09,   1.82564e-16,   7.77673e-24,
     1.32481e-01,   5.34093e-09,   2.21471e-16,   9.43407e-24,
     1.58907e-01,   6.40629e-09,   2.65648e-16,   1.13159e-23,
     1.91565e-01,   7.72285e-09,   3.20241e-16,   1.36414e-23,
     2.30490e-01,   9.29212e-09,   3.85313e-16,   1.64133e-23,
     2.76795e-01,   1.11589e-08,   4.62722e-16,   1.97107e-23,
     3.33074e-01,   1.34277e-08,   5.56804e-16,   2.37184e-23,
     4.05325e-01,   1.63405e-08,   6.77586e-16,   2.88634e-23,
     5.08160e-01,   2.04863e-08,   8.49498e-16,   3.61863e-23,
     6.72472e-01,   2.71105e-08,   1.12418e-15,   4.78871e-23,
     9.48549e-01,   3.82404e-08,   1.58570e-15,   6.75467e-23,
     1.41787e+00,   5.71610e-08,   2.37028e-15,   1.00968e-22,
     2.19502e+00,   8.84912e-08,   3.66944e-15,   1.56308e-22,
     3.46719e+00,   1.39778e-07,   5.79615e-15,   2.46901e-22,
     5.76852e+00,   2.32556e-07,   9.64331e-15,   4.10780e-22,
     1.17194e+01,   4.72463e-07,   1.95915e-14,   8.34544e-22,
     3.16449e+01,   1.27575e-06,   5.29012e-14,   2.25345e-21,
     8.68296e+01,   3.50050e-06,   1.45154e-13,   6.18318e-21,
     1.92300e+02,   7.75250e-06,   3.21470e-13,   1.36938e-20,
     3.36231e+02,   1.35550e-05,   5.62081e-13,   2.39432e-20,
     5.05825e+02,   2.03921e-05,   8.45594e-13,   3.60201e-20,
     7.20624e+02,   2.90517e-05,   1.20468e-12,   5.13160e-20,
     9.77376e+02,   3.94025e-05,   1.63389e-12,   6.95995e-20,
     1.18646e+03,   4.78315e-05,   1.98341e-12,   8.44881e-20,
     1.23845e+03,   4.99275e-05,   2.07033e-12,   8.81904e-20,
     1.11188e+03,   4.48251e-05,   1.85875e-12,   7.91777e-20,
     8.76396e+02,   3.53316e-05,   1.46508e-12,   6.24086e-20,
     6.22207e+02,   2.50840e-05,   1.04015e-12,   4.43077e-20,
     4.07290e+02,   1.64197e-05,   6.80872e-13,   2.90033e-20,
     2.50573e+02,   1.01017e-05,   4.18886e-13,   1.78434e-20,
     1.47073e+02,   5.92921e-06,   2.45865e-13,   1.04732e-20,
     8.33122e+01,   3.35870e-06,   1.39274e-13,   5.93271e-21,
     4.59585e+01,   1.85280e-06,   7.68295e-14,   3.27273e-21  };

  double F13_kpAl2O3[] = 
  {  9.93250e-04,   4.00425e-11,   1.66043e-18,   7.07299e-26,
     1.81240e-03,   7.30662e-11,   3.02981e-18,   1.29062e-25,
     2.84365e-03,   1.14641e-10,   4.75376e-18,   2.02498e-25,
     4.14191e-03,   1.66980e-10,   6.92409e-18,   2.94948e-25,
     7.18271e-03,   2.89568e-10,   1.20074e-17,   5.11485e-25,
     1.13364e-02,   4.57021e-10,   1.89512e-17,   8.07269e-25,
     1.77361e-02,   7.15022e-10,   2.96496e-17,   1.26300e-24,
     2.59477e-02,   1.04607e-09,   4.33772e-17,   1.84775e-24,
     3.45425e-02,   1.39257e-09,   5.77452e-17,   2.45979e-24,
     4.22006e-02,   1.70130e-09,   7.05474e-17,   3.00513e-24,
     4.71420e-02,   1.90051e-09,   7.88079e-17,   3.35701e-24,
     4.91934e-02,   1.98321e-09,   8.22373e-17,   3.50309e-24,
     5.05162e-02,   2.03654e-09,   8.44486e-17,   3.59729e-24,
     5.78201e-02,   2.33100e-09,   9.66587e-17,   4.11740e-24,
     8.84237e-02,   3.56477e-09,   1.47819e-16,   6.29670e-24,
     1.78786e-01,   7.20769e-09,   2.98879e-16,   1.27315e-23,
     4.36404e-01,   1.75935e-08,   7.29542e-16,   3.10766e-23,
     1.63796e+00,   6.60337e-08,   2.73820e-15,   1.16640e-22,
     8.50817e+00,   3.43004e-07,   1.42232e-14,   6.05872e-22,
     3.92751e+01,   1.58336e-06,   6.56567e-14,   2.79680e-21,
     1.41436e+02,   5.70194e-06,   2.36441e-13,   1.00718e-20,
     3.83709e+02,   1.54691e-05,   6.41451e-13,   2.73241e-20,
     7.70411e+02,   3.10588e-05,   1.28791e-12,   5.48614e-20,
     1.16399e+03,   4.69258e-05,   1.94586e-12,   8.28883e-20,
     1.37566e+03,   5.54594e-05,   2.29972e-12,   9.79619e-20,
     1.33070e+03,   5.36466e-05,   2.22455e-12,   9.47599e-20,
     1.09978e+03,   4.43371e-05,   1.83851e-12,   7.83159e-20,
     8.05638e+02,   3.24790e-05,   1.34680e-12,   5.73700e-20,
     5.38690e+02,   2.17171e-05,   9.00535e-13,   3.83604e-20,
     3.36338e+02,   1.35593e-05,   5.62261e-13,   2.39508e-20,
     1.99460e+02,   8.04115e-06,   3.33440e-13,   1.42037e-20,
     1.13787e+02,   4.58728e-06,   1.90220e-13,   8.10285e-21,
     6.30411e+01,   2.54148e-06,   1.05387e-13,   4.48919e-21,
     3.41529e+01,   1.37686e-06,   5.70939e-14,   2.43205e-21,
     1.81893e+01,   7.33293e-07,   3.04073e-14,   1.29527e-21  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpFeM     [itab0] = F13_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = F13_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = F13_kpMgSiO3  [itab];
      my_rates->SN0_kpFe3O4   [itab0] = F13_kpFe3O4   [itab];
      my_rates->SN0_kpAC      [itab0] = F13_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = F13_kpSiO2D   [itab];
      my_rates->SN0_kpAl2O3   [itab0] = F13_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_F15(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   6.46299e-01;
  my_rates->SN0_XO [iSN] =   3.53548e-01;
  my_rates->SN0_XMg[iSN] =   1.29204e-04;
  my_rates->SN0_XAl[iSN] =   2.22729e-07;
  my_rates->SN0_XSi[iSN] =   1.32242e-05;
  my_rates->SN0_XS [iSN] =   0.00000e+00;
  my_rates->SN0_XFe[iSN] =   9.66658e-06;

  my_rates->SN0_fC [iSN] =   4.57071e-01;
  my_rates->SN0_fO [iSN] =   3.53548e-01;
  my_rates->SN0_fMg[iSN] =   1.29204e-04;
  my_rates->SN0_fAl[iSN] =   2.22729e-07;
  my_rates->SN0_fSi[iSN] =   1.32242e-05;
  my_rates->SN0_fS [iSN] =   0.00000e+00;
  my_rates->SN0_fFe[iSN] =   9.66658e-06;

  my_rates->SN0_fFeM     [iSN] =   1.53361e-25;
  my_rates->SN0_fMg2SiO4 [iSN] =   1.56864e-15;
  my_rates->SN0_fMgSiO3  [iSN] =   2.13810e-14;
  my_rates->SN0_fFe3O4   [iSN] =   1.22287e-14;
  my_rates->SN0_fAC      [iSN] =   1.89229e-01;
  my_rates->SN0_fSiO2D   [iSN] =   1.47463e-15;
  my_rates->SN0_fAl2O3   [iSN] =   2.15191e-16;

  itab0 = 3 * iSN;
  my_rates->SN0_r0FeM     [itab0 + 0] =   4.02634e-08;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   4.03318e-08;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   4.03159e-08;
  my_rates->SN0_r0Fe3O4   [itab0 + 0] =   4.03301e-08;
  my_rates->SN0_r0AC      [itab0 + 0] =   1.14540e-05;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   4.03146e-08;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   4.03146e-08;

  my_rates->SN0_r0FeM     [itab0 + 1] =   1.66860e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   1.67341e-15;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.67184e-15;
  my_rates->SN0_r0Fe3O4   [itab0 + 1] =   1.67324e-15;
  my_rates->SN0_r0AC      [itab0 + 1] =   1.60512e-10;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   1.67171e-15;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   1.67171e-15;

  my_rates->SN0_r0FeM     [itab0 + 2] =   7.10566e-23;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   7.13397e-23;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   7.12201e-23;
  my_rates->SN0_r0Fe3O4   [itab0 + 2] =   7.13269e-23;
  my_rates->SN0_r0AC      [itab0 + 2] =   2.55303e-15;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   7.12105e-23;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   7.12106e-23;

  NTd =            35;
 Nmom =             4;

  double F15_kpFeM[] = 
  {  1.23614e-05,   4.98551e-13,   2.06942e-20,   8.82572e-28,
     2.19525e-05,   8.85374e-13,   3.67505e-20,   1.56735e-27,
     3.40270e-05,   1.37235e-12,   5.69643e-20,   2.42942e-27,
     4.92280e-05,   1.98543e-12,   8.24119e-20,   3.51472e-27,
     8.08465e-05,   3.26062e-12,   1.35343e-19,   5.77209e-27,
     1.25012e-04,   5.04181e-12,   2.09274e-19,   8.92500e-27,
     1.96574e-04,   7.92769e-12,   3.29050e-19,   1.40328e-26,
     3.14473e-04,   1.26818e-11,   5.26348e-19,   2.24457e-26,
     5.06822e-04,   2.04373e-11,   8.48182e-19,   3.61680e-26,
     8.07243e-04,   3.25496e-11,   1.35079e-18,   5.75968e-26,
     1.28661e-03,   5.18753e-11,   2.15265e-18,   9.17827e-26,
     2.05232e-03,   8.27416e-11,   3.43325e-18,   1.46374e-25,
     3.27011e-03,   1.31826e-10,   5.46945e-18,   2.33165e-25,
     5.23877e-03,   2.11162e-10,   8.76012e-18,   3.73409e-25,
     8.44994e-03,   3.40550e-10,   1.41259e-17,   6.02059e-25,
     1.37154e-02,   5.52677e-10,   2.29218e-17,   9.76820e-25,
     2.24094e-02,   9.02883e-10,   3.74411e-17,   1.59536e-24,
     3.70035e-02,   1.49068e-09,   6.18079e-17,   2.63331e-24,
     6.21576e-02,   2.50369e-09,   1.03798e-16,   4.42178e-24,
     1.07032e-01,   4.31074e-09,   1.78696e-16,   7.61171e-24,
     1.90087e-01,   7.65516e-09,   3.17308e-16,   1.35149e-23,
     3.49468e-01,   1.40728e-08,   5.83283e-16,   2.48420e-23,
     6.64945e-01,   2.67755e-08,   1.10973e-15,   4.72613e-23,
     1.30413e+00,   5.25119e-08,   2.17633e-15,   9.26832e-23,
     2.61639e+00,   1.05349e-07,   4.36605e-15,   1.85933e-22,
     5.31790e+00,   2.14123e-07,   8.87391e-15,   3.77901e-22,
     1.08366e+01,   4.36326e-07,   1.80825e-14,   7.70051e-22,
     2.19131e+01,   8.82310e-07,   3.65652e-14,   1.55714e-21,
     4.35353e+01,   1.75290e-06,   7.26444e-14,   3.09356e-21,
     8.42362e+01,   3.39166e-06,   1.40558e-13,   5.98567e-21,
     1.57704e+02,   6.34974e-06,   2.63148e-13,   1.12061e-20,
     2.84822e+02,   1.14680e-05,   4.75258e-13,   2.02388e-20,
     4.96653e+02,   1.99971e-05,   8.28723e-13,   3.52910e-20,
     8.39966e+02,   3.38201e-05,   1.40158e-12,   5.96860e-20,
     1.38932e+03,   5.59392e-05,   2.31824e-12,   9.87221e-20  };

  double F15_kpMg2SiO4[] = 
  {  1.05240e-01,   4.24452e-09,   1.76110e-16,   7.50779e-24,
     1.32588e-01,   5.34750e-09,   2.21874e-16,   9.45877e-24,
     1.67016e-01,   6.73607e-09,   2.79487e-16,   1.19149e-23,
     2.10360e-01,   8.48418e-09,   3.52018e-16,   1.50070e-23,
     2.71887e-01,   1.09657e-08,   4.54979e-16,   1.93963e-23,
     3.55694e-01,   1.43458e-08,   5.95222e-16,   2.53751e-23,
     4.84932e-01,   1.95582e-08,   8.11491e-16,   3.45949e-23,
     6.99767e-01,   2.82229e-08,   1.17100e-15,   4.99212e-23,
     1.05860e+00,   4.26950e-08,   1.77146e-15,   7.55199e-23,
     1.62902e+00,   6.57012e-08,   2.72602e-15,   1.16214e-22,
     2.54260e+00,   1.02548e-07,   4.25481e-15,   1.81388e-22,
     3.96488e+00,   1.59911e-07,   6.63487e-15,   2.82853e-22,
     6.10630e+00,   2.46278e-07,   1.02183e-14,   4.35621e-22,
     9.28766e+00,   3.74588e-07,   1.55421e-14,   6.62579e-22,
     1.39265e+01,   5.61679e-07,   2.33047e-14,   9.93510e-22,
     2.05382e+01,   8.28344e-07,   3.43689e-14,   1.46519e-21,
     3.00649e+01,   1.21257e-06,   5.03110e-14,   2.14482e-21,
     4.55102e+01,   1.83551e-06,   7.61572e-14,   3.24668e-21,
     7.47839e+01,   3.01617e-06,   1.25144e-13,   5.33506e-21,
     1.29621e+02,   5.22783e-06,   2.16909e-13,   9.24710e-21,
     2.14820e+02,   8.66407e-06,   3.59482e-13,   1.53252e-20,
     3.20002e+02,   1.29062e-05,   5.35494e-13,   2.28288e-20,
     4.29768e+02,   1.73333e-05,   7.19179e-13,   3.06596e-20,
     5.30827e+02,   2.14092e-05,   8.88291e-13,   3.78690e-20,
     5.99694e+02,   2.41867e-05,   1.00353e-12,   4.27820e-20,
     6.06537e+02,   2.44627e-05,   1.01498e-12,   4.32701e-20,
     5.43262e+02,   2.19107e-05,   9.09100e-13,   3.87561e-20,
     4.33545e+02,   1.74857e-05,   7.25499e-13,   3.09290e-20,
     3.13324e+02,   1.26369e-05,   5.24319e-13,   2.23524e-20,
     2.09006e+02,   8.42958e-06,   3.49753e-13,   1.49104e-20,
     1.31150e+02,   5.28951e-06,   2.19468e-13,   9.35620e-21,
     7.90681e+01,   3.18896e-06,   1.32313e-13,   5.64070e-21,
     4.73389e+01,   1.90926e-06,   7.92175e-14,   3.37715e-21,
     2.96409e+01,   1.19547e-06,   4.96014e-14,   2.11457e-21,
     2.04708e+01,   8.25623e-07,   3.42560e-14,   1.46038e-21  };

  double F15_kpMgSiO3[] = 
  {  2.19890e-02,   8.86506e-10,   3.67621e-17,   1.56606e-24,
     3.90612e-02,   1.57479e-09,   6.53041e-17,   2.78195e-24,
     6.05539e-02,   2.44129e-09,   1.01236e-16,   4.31266e-24,
     8.76116e-02,   3.53214e-09,   1.46473e-16,   6.23971e-24,
     1.43288e-01,   5.77677e-09,   2.39554e-16,   1.02050e-23,
     2.19266e-01,   8.83991e-09,   3.66578e-16,   1.56162e-23,
     3.36256e-01,   1.35564e-08,   5.62165e-16,   2.39482e-23,
     5.14336e-01,   2.07359e-08,   8.59887e-16,   3.66311e-23,
     7.97216e-01,   3.21405e-08,   1.33282e-15,   5.67779e-23,
     1.25414e+00,   5.05618e-08,   2.09672e-15,   8.93200e-23,
     2.03450e+00,   8.20227e-08,   3.40136e-15,   1.44897e-22,
     3.34648e+00,   1.34916e-07,   5.59477e-15,   2.38337e-22,
     5.45893e+00,   2.20082e-07,   9.12646e-15,   3.88786e-22,
     8.82117e+00,   3.55633e-07,   1.47476e-14,   6.28245e-22,
     1.41825e+01,   5.71780e-07,   2.37109e-14,   1.01008e-21,
     2.28419e+01,   9.20892e-07,   3.81880e-14,   1.62680e-21,
     3.71178e+01,   1.49644e-06,   6.20550e-14,   2.64353e-21,
     6.14272e+01,   2.47649e-06,   1.02696e-13,   4.37486e-21,
     1.03847e+02,   4.18667e-06,   1.73615e-13,   7.39597e-21,
     1.75507e+02,   7.07570e-06,   2.93419e-13,   1.24996e-20,
     2.82060e+02,   1.13715e-05,   4.71559e-13,   2.00884e-20,
     4.14519e+02,   1.67117e-05,   6.93009e-13,   2.95221e-20,
     5.59961e+02,   2.25753e-05,   9.36165e-13,   3.98805e-20,
     7.11024e+02,   2.86656e-05,   1.18872e-12,   5.06393e-20,
     8.40805e+02,   3.38978e-05,   1.40569e-12,   5.98823e-20,
     8.95312e+02,   3.60953e-05,   1.49682e-12,   6.37643e-20,
     8.40415e+02,   3.38821e-05,   1.40504e-12,   5.98544e-20,
     6.96693e+02,   2.80878e-05,   1.16476e-12,   4.96186e-20,
     5.18202e+02,   2.08918e-05,   8.66351e-13,   3.69065e-20,
     3.52864e+02,   1.42260e-05,   5.89932e-13,   2.51310e-20,
     2.24210e+02,   9.03922e-06,   3.74843e-13,   1.59683e-20,
     1.35138e+02,   5.44820e-06,   2.25929e-13,   9.62454e-21,
     7.83119e+01,   3.15721e-06,   1.30925e-13,   5.57739e-21,
     4.41553e+01,   1.78016e-06,   7.38206e-14,   3.14475e-21,
     2.45888e+01,   9.91321e-07,   4.11086e-14,   1.75122e-21  };

  double F15_kpFe3O4[] = 
  {  1.47700e-02,   5.95675e-10,   2.47138e-17,   1.05350e-24,
     2.47694e-02,   9.98953e-10,   4.14453e-17,   1.76673e-24,
     3.73580e-02,   1.50665e-09,   6.25090e-17,   2.66463e-24,
     5.32060e-02,   2.14580e-09,   8.90266e-17,   3.79502e-24,
     8.50036e-02,   3.42820e-09,   1.42232e-16,   6.06304e-24,
     1.29213e-01,   5.21117e-09,   2.16205e-16,   9.21636e-24,
     2.00170e-01,   8.07286e-09,   3.34932e-16,   1.42775e-23,
     3.15560e-01,   1.27266e-08,   5.28008e-16,   2.25079e-23,
     5.01384e-01,   2.02209e-08,   8.38937e-16,   3.57622e-23,
     7.88907e-01,   3.18167e-08,   1.32003e-15,   5.62703e-23,
     1.24250e+00,   5.01101e-08,   2.07900e-15,   8.86237e-23,
     1.95225e+00,   7.87342e-08,   3.26658e-15,   1.39248e-22,
     3.04002e+00,   1.22604e-07,   5.08669e-15,   2.16835e-22,
     4.68918e+00,   1.89115e-07,   7.84613e-15,   3.34465e-22,
     7.12599e+00,   2.87392e-07,   1.19235e-14,   5.08275e-22,
     1.05834e+01,   4.26829e-07,   1.77086e-14,   7.54881e-22,
     1.52356e+01,   6.14453e-07,   2.54928e-14,   1.08671e-21,
     2.13345e+01,   8.60424e-07,   3.56978e-14,   1.52173e-21,
     2.98061e+01,   1.20208e-06,   4.98728e-14,   2.12597e-21,
     4.27642e+01,   1.72468e-06,   7.15548e-14,   3.05024e-21,
     6.30370e+01,   2.54229e-06,   1.05476e-13,   4.49623e-21,
     9.29361e+01,   3.74812e-06,   1.55505e-13,   6.62884e-21,
     1.32987e+02,   5.36337e-06,   2.22519e-13,   9.48555e-21,
     1.82150e+02,   7.34614e-06,   3.04782e-13,   1.29922e-20,
     2.40388e+02,   9.69485e-06,   4.02226e-13,   1.71461e-20,
     3.12065e+02,   1.25856e-05,   5.22160e-13,   2.22586e-20,
     4.08414e+02,   1.64714e-05,   6.83375e-13,   2.91309e-20,
     5.49591e+02,   2.21650e-05,   9.19598e-13,   3.92006e-20,
     7.67451e+02,   3.09514e-05,   1.28413e-12,   5.47399e-20,
     1.10725e+03,   4.46557e-05,   1.85271e-12,   7.89771e-20,
     1.62060e+03,   6.53589e-05,   2.71166e-12,   1.15592e-19,
     2.33999e+03,   9.43719e-05,   3.91537e-12,   1.66904e-19,
     3.24367e+03,   1.30818e-04,   5.42745e-12,   2.31361e-19,
     4.25716e+03,   1.71692e-04,   7.12327e-12,   3.03651e-19,
     5.34010e+03,   2.15367e-04,   8.93529e-12,   3.80893e-19  };

  double F15_kpAC[] = 
  {  3.27956e-01,   3.75639e-06,   5.26403e-11,   8.37270e-16,
     4.38770e-01,   5.02579e-06,   7.04309e-11,   1.12026e-15,
     5.78277e-01,   6.62387e-06,   9.28278e-11,   1.47653e-15,
     7.53906e-01,   8.63573e-06,   1.21024e-10,   1.92505e-15,
     1.04036e+00,   1.19174e-05,   1.67020e-10,   2.65676e-15,
     1.41778e+00,   1.62414e-05,   2.27626e-10,   3.62091e-15,
     1.95366e+00,   2.23811e-05,   3.13688e-10,   4.99012e-15,
     2.71664e+00,   3.11237e-05,   4.36248e-10,   6.94018e-15,
     3.79935e+00,   4.35320e-05,   6.10219e-10,   9.70860e-15,
     5.30234e+00,   6.07596e-05,   8.51798e-10,   1.35534e-14,
     7.38743e+00,   8.46654e-05,   1.18710e-09,   1.88910e-14,
     1.02340e+01,   1.17314e-04,   1.64516e-09,   2.61851e-14,
     1.40739e+01,   1.61376e-04,   2.26365e-09,   3.60378e-14,
     1.92619e+01,   2.20951e-04,   3.10040e-09,   4.93757e-14,
     2.62735e+01,   3.01546e-04,   4.23341e-09,   6.74511e-14,
     3.57407e+01,   4.10519e-04,   5.76727e-09,   9.19507e-14,
     4.85587e+01,   5.58361e-04,   7.85200e-09,   1.25306e-13,
     6.60737e+01,   7.60948e-04,   1.07159e-08,   1.71237e-13,
     9.02145e+01,   1.04125e-03,   1.46920e-08,   2.35213e-13,
     1.23479e+02,   1.42959e-03,   2.02271e-08,   3.24675e-13,
     1.68976e+02,   1.96503e-03,   2.79128e-08,   4.49723e-13,
     2.30962e+02,   2.70317e-03,   3.86182e-08,   6.25581e-13,
     3.16043e+02,   3.73154e-03,   5.37252e-08,   8.76701e-13,
     4.33583e+02,   5.17180e-03,   7.51313e-08,   1.23638e-12,
     5.93567e+02,   7.14700e-03,   1.04677e-07,   1.73582e-12,
     8.03334e+02,   9.73553e-03,   1.43375e-07,   2.38969e-12,
     1.06890e+03,   1.29850e-02,   1.91575e-07,   3.19823e-12,
     1.40232e+03,   1.70036e-02,   2.50333e-07,   4.17041e-12,
     1.82952e+03,   2.20557e-02,   3.22829e-07,   5.34792e-12,
     2.39537e+03,   2.86193e-02,   4.15124e-07,   6.81663e-12,
     3.17137e+03,   3.74648e-02,   5.37088e-07,   8.71845e-12,
     4.27121e+03,   4.98190e-02,   7.04376e-07,   1.12776e-11,
     5.87303e+03,   6.75884e-02,   9.40996e-07,   1.48329e-11,
     8.22455e+03,   9.33574e-02,   1.27863e-06,   1.98186e-11,
     1.15719e+04,   1.29494e-01,   1.74429e-06,   2.65764e-11  };

  double F15_kpSiO2D[] = 
  {  7.60360e-02,   3.06536e-09,   1.27110e-16,   5.41456e-24,
     9.07207e-02,   3.65737e-09,   1.51659e-16,   6.46027e-24,
     1.09208e-01,   4.40266e-09,   1.82564e-16,   7.77673e-24,
     1.32481e-01,   5.34093e-09,   2.21471e-16,   9.43407e-24,
     1.58907e-01,   6.40629e-09,   2.65648e-16,   1.13159e-23,
     1.91565e-01,   7.72285e-09,   3.20241e-16,   1.36414e-23,
     2.30490e-01,   9.29212e-09,   3.85313e-16,   1.64133e-23,
     2.76795e-01,   1.11589e-08,   4.62722e-16,   1.97107e-23,
     3.33074e-01,   1.34277e-08,   5.56804e-16,   2.37184e-23,
     4.05325e-01,   1.63405e-08,   6.77586e-16,   2.88634e-23,
     5.08160e-01,   2.04863e-08,   8.49498e-16,   3.61864e-23,
     6.72472e-01,   2.71105e-08,   1.12418e-15,   4.78871e-23,
     9.48549e-01,   3.82404e-08,   1.58570e-15,   6.75467e-23,
     1.41787e+00,   5.71610e-08,   2.37028e-15,   1.00968e-22,
     2.19502e+00,   8.84912e-08,   3.66944e-15,   1.56308e-22,
     3.46719e+00,   1.39778e-07,   5.79615e-15,   2.46901e-22,
     5.76852e+00,   2.32556e-07,   9.64331e-15,   4.10780e-22,
     1.17194e+01,   4.72463e-07,   1.95915e-14,   8.34544e-22,
     3.16449e+01,   1.27575e-06,   5.29013e-14,   2.25345e-21,
     8.68296e+01,   3.50050e-06,   1.45154e-13,   6.18318e-21,
     1.92300e+02,   7.75250e-06,   3.21470e-13,   1.36938e-20,
     3.36231e+02,   1.35550e-05,   5.62081e-13,   2.39432e-20,
     5.05825e+02,   2.03921e-05,   8.45594e-13,   3.60201e-20,
     7.20624e+02,   2.90517e-05,   1.20468e-12,   5.13160e-20,
     9.77376e+02,   3.94025e-05,   1.63389e-12,   6.95995e-20,
     1.18646e+03,   4.78315e-05,   1.98341e-12,   8.44881e-20,
     1.23845e+03,   4.99275e-05,   2.07033e-12,   8.81904e-20,
     1.11188e+03,   4.48251e-05,   1.85875e-12,   7.91777e-20,
     8.76396e+02,   3.53316e-05,   1.46508e-12,   6.24086e-20,
     6.22207e+02,   2.50840e-05,   1.04015e-12,   4.43077e-20,
     4.07290e+02,   1.64197e-05,   6.80872e-13,   2.90033e-20,
     2.50573e+02,   1.01017e-05,   4.18886e-13,   1.78434e-20,
     1.47073e+02,   5.92921e-06,   2.45865e-13,   1.04732e-20,
     8.33122e+01,   3.35870e-06,   1.39274e-13,   5.93271e-21,
     4.59585e+01,   1.85280e-06,   7.68295e-14,   3.27273e-21  };

  double F15_kpAl2O3[] = 
  {  9.93250e-04,   4.00425e-11,   1.66043e-18,   7.07299e-26,
     1.81240e-03,   7.30662e-11,   3.02981e-18,   1.29062e-25,
     2.84365e-03,   1.14641e-10,   4.75376e-18,   2.02498e-25,
     4.14191e-03,   1.66980e-10,   6.92409e-18,   2.94948e-25,
     7.18271e-03,   2.89568e-10,   1.20074e-17,   5.11485e-25,
     1.13364e-02,   4.57021e-10,   1.89512e-17,   8.07269e-25,
     1.77361e-02,   7.15022e-10,   2.96496e-17,   1.26300e-24,
     2.59477e-02,   1.04607e-09,   4.33772e-17,   1.84775e-24,
     3.45425e-02,   1.39257e-09,   5.77452e-17,   2.45979e-24,
     4.22006e-02,   1.70130e-09,   7.05474e-17,   3.00513e-24,
     4.71420e-02,   1.90051e-09,   7.88079e-17,   3.35701e-24,
     4.91934e-02,   1.98321e-09,   8.22373e-17,   3.50309e-24,
     5.05162e-02,   2.03654e-09,   8.44486e-17,   3.59729e-24,
     5.78201e-02,   2.33100e-09,   9.66587e-17,   4.11741e-24,
     8.84237e-02,   3.56477e-09,   1.47819e-16,   6.29670e-24,
     1.78786e-01,   7.20769e-09,   2.98879e-16,   1.27315e-23,
     4.36404e-01,   1.75935e-08,   7.29542e-16,   3.10766e-23,
     1.63796e+00,   6.60337e-08,   2.73820e-15,   1.16640e-22,
     8.50817e+00,   3.43004e-07,   1.42232e-14,   6.05872e-22,
     3.92751e+01,   1.58336e-06,   6.56567e-14,   2.79680e-21,
     1.41436e+02,   5.70194e-06,   2.36441e-13,   1.00718e-20,
     3.83709e+02,   1.54691e-05,   6.41451e-13,   2.73241e-20,
     7.70411e+02,   3.10588e-05,   1.28791e-12,   5.48614e-20,
     1.16399e+03,   4.69258e-05,   1.94586e-12,   8.28883e-20,
     1.37566e+03,   5.54594e-05,   2.29972e-12,   9.79619e-20,
     1.33070e+03,   5.36466e-05,   2.22455e-12,   9.47599e-20,
     1.09978e+03,   4.43371e-05,   1.83851e-12,   7.83159e-20,
     8.05638e+02,   3.24790e-05,   1.34680e-12,   5.73700e-20,
     5.38690e+02,   2.17171e-05,   9.00535e-13,   3.83604e-20,
     3.36338e+02,   1.35593e-05,   5.62261e-13,   2.39508e-20,
     1.99460e+02,   8.04115e-06,   3.33440e-13,   1.42037e-20,
     1.13787e+02,   4.58728e-06,   1.90220e-13,   8.10285e-21,
     6.30411e+01,   2.54148e-06,   1.05387e-13,   4.48919e-21,
     3.41529e+01,   1.37686e-06,   5.70939e-14,   2.43205e-21,
     1.81893e+01,   7.33293e-07,   3.04073e-14,   1.29527e-21  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpFeM     [itab0] = F15_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = F15_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = F15_kpMgSiO3  [itab];
      my_rates->SN0_kpFe3O4   [itab0] = F15_kpFe3O4   [itab];
      my_rates->SN0_kpAC      [itab0] = F15_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = F15_kpSiO2D   [itab];
      my_rates->SN0_kpAl2O3   [itab0] = F15_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_F50(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   2.79167e-01;
  my_rates->SN0_XO [iSN] =   7.20575e-01;
  my_rates->SN0_XMg[iSN] =   2.49794e-04;
  my_rates->SN0_XAl[iSN] =   1.66468e-08;
  my_rates->SN0_XSi[iSN] =   4.01099e-06;
  my_rates->SN0_XS [iSN] =   0.00000e+00;
  my_rates->SN0_XFe[iSN] =   4.15804e-06;

  my_rates->SN0_fC [iSN] =   2.79057e-01;
  my_rates->SN0_fO [iSN] =   7.20575e-01;
  my_rates->SN0_fMg[iSN] =   2.49793e-04;
  my_rates->SN0_fAl[iSN] =   1.66468e-08;
  my_rates->SN0_fSi[iSN] =   4.01058e-06;
  my_rates->SN0_fS [iSN] =   0.00000e+00;
  my_rates->SN0_fFe[iSN] =   4.15804e-06;

  my_rates->SN0_fFeM     [iSN] =   2.33171e-24;
  my_rates->SN0_fMg2SiO4 [iSN] =   2.62486e-10;
  my_rates->SN0_fMgSiO3  [iSN] =   1.21446e-09;
  my_rates->SN0_fFe3O4   [iSN] =   2.41799e-13;
  my_rates->SN0_fAC      [iSN] =   1.09849e-04;
  my_rates->SN0_fSiO2D   [iSN] =   3.41863e-11;
  my_rates->SN0_fAl2O3   [iSN] =   2.53950e-17;

  itab0 = 3 * iSN;
  my_rates->SN0_r0FeM     [itab0 + 0] =   4.02891e-08;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   1.68491e-07;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   1.33003e-07;
  my_rates->SN0_r0Fe3O4   [itab0 + 0] =   5.89806e-08;
  my_rates->SN0_r0AC      [itab0 + 0] =   6.81790e-07;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   9.81613e-08;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   4.03146e-08;

  my_rates->SN0_r0FeM     [itab0 + 1] =   1.67016e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   3.02634e-14;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.84568e-14;
  my_rates->SN0_r0Fe3O4   [itab0 + 1] =   3.51732e-15;
  my_rates->SN0_r0AC      [itab0 + 1] =   6.53175e-13;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   9.72845e-15;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   1.67172e-15;

  my_rates->SN0_r0FeM     [itab0 + 2] =   7.11339e-23;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   5.60369e-21;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   2.62630e-21;
  my_rates->SN0_r0Fe3O4   [itab0 + 2] =   2.11807e-22;
  my_rates->SN0_r0AC      [itab0 + 2] =   7.65748e-19;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   9.68327e-22;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   7.12107e-23;

  NTd =            35;
 Nmom =             4;

  double F50_kpFeM[] = 
  {  1.23620e-05,   4.98882e-13,   2.07138e-20,   8.83538e-28,
     2.19537e-05,   8.85960e-13,   3.67854e-20,   1.56906e-27,
     3.40288e-05,   1.37326e-12,   5.70184e-20,   2.43208e-27,
     4.92305e-05,   1.98674e-12,   8.24901e-20,   3.51856e-27,
     8.08507e-05,   3.26278e-12,   1.35471e-19,   5.77841e-27,
     1.25019e-04,   5.04514e-12,   2.09472e-19,   8.93477e-27,
     1.96584e-04,   7.93293e-12,   3.29363e-19,   1.40481e-26,
     3.14489e-04,   1.26902e-11,   5.26848e-19,   2.24703e-26,
     5.06846e-04,   2.04508e-11,   8.48987e-19,   3.62076e-26,
     8.07280e-04,   3.25711e-11,   1.35207e-18,   5.76597e-26,
     1.28667e-03,   5.19094e-11,   2.15469e-18,   9.18830e-26,
     2.05240e-03,   8.27959e-11,   3.43651e-18,   1.46534e-25,
     3.27023e-03,   1.31912e-10,   5.47463e-18,   2.33420e-25,
     5.23895e-03,   2.11300e-10,   8.76840e-18,   3.73817e-25,
     8.45019e-03,   3.40772e-10,   1.41393e-17,   6.02717e-25,
     1.37157e-02,   5.53037e-10,   2.29434e-17,   9.77885e-25,
     2.24099e-02,   9.03468e-10,   3.74764e-17,   1.59710e-24,
     3.70041e-02,   1.49164e-09,   6.18661e-17,   2.63618e-24,
     6.21584e-02,   2.50530e-09,   1.03896e-16,   4.42660e-24,
     1.07033e-01,   4.31351e-09,   1.78864e-16,   7.62000e-24,
     1.90089e-01,   7.66007e-09,   3.17606e-16,   1.35297e-23,
     3.49470e-01,   1.40818e-08,   5.83831e-16,   2.48691e-23,
     6.64947e-01,   2.67926e-08,   1.11077e-15,   4.73128e-23,
     1.30413e+00,   5.25455e-08,   2.17837e-15,   9.27841e-23,
     2.61640e+00,   1.05417e-07,   4.37015e-15,   1.86136e-22,
     5.31790e+00,   2.14259e-07,   8.88223e-15,   3.78313e-22,
     1.08366e+01,   4.36604e-07,   1.80995e-14,   7.70889e-22,
     2.19132e+01,   8.82874e-07,   3.65995e-14,   1.55883e-21,
     4.35354e+01,   1.75402e-06,   7.27125e-14,   3.09693e-21,
     8.42362e+01,   3.39383e-06,   1.40690e-13,   5.99219e-21,
     1.57704e+02,   6.35380e-06,   2.63395e-13,   1.12183e-20,
     2.84822e+02,   1.14753e-05,   4.75704e-13,   2.02608e-20,
     4.96653e+02,   2.00098e-05,   8.29500e-13,   3.53294e-20,
     8.39966e+02,   3.38417e-05,   1.40289e-12,   5.97510e-20,
     1.38932e+03,   5.59749e-05,   2.32042e-12,   9.88295e-20  };

  double F50_kpMg2SiO4[] = 
  {  1.05240e-01,   1.77320e-08,   3.18492e-15,   5.89732e-22,
     1.32588e-01,   2.23399e-08,   4.01256e-15,   7.42980e-22,
     1.67016e-01,   2.81408e-08,   5.05449e-15,   9.35908e-22,
     2.10360e-01,   3.54438e-08,   6.36620e-15,   1.17879e-21,
     2.71887e-01,   4.58106e-08,   8.22824e-15,   1.52357e-21,
     3.55694e-01,   5.99313e-08,   1.07645e-14,   1.99320e-21,
     4.84932e-01,   8.17069e-08,   1.46757e-14,   2.71741e-21,
     6.99767e-01,   1.17905e-07,   2.11774e-14,   3.92128e-21,
     1.05860e+00,   1.78364e-07,   3.20367e-14,   5.93204e-21,
     1.62902e+00,   2.74475e-07,   4.92997e-14,   9.12851e-21,
     2.54260e+00,   4.28406e-07,   7.69478e-14,   1.42479e-20,
     3.96488e+00,   6.68047e-07,   1.19991e-13,   2.22179e-20,
     6.10630e+00,   1.02886e-06,   1.84797e-13,   3.42178e-20,
     9.28766e+00,   1.56489e-06,   2.81076e-13,   5.20451e-20,
     1.39265e+01,   2.34649e-06,   4.21463e-13,   7.80396e-20,
     2.05382e+01,   3.46051e-06,   6.21558e-13,   1.15090e-19,
     3.00649e+01,   5.06568e-06,   9.09868e-13,   1.68475e-19,
     4.55102e+01,   7.66807e-06,   1.37729e-12,   2.55025e-19,
     7.47839e+01,   1.26004e-05,   2.26322e-12,   4.19065e-19,
     1.29621e+02,   2.18400e-05,   3.92277e-12,   7.26354e-19,
     2.14820e+02,   3.61953e-05,   6.50119e-12,   1.20378e-18,
     3.20002e+02,   5.39175e-05,   9.68435e-12,   1.79319e-18,
     4.29768e+02,   7.24122e-05,   1.30063e-11,   2.40829e-18,
     5.30829e+02,   8.94401e-05,   1.60647e-11,   2.97460e-18,
     5.99696e+02,   1.01044e-04,   1.81489e-11,   3.36051e-18,
     6.06539e+02,   1.02197e-04,   1.83560e-11,   3.39886e-18,
     5.43264e+02,   9.15353e-05,   1.64410e-11,   3.04428e-18,
     4.33547e+02,   7.30489e-05,   1.31206e-11,   2.42946e-18,
     3.13324e+02,   5.27923e-05,   9.48225e-12,   1.75577e-18,
     2.09008e+02,   3.52161e-05,   6.32532e-12,   1.17122e-18,
     1.31150e+02,   2.20976e-05,   3.96905e-12,   7.34924e-19,
     7.90692e+01,   1.33225e-05,   2.39291e-12,   4.43080e-19,
     4.73403e+01,   7.97645e-06,   1.43269e-12,   2.65282e-19,
     2.96433e+01,   4.99469e-06,   8.97122e-13,   1.66115e-19,
     2.04754e+01,   3.45001e-06,   6.19680e-13,   1.14743e-19  };

  double F50_kpMgSiO3[] = 
  {  2.19890e-02,   2.92460e-09,   4.05846e-16,   5.77498e-23,
     3.90612e-02,   5.19526e-09,   7.20944e-16,   1.02587e-22,
     6.05539e-02,   8.05385e-09,   1.11763e-15,   1.59033e-22,
     8.76116e-02,   1.16526e-08,   1.61703e-15,   2.30095e-22,
     1.43288e-01,   1.90577e-08,   2.64462e-15,   3.76317e-22,
     2.19266e-01,   2.91631e-08,   4.04694e-15,   5.75860e-22,
     3.36256e-01,   4.47230e-08,   6.20619e-15,   8.83109e-22,
     5.14336e-01,   6.84082e-08,   9.49297e-15,   1.35080e-21,
     7.97216e-01,   1.06032e-07,   1.47140e-14,   2.09373e-21,
     1.25414e+00,   1.66804e-07,   2.31474e-14,   3.29375e-21,
     2.03450e+00,   2.70594e-07,   3.75503e-14,   5.34322e-21,
     3.34648e+00,   4.45091e-07,   6.17651e-14,   8.78886e-21,
     5.45893e+00,   7.26054e-07,   1.00754e-13,   1.43368e-20,
     8.82117e+00,   1.17324e-06,   1.62810e-13,   2.31671e-20,
     1.41825e+01,   1.88631e-06,   2.61763e-13,   3.72475e-20,
     2.28419e+01,   3.03804e-06,   4.21588e-13,   5.99898e-20,
     3.71178e+01,   4.93677e-06,   6.85074e-13,   9.74826e-20,
     6.14273e+01,   8.17002e-06,   1.13375e-12,   1.61327e-19,
     1.03847e+02,   1.38119e-05,   1.91667e-12,   2.72733e-19,
     1.75507e+02,   2.33429e-05,   3.23928e-12,   4.60933e-19,
     2.82060e+02,   3.75148e-05,   5.20591e-12,   7.40775e-19,
     4.14519e+02,   5.51322e-05,   7.65067e-12,   1.08865e-18,
     5.59962e+02,   7.44767e-05,   1.03351e-11,   1.47063e-18,
     7.11026e+02,   9.45685e-05,   1.31232e-11,   1.86737e-18,
     8.40809e+02,   1.11830e-04,   1.55186e-11,   2.20822e-18,
     8.95315e+02,   1.19080e-04,   1.65246e-11,   2.35137e-18,
     8.40416e+02,   1.11778e-04,   1.55114e-11,   2.20719e-18,
     6.96694e+02,   9.26624e-05,   1.28587e-11,   1.82973e-18,
     5.18204e+02,   6.89226e-05,   9.56436e-12,   1.36096e-18,
     3.52865e+02,   4.69321e-05,   6.51275e-12,   9.26731e-19,
     2.24211e+02,   2.98208e-05,   4.13822e-12,   5.88848e-19,
     1.35138e+02,   1.79738e-05,   2.49421e-12,   3.54913e-19,
     7.83122e+01,   1.04158e-05,   1.44539e-12,   2.05672e-19,
     4.41556e+01,   5.87282e-06,   8.14969e-13,   1.15966e-19,
     2.45896e+01,   3.27050e-06,   4.53847e-13,   6.45802e-20  };

  double F50_kpFe3O4[] = 
  {  1.47700e-02,   8.71144e-10,   5.19508e-17,   3.12839e-24,
     2.47694e-02,   1.46092e-09,   8.71220e-17,   5.24635e-24,
     3.73580e-02,   2.20340e-09,   1.31400e-16,   7.91270e-24,
     5.32060e-02,   3.13813e-09,   1.87143e-16,   1.12694e-23,
     8.50036e-02,   5.01357e-09,   2.98985e-16,   1.80044e-23,
     1.29213e-01,   7.62106e-09,   4.54483e-16,   2.73683e-23,
     2.00170e-01,   1.18061e-08,   7.04061e-16,   4.23974e-23,
     3.15560e-01,   1.86119e-08,   1.10993e-15,   6.68379e-23,
     5.01384e-01,   2.95719e-08,   1.76353e-15,   1.06197e-22,
     7.88907e-01,   4.65303e-08,   2.77484e-15,   1.67096e-22,
     1.24250e+00,   7.32834e-08,   4.37027e-15,   2.63171e-22,
     1.95225e+00,   1.15145e-07,   6.86667e-15,   4.13500e-22,
     3.04002e+00,   1.79302e-07,   1.06927e-14,   6.43899e-22,
     4.68918e+00,   2.76571e-07,   1.64933e-14,   9.93203e-22,
     7.12599e+00,   4.20295e-07,   2.50644e-14,   1.50934e-21,
     1.05834e+01,   6.24215e-07,   3.72252e-14,   2.24164e-21,
     1.52356e+01,   8.98605e-07,   5.35885e-14,   3.22701e-21,
     2.13345e+01,   1.25832e-06,   7.50404e-14,   4.51881e-21,
     2.98061e+01,   1.75798e-06,   1.04837e-13,   6.31314e-21,
     4.27642e+01,   2.52226e-06,   1.50415e-13,   9.05777e-21,
     6.30370e+01,   3.71796e-06,   2.21721e-13,   1.33517e-20,
     9.29361e+01,   5.48143e-06,   3.26886e-13,   1.96845e-20,
     1.32987e+02,   7.84365e-06,   4.67758e-13,   2.81676e-20,
     1.82150e+02,   1.07433e-05,   6.40681e-13,   3.85808e-20,
     2.40388e+02,   1.41782e-05,   8.45520e-13,   5.09158e-20,
     3.12065e+02,   1.84058e-05,   1.09763e-12,   6.60977e-20,
     4.08414e+02,   2.40885e-05,   1.43652e-12,   8.65051e-20,
     5.49591e+02,   3.24152e-05,   1.93309e-12,   1.16407e-19,
     7.67451e+02,   4.52647e-05,   2.69937e-12,   1.62552e-19,
     1.10725e+03,   6.53066e-05,   3.89457e-12,   2.34525e-19,
     1.62060e+03,   9.55840e-05,   5.70017e-12,   3.43255e-19,
     2.33999e+03,   1.38014e-04,   8.23049e-12,   4.95627e-19,
     3.24369e+03,   1.91315e-04,   1.14091e-11,   6.87037e-19,
     4.25718e+03,   2.51091e-04,   1.49739e-11,   9.01702e-19,
     5.34014e+03,   3.14965e-04,   1.87830e-11,   1.13108e-18  };

  double F50_kpAC[] = 
  {  3.27960e-01,   2.23600e-07,   2.14215e-13,   2.51135e-19,
     4.38752e-01,   2.99136e-07,   2.86582e-13,   3.35973e-19,
     5.78230e-01,   3.94231e-07,   3.77685e-13,   4.42778e-19,
     7.53823e-01,   5.13949e-07,   4.92378e-13,   5.77238e-19,
     1.04013e+00,   7.09149e-07,   6.79388e-13,   7.96479e-19,
     1.41735e+00,   9.66336e-07,   9.25781e-13,   1.08534e-18,
     1.95292e+00,   1.33148e-06,   1.27560e-12,   1.49545e-18,
     2.71530e+00,   1.85127e-06,   1.77357e-12,   2.07925e-18,
     3.79675e+00,   2.58859e-06,   2.47995e-12,   2.90737e-18,
     5.29742e+00,   3.61173e-06,   3.46016e-12,   4.05652e-18,
     7.37832e+00,   5.03048e-06,   4.81937e-12,   5.65000e-18,
     1.02167e+01,   6.96570e-06,   6.67339e-12,   7.82359e-18,
     1.40420e+01,   9.57377e-06,   9.17203e-12,   1.07529e-17,
     1.92020e+01,   1.30918e-05,   1.25425e-11,   1.47044e-17,
     2.61614e+01,   1.78368e-05,   1.70885e-11,   2.00340e-17,
     3.55303e+01,   2.42246e-05,   2.32084e-11,   2.72089e-17,
     4.81604e+01,   3.28360e-05,   3.14588e-11,   3.68819e-17,
     6.53156e+01,   4.45330e-05,   4.26657e-11,   5.00213e-17,
     8.87641e+01,   6.05213e-05,   5.79846e-11,   6.79824e-17,
     1.20700e+02,   8.22973e-05,   7.88495e-11,   9.24471e-17,
     1.63611e+02,   1.11558e-04,   1.06888e-10,   1.25326e-16,
     2.20556e+02,   1.50393e-04,   1.44104e-10,   1.68969e-16,
     2.96069e+02,   2.01894e-04,   1.93464e-10,   2.26864e-16,
     3.96959e+02,   2.70713e-04,   2.59434e-10,   3.04254e-16,
     5.31398e+02,   3.62437e-04,   3.47381e-10,   4.07456e-16,
     7.06744e+02,   4.82105e-04,   4.62166e-10,   5.42204e-16,
     9.30503e+02,   6.34887e-04,   6.08797e-10,   7.14450e-16,
     1.21574e+03,   8.29802e-04,   7.96038e-10,   9.34622e-16,
     1.58603e+03,   1.08311e-03,   1.03970e-09,   1.22157e-15,
     2.07753e+03,   1.41988e-03,   1.36429e-09,   1.60466e-15,
     2.74067e+03,   1.87531e-03,   1.80443e-09,   2.12568e-15,
     3.64502e+03,   2.49832e-03,   2.40876e-09,   2.84400e-15,
     4.89065e+03,   3.35972e-03,   3.24810e-09,   3.84654e-15,
     6.62881e+03,   4.56615e-03,   4.42856e-09,   5.26292e-15,
     9.09708e+03,   6.28233e-03,   6.11057e-09,   7.28448e-15  };

  double F50_kpSiO2D[] = 
  {  7.60360e-02,   7.46380e-09,   7.39712e-16,   7.36277e-23,
     9.07207e-02,   8.90526e-09,   8.82572e-16,   8.78473e-23,
     1.09208e-01,   1.07200e-08,   1.06242e-15,   1.05749e-22,
     1.32481e-01,   1.30045e-08,   1.28884e-15,   1.28285e-22,
     1.58907e-01,   1.55986e-08,   1.54592e-15,   1.53874e-22,
     1.91565e-01,   1.88042e-08,   1.86363e-15,   1.85497e-22,
     2.30490e-01,   2.26252e-08,   2.24231e-15,   2.23190e-22,
     2.76795e-01,   2.71706e-08,   2.69279e-15,   2.68028e-22,
     3.33074e-01,   3.26950e-08,   3.24029e-15,   3.22524e-22,
     4.05325e-01,   3.97872e-08,   3.94318e-15,   3.92487e-22,
     5.08160e-01,   4.98817e-08,   4.94361e-15,   4.92065e-22,
     6.72472e-01,   6.60108e-08,   6.54211e-15,   6.51173e-22,
     9.48549e-01,   9.31109e-08,   9.22792e-15,   9.18506e-22,
     1.41787e+00,   1.39180e-07,   1.37937e-14,   1.37297e-21,
     2.19502e+00,   2.15466e-07,   2.13541e-14,   2.12549e-21,
     3.46719e+00,   3.40344e-07,   3.37304e-14,   3.35738e-21,
     5.76852e+00,   5.66246e-07,   5.61188e-14,   5.58582e-21,
     1.17194e+01,   1.15039e-06,   1.14012e-13,   1.13482e-20,
     3.16449e+01,   3.10631e-06,   3.07856e-13,   3.06426e-20,
     8.68296e+01,   8.52331e-06,   8.44717e-13,   8.40795e-20,
     1.92300e+02,   1.88764e-05,   1.87078e-12,   1.86209e-19,
     3.36231e+02,   3.30049e-05,   3.27100e-12,   3.25581e-19,
     5.05825e+02,   4.96525e-05,   4.92089e-12,   4.89804e-19,
     7.20624e+02,   7.07374e-05,   7.01055e-12,   6.97800e-19,
     9.77376e+02,   9.59406e-05,   9.50836e-12,   9.46420e-19,
     1.18646e+03,   1.16464e-04,   1.15424e-11,   1.14888e-18,
     1.23845e+03,   1.21568e-04,   1.20482e-11,   1.19922e-18,
     1.11188e+03,   1.09144e-04,   1.08169e-11,   1.07667e-18,
     8.76396e+02,   8.60282e-05,   8.52597e-12,   8.48638e-19,
     6.22207e+02,   6.10766e-05,   6.05311e-12,   6.02500e-19,
     4.07290e+02,   3.99801e-05,   3.96230e-12,   3.94390e-19,
     2.50573e+02,   2.45966e-05,   2.43769e-12,   2.42636e-19,
     1.47073e+02,   1.44369e-05,   1.43080e-12,   1.42415e-19,
     8.33122e+01,   8.17804e-06,   8.10499e-13,   8.06735e-20,
     4.59586e+01,   4.51135e-06,   4.47106e-13,   4.45029e-20  };

  double F50_kpAl2O3[] = 
  {  9.93250e-04,   4.00425e-11,   1.66043e-18,   7.07300e-26,
     1.81240e-03,   7.30662e-11,   3.02982e-18,   1.29062e-25,
     2.84365e-03,   1.14641e-10,   4.75377e-18,   2.02498e-25,
     4.14191e-03,   1.66980e-10,   6.92410e-18,   2.94948e-25,
     7.18271e-03,   2.89568e-10,   1.20074e-17,   5.11486e-25,
     1.13364e-02,   4.57021e-10,   1.89512e-17,   8.07270e-25,
     1.77361e-02,   7.15022e-10,   2.96496e-17,   1.26300e-24,
     2.59477e-02,   1.04607e-09,   4.33772e-17,   1.84776e-24,
     3.45425e-02,   1.39257e-09,   5.77453e-17,   2.45980e-24,
     4.22006e-02,   1.70130e-09,   7.05475e-17,   3.00514e-24,
     4.71420e-02,   1.90051e-09,   7.88080e-17,   3.35701e-24,
     4.91934e-02,   1.98321e-09,   8.22374e-17,   3.50310e-24,
     5.05162e-02,   2.03654e-09,   8.44487e-17,   3.59729e-24,
     5.78201e-02,   2.33100e-09,   9.66588e-17,   4.11741e-24,
     8.84237e-02,   3.56477e-09,   1.47819e-16,   6.29671e-24,
     1.78786e-01,   7.20770e-09,   2.98880e-16,   1.27315e-23,
     4.36404e-01,   1.75935e-08,   7.29543e-16,   3.10766e-23,
     1.63796e+00,   6.60337e-08,   2.73820e-15,   1.16640e-22,
     8.50817e+00,   3.43004e-07,   1.42232e-14,   6.05873e-22,
     3.92751e+01,   1.58336e-06,   6.56568e-14,   2.79681e-21,
     1.41436e+02,   5.70194e-06,   2.36441e-13,   1.00718e-20,
     3.83709e+02,   1.54691e-05,   6.41452e-13,   2.73242e-20,
     7.70411e+02,   3.10588e-05,   1.28791e-12,   5.48615e-20,
     1.16399e+03,   4.69258e-05,   1.94586e-12,   8.28885e-20,
     1.37566e+03,   5.54594e-05,   2.29972e-12,   9.79620e-20,
     1.33070e+03,   5.36466e-05,   2.22455e-12,   9.47600e-20,
     1.09978e+03,   4.43372e-05,   1.83852e-12,   7.83160e-20,
     8.05638e+02,   3.24790e-05,   1.34680e-12,   5.73701e-20,
     5.38690e+02,   2.17171e-05,   9.00536e-13,   3.83605e-20,
     3.36338e+02,   1.35593e-05,   5.62261e-13,   2.39509e-20,
     1.99460e+02,   8.04115e-06,   3.33440e-13,   1.42037e-20,
     1.13787e+02,   4.58729e-06,   1.90220e-13,   8.10286e-21,
     6.30411e+01,   2.54148e-06,   1.05387e-13,   4.48920e-21,
     3.41529e+01,   1.37686e-06,   5.70940e-14,   2.43205e-21,
     1.81893e+01,   7.33294e-07,   3.04073e-14,   1.29527e-21  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpFeM     [itab0] = F50_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = F50_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = F50_kpMgSiO3  [itab];
      my_rates->SN0_kpFe3O4   [itab0] = F50_kpFe3O4   [itab];
      my_rates->SN0_kpAC      [itab0] = F50_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = F50_kpSiO2D   [itab];
      my_rates->SN0_kpAl2O3   [itab0] = F50_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_F80(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   2.52563e-01;
  my_rates->SN0_XO [iSN] =   7.46061e-01;
  my_rates->SN0_XMg[iSN] =   1.36917e-03;
  my_rates->SN0_XAl[iSN] =   1.55602e-08;
  my_rates->SN0_XSi[iSN] =   3.63906e-06;
  my_rates->SN0_XS [iSN] =   0.00000e+00;
  my_rates->SN0_XFe[iSN] =   2.43915e-06;

  my_rates->SN0_fC [iSN] =   2.43883e-01;
  my_rates->SN0_fO [iSN] =   7.46061e-01;
  my_rates->SN0_fMg[iSN] =   1.36917e-03;
  my_rates->SN0_fAl[iSN] =   1.55602e-08;
  my_rates->SN0_fSi[iSN] =   3.63906e-06;
  my_rates->SN0_fS [iSN] =   0.00000e+00;
  my_rates->SN0_fFe[iSN] =   2.43915e-06;

  my_rates->SN0_fFeM     [iSN] =   3.87590e-26;
  my_rates->SN0_fMg2SiO4 [iSN] =   2.36180e-13;
  my_rates->SN0_fMgSiO3  [iSN] =   2.48190e-12;
  my_rates->SN0_fFe3O4   [iSN] =   3.01120e-15;
  my_rates->SN0_fAC      [iSN] =   8.68025e-03;
  my_rates->SN0_fSiO2D   [iSN] =   3.70132e-14;
  my_rates->SN0_fAl2O3   [iSN] =   3.77811e-18;

  itab0 = 3 * iSN;
  my_rates->SN0_r0FeM     [itab0 + 0] =   4.02891e-08;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   5.88698e-08;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   5.87709e-08;
  my_rates->SN0_r0Fe3O4   [itab0 + 0] =   4.03342e-08;
  my_rates->SN0_r0AC      [itab0 + 0] =   4.22607e-06;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   4.03439e-08;
  my_rates->SN0_r0Al2O3   [itab0 + 0] =   4.03146e-08;

  my_rates->SN0_r0FeM     [itab0 + 1] =   1.67016e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   3.50624e-15;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   3.49547e-15;
  my_rates->SN0_r0Fe3O4   [itab0 + 1] =   1.67365e-15;
  my_rates->SN0_r0AC      [itab0 + 1] =   2.30435e-11;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   1.67461e-15;
  my_rates->SN0_r0Al2O3   [itab0 + 1] =   1.67171e-15;

  my_rates->SN0_r0FeM     [itab0 + 2] =   7.11339e-23;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   2.10950e-22;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   2.10029e-22;
  my_rates->SN0_r0Fe3O4   [itab0 + 2] =   7.13577e-23;
  my_rates->SN0_r0AC      [itab0 + 2] =   1.46801e-16;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   7.14309e-23;
  my_rates->SN0_r0Al2O3   [itab0 + 2] =   7.12106e-23;

  NTd =            35;
 Nmom =             4;

  double F80_kpFeM[] = 
  {  1.23620e-05,   4.98882e-13,   2.07138e-20,   8.83538e-28,
     2.19537e-05,   8.85960e-13,   3.67854e-20,   1.56906e-27,
     3.40288e-05,   1.37326e-12,   5.70184e-20,   2.43208e-27,
     4.92305e-05,   1.98674e-12,   8.24901e-20,   3.51856e-27,
     8.08507e-05,   3.26278e-12,   1.35471e-19,   5.77841e-27,
     1.25019e-04,   5.04514e-12,   2.09472e-19,   8.93477e-27,
     1.96584e-04,   7.93293e-12,   3.29363e-19,   1.40481e-26,
     3.14489e-04,   1.26902e-11,   5.26848e-19,   2.24703e-26,
     5.06846e-04,   2.04508e-11,   8.48987e-19,   3.62076e-26,
     8.07280e-04,   3.25711e-11,   1.35207e-18,   5.76597e-26,
     1.28667e-03,   5.19094e-11,   2.15469e-18,   9.18830e-26,
     2.05240e-03,   8.27959e-11,   3.43651e-18,   1.46534e-25,
     3.27023e-03,   1.31912e-10,   5.47463e-18,   2.33420e-25,
     5.23895e-03,   2.11300e-10,   8.76840e-18,   3.73817e-25,
     8.45019e-03,   3.40772e-10,   1.41393e-17,   6.02717e-25,
     1.37157e-02,   5.53037e-10,   2.29434e-17,   9.77885e-25,
     2.24099e-02,   9.03468e-10,   3.74764e-17,   1.59710e-24,
     3.70041e-02,   1.49164e-09,   6.18661e-17,   2.63618e-24,
     6.21584e-02,   2.50530e-09,   1.03896e-16,   4.42660e-24,
     1.07033e-01,   4.31351e-09,   1.78864e-16,   7.62000e-24,
     1.90089e-01,   7.66007e-09,   3.17606e-16,   1.35297e-23,
     3.49470e-01,   1.40818e-08,   5.83831e-16,   2.48691e-23,
     6.64947e-01,   2.67926e-08,   1.11077e-15,   4.73128e-23,
     1.30413e+00,   5.25455e-08,   2.17837e-15,   9.27841e-23,
     2.61640e+00,   1.05417e-07,   4.37015e-15,   1.86136e-22,
     5.31790e+00,   2.14259e-07,   8.88223e-15,   3.78313e-22,
     1.08366e+01,   4.36604e-07,   1.80995e-14,   7.70889e-22,
     2.19132e+01,   8.82874e-07,   3.65995e-14,   1.55883e-21,
     4.35354e+01,   1.75402e-06,   7.27125e-14,   3.09693e-21,
     8.42362e+01,   3.39383e-06,   1.40690e-13,   5.99219e-21,
     1.57704e+02,   6.35380e-06,   2.63395e-13,   1.12183e-20,
     2.84822e+02,   1.14753e-05,   4.75704e-13,   2.02608e-20,
     4.96653e+02,   2.00098e-05,   8.29500e-13,   3.53294e-20,
     8.39966e+02,   3.38417e-05,   1.40289e-12,   5.97510e-20,
     1.38932e+03,   5.59749e-05,   2.32042e-12,   9.88295e-20  };

  double F80_kpMg2SiO4[] = 
  {  1.05240e-01,   6.19546e-09,   3.68996e-16,   2.22004e-23,
     1.32588e-01,   7.80541e-09,   4.64884e-16,   2.79694e-23,
     1.67016e-01,   9.83223e-09,   5.85599e-16,   3.52321e-23,
     2.10360e-01,   1.23838e-08,   7.37570e-16,   4.43753e-23,
     2.71887e-01,   1.60059e-08,   9.53300e-16,   5.73545e-23,
     3.55694e-01,   2.09396e-08,   1.24715e-15,   7.50336e-23,
     4.84932e-01,   2.85479e-08,   1.70029e-15,   1.02296e-22,
     6.99767e-01,   4.11952e-08,   2.45355e-15,   1.47616e-22,
     1.05860e+00,   6.23193e-08,   3.71169e-15,   2.23311e-22,
     1.62902e+00,   9.59000e-08,   5.71172e-15,   3.43641e-22,
     2.54260e+00,   1.49682e-07,   8.91495e-15,   5.36361e-22,
     3.96488e+00,   2.33412e-07,   1.39018e-14,   8.36390e-22,
     6.10630e+00,   3.59476e-07,   2.14101e-14,   1.28812e-21,
     9.28766e+00,   5.46763e-07,   3.25647e-14,   1.95923e-21,
     1.39265e+01,   8.19849e-07,   4.88295e-14,   2.93779e-21,
     2.05382e+01,   1.20908e-06,   7.20119e-14,   4.33254e-21,
     3.00649e+01,   1.76992e-06,   1.05415e-13,   6.34219e-21,
     4.55102e+01,   2.67918e-06,   1.59569e-13,   9.60036e-21,
     7.47839e+01,   4.40251e-06,   2.62210e-13,   1.57756e-20,
     1.29621e+02,   7.63075e-06,   4.54481e-13,   2.73435e-20,
     2.14820e+02,   1.26464e-05,   7.53210e-13,   4.53162e-20,
     3.20002e+02,   1.88384e-05,   1.12200e-12,   6.75043e-20,
     4.29768e+02,   2.53004e-05,   1.50687e-12,   9.06596e-20,
     5.30827e+02,   3.12497e-05,   1.86120e-12,   1.11978e-19,
     5.99694e+02,   3.53039e-05,   2.10267e-12,   1.26505e-19,
     6.06537e+02,   3.57067e-05,   2.12666e-12,   1.27949e-19,
     5.43262e+02,   3.19817e-05,   1.90480e-12,   1.14601e-19,
     4.33545e+02,   2.55227e-05,   1.52011e-12,   9.14563e-20,
     3.13324e+02,   1.84453e-05,   1.09859e-12,   6.60956e-20,
     2.09006e+02,   1.23041e-05,   7.32824e-13,   4.40898e-20,
     1.31150e+02,   7.72078e-06,   4.59843e-13,   2.76661e-20,
     7.90683e+01,   4.65474e-06,   2.77232e-13,   1.66794e-20,
     4.73389e+01,   2.78683e-06,   1.65981e-13,   9.98614e-21,
     2.96409e+01,   1.74496e-06,   1.03928e-13,   6.25275e-21,
     2.04709e+01,   1.20512e-06,   7.17759e-14,   4.31834e-21  };

  double F80_kpMgSiO3[] = 
  {  2.19890e-02,   1.29231e-09,   7.68620e-17,   4.61833e-24,
     3.90612e-02,   2.29567e-09,   1.36538e-16,   8.20401e-24,
     6.05539e-02,   3.55881e-09,   2.11665e-16,   1.27181e-23,
     8.76116e-02,   5.14902e-09,   3.06244e-16,   1.84010e-23,
     1.43288e-01,   8.42115e-09,   5.00858e-16,   3.00946e-23,
     2.19266e-01,   1.28865e-08,   7.66440e-16,   4.60523e-23,
     3.36256e-01,   1.97621e-08,   1.17537e-15,   7.06235e-23,
     5.14336e-01,   3.02280e-08,   1.79785e-15,   1.08026e-22,
     7.97216e-01,   4.68532e-08,   2.78665e-15,   1.67439e-22,
     1.25414e+00,   7.37070e-08,   4.38381e-15,   2.63406e-22,
     2.03450e+00,   1.19569e-07,   7.11154e-15,   4.27305e-22,
     3.34648e+00,   1.96676e-07,   1.16975e-14,   7.02858e-22,
     5.45893e+00,   3.20827e-07,   1.90816e-14,   1.14654e-21,
     8.82117e+00,   5.18428e-07,   3.08342e-14,   1.85270e-21,
     1.41825e+01,   8.33519e-07,   4.95746e-14,   2.97874e-21,
     2.28419e+01,   1.34244e-06,   7.98433e-14,   4.79747e-21,
     3.71178e+01,   2.18145e-06,   1.29744e-13,   7.79582e-21,
     6.14272e+01,   3.61014e-06,   2.14717e-13,   1.29015e-20,
     1.03847e+02,   6.10317e-06,   3.62993e-13,   2.18108e-20,
     1.75507e+02,   1.03147e-05,   6.13479e-13,   3.68615e-20,
     2.82060e+02,   1.65769e-05,   9.85934e-13,   5.92409e-20,
     4.14519e+02,   2.43617e-05,   1.44894e-12,   8.70611e-20,
     5.59961e+02,   3.29094e-05,   1.95733e-12,   1.17608e-19,
     7.11024e+02,   4.17876e-05,   2.48537e-12,   1.49336e-19,
     8.40806e+02,   4.94150e-05,   2.93902e-12,   1.76594e-19,
     8.95312e+02,   5.26183e-05,   3.12954e-12,   1.88042e-19,
     8.40415e+02,   4.93920e-05,   2.93765e-12,   1.76512e-19,
     6.96693e+02,   4.09453e-05,   2.43527e-12,   1.46326e-19,
     5.18202e+02,   3.04552e-05,   1.81136e-12,   1.08838e-19,
     3.52864e+02,   2.07382e-05,   1.23343e-12,   7.41118e-20,
     2.24210e+02,   1.31770e-05,   7.83720e-13,   4.70907e-20,
     1.35138e+02,   7.94218e-06,   4.72371e-13,   2.83829e-20,
     7.83119e+01,   4.60247e-06,   2.73737e-13,   1.64478e-20,
     4.41553e+01,   2.59505e-06,   1.54344e-13,   9.27391e-21,
     2.45888e+01,   1.44511e-06,   8.59497e-14,   5.16438e-21  };

  double F80_kpFe3O4[] = 
  {  1.47700e-02,   5.95736e-10,   2.47198e-17,   1.05395e-24,
     2.47694e-02,   9.99055e-10,   4.14553e-17,   1.76749e-24,
     3.73580e-02,   1.50680e-09,   6.25240e-17,   2.66578e-24,
     5.32060e-02,   2.14602e-09,   8.90481e-17,   3.79666e-24,
     8.50036e-02,   3.42855e-09,   1.42266e-16,   6.06566e-24,
     1.29213e-01,   5.21170e-09,   2.16257e-16,   9.22034e-24,
     2.00170e-01,   8.07368e-09,   3.35013e-16,   1.42836e-23,
     3.15560e-01,   1.27279e-08,   5.28136e-16,   2.25176e-23,
     5.01384e-01,   2.02229e-08,   8.39139e-16,   3.57776e-23,
     7.88907e-01,   3.18199e-08,   1.32035e-15,   5.62946e-23,
     1.24250e+00,   5.01152e-08,   2.07951e-15,   8.86619e-23,
     1.95225e+00,   7.87422e-08,   3.26737e-15,   1.39308e-22,
     3.04002e+00,   1.22617e-07,   5.08792e-15,   2.16929e-22,
     4.68918e+00,   1.89134e-07,   7.84803e-15,   3.34609e-22,
     7.12599e+00,   2.87421e-07,   1.19264e-14,   5.08494e-22,
     1.05834e+01,   4.26872e-07,   1.77129e-14,   7.55206e-22,
     1.52356e+01,   6.14515e-07,   2.54990e-14,   1.08718e-21,
     2.13345e+01,   8.60511e-07,   3.57065e-14,   1.52238e-21,
     2.98061e+01,   1.20220e-06,   4.98848e-14,   2.12689e-21,
     4.27642e+01,   1.72486e-06,   7.15721e-14,   3.05155e-21,
     6.30370e+01,   2.54254e-06,   1.05502e-13,   4.49817e-21,
     9.29361e+01,   3.74850e-06,   1.55542e-13,   6.63170e-21,
     1.32987e+02,   5.36392e-06,   2.22573e-13,   9.48964e-21,
     1.82150e+02,   7.34689e-06,   3.04855e-13,   1.29978e-20,
     2.40388e+02,   9.69583e-06,   4.02324e-13,   1.71535e-20,
     3.12065e+02,   1.25869e-05,   5.22286e-13,   2.22682e-20,
     4.08414e+02,   1.64730e-05,   6.83540e-13,   2.91435e-20,
     5.49591e+02,   2.21673e-05,   9.19821e-13,   3.92175e-20,
     7.67451e+02,   3.09545e-05,   1.28444e-12,   5.47635e-20,
     1.10725e+03,   4.46602e-05,   1.85315e-12,   7.90111e-20,
     1.62060e+03,   6.53655e-05,   2.71231e-12,   1.15642e-19,
     2.33999e+03,   9.43815e-05,   3.91631e-12,   1.66976e-19,
     3.24367e+03,   1.30831e-04,   5.42876e-12,   2.31461e-19,
     4.25716e+03,   1.71709e-04,   7.12499e-12,   3.03781e-19,
     5.34010e+03,   2.15389e-04,   8.93744e-12,   3.81058e-19  };

  double F80_kpAC[] = 
  {  3.27960e-01,   1.38598e-06,   7.55735e-12,   4.81450e-17,
     4.38752e-01,   1.85420e-06,   1.01104e-11,   6.44096e-17,
     5.78230e-01,   2.44365e-06,   1.33245e-11,   8.48855e-17,
     7.53823e-01,   3.18572e-06,   1.73708e-11,   1.10663e-16,
     1.04016e+00,   4.39581e-06,   2.39692e-11,   1.52700e-16,
     1.41740e+00,   5.99012e-06,   3.26627e-11,   2.08084e-16,
     1.95298e+00,   8.25351e-06,   4.50044e-11,   2.86710e-16,
     2.71536e+00,   1.14754e-05,   6.25729e-11,   3.98636e-16,
     3.79688e+00,   1.60461e-05,   8.74963e-11,   5.57421e-16,
     5.29770e+00,   2.23889e-05,   1.22083e-10,   7.77774e-16,
     7.37879e+00,   3.11841e-05,   1.70043e-10,   1.08333e-15,
     1.02177e+01,   4.31826e-05,   2.35472e-10,   1.50021e-15,
     1.40437e+01,   5.93525e-05,   3.23651e-10,   2.06204e-15,
     1.92054e+01,   8.11692e-05,   4.42628e-10,   2.82017e-15,
     2.61678e+01,   1.10598e-04,   6.03124e-10,   3.84293e-15,
     3.55420e+01,   1.50223e-04,   8.19247e-10,   5.22035e-15,
     4.81826e+01,   2.03661e-04,   1.11074e-09,   7.07845e-15,
     6.53571e+01,   2.76274e-04,   1.50689e-09,   9.60430e-15,
     8.88432e+01,   3.75591e-04,   2.04885e-09,   1.30610e-14,
     1.20848e+02,   5.10965e-04,   2.78778e-09,   1.77762e-14,
     1.63896e+02,   6.93122e-04,   3.78255e-09,   2.41285e-14,
     2.21094e+02,   9.35280e-04,   5.10587e-09,   3.25877e-14,
     2.97091e+02,   1.25728e-03,   6.86710e-09,   4.38623e-14,
     3.98900e+02,   1.68904e-03,   9.23101e-09,   5.90159e-14,
     5.35045e+02,   2.26701e-03,   1.23979e-08,   7.93362e-14,
     7.13531e+02,   3.02568e-03,   1.65576e-08,   1.06036e-13,
     9.43152e+02,   4.00348e-03,   2.19225e-08,   1.40472e-13,
     1.23958e+03,   5.26932e-03,   2.88749e-08,   1.85093e-13,
     1.63174e+03,   6.95105e-03,   3.81274e-08,   2.44493e-13,
     2.16646e+03,   9.25803e-03,   5.08543e-08,   3.26274e-13,
     2.91544e+03,   1.25164e-02,   6.89019e-08,   4.42468e-13,
     3.98965e+03,   1.72412e-02,   9.52155e-08,   6.12399e-13,
     5.56191e+03,   2.42454e-02,   1.34481e-07,   8.66957e-13,
     7.88390e+03,   3.47095e-02,   1.93494e-07,   1.25089e-12,
     1.12585e+04,   4.99889e-02,   2.79860e-07,   1.81329e-12  };

  double F80_kpSiO2D[] = 
  {  7.60360e-02,   3.06759e-09,   1.27330e-16,   5.43132e-24,
     9.07207e-02,   3.66003e-09,   1.51922e-16,   6.48026e-24,
     1.09208e-01,   4.40586e-09,   1.82880e-16,   7.80079e-24,
     1.32481e-01,   5.34481e-09,   2.21854e-16,   9.46325e-24,
     1.58907e-01,   6.41094e-09,   2.66108e-16,   1.13509e-23,
     1.91565e-01,   7.72847e-09,   3.20796e-16,   1.36836e-23,
     2.30490e-01,   9.29887e-09,   3.85981e-16,   1.64641e-23,
     2.76795e-01,   1.11670e-08,   4.63523e-16,   1.97717e-23,
     3.33074e-01,   1.34375e-08,   5.57768e-16,   2.37917e-23,
     4.05325e-01,   1.63524e-08,   6.78759e-16,   2.89527e-23,
     5.08160e-01,   2.05012e-08,   8.50969e-16,   3.62983e-23,
     6.72472e-01,   2.71301e-08,   1.12613e-15,   4.80353e-23,
     9.48549e-01,   3.82682e-08,   1.58845e-15,   6.77557e-23,
     1.41787e+00,   5.72025e-08,   2.37438e-15,   1.01280e-22,
     2.19502e+00,   8.85555e-08,   3.67579e-15,   1.56792e-22,
     3.46719e+00,   1.39880e-07,   5.80619e-15,   2.47665e-22,
     5.76852e+00,   2.32725e-07,   9.66001e-15,   4.12050e-22,
     1.17194e+01,   4.72806e-07,   1.96254e-14,   8.37126e-22,
     3.16449e+01,   1.27668e-06,   5.29928e-14,   2.26042e-21,
     8.68296e+01,   3.50304e-06,   1.45406e-13,   6.20231e-21,
     1.92300e+02,   7.75813e-06,   3.22027e-13,   1.37362e-20,
     3.36231e+02,   1.35649e-05,   5.63055e-13,   2.40173e-20,
     5.05825e+02,   2.04070e-05,   8.47058e-13,   3.61315e-20,
     7.20624e+02,   2.90728e-05,   1.20676e-12,   5.14748e-20,
     9.77376e+02,   3.94312e-05,   1.63672e-12,   6.98148e-20,
     1.18646e+03,   4.78662e-05,   1.98685e-12,   8.47495e-20,
     1.23845e+03,   4.99638e-05,   2.07391e-12,   8.84633e-20,
     1.11188e+03,   4.48576e-05,   1.86197e-12,   7.94227e-20,
     8.76396e+02,   3.53572e-05,   1.46762e-12,   6.26017e-20,
     6.22207e+02,   2.51022e-05,   1.04195e-12,   4.44448e-20,
     4.07290e+02,   1.64317e-05,   6.82051e-13,   2.90931e-20,
     2.50573e+02,   1.01091e-05,   4.19611e-13,   1.78986e-20,
     1.47073e+02,   5.93352e-06,   2.46290e-13,   1.05056e-20,
     8.33122e+01,   3.36114e-06,   1.39515e-13,   5.95106e-21,
     4.59585e+01,   1.85415e-06,   7.69625e-14,   3.28286e-21  };

  double F80_kpAl2O3[] = 
  {  9.93250e-04,   4.00425e-11,   1.66043e-18,   7.07300e-26,
     1.81240e-03,   7.30662e-11,   3.02981e-18,   1.29062e-25,
     2.84365e-03,   1.14641e-10,   4.75377e-18,   2.02498e-25,
     4.14191e-03,   1.66980e-10,   6.92409e-18,   2.94948e-25,
     7.18271e-03,   2.89568e-10,   1.20074e-17,   5.11485e-25,
     1.13364e-02,   4.57021e-10,   1.89512e-17,   8.07269e-25,
     1.77361e-02,   7.15022e-10,   2.96496e-17,   1.26300e-24,
     2.59477e-02,   1.04607e-09,   4.33772e-17,   1.84775e-24,
     3.45425e-02,   1.39257e-09,   5.77452e-17,   2.45980e-24,
     4.22006e-02,   1.70130e-09,   7.05474e-17,   3.00513e-24,
     4.71420e-02,   1.90051e-09,   7.88080e-17,   3.35701e-24,
     4.91934e-02,   1.98321e-09,   8.22373e-17,   3.50309e-24,
     5.05162e-02,   2.03654e-09,   8.44486e-17,   3.59729e-24,
     5.78201e-02,   2.33100e-09,   9.66587e-17,   4.11741e-24,
     8.84237e-02,   3.56477e-09,   1.47819e-16,   6.29671e-24,
     1.78786e-01,   7.20770e-09,   2.98879e-16,   1.27315e-23,
     4.36404e-01,   1.75935e-08,   7.29543e-16,   3.10766e-23,
     1.63796e+00,   6.60337e-08,   2.73820e-15,   1.16640e-22,
     8.50817e+00,   3.43004e-07,   1.42232e-14,   6.05872e-22,
     3.92751e+01,   1.58336e-06,   6.56567e-14,   2.79680e-21,
     1.41436e+02,   5.70194e-06,   2.36441e-13,   1.00718e-20,
     3.83709e+02,   1.54691e-05,   6.41451e-13,   2.73241e-20,
     7.70411e+02,   3.10588e-05,   1.28791e-12,   5.48615e-20,
     1.16399e+03,   4.69258e-05,   1.94586e-12,   8.28884e-20,
     1.37566e+03,   5.54594e-05,   2.29972e-12,   9.79619e-20,
     1.33070e+03,   5.36466e-05,   2.22455e-12,   9.47599e-20,
     1.09978e+03,   4.43371e-05,   1.83852e-12,   7.83159e-20,
     8.05638e+02,   3.24790e-05,   1.34680e-12,   5.73700e-20,
     5.38690e+02,   2.17171e-05,   9.00536e-13,   3.83605e-20,
     3.36338e+02,   1.35593e-05,   5.62261e-13,   2.39508e-20,
     1.99460e+02,   8.04115e-06,   3.33440e-13,   1.42037e-20,
     1.13787e+02,   4.58729e-06,   1.90220e-13,   8.10285e-21,
     6.30411e+01,   2.54148e-06,   1.05387e-13,   4.48920e-21,
     3.41529e+01,   1.37686e-06,   5.70939e-14,   2.43205e-21,
     1.81893e+01,   7.33293e-07,   3.04073e-14,   1.29527e-21  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpFeM     [itab0] = F80_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = F80_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = F80_kpMgSiO3  [itab];
      my_rates->SN0_kpFe3O4   [itab0] = F80_kpFe3O4   [itab];
      my_rates->SN0_kpAC      [itab0] = F80_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = F80_kpSiO2D   [itab];
      my_rates->SN0_kpAl2O3   [itab0] = F80_kpAl2O3   [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_P170(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   5.29975e-02;
  my_rates->SN0_XO [iSN] =   5.60864e-01;
  my_rates->SN0_XMg[iSN] =   3.58367e-02;
  my_rates->SN0_XAl[iSN] =   3.27680e-04;
  my_rates->SN0_XSi[iSN] =   1.52750e-01;
  my_rates->SN0_XS [iSN] =   8.06035e-02;
  my_rates->SN0_XFe[iSN] =   5.29729e-02;

  my_rates->SN0_fC [iSN] =   5.29528e-02;
  my_rates->SN0_fO [iSN] =   5.60799e-01;
  my_rates->SN0_fMg[iSN] =   3.58366e-02;
  my_rates->SN0_fAl[iSN] =   3.27680e-04;
  my_rates->SN0_fSi[iSN] =   1.39585e-01;
  my_rates->SN0_fS [iSN] =   8.06035e-02;
  my_rates->SN0_fFe[iSN] =   5.29394e-02;

  my_rates->SN0_fSiM     [iSN] =   1.31079e-02;
  my_rates->SN0_fFeM     [iSN] =   3.34688e-05;
  my_rates->SN0_fMg2SiO4 [iSN] =   2.84952e-13;
  my_rates->SN0_fMgSiO3  [iSN] =   7.72302e-25;
  my_rates->SN0_fAC      [iSN] =   4.47758e-05;
  my_rates->SN0_fSiO2D   [iSN] =   1.23405e-04;
  my_rates->SN0_fMgO     [iSN] =   1.41247e-07;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   2.72050e-06;
  my_rates->SN0_r0FeM     [itab0 + 0] =   1.08069e-05;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   1.79010e-05;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   2.51189e-05;
  my_rates->SN0_r0AC      [itab0 + 0] =   8.32266e-07;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   2.12560e-05;
  my_rates->SN0_r0MgO     [itab0 + 0] =   1.60812e-05;

  my_rates->SN0_r0SiM     [itab0 + 1] =   2.87427e-11;
  my_rates->SN0_r0FeM     [itab0 + 1] =   1.19634e-10;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   3.24658e-10;
  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   6.30957e-10;
  my_rates->SN0_r0AC      [itab0 + 1] =   1.33383e-12;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   4.59721e-10;
  my_rates->SN0_r0MgO     [itab0 + 1] =   2.65603e-10;

  my_rates->SN0_r0SiM     [itab0 + 2] =   7.09270e-16;
  my_rates->SN0_r0FeM     [itab0 + 2] =   1.36724e-15;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   5.96244e-15;
  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   1.58489e-14;
  my_rates->SN0_r0AC      [itab0 + 2] =   4.37739e-18;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   1.01590e-14;
  my_rates->SN0_r0MgO     [itab0 + 2] =   4.50188e-15;

  NTd =            35;
 Nmom =             4;

  double P170_kpSiM[] = 
  {  1.54566e-01,   4.18815e-07,   4.38802e-12,   1.07634e-16,
     1.94597e-01,   5.27516e-07,   5.53200e-12,   1.35785e-16,
     2.44992e-01,   6.64362e-07,   6.97216e-12,   1.71223e-16,
     3.08437e-01,   8.36640e-07,   8.78521e-12,   2.15837e-16,
     3.88320e-01,   1.05373e-06,   1.10736e-11,   2.72217e-16,
     4.88890e-01,   1.32707e-06,   1.39559e-11,   3.43241e-16,
     6.15507e-01,   1.67134e-06,   1.75891e-11,   4.32825e-16,
     7.74906e-01,   2.10492e-06,   2.21686e-11,   5.45804e-16,
     9.75270e-01,   2.65013e-06,   2.79315e-11,   6.88059e-16,
     1.22485e+00,   3.32948e-06,   3.51169e-11,   8.65513e-16,
     1.52110e+00,   4.13615e-06,   4.36560e-11,   1.07652e-15,
     1.83679e+00,   4.99624e-06,   5.27698e-11,   1.30190e-15,
     2.15666e+00,   5.86838e-06,   6.20278e-11,   1.53115e-15,
     2.55518e+00,   6.95661e-06,   7.36137e-11,   1.81866e-15,
     3.22834e+00,   8.79676e-06,   9.32496e-11,   2.30667e-15,
     4.33225e+00,   1.18177e-05,   1.25559e-10,   3.11092e-15,
     5.81697e+00,   1.58892e-05,   1.69288e-10,   4.20258e-15,
     7.48671e+00,   2.04890e-05,   2.19146e-10,   5.45508e-15,
     9.22042e+00,   2.53165e-05,   2.72602e-10,   6.81745e-15,
     1.12094e+01,   3.09883e-05,   3.38326e-10,   8.54226e-15,
     1.40327e+01,   3.92793e-05,   4.39589e-10,   1.12871e-14,
     1.79556e+01,   5.11182e-05,   5.91090e-10,   1.55105e-14,
     2.19076e+01,   6.38759e-05,   7.72200e-10,   2.08556e-14,
     2.40396e+01,   7.32505e-05,   9.56570e-10,   2.71081e-14,
     2.35050e+01,   7.74606e-05,   1.14341e-09,   3.48207e-14,
     2.10698e+01,   7.75334e-05,   1.32500e-09,   4.37333e-14,
     1.81145e+01,   7.52350e-05,   1.46962e-09,   5.20115e-14,
     1.55283e+01,   7.16932e-05,   1.54535e-09,   5.74821e-14,
     1.35973e+01,   6.77513e-05,   1.55041e-09,   5.93660e-14,
     1.28568e+01,   6.72854e-05,   1.57294e-09,   6.05730e-14,
     1.62414e+01,   9.03391e-05,   2.02876e-09,   7.45675e-14,
     3.62078e+01,   2.25165e-04,   4.46960e-09,   1.44123e-13,
     1.13353e+02,   7.52538e-04,   1.29670e-08,   3.61005e-13,
     3.52139e+02,   2.31091e-03,   3.53205e-08,   8.74277e-13,
     9.61671e+02,   5.84143e-03,   8.01724e-08,   1.81153e-12  };

  double P170_kpFeM[] = 
  {  1.89038e-02,   2.14190e-07,   2.52372e-12,   3.12761e-17,
     3.29962e-02,   3.72180e-07,   4.35946e-12,   5.36210e-17,
     5.06585e-02,   5.69928e-07,   6.65321e-12,   8.14803e-17,
     7.28477e-02,   8.18206e-07,   9.53077e-12,   1.16395e-16,
     1.17388e-01,   1.31311e-06,   1.52141e-11,   1.84531e-16,
     1.77851e-01,   1.98324e-06,   2.28841e-11,   2.76089e-16,
     2.70459e-01,   3.00576e-06,   3.45287e-11,   4.14188e-16,
     4.11194e-01,   4.55413e-06,   5.20790e-11,   6.21046e-16,
     6.23355e-01,   6.88052e-06,   7.83314e-11,   9.28691e-16,
     9.33206e-01,   1.02680e-05,   1.16407e-10,   1.37259e-15,
     1.38665e+00,   1.52105e-05,   1.71743e-10,   2.01438e-15,
     2.02921e+00,   2.21948e-05,   2.49647e-10,   2.91351e-15,
     2.90114e+00,   3.16481e-05,   3.54739e-10,   4.12109e-15,
     4.04178e+00,   4.39869e-05,   4.91497e-10,   5.68634e-15,
     5.48102e+00,   5.95254e-05,   6.63275e-10,   7.64576e-15,
     7.23953e+00,   7.84801e-05,   8.72371e-10,   1.00241e-14,
     9.32332e+00,   1.00911e-04,   1.11939e-09,   1.28274e-14,
     1.17250e+01,   1.26738e-04,   1.40342e-09,   1.60451e-14,
     1.44493e+01,   1.56015e-04,   1.72511e-09,   1.96853e-14,
     1.75655e+01,   1.89497e-04,   2.09290e-09,   2.38454e-14,
     2.12689e+01,   2.29295e-04,   2.53020e-09,   2.87941e-14,
     2.59441e+01,   2.79574e-04,   3.08321e-09,   3.50603e-14,
     3.22597e+01,   3.47571e-04,   3.83222e-09,   4.35643e-14,
     4.13113e+01,   4.45148e-04,   4.90883e-09,   5.58148e-14,
     5.48617e+01,   5.91330e-04,   6.52336e-09,   7.42101e-14,
     7.57474e+01,   8.16590e-04,   9.01042e-09,   1.02534e-13,
     1.08611e+02,   1.17059e-03,   1.29124e-08,   1.46875e-13,
     1.61171e+02,   1.73570e-03,   1.91260e-08,   2.17252e-13,
     2.46323e+02,   2.64930e-03,   2.91428e-08,   3.30280e-13,
     3.85260e+02,   4.13661e-03,   4.54013e-08,   5.13015e-13,
     6.12075e+02,   6.55899e-03,   7.17997e-08,   8.08495e-13,
     9.79548e+02,   1.04742e-02,   1.14330e-07,   1.28253e-12,
     1.56260e+03,   1.66712e-02,   1.81431e-07,   2.02722e-12,
     2.44258e+03,   2.60024e-02,   2.82158e-07,   3.14054e-12,
     3.65611e+03,   3.88442e-02,   4.20404e-07,   4.66303e-12  };

  double P170_kpMg2SiO4[] = 
  {  1.05240e-01,   1.88391e-06,   3.41670e-11,   6.27487e-16,
     1.32588e-01,   2.37346e-06,   4.30456e-11,   7.90546e-16,
     1.67016e-01,   2.98977e-06,   5.42231e-11,   9.95825e-16,
     2.10360e-01,   3.76566e-06,   6.82948e-11,   1.25426e-15,
     2.71891e-01,   4.86713e-06,   8.82715e-11,   1.62113e-15,
     3.55703e-01,   6.36747e-06,   1.15482e-10,   2.12086e-15,
     4.84952e-01,   8.68116e-06,   1.57444e-10,   2.89150e-15,
     6.99815e-01,   1.25274e-05,   2.27201e-10,   4.17262e-15,
     1.05872e+00,   1.89523e-05,   3.43724e-10,   6.31261e-15,
     1.62931e+00,   2.91665e-05,   5.28972e-10,   9.71478e-15,
     2.54332e+00,   4.55285e-05,   8.25720e-10,   1.51647e-14,
     3.96674e+00,   7.10095e-05,   1.28786e-09,   2.36521e-14,
     6.11063e+00,   1.09388e-04,   1.98392e-09,   3.64359e-14,
     9.29771e+00,   1.66443e-04,   3.01871e-09,   5.54408e-14,
     1.39496e+01,   2.49721e-04,   4.52916e-09,   8.31823e-14,
     2.05907e+01,   3.68616e-04,   6.68568e-09,   1.22791e-13,
     3.01907e+01,   5.40495e-04,   9.80342e-09,   1.80058e-13,
     4.58362e+01,   8.20646e-04,   1.48857e-08,   2.73420e-13,
     7.56329e+01,   1.35424e-03,   2.45668e-08,   4.51281e-13,
     1.31563e+02,   2.35588e-03,   4.27405e-08,   7.85179e-13,
     2.18570e+02,   3.91412e-03,   7.10137e-08,   1.30465e-12,
     3.26293e+02,   5.84347e-03,   1.06022e-07,   1.94790e-12,
     4.39310e+02,   7.86786e-03,   1.42760e-07,   2.62299e-12,
     5.43909e+02,   9.74168e-03,   1.76768e-07,   3.24798e-12,
     6.15422e+02,   1.10228e-02,   2.00021e-07,   3.67534e-12,
     6.22866e+02,   1.11563e-02,   2.02446e-07,   3.71993e-12,
     5.58030e+02,   9.99504e-03,   1.81374e-07,   3.33274e-12,
     4.45434e+02,   7.97833e-03,   1.44778e-07,   2.66031e-12,
     3.22068e+02,   5.76874e-03,   1.04683e-07,   1.92357e-12,
     2.15071e+02,   3.85233e-03,   6.99082e-08,   1.28461e-12,
     1.35309e+02,   2.42376e-03,   4.39862e-08,   8.08311e-13,
     8.22952e+01,   1.47437e-03,   2.67608e-08,   4.91842e-13,
     5.06858e+01,   9.08497e-04,   1.64974e-08,   3.03339e-13,
     3.39204e+01,   6.08602e-04,   1.10623e-08,   2.03592e-13,
     2.60778e+01,   4.68473e-04,   8.52552e-09,   1.57083e-13  };

  double P170_kpMgSiO3[] = 
  {  2.19890e-02,   5.52339e-07,   1.38741e-11,   3.48502e-16,
     3.90618e-02,   9.81187e-07,   2.46463e-11,   6.19087e-16,
     6.05551e-02,   1.52108e-06,   3.82077e-11,   9.59734e-16,
     8.76136e-02,   2.20075e-06,   5.52805e-11,   1.38858e-15,
     1.43295e-01,   3.59940e-06,   9.04128e-11,   2.27107e-15,
     2.19281e-01,   5.50809e-06,   1.38357e-10,   3.47537e-15,
     3.36283e-01,   8.44705e-06,   2.12180e-10,   5.32973e-15,
     5.14405e-01,   1.29213e-05,   3.24567e-10,   8.15276e-15,
     7.97396e-01,   2.00297e-05,   5.03123e-10,   1.26379e-14,
     1.25461e+00,   3.15143e-05,   7.91604e-10,   1.98842e-14,
     2.03570e+00,   5.11345e-05,   1.28444e-09,   3.22637e-14,
     3.34964e+00,   8.41390e-05,   2.11348e-09,   5.30881e-14,
     5.46697e+00,   1.37324e-04,   3.44942e-09,   8.66456e-14,
     8.84149e+00,   2.22088e-04,   5.57860e-09,   1.40128e-13,
     1.42340e+01,   3.57543e-04,   8.98107e-09,   2.25594e-13,
     2.29730e+01,   5.77056e-04,   1.44950e-08,   3.64097e-13,
     3.74525e+01,   9.40765e-04,   2.36310e-08,   5.93583e-13,
     6.22721e+01,   1.56420e-03,   3.92910e-08,   9.86946e-13,
     1.05860e+02,   2.65909e-03,   6.67934e-08,   1.67777e-12,
     1.79797e+02,   4.51630e-03,   1.13444e-07,   2.84959e-12,
     2.89938e+02,   7.28291e-03,   1.82938e-07,   4.59521e-12,
     4.27133e+02,   1.07291e-02,   2.69503e-07,   6.76961e-12,
     5.77898e+02,   1.45162e-02,   3.64629e-07,   9.15907e-12,
     7.33244e+02,   1.84183e-02,   4.62646e-07,   1.16211e-11,
     8.63989e+02,   2.17024e-02,   5.45140e-07,   1.36933e-11,
     9.15488e+02,   2.29960e-02,   5.77634e-07,   1.45095e-11,
     8.55336e+02,   2.14851e-02,   5.39680e-07,   1.35562e-11,
     7.06445e+02,   1.77451e-02,   4.45737e-07,   1.11964e-11,
     5.24103e+02,   1.31649e-02,   3.30686e-07,   8.30647e-12,
     3.56328e+02,   8.95056e-03,   2.24828e-07,   5.64742e-12,
     2.26296e+02,   5.68431e-03,   1.42783e-07,   3.58656e-12,
     1.36494e+02,   3.42857e-03,   8.61218e-08,   2.16328e-12,
     7.93347e+01,   1.99280e-03,   5.00568e-08,   1.25737e-12,
     4.50672e+01,   1.13204e-03,   2.84355e-08,   7.14267e-13,
     2.56082e+01,   6.43249e-04,   1.61577e-08,   4.05863e-13  };

  double P170_kpAC[] = 
  {  3.27960e-01,   2.72950e-07,   4.37442e-13,   1.43560e-18,
     4.38752e-01,   3.65158e-07,   5.85222e-13,   1.92064e-18,
     5.78230e-01,   4.81242e-07,   7.71265e-13,   2.53127e-18,
     7.53823e-01,   6.27382e-07,   1.00548e-12,   3.30001e-18,
     1.04013e+00,   8.65671e-07,   1.38740e-12,   4.55369e-18,
     1.41735e+00,   1.17963e-06,   1.89059e-12,   6.20550e-18,
     1.95292e+00,   1.62537e-06,   2.60499e-12,   8.55074e-18,
     2.71531e+00,   2.25988e-06,   3.62196e-12,   1.18896e-17,
     3.79676e+00,   3.15995e-06,   5.06464e-12,   1.66270e-17,
     5.29743e+00,   4.40895e-06,   7.06667e-12,   2.32024e-17,
     7.37835e+00,   6.14088e-06,   9.84291e-12,   3.23229e-17,
     1.02168e+01,   8.50337e-06,   1.36303e-11,   4.47703e-17,
     1.40421e+01,   1.16873e-05,   1.87349e-11,   6.15552e-17,
     1.92022e+01,   1.59823e-05,   2.56223e-11,   8.42199e-17,
     2.61618e+01,   2.17754e-05,   3.49137e-11,   1.14828e-16,
     3.55310e+01,   2.95745e-05,   4.74267e-11,   1.56111e-16,
     4.81618e+01,   4.00896e-05,   6.43047e-11,   2.11918e-16,
     6.53182e+01,   5.43738e-05,   8.72469e-11,   2.88010e-16,
     8.87691e+01,   7.39018e-05,   1.18639e-10,   3.92576e-16,
     1.20709e+02,   1.00504e-04,   1.61458e-10,   5.36074e-16,
     1.63629e+02,   1.36264e-04,   2.19128e-10,   7.31160e-16,
     2.20590e+02,   1.83746e-04,   2.95924e-10,   9.94659e-16,
     2.96134e+02,   2.46760e-04,   3.98255e-10,   1.35250e-15,
     3.97082e+02,   3.31036e-04,   5.35697e-10,   1.84221e-15,
     5.31629e+02,   4.43468e-04,   7.19648e-10,   2.50526e-15,
     7.07172e+02,   5.90324e-04,   9.60341e-10,   3.37390e-15,
     9.31294e+02,   7.78129e-04,   1.26840e-09,   4.47614e-15,
     1.21722e+03,   1.01834e-03,   1.66280e-09,   5.86569e-15,
     1.58885e+03,   1.33177e-03,   2.17875e-09,   7.65284e-15,
     2.08299e+03,   1.75095e-03,   2.87268e-09,   1.00241e-14,
     2.75132e+03,   2.32265e-03,   3.82862e-09,   1.32680e-14,
     3.66587e+03,   3.11403e-03,   5.17218e-09,   1.78307e-14,
     4.93093e+03,   4.22440e-03,   7.09398e-09,   2.43990e-14,
     6.70346e+03,   5.80132e-03,   9.87242e-09,   3.39492e-14,
     9.22475e+03,   8.05642e-03,   1.38679e-08,   4.76003e-14  };

  double P170_kpSiO2D[] = 
  {  7.60055e-02,   1.61556e-06,   3.49407e-11,   7.72118e-16,
     9.06897e-02,   1.92769e-06,   4.16913e-11,   9.21295e-16,
     1.09176e-01,   2.32063e-06,   5.01898e-11,   1.10910e-15,
     1.32449e-01,   2.81532e-06,   6.08888e-11,   1.34553e-15,
     1.58880e-01,   3.37714e-06,   7.30398e-11,   1.61404e-15,
     1.91543e-01,   4.07143e-06,   8.80558e-11,   1.94587e-15,
     2.30474e-01,   4.89894e-06,   1.05953e-10,   2.34138e-15,
     2.76789e-01,   5.88341e-06,   1.27245e-10,   2.81190e-15,
     3.33087e-01,   7.08010e-06,   1.53127e-10,   3.38385e-15,
     4.05372e-01,   8.61660e-06,   1.86359e-10,   4.11822e-15,
     5.08264e-01,   1.08037e-05,   2.33661e-10,   5.16354e-15,
     6.72714e-01,   1.42993e-05,   3.09265e-10,   6.83430e-15,
     9.49132e-01,   2.01750e-05,   4.36348e-10,   9.64271e-15,
     1.41935e+00,   3.01704e-05,   6.52536e-10,   1.44204e-14,
     2.19866e+00,   4.67362e-05,   1.01084e-09,   2.23390e-14,
     3.47731e+00,   7.39180e-05,   1.59880e-09,   3.53339e-14,
     5.80454e+00,   1.23397e-04,   2.66922e-09,   5.89965e-14,
     1.18802e+01,   2.52595e-04,   5.46495e-09,   1.20817e-13,
     3.22814e+01,   6.86446e-04,   1.48537e-08,   3.28442e-13,
     8.86792e+01,   1.88575e-03,   4.08062e-08,   9.02327e-13,
     1.96255e+02,   4.17327e-03,   9.03045e-08,   1.99681e-12,
     3.42549e+02,   7.28388e-03,   1.57607e-07,   3.48479e-12,
     5.11736e+02,   1.08798e-02,   2.35370e-07,   5.20299e-12,
     7.16479e+02,   1.52272e-02,   3.29269e-07,   7.27459e-12,
     9.48899e+02,   2.01568e-02,   4.35598e-07,   9.61645e-12,
     1.12707e+03,   2.39309e-02,   5.16871e-07,   1.14029e-11,
     1.15762e+03,   2.45716e-02,   5.30490e-07,   1.16974e-11,
     1.02786e+03,   2.18125e-02,   4.70793e-07,   1.03775e-11,
     8.04358e+02,   1.70671e-02,   3.68303e-07,   8.11658e-12,
     5.68400e+02,   1.20593e-02,   2.60207e-07,   5.73357e-12,
     3.70952e+02,   7.86976e-03,   1.69795e-07,   3.74103e-12,
     2.27785e+02,   4.83229e-03,   1.04255e-07,   2.29687e-12,
     1.33544e+02,   2.83297e-03,   6.11184e-08,   1.34647e-12,
     7.56194e+01,   1.60416e-03,   3.46077e-08,   7.62416e-13,
     4.17402e+01,   8.85464e-04,   1.91029e-08,   4.20847e-13  };

  double P170_kpMgO[] = 
  {  2.25358e-04,   3.62400e-09,   5.98548e-14,   1.01451e-18,
     4.04933e-04,   6.51178e-09,   1.07551e-13,   1.82293e-18,
     6.31005e-04,   1.01473e-08,   1.67596e-13,   2.84068e-18,
     9.15612e-04,   1.47241e-08,   2.43188e-13,   4.12194e-18,
     1.52196e-03,   2.44750e-08,   4.04238e-13,   6.85169e-18,
     2.37413e-03,   3.81788e-08,   6.30577e-13,   1.06881e-17,
     3.77230e-03,   6.06633e-08,   1.00194e-12,   1.69826e-17,
     6.14410e-03,   9.88050e-08,   1.63191e-12,   2.76605e-17,
     1.01925e-02,   1.63910e-07,   2.70723e-12,   4.58870e-17,
     1.68946e-02,   2.71691e-07,   4.48743e-12,   7.60616e-17,
     2.96294e-02,   4.76492e-07,   7.87018e-12,   1.33401e-16,
     6.11307e-02,   9.83118e-07,   1.62386e-11,   2.75255e-16,
     1.43660e-01,   2.31047e-06,   3.81646e-11,   6.46939e-16,
     3.28203e-01,   5.27865e-06,   8.71964e-11,   1.47815e-15,
     6.41613e-01,   1.03199e-05,   1.70479e-10,   2.89009e-15,
     1.05901e+00,   1.70371e-05,   2.81505e-10,   4.77324e-15,
     1.60722e+00,   2.58913e-05,   4.28363e-10,   7.27252e-15,
     3.19739e+00,   5.16560e-05,   8.57040e-10,   1.45895e-14,
     1.40905e+01,   2.27322e-04,   3.76636e-09,   6.40310e-14,
     7.15093e+01,   1.15025e-03,   1.90029e-08,   3.22170e-13,
     2.48952e+02,   3.99764e-03,   6.59323e-08,   1.11600e-12,
     5.77048e+02,   9.25757e-03,   1.52545e-07,   2.57979e-12,
     9.46693e+02,   1.51796e-02,   2.49995e-07,   4.22569e-12,
     1.17600e+03,   1.88501e-02,   3.10343e-07,   5.24410e-12,
     1.17447e+03,   1.88215e-02,   3.09806e-07,   5.23395e-12,
     9.91612e+02,   1.58888e-02,   2.61497e-07,   4.41720e-12,
     7.37095e+02,   1.18095e-02,   1.94341e-07,   3.28252e-12,
     4.97965e+02,   7.97772e-03,   1.31276e-07,   2.21717e-12,
     3.13212e+02,   5.01764e-03,   8.25631e-08,   1.39438e-12,
     1.86742e+02,   2.99150e-03,   4.92223e-08,   8.31277e-13,
     1.06960e+02,   1.71340e-03,   2.81918e-08,   4.76098e-13,
     5.94287e+01,   9.51978e-04,   1.56633e-08,   2.64516e-13,
     3.22675e+01,   5.16881e-04,   8.50437e-09,   1.43617e-13,
     1.72142e+01,   2.75746e-04,   4.53688e-09,   7.66156e-14,
     9.06015e+00,   1.45130e-04,   2.38782e-09,   4.03237e-14  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = P170_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = P170_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = P170_kpMg2SiO4 [itab];
      my_rates->SN0_kpMgSiO3  [itab0] = P170_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = P170_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = P170_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = P170_kpMgO     [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_P200(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   3.65050e-02;
  my_rates->SN0_XO [iSN] =   4.88552e-01;
  my_rates->SN0_XMg[iSN] =   2.69665e-02;
  my_rates->SN0_XAl[iSN] =   1.36872e-04;
  my_rates->SN0_XSi[iSN] =   1.87324e-01;
  my_rates->SN0_XS [iSN] =   1.15582e-01;
  my_rates->SN0_XFe[iSN] =   6.79294e-02;

  my_rates->SN0_fC [iSN] =   3.64677e-02;
  my_rates->SN0_fO [iSN] =   4.88307e-01;
  my_rates->SN0_fMg[iSN] =   2.69665e-02;
  my_rates->SN0_fAl[iSN] =   1.36872e-04;
  my_rates->SN0_fSi[iSN] =   1.87051e-01;
  my_rates->SN0_fS [iSN] =   1.15582e-01;
  my_rates->SN0_fFe[iSN] =   6.75026e-02;

  my_rates->SN0_fSiM     [iSN] =   5.90622e-05;
  my_rates->SN0_fFeM     [iSN] =   4.26809e-04;
  my_rates->SN0_fMg2SiO4 [iSN] =   4.08246e-15;
  my_rates->SN0_fAC      [iSN] =   3.72287e-05;
  my_rates->SN0_fSiO2D   [iSN] =   4.59330e-04;
  my_rates->SN0_fMgO     [iSN] =   5.38389e-09;

  itab0 = 3 * iSN;
  my_rates->SN0_r0SiM     [itab0 + 0] =   8.86269e-07;
  my_rates->SN0_r0FeM     [itab0 + 0] =   2.02272e-06;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 0] =   1.42189e-05;
  my_rates->SN0_r0AC      [itab0 + 0] =   7.46096e-07;
  my_rates->SN0_r0SiO2D   [itab0 + 0] =   1.73471e-05;
  my_rates->SN0_r0MgO     [itab0 + 0] =   1.26307e-05;

  my_rates->SN0_r0SiM     [itab0 + 1] =   1.71166e-12;
  my_rates->SN0_r0FeM     [itab0 + 1] =   5.41308e-12;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 1] =   2.04834e-10;
  my_rates->SN0_r0AC      [itab0 + 1] =   9.32091e-13;
  my_rates->SN0_r0SiO2D   [itab0 + 1] =   3.08556e-10;
  my_rates->SN0_r0MgO     [itab0 + 1] =   1.59673e-10;

  my_rates->SN0_r0SiM     [itab0 + 2] =   5.46663e-18;
  my_rates->SN0_r0FeM     [itab0 + 2] =   2.06248e-17;
  my_rates->SN0_r0Mg2SiO4 [itab0 + 2] =   2.98805e-15;
  my_rates->SN0_r0AC      [itab0 + 2] =   1.99556e-18;
  my_rates->SN0_r0SiO2D   [itab0 + 2] =   5.66409e-15;
  my_rates->SN0_r0MgO     [itab0 + 2] =   2.02075e-15;

  NTd =            35;
 Nmom =             4;

  double P200_kpSiM[] = 
  {  1.54645e-01,   1.37048e-07,   2.64662e-13,   8.45209e-19,
     1.94685e-01,   1.72534e-07,   3.33193e-13,   1.06408e-18,
     2.45092e-01,   2.17207e-07,   4.19469e-13,   1.33961e-18,
     3.08551e-01,   2.73447e-07,   5.28084e-13,   1.68649e-18,
     3.88446e-01,   3.44254e-07,   6.64830e-13,   2.12322e-18,
     4.89028e-01,   4.33396e-07,   8.36986e-13,   2.67303e-18,
     6.15654e-01,   5.45619e-07,   1.05372e-12,   3.36522e-18,
     7.75056e-01,   6.86892e-07,   1.32656e-12,   4.23661e-18,
     9.75414e-01,   8.64464e-07,   1.66951e-12,   5.33191e-18,
     1.22498e+00,   1.08564e-06,   2.09668e-12,   6.69620e-18,
     1.52119e+00,   1.34817e-06,   2.60370e-12,   8.31550e-18,
     1.83682e+00,   1.62792e-06,   3.14401e-12,   1.00412e-17,
     2.15660e+00,   1.91133e-06,   3.69138e-12,   1.17894e-17,
     2.55494e+00,   2.26439e-06,   4.37331e-12,   1.39675e-17,
     3.22768e+00,   2.86067e-06,   5.52501e-12,   1.76460e-17,
     4.33077e+00,   3.83839e-06,   7.41349e-12,   2.36779e-17,
     5.81399e+00,   5.15308e-06,   9.95292e-12,   3.17892e-17,
     7.48107e+00,   6.63084e-06,   1.28076e-11,   4.09082e-17,
     9.20965e+00,   8.16336e-06,   1.57686e-11,   5.03682e-17,
     1.11868e+01,   9.91669e-06,   1.91573e-11,   6.11981e-17,
     1.39824e+01,   1.23969e-05,   2.39531e-11,   7.65315e-17,
     1.78531e+01,   1.58317e-05,   3.05972e-11,   9.77803e-17,
     2.17179e+01,   1.92624e-05,   3.72355e-11,   1.19020e-16,
     2.37018e+01,   2.10255e-05,   4.06526e-11,   1.29977e-16,
     2.29438e+01,   2.03578e-05,   3.93739e-11,   1.25940e-16,
     2.02577e+01,   1.79829e-05,   3.48015e-11,   1.11389e-16,
     1.70972e+01,   1.51927e-05,   2.94374e-11,   9.43240e-17,
     1.43876e+01,   1.28090e-05,   2.48736e-11,   7.98414e-17,
     1.23927e+01,   1.10758e-05,   2.16034e-11,   6.95714e-17,
     1.14807e+01,   1.03984e-05,   2.05890e-11,   6.70125e-17,
     1.34941e+01,   1.30878e-05,   2.78918e-11,   9.54730e-17,
     2.46530e+01,   2.87907e-05,   7.29927e-11,   2.79411e-16,
     6.22628e+01,   8.83029e-05,   2.60128e-10,   1.09021e-15,
     1.73608e+02,   2.73721e-04,   8.65552e-10,   3.78030e-15,
     4.81453e+02,   7.67244e-04,   2.43830e-09,   1.06714e-14  };

  double P200_kpFeM[] = 
  {  9.09085e-04,   3.41941e-09,   1.80920e-14,   1.27369e-19,
     1.61039e-03,   6.04845e-09,   3.18823e-14,   2.22622e-19,
     2.49320e-03,   9.35673e-09,   4.92187e-14,   3.42113e-19,
     3.60452e-03,   1.35207e-08,   7.10297e-14,   4.92299e-19,
     5.89215e-03,   2.20677e-08,   1.15535e-13,   7.95087e-19,
     9.02912e-03,   3.37737e-08,   1.76355e-13,   1.20713e-18,
     1.39130e-02,   5.19580e-08,   2.70494e-13,   1.84076e-18,
     2.14719e-02,   8.00206e-08,   4.15192e-13,   2.80850e-18,
     3.30912e-02,   1.23012e-07,   6.35936e-13,   4.27570e-18,
     5.03917e-02,   1.86799e-07,   9.62098e-13,   6.43127e-18,
     7.63291e-02,   2.81960e-07,   1.44615e-12,   9.61004e-18,
     1.14179e-01,   4.19929e-07,   2.14346e-12,   1.41576e-17,
     1.67303e-01,   6.12055e-07,   3.10728e-12,   2.03974e-17,
     2.39689e-01,   8.71248e-07,   4.39580e-12,   2.86715e-17,
     3.35470e-01,   1.21013e-06,   6.06256e-12,   3.92787e-17,
     4.59119e-01,   1.64148e-06,   8.15789e-12,   5.24812e-17,
     6.15407e-01,   2.17768e-06,   1.07246e-11,   6.84719e-17,
     8.10332e-01,   2.83301e-06,   1.38055e-11,   8.74075e-17,
     1.05438e+00,   3.63322e-06,   1.74835e-11,   1.09639e-16,
     1.37021e+00,   4.63773e-06,   2.19726e-11,   1.36212e-16,
     1.80496e+00,   5.97282e-06,   2.77433e-11,   1.69529e-16,
     2.45118e+00,   7.88356e-06,   3.57016e-11,   2.14206e-16,
     3.48432e+00,   1.08228e-05,   4.74799e-11,   2.78412e-16,
     5.23381e+00,   1.56143e-05,   6.59474e-11,   3.76121e-16,
     8.31585e+00,   2.37506e-05,   9.61230e-11,   5.31044e-16,
     1.38802e+01,   3.79407e-05,   1.46823e-10,   7.83600e-16,
     2.40529e+01,   6.30826e-05,   2.33568e-10,   1.20322e-15,
     4.26908e+01,   1.07918e-04,   3.83527e-10,   1.90929e-15,
     7.65687e+01,   1.87651e-04,   6.43418e-10,   3.10502e-15,
     1.37077e+02,   3.27739e-04,   1.09108e-09,   5.12746e-15,
     2.42517e+02,   5.69117e-04,   1.85185e-09,   8.51948e-15,
     4.21112e+02,   9.75067e-04,   3.11997e-09,   1.41236e-14,
     7.14738e+02,   1.63904e-03,   5.18038e-09,   2.31656e-14,
     1.18235e+03,   2.68818e-03,   8.40281e-09,   3.71602e-14,
     1.90342e+03,   4.27908e-03,   1.31829e-08,   5.74720e-14  };

  double P200_kpMg2SiO4[] = 
  {  1.05240e-01,   1.49640e-06,   2.15567e-11,   3.14463e-16,
     1.32588e-01,   1.88525e-06,   2.71584e-11,   3.96179e-16,
     1.67016e-01,   2.37479e-06,   3.42106e-11,   4.99054e-16,
     2.10360e-01,   2.99108e-06,   4.30887e-11,   6.28566e-16,
     2.71889e-01,   3.86597e-06,   5.56922e-11,   8.12421e-16,
     3.55700e-01,   5.05766e-06,   7.28594e-11,   1.06285e-15,
     4.84944e-01,   6.89538e-06,   9.93331e-11,   1.44904e-15,
     6.99796e-01,   9.95035e-06,   1.43342e-10,   2.09104e-15,
     1.05867e+00,   1.50532e-05,   2.16853e-10,   3.16340e-15,
     1.62919e+00,   2.31655e-05,   3.33717e-10,   4.86818e-15,
     2.54303e+00,   3.61594e-05,   5.20906e-10,   7.59887e-15,
     3.96603e+00,   5.63932e-05,   8.12393e-10,   1.18511e-14,
     6.10897e+00,   8.68643e-05,   1.25136e-09,   1.82548e-14,
     9.29384e+00,   1.32151e-04,   1.90378e-09,   2.77724e-14,
     1.39407e+01,   1.98229e-04,   2.85574e-09,   4.16600e-14,
     2.05705e+01,   2.92507e-04,   4.21400e-09,   6.14759e-14,
     3.01423e+01,   4.28630e-04,   6.17529e-09,   9.00909e-14,
     4.57104e+01,   6.50055e-04,   9.36595e-09,   1.36648e-13,
     7.53064e+01,   1.07104e-03,   1.54328e-08,   2.25181e-13,
     1.30815e+02,   1.86065e-03,   2.68126e-08,   3.91253e-13,
     2.17127e+02,   3.08849e-03,   4.45083e-08,   6.49502e-13,
     3.23880e+02,   4.60717e-03,   6.63972e-08,   9.68963e-13,
     4.35662e+02,   6.19759e-03,   8.93224e-08,   1.30358e-12,
     5.38932e+02,   7.66706e-03,   1.10506e-07,   1.61282e-12,
     6.09467e+02,   8.67078e-03,   1.24977e-07,   1.82406e-12,
     6.16707e+02,   8.77390e-03,   1.26464e-07,   1.84580e-12,
     5.52479e+02,   7.86016e-03,   1.13294e-07,   1.65358e-12,
     4.40975e+02,   6.27379e-03,   9.04293e-08,   1.31986e-12,
     3.18791e+02,   4.53552e-03,   6.53747e-08,   9.54183e-13,
     2.12803e+02,   3.02766e-03,   4.36414e-08,   6.36986e-13,
     1.33761e+02,   1.90318e-03,   2.74344e-08,   4.00448e-13,
     8.11305e+01,   1.15453e-03,   1.66451e-08,   2.42996e-13,
     4.95598e+01,   7.05594e-04,   1.01773e-08,   1.48640e-13,
     3.25906e+01,   4.64465e-04,   6.70585e-09,   9.80294e-14,
     2.45153e+01,   3.49814e-04,   5.05664e-09,   7.40047e-14  };

  double P200_kpAC[] = 
  {  3.27960e-01,   2.44690e-07,   3.05689e-13,   6.54462e-19,
     4.38752e-01,   3.27351e-07,   4.08957e-13,   8.75569e-19,
     5.78230e-01,   4.31415e-07,   5.38965e-13,   1.15393e-18,
     7.53823e-01,   5.62424e-07,   7.02635e-13,   1.50436e-18,
     1.04013e+00,   7.76039e-07,   9.69511e-13,   2.07581e-18,
     1.41735e+00,   1.05749e-06,   1.32113e-12,   2.82873e-18,
     1.95292e+00,   1.45707e-06,   1.82034e-12,   3.89771e-18,
     2.71531e+00,   2.02589e-06,   2.53098e-12,   5.41951e-18,
     3.79675e+00,   2.83276e-06,   3.53905e-12,   7.57851e-18,
     5.29742e+00,   3.95242e-06,   4.93793e-12,   1.05749e-17,
     7.37833e+00,   5.50501e-06,   6.87773e-12,   1.47304e-17,
     1.02168e+01,   7.62282e-06,   9.52385e-12,   2.04006e-17,
     1.40421e+01,   1.04770e-05,   1.30901e-11,   2.80446e-17,
     1.92021e+01,   1.43270e-05,   1.79012e-11,   3.83618e-17,
     2.61616e+01,   1.95199e-05,   2.43907e-11,   5.22874e-17,
     3.55307e+01,   2.65108e-05,   3.31285e-11,   7.10543e-17,
     4.81611e+01,   3.59357e-05,   4.49107e-11,   9.63939e-17,
     6.53170e+01,   4.87382e-05,   6.09196e-11,   1.30888e-16,
     8.87668e+01,   6.62390e-05,   8.28118e-11,   1.78181e-16,
     1.20705e+02,   9.00772e-05,   1.12647e-10,   2.42873e-16,
     1.63621e+02,   1.22115e-04,   1.52778e-10,   3.30388e-16,
     2.20575e+02,   1.64642e-04,   2.06114e-10,   4.47715e-16,
     2.96104e+02,   2.21060e-04,   2.76991e-10,   6.05473e-16,
     3.97026e+02,   2.96479e-04,   3.71917e-10,   8.19263e-16,
     5.31526e+02,   3.97049e-04,   4.98696e-10,   1.10695e-15,
     7.06983e+02,   5.28345e-04,   6.64405e-10,   1.48343e-15,
     9.30949e+02,   6.96139e-04,   8.76417e-10,   1.96294e-15,
     1.21659e+03,   9.10536e-04,   1.14776e-09,   2.57160e-15,
     1.58765e+03,   1.18982e-03,   1.50225e-09,   3.36001e-15,
     2.08069e+03,   1.56242e-03,   1.97748e-09,   4.41103e-15,
     2.74685e+03,   2.06876e-03,   2.62823e-09,   5.84923e-15,
     3.65712e+03,   2.76611e-03,   3.53431e-09,   7.86205e-15,
     4.91402e+03,   3.73838e-03,   4.81486e-09,   1.07337e-14,
     6.67206e+03,   5.11081e-03,   6.64539e-09,   1.48739e-14,
     9.17096e+03,   7.06881e-03,   9.26820e-09,   2.07940e-14  };

  double P200_kpSiO2D[] = 
  {  7.60136e-02,   1.31860e-06,   2.34539e-11,   4.30529e-16,
     9.06981e-02,   1.57333e-06,   2.79848e-11,   5.13702e-16,
     1.09185e-01,   1.89402e-06,   3.36890e-11,   6.18411e-16,
     1.32458e-01,   2.29774e-06,   4.08700e-11,   7.50232e-16,
     1.58888e-01,   2.75623e-06,   4.90252e-11,   8.99935e-16,
     1.91549e-01,   3.32281e-06,   5.91032e-11,   1.08494e-15,
     2.30478e-01,   3.99813e-06,   7.11152e-11,   1.30544e-15,
     2.76790e-01,   4.80151e-06,   8.54053e-11,   1.56776e-15,
     3.33085e-01,   5.77806e-06,   1.02776e-10,   1.88663e-15,
     4.05362e-01,   7.03188e-06,   1.25078e-10,   2.29603e-15,
     5.08238e-01,   8.81651e-06,   1.56822e-10,   2.87877e-15,
     6.72650e-01,   1.16686e-05,   2.07555e-10,   3.81010e-15,
     9.48976e-01,   1.64623e-05,   2.92824e-10,   5.37545e-15,
     1.41897e+00,   2.46158e-05,   4.37862e-10,   8.03808e-15,
     2.19768e+00,   3.81252e-05,   6.78181e-10,   1.24501e-14,
     3.47459e+00,   6.02792e-05,   1.07231e-09,   1.96868e-14,
     5.79486e+00,   1.00542e-04,   1.78877e-09,   3.28453e-14,
     1.18370e+01,   2.05418e-04,   3.65563e-09,   6.71473e-14,
     3.21108e+01,   5.57348e-04,   9.92083e-09,   1.82279e-13,
     8.81851e+01,   1.53068e-03,   2.72473e-08,   5.00649e-13,
     1.95201e+02,   3.38815e-03,   6.03099e-08,   1.10811e-12,
     3.40883e+02,   5.91645e-03,   1.05307e-07,   1.93471e-12,
     5.10268e+02,   8.85442e-03,   1.57557e-07,   2.89365e-12,
     7.17917e+02,   1.24511e-02,   2.21410e-07,   4.06297e-12,
     9.57120e+02,   1.65879e-02,   2.94709e-07,   5.40198e-12,
     1.14369e+03,   1.98088e-02,   3.51649e-07,   6.43917e-12,
     1.17995e+03,   2.04271e-02,   3.62410e-07,   6.63127e-12,
     1.05086e+03,   1.81866e-02,   3.22531e-07,   5.89860e-12,
     8.23986e+02,   1.42573e-02,   2.52780e-07,   4.62145e-12,
     5.83014e+02,   1.00864e-02,   1.78801e-07,   3.26825e-12,
     3.80801e+02,   6.58750e-03,   1.16763e-07,   2.13399e-12,
     2.33954e+02,   4.04697e-03,   7.17275e-08,   1.31079e-12,
     1.37203e+02,   2.37328e-03,   4.20617e-08,   7.68622e-13,
     7.77017e+01,   1.34403e-03,   2.38199e-08,   4.35268e-13,
     4.28856e+01,   7.41814e-04,   1.31471e-08,   2.40245e-13  };

  double P200_kpMgO[] = 
  {  2.25374e-04,   2.84663e-09,   3.59861e-14,   4.55424e-19,
     4.04950e-04,   5.11481e-09,   6.46597e-14,   8.18303e-19,
     6.31024e-04,   7.97028e-09,   1.00758e-13,   1.27514e-18,
     9.15633e-04,   1.15651e-08,   1.46202e-13,   1.85026e-18,
     1.52197e-03,   1.92235e-08,   2.43017e-13,   3.07552e-18,
     2.37410e-03,   2.99866e-08,   3.79080e-13,   4.79746e-18,
     3.77219e-03,   4.76455e-08,   6.02318e-13,   7.62266e-18,
     6.14378e-03,   7.76004e-08,   9.80998e-13,   1.24151e-17,
     1.01916e-02,   1.28728e-07,   1.62733e-12,   2.05948e-17,
     1.68921e-02,   2.13360e-07,   2.69722e-12,   3.41348e-17,
     2.96213e-02,   3.74139e-07,   4.72974e-12,   5.98574e-17,
     6.10970e-02,   7.71701e-07,   9.75560e-12,   1.23463e-16,
     1.43533e-01,   1.81293e-06,   2.29185e-11,   2.90047e-16,
     3.27806e-01,   4.14044e-06,   5.23423e-11,   6.62424e-16,
     6.40565e-01,   8.09084e-06,   1.02282e-10,   1.29445e-15,
     1.05534e+00,   1.33298e-05,   1.68513e-10,   2.13267e-15,
     1.58348e+00,   2.00015e-05,   2.52868e-10,   3.20044e-15,
     3.07310e+00,   3.88206e-05,   4.90841e-10,   6.21320e-15,
     1.37089e+01,   1.73169e-04,   2.18940e-09,   2.77122e-14,
     7.13521e+01,   9.01235e-04,   1.13932e-08,   1.44189e-13,
     2.52035e+02,   3.18325e-03,   4.02394e-08,   5.09217e-13,
     5.88759e+02,   7.43595e-03,   9.39947e-08,   1.18943e-12,
     9.70270e+02,   1.22542e-02,   1.54897e-07,   1.96005e-12,
     1.20866e+03,   1.52648e-02,   1.92951e-07,   2.44153e-12,
     1.20927e+03,   1.52725e-02,   1.93046e-07,   2.44272e-12,
     1.02223e+03,   1.29102e-02,   1.63185e-07,   2.06485e-12,
     7.60460e+02,   9.60416e-03,   1.21396e-07,   1.53608e-12,
     5.14029e+02,   6.49187e-03,   8.20570e-08,   1.03830e-12,
     3.23438e+02,   4.08481e-03,   5.16318e-08,   6.53316e-13,
     1.92889e+02,   2.43607e-03,   3.07917e-08,   3.89618e-13,
     1.10502e+02,   1.39557e-03,   1.76399e-08,   2.23204e-13,
     6.14046e+01,   7.75500e-04,   9.80225e-09,   1.24031e-13,
     3.33434e+01,   4.21106e-04,   5.32274e-09,   6.73503e-14,
     1.77894e+01,   2.24668e-04,   2.83978e-09,   3.59327e-14,
     9.36317e+00,   1.18251e-04,   1.49468e-09,   1.89126e-14  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpSiM     [itab0] = P200_kpSiM     [itab];
      my_rates->SN0_kpFeM     [itab0] = P200_kpFeM     [itab];
      my_rates->SN0_kpMg2SiO4 [itab0] = P200_kpMg2SiO4 [itab];
      my_rates->SN0_kpAC      [itab0] = P200_kpAC      [itab];
      my_rates->SN0_kpSiO2D   [itab0] = P200_kpSiO2D   [itab];
      my_rates->SN0_kpMgO     [itab0] = P200_kpMgO     [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}

int calc_rates_dust_Y19(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{

  int NTd, Nmom;
  int iTd, imom, itab0, itab;

  my_rates->SN0_XC [iSN] =   2.50000e-01;
  my_rates->SN0_XO [iSN] =   2.93867e-01;
  my_rates->SN0_XMg[iSN] =   6.00000e-02;
  my_rates->SN0_XAl[iSN] =   2.85361e-03;
  my_rates->SN0_XSi[iSN] =   7.00000e-02;
  my_rates->SN0_XS [iSN] =   1.58191e-02;
  my_rates->SN0_XFe[iSN] =   6.64078e-02;

  my_rates->SN0_fC [iSN] =   0.00000e+00;
  my_rates->SN0_fO [iSN] =   1.73867e-01;
  my_rates->SN0_fMg[iSN] =   0.00000e+00;
  my_rates->SN0_fAl[iSN] =   2.85361e-03;
  my_rates->SN0_fSi[iSN] =   0.00000e+00;
  my_rates->SN0_fS [iSN] =   1.58191e-02;
  my_rates->SN0_fFe[iSN] =   6.64078e-02;

  my_rates->SN0_fMgSiO3  [iSN] =   2.50000e-01;
  my_rates->SN0_fAC      [iSN] =   2.50000e-01;

  itab0 = 3 * iSN;
  my_rates->SN0_r0MgSiO3  [itab0 + 0] =   1.00000e-05;
  my_rates->SN0_r0AC      [itab0 + 0] =   1.00000e-05;

  my_rates->SN0_r0MgSiO3  [itab0 + 1] =   1.00000e-10;
  my_rates->SN0_r0AC      [itab0 + 1] =   1.00000e-10;

  my_rates->SN0_r0MgSiO3  [itab0 + 2] =   1.00000e-15;
  my_rates->SN0_r0AC      [itab0 + 2] =   1.00000e-15;

  NTd =            35;
 Nmom =             4;

  double Y19_kpMgSiO3[] = 
  {  2.19890e-02,   2.19890e-07,   2.19890e-12,   2.19890e-17,
     3.90612e-02,   3.90612e-07,   3.90612e-12,   3.90612e-17,
     6.05539e-02,   6.05539e-07,   6.05539e-12,   6.05539e-17,
     8.76116e-02,   8.76116e-07,   8.76116e-12,   8.76116e-17,
     1.43288e-01,   1.43288e-06,   1.43288e-11,   1.43288e-16,
     2.19266e-01,   2.19266e-06,   2.19266e-11,   2.19266e-16,
     3.36256e-01,   3.36256e-06,   3.36256e-11,   3.36256e-16,
     5.14336e-01,   5.14336e-06,   5.14336e-11,   5.14336e-16,
     7.97216e-01,   7.97216e-06,   7.97216e-11,   7.97216e-16,
     1.25414e+00,   1.25414e-05,   1.25414e-10,   1.25414e-15,
     2.03450e+00,   2.03450e-05,   2.03450e-10,   2.03450e-15,
     3.34654e+00,   3.34654e-05,   3.34654e-10,   3.34654e-15,
     5.45913e+00,   5.45913e-05,   5.45913e-10,   5.45913e-15,
     8.82166e+00,   8.82166e-05,   8.82166e-10,   8.82166e-15,
     1.41836e+01,   1.41836e-04,   1.41836e-09,   1.41836e-14,
     2.28449e+01,   2.28449e-04,   2.28449e-09,   2.28449e-14,
     3.71258e+01,   3.71258e-04,   3.71258e-09,   3.71258e-14,
     6.14485e+01,   6.14485e-04,   6.14485e-09,   6.14485e-14,
     1.03898e+02,   1.03898e-03,   1.03898e-08,   1.03898e-13,
     1.75627e+02,   1.75627e-03,   1.75627e-08,   1.75627e-13,
     2.82290e+02,   2.82290e-03,   2.82290e-08,   2.82290e-13,
     4.14908e+02,   4.14908e-03,   4.14908e-08,   4.14908e-13,
     5.60606e+02,   5.60606e-03,   5.60606e-08,   5.60606e-13,
     7.12020e+02,   7.12020e-03,   7.12020e-08,   7.12020e-13,
     8.42130e+02,   8.42130e-03,   8.42130e-08,   8.42130e-13,
     8.96812e+02,   8.96812e-03,   8.96812e-08,   8.96812e-13,
     8.41845e+02,   8.41845e-03,   8.41845e-08,   8.41845e-13,
     6.97883e+02,   6.97883e-03,   6.97883e-08,   6.97883e-13,
     5.19082e+02,   5.19082e-03,   5.19082e-08,   5.19082e-13,
     3.53464e+02,   3.53464e-03,   3.53464e-08,   3.53464e-13,
     2.24610e+02,   2.24610e-03,   2.24610e-08,   2.24610e-13,
     1.35389e+02,   1.35389e-03,   1.35389e-08,   1.35389e-13,
     7.84898e+01,   7.84898e-04,   7.84898e-09,   7.84898e-14,
     4.43113e+01,   4.43113e-04,   4.43113e-09,   4.43113e-14,
     2.49396e+01,   2.49396e-04,   2.49396e-09,   2.49396e-14  };

  double Y19_kpAC[] = 
  {  6.76020e-02,   6.76020e-07,   6.76020e-12,   6.76020e-17,
     1.20181e-01,   1.20181e-06,   1.20181e-11,   1.20181e-16,
     1.86375e-01,   1.86375e-06,   1.86375e-11,   1.86375e-16,
     2.69708e-01,   2.69708e-06,   2.69708e-11,   2.69708e-16,
     4.44368e-01,   4.44368e-06,   4.44368e-11,   4.44368e-16,
     6.87406e-01,   6.87406e-06,   6.87406e-11,   6.87406e-16,
     1.07797e+00,   1.07797e-05,   1.07797e-10,   1.07797e-15,
     1.71241e+00,   1.71241e-05,   1.71241e-10,   1.71241e-15,
     2.74163e+00,   2.74163e-05,   2.74163e-10,   2.74163e-15,
     4.35812e+00,   4.35812e-05,   4.35812e-10,   4.35812e-15,
     6.98720e+00,   6.98720e-05,   6.98720e-10,   6.98720e-15,
     1.13206e+01,   1.13206e-04,   1.13206e-09,   1.13206e-14,
     1.85159e+01,   1.85159e-04,   1.85159e-09,   1.85159e-14,
     3.09414e+01,   3.09414e-04,   3.09414e-09,   3.09414e-14,
     5.34575e+01,   5.34575e-04,   5.34575e-09,   5.34575e-14,
     9.60912e+01,   9.60912e-04,   9.60912e-09,   9.60912e-14,
     1.76000e+02,   1.76000e-03,   1.76000e-08,   1.76000e-13,
     3.10598e+02,   3.10598e-03,   3.10598e-08,   3.10598e-13,
     4.95502e+02,   4.95502e-03,   4.95502e-08,   4.95502e-13,
     6.87389e+02,   6.87389e-03,   6.87389e-08,   6.87389e-13,
     8.21760e+02,   8.21760e-03,   8.21760e-08,   8.21760e-13,
     8.55209e+02,   8.55209e-03,   8.55209e-08,   8.55209e-13,
     7.90315e+02,   7.90315e-03,   7.90315e-08,   7.90315e-13,
     6.63814e+02,   6.63814e-03,   6.63814e-08,   6.63814e-13,
     5.19410e+02,   5.19410e-03,   5.19410e-08,   5.19410e-13,
     3.88956e+02,   3.88956e-03,   3.88956e-08,   3.88956e-13,
     2.88141e+02,   2.88141e-03,   2.88141e-08,   2.88141e-13,
     2.20698e+02,   2.20698e-03,   2.20698e-08,   2.20698e-13,
     1.84716e+02,   1.84716e-03,   1.84716e-08,   1.84716e-13,
     1.78316e+02,   1.78316e-03,   1.78316e-08,   1.78316e-13,
     2.05010e+02,   2.05010e-03,   2.05010e-08,   2.05010e-13,
     2.82760e+02,   2.82760e-03,   2.82760e-08,   2.82760e-13,
     4.70437e+02,   4.70437e-03,   4.70437e-08,   4.70437e-13,
     9.49808e+02,   9.49808e-03,   9.49808e-08,   9.49808e-13,
     2.18634e+03,   2.18634e-02,   2.18634e-07,   2.18634e-12  };


  itab0 = Nmom * NTd * iSN;
  itab  = 0;
  for(imom = 0; imom < Nmom; imom++) {
    for(iTd = 0; iTd < NTd; iTd++) {
      my_rates->SN0_kpMgSiO3  [itab0] = Y19_kpMgSiO3  [itab];
      my_rates->SN0_kpAC      [itab0] = Y19_kpAC      [itab];
      itab0++;
      itab ++;
    }
  }

  return SUCCESS;
}
