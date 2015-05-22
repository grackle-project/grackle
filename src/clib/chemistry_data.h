/***********************************************************************
/
/ Chemistry data structure
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __CHEMISTRY_DATA_H__
#define __CHEMISTRY_DATA_H__

typedef struct
{

  // Rank of dataset.
  long long grid_rank;

  // Dimension of dataset.
  long long *grid_dimension;

  // Dataset parameter values.
  double **grid_parameters;

  // Heating values
  double *heating_data;

  // Cooling values
  double *cooling_data;

  // Mean Molecular Weight values
  double *mmw_data;

  // Length of 1D flattened data
  long long data_size;

} cloudy_data;

typedef struct 
{

  // adiabatic index
  double Gamma;

  // HDF5 file containing Cloudy cooling table and UV background
  char *grackle_data_file;


  /****************************************
   *** chemistry and cooling parameters ***
   ****************************************/

  int use_grackle;  // 0) grackle off, 1) grackle on
  int with_radiative_cooling; // include cooling in chemistry solver
                                 // 0) no, 1) yes
  int primordial_chemistry; // 1) HI, HII, HeI, HeII, HeIII, e
                               // 2) + H2, H2I+, H-
                               // 3) + D, D+, HD
  int metal_cooling;        // 0) off, 1) on using Cloudy tables
  int h2_on_dust;  // H2 formation on dust
                      // not well tested, should be left off for now

  // Use a CMB temperature floor.
  int cmb_temperature_floor;

  /* additional H2 chemistry parameters
     best left unchanged. */
  
  int three_body_rate;
  int cie_cooling;
  int h2_optical_depth_approximation;

  /* photo-electric heating from irradiated dust */

  int photoelectric_heating;
  double photoelectric_heating_rate; // in CGS

  /***************************************
   *** radiation background parameters ***
   ***************************************/

  int UVbackground;

  struct UVBtable {
    long long Nz;

    double zmin, zmax;    
    double *z;

    double *k24;
    double *k25;
    double *k26;
    double *k27;
    double *k28;
    double *k29;
    double *k30;
    double *k31;

    double *piHI;
    double *piHeII;
    double *piHeI;
  } UVbackground_table;

  double UVbackground_redshift_on;
  double UVbackground_redshift_off;
  double UVbackground_redshift_fullon;
  double UVbackground_redshift_drop;

  int Compton_xray_heating;

  double LWbackground_intensity;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  int LWbackground_sawtooth_suppression;

  /**************************************
   *** primordial chemistry rate data ***
   **************************************/

  int NumberOfTemperatureBins;   
  int CaseBRecombination;
  double TemperatureStart;        // range of temperature in K
  double TemperatureEnd;

  /* 6 species rates */

  double *k1;
  double *k2;
  double *k3;
  double *k4;
  double *k5;
  double *k6;

  /* 9 species rates (including H2) */

  double *k7;
  double *k8;
  double *k9;
  double *k10;
  double *k11;
  double *k12;
  double *k13;
  double *k14;
  double *k15;
  double *k16;
  double *k17;
  double *k18;
  double *k19;
  double *k20;  /* currently not used */
  double *k21;  /* currently not used */
  double *k22;  /* 3-body H2 formation */
  double *k23;  /* H2-H2 dissociation */

  double *k13dd;  /* density dependent version of k13 (collisional H2
                    dissociation); actually 7 functions instead of 1. */

  /* Radiative rates for 6-species (for external field). */

  double k24;
  double k25;
  double k26;

  /* Radiative rates for 9-species (for external field). */

  double k27;
  double k28;
  double k29;
  double k30;
  double k31;

  /* 12 species rates (with Deuterium). */

  double *k50;
  double *k51;
  double *k52;
  double *k53;
  double *k54;
  double *k55;
  double *k56;

  /* H2 formation on dust. */

  int NumberOfDustTemperatureBins;   
  double DustTemperatureStart;        // range of temperature in K
  double DustTemperatureEnd;
  double *h2dust;                     // function of Tgas and Tdust

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */

  double *n_cr_n;
  double *n_cr_d1;
  double *n_cr_d2;

  /********************
   *** cooling data ***
   ********************/

  int ih2co;                     // flag for H2 cooling (0-off/1-on)
  int ipiht;                     // flag for photoionization cooling

  double HydrogenFractionByMass;
  double DeuteriumToHydrogenRatio;
  double SolarMetalFractionByMass;

  /* 6 species rates */

  double *ceHI;                   // collisional excitation rates
  double *ceHeI;
  double *ceHeII;
  double *ciHI;                   // collisional ionization
  double *ciHeI;
  double *ciHeIS;
  double *ciHeII;
  double *reHII;                  // recombination
  double *reHeII1;
  double *reHeII2;
  double *reHeIII;
  double *brem;                   // free-free (Bremsstrahlung)
  double comp;                    // Compton cooling
  double comp_xray;               // X-ray compton heating coefficient
  double temp_xray;               // X-ray compton heating temperature (K)
  double gammah;                  // Photoelectric heating (code units)

  /* radiative rates (external field). */

  double piHI;                    // photo-ionization cooling
  double piHeI;                   //    (no temperature dependance)
  double piHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi_g.F. */

  double *hyd01k;
  double *h2k01;
  double *vibh;
  double *roth;
  double *rotl;

  double *GP99LowDensityLimit;
  double *GP99HighDensityLimit;

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  double *GAHI;
  double *GAH2;
  double *GAHe;
  double *GAHp;
  double *GAel;

  /* 12 species rates (including HD) */

  double *HDlte;
  double *HDlow;

  /* CIE cooling */
  double *cieco;

  /* Gas/grain energy transfer. */
  double *gas_grain;

  // Primordial cooling data

  cloudy_data cloudy_primordial;

  // Metal cooling data

  cloudy_data cloudy_metal;

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  double cloudy_electron_fraction_factor;

} chemistry_data;

#endif
