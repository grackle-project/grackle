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

/**********************************
 *** Grackle runtime parameters ***
 **********************************/
typedef struct
{

  /* grackle on/off flag
     0) off, 1) on */
  int use_grackle;

  /* include cooling in chemistry solver
     0) no, 1) yes */
  int with_radiative_cooling;

  /* chemistry network
     0) tabulated cooling (no chemistry)
     1) HI, HII, HeI, HeII, HeIII, e
     2) + H2, H2I+, H-
     3) + D, D+, HD */
  int primordial_chemistry;

  /* metal cooling on/off
     0) off, 1) on */
  int metal_cooling;

  /* add heating from UV background model
     0) off, 1) on */
  int UVbackground;

  /* data file containing cooling and UV background tables */
  char *grackle_data_file;

  /* Use a CMB temperature floor
     0) no, 1) yes */
  int cmb_temperature_floor;

  /* adiabatic index */
  double Gamma;

  /* H2 formation on dust grains and dust cooling
     0) off, 1) on */
  int h2_on_dust;

  /* photo-electric heating from irradiated dust */

  int photoelectric_heating;
  double photoelectric_heating_rate; // in CGS

  /* flags to signal that arrays of volumetric or
     specific heating rates are being provided */

  int use_volumetric_heating_rate;
  int use_specific_heating_rate;

  /* additional chemistry solver parameters */
  int three_body_rate;
  int cie_cooling;
  int h2_optical_depth_approximation;
  int ih2co; // flag for H2 cooling (0-off/1-on)
  int ipiht; // flag for photoionization cooling
  double HydrogenFractionByMass;
  double DeuteriumToHydrogenRatio;
  double SolarMetalFractionByMass;
  int NumberOfTemperatureBins;
  int CaseBRecombination;
  double TemperatureStart;
  double TemperatureEnd;
  int NumberOfDustTemperatureBins;
  double DustTemperatureStart;
  double DustTemperatureEnd;

  /* additional radiation background parameters */
  int Compton_xray_heating;
  int LWbackground_sawtooth_suppression;
  double LWbackground_intensity;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  double UVbackground_redshift_on;
  double UVbackground_redshift_off;
  double UVbackground_redshift_fullon;
  double UVbackground_redshift_drop;

  /* Factor to account for extra electrons from metals.
     Only for old-style Cloudy tables.
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3. */
  double cloudy_electron_fraction_factor;

  /* flags and parameters to signal that RT
     is being used, and appropriate parameters
     for setting RT solvers */
  int use_radiative_transfer;
  int radiative_transfer_coupled_rate_solver;
  int radiative_transfer_intermediate_step;
  int radiative_transfer_hydrogen_only;

  /* flag for approximiate self-shielding as well
     as spectrum averaged photo heating and
     photo ionization shielding factors */
  int self_shielding_method;

  /* flag for Wolcott-Green+ 2011 H2 self-shielding */
  int H2_self_shielding;

  /* number of OpenMP threads, if supported */
# ifdef _OPENMP
  int omp_nthreads;
# endif

} chemistry_data;

/*****************************
 *** Cooling table storage ***
 *****************************/
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

/**********************************
 *** UVbackground table storage ***
 **********************************/

typedef struct{

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

    // spectrum averaged absorption cross sections
    double *crsHI;
    double *crsHeII;
    double *crsHeI;

} UVBtable;

/******************************************
 *** Chemistry and cooling data storage ***
 ******************************************/
typedef struct
{

  /**********************************
   * primordial chemistry rate data *
   **********************************/

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

  /* New H-ionizing reactions, used for 6, 9 & 12 species chemistry */
  double *k57;
  double *k58;

  /* H2 formation on dust grains */
  double *h2dust;

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */
  double *n_cr_n;
  double *n_cr_d1;
  double *n_cr_d2;

  /********************************
   * primordial cooling rate data *
   ********************************/

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

  /* Updated H2 LTE rate from Glover (2015, MNRAS, 451, 2082) */
  double *H2LTE;

  /* 12 species rates (including HD) */
  double *HDlte;
  double *HDlow;

  /* CIE cooling */
  double *cieco;

  /* Gas/grain energy transfer. */
  double *gas_grain;

  /* UV background data */
  UVBtable UVbackground_table;

  /* Primordial cooling table */
  cloudy_data cloudy_primordial;

  /* Metal cooling table */
  cloudy_data cloudy_metal;

  /* New/old cloudy data flag */
  int cloudy_data_new;

  /* gray averaged cross sections for UVB self shielding */
  double hi_avg_crs;
  double heii_avg_crs;
  double hei_avg_crs;

} chemistry_data_storage;

#endif
