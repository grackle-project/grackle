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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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

  /* dust chemistry and cooling
     0) no dust chemistry or cooling
     1) basic dust treatment
        - photo-electric heating
        - recombination cooling
        - with primordial_chemistry > 1:
          - H2 formation on dust grains
          - gas/grain heat transfer */
  int dust_chemistry;

  /* metal cooling on/off
     0) off, 1) on */
  int metal_cooling;

  /* add heating from UV background model
     0) off, 1) on */
  int UVbackground;

  /* data file containing cooling and UV background tables */
  const char *grackle_data_file;

  /* Use a CMB temperature floor
     0) no, 1) yes */
  int cmb_temperature_floor;

  /* adiabatic index */
  double Gamma;

  /* H2 formation on dust grains and dust cooling
     0) off, 1) on */
  int h2_on_dust;

  /* Flag to supply a dust density field */
  int use_dust_density_field;

  /* Flag for dust recombination cooling */
  int dust_recombination_cooling;

  /* Flag to solve metal chemistry */
  int metal_chemistry;

  /* Flag to solve grain growth reactions */
  int grain_growth;

  /* Flag to enable tracking of multiple metal sources for dust evolution */
  int multi_metals;

  /* if metal_chemistry>0 and multi_metals=0,
     this flag selects a single metal source for dust evolution from the
     following options:
     0. metal/dust abundances of local ISM (Pollack et al. 1994)
     1-4. Pop III normal core-collapse supernovae (Nozawa et al. 2007)
     with progenitor masses 13, 20, 25 and 30 Msun
     5-8. Pop III faint supernovae (Marassi et al. 2014)
     with progenitor masses 13, 50 and 80 Msun
     9-10. Pop III pair-instability supernovae (Nozawa et al. 2007)
     with progenitor masses 170 and 200 Msun
     11. simple dust model (only include silicate and graphite; Yajima et al. 2019)
  */
  int metal_abundances;

  /* Flag to solve multiple grain species
     1. enstatite + amorphous carbon (also follow Mg metal density)
     2. + metallic silicon + metallic iron + forsterite + magnetite
        + silica + magnesia + troilite + alumina
        (also follow Al, S, Fe metal densities)
     3. + water ice + volatile organics + refractory organics
  */
  int dust_species;

  /* Flag to solve temperatures of multiple grain species */
  int use_multiple_dust_temperatures;

  /* Flag to supply dust sublimation */
  int dust_sublimation;

  /* photo-electric heating from irradiated dust */
  int photoelectric_heating;
  double photoelectric_heating_rate;

  /* Flag to supply a field for the interstellar radiation field */
  int use_isrf_field;

  // local FUV interstellar radiation field in Habing units
  double interstellar_radiation_field;

  /* flags to signal that arrays of volumetric or
     specific heating rates are being provided */

  int use_volumetric_heating_rate;
  int use_specific_heating_rate;

  /* flag to use a field specifying a temperature floor */

  int use_temperature_floor;

  /* use to specify scalar temperature floor */

  double temperature_floor_scalar;

  /* H2 cooling rate
   * 0: Lepp & Shull (1983)
   * 1: Galli & Palla (1998)
   * 2: Glover & Abel (2008)
   * 3: Chiaki & Wise (2019)
   */
  int h2_cooling_rate;

  /* HD cooling rate
   * 0: Coppola et al (2011) and Wrathmall, Gusdorf, & Flower (2007)
   * 1: Chiaki & Wise (2019)
   */
  int hd_cooling_rate;

  /* H2 formation from 3-body reactions
   * 0: Abel, Bryan & Norman (2002)
   * 1: Palla, Salpeter & Stahler (1983)
   * 2: Cohen & Westberg (1983)
   * 3: Flower & Harris (2007)
   * 4: Glover (2008)
   * 5: Forrey (2013)
   */
  int three_body_rate;

  /* H2 collisionally-induced emission cooling
   * 0: off
   * 1: Ripamonti & Abel (2003)
   * 2: Yoshida et al. (2006)
   */
  int cie_cooling;

  /* H2 cooling attenuation from Ripamonti & Abel (2004) */
  int h2_optical_depth_approximation;

  int ih2co; // flag for H2 cooling (0-off/1-on)
  int ipiht; // flag for photoionization cooling
  double HydrogenFractionByMass;
  double DeuteriumToHydrogenRatio;
  double SolarMetalFractionByMass;
  double local_dust_to_gas_ratio;
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
  int radiative_transfer_HDI_dissociation;
  int radiative_transfer_metal_ionization;
  int radiative_transfer_metal_dissociation;

  /* flag to signal H2 self-shielding is being done in hydro code */
  int radiative_transfer_use_H2_shielding;

  /* flag for approximiate self-shielding as well
     as spectrum averaged photo heating and
     photo ionization shielding factors */
  int self_shielding_method;

  /* flag for Wolcott-Green+ 2011 H2 self-shielding. Can be set to either 0,1,2 or 3.
     These determine the length scale used in the calculation of the H2 column density.
     Please refer to the grackle documentation for specifics. */
  int H2_self_shielding;

  /* flag for custom H2-shielding factor. The factor is provided as an additional field 
     by the user and is multiplied to the rate for radiative H2 dissocitation */
  int H2_custom_shielding;

  /* flag to select which formula for calculating k11 you want to use. 
     Setting to 1 will use Savin 2004, 2 will use Abel et al. 1996  */
  int h2_charge_exchange_rate;

  /* flag to select which formula for calculating h2dust you want to use. 
     Setting to 1 will use Omukai 2000, 2 will use Hollenbach & McKee (1979) */ 
  int h2_dust_rate;

  /* flag to select which formula for calculating low density H2 cooling rate 
     due to H collisions. Setting to 1 will use Lique 2015, 2 will use Glover and Abel 2008 */
  int h2_h_cooling_rate;
  
  /* flags specific to calc_rates_g (1 is on, 0 is off) -- will be set to default values 
     if unspecified  */
  int collisional_excitation_rates; //Collisional excitation
  int collisional_ionisation_rates; //Collisional ionisation
  int recombination_cooling_rates; //Recombination cooling
  int bremsstrahlung_cooling_rates; //Bremsstrahlung cooling

  /* flag to add primordial continuum opacity */
  int use_primordial_continuum_opacity;

  /* Alternative rates for HD-related reactions (k50-k56)
   * 0: multiple sources (see rate_functions.c for details)
   * 1: Stancil, Lepp & Dalgarno (1998)
   */
  int hd_reaction_rates;

  /* Alternative gas-grain heat transfer rate.
   * 0: Hollenbach & McKee (1989)
   * 1: Omukai (2000) - see rate_functions.c for more details
   */
  int gas_grain_cooling_rate;

  /* flags for alternative rate calculations */
  int use_uniform_grain_dist_gamma_isrf; //Alternative calculation scheme for gamma_isrf

  /* maximum number of subcycle iterations for solve_chemistry */
  int max_iterations;

  /* flag to exit if max iterations exceeded */
  int exit_after_iterations_exceeded;

  /* number of OpenMP threads, if supported */
  int omp_nthreads;

} chemistry_data;

/*****************************
 ***** code units struct *****
 *****************************/

typedef struct
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;

} code_units;

/*****************************
 *** Cooling table storage ***
 *****************************/

// this is an implementation detail that may change in the future. Always rely
// upon the grid_rank value stored in cloudy_data
#define GRACKLE_CLOUDY_TABLE_MAX_DIMENSION 5

typedef struct
{

  // Rank of dataset.
  long long grid_rank;

  // Dimension of dataset.
  long long grid_dimension[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  // Dataset parameter values.
  double *grid_parameters[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

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

typedef struct
{

    long long Nz;

    double zmin, zmax;    
    double *z;

    /* Radiative rates for 6-species. */
    double *k24;
    double *k25;
    double *k26;
    /* Radiative rates for 9-species. */
    double *k27;
    double *k28;
    double *k29;
    double *k30;
    double *k31;

    double *piHI;
    double *piHeI;
    double *piHeII;

    // spectrum averaged absorption cross sections
    double *crsHI;
    double *crsHeI;
    double *crsHeII;

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

  /* 15 species rates (with DM, HDII, HeHII) */
  double *k125;
  double *k129;
  double *k130;
  double *k131;
  double *k132;
  double *k133;
  double *k134;
  double *k135;
  double *k136;
  double *k137;
  double *k148;
  double *k149;
  double *k150;
  double *k151;
  double *k152;
  double *k153;

  /* Metal species */
  double *kz15;
  double *kz16;
  double *kz17;
  double *kz18;
  double *kz19;
  double *kz20;
  double *kz21;
  double *kz22;
  double *kz23;
  double *kz24;
  double *kz25;
  double *kz26;
  double *kz27;
  double *kz28;
  double *kz29;
  double *kz30;
  double *kz31;
  double *kz32;
  double *kz33;
  double *kz34;
  double *kz35;
  double *kz36;
  double *kz37;
  double *kz38;
  double *kz39;
  double *kz40;
  double *kz41;
  double *kz42;
  double *kz43;
  double *kz44;
  double *kz45;
  double *kz46;
  double *kz47;
  double *kz48;
  double *kz49;
  double *kz50;
  double *kz51;
  double *kz52;
  double *kz53;
  double *kz54;

  /* H2 formation on dust grains */
  double *h2dust;
  double *h2dustS;
  double *h2dustC;

  /* Grain growth rate */
  double *grain_growth_rate;

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

  /* radiative rates (external field). */
  double piHI;                    // photo-ionization cooling
  double piHeI;                   //    (no temperature dependance)
  double piHeII;

  // spectrum averaged absorption cross sections
  double crsHI;
  double crsHeI;
  double crsHeII;

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

  /*******************************
   * dust chemistry/cooling data *
   *******************************/

  // Photo-electric heating (code units)
  double gammah;

  // Electron recombination onto dust grains
  double *regr;

  // Heating of dust by interstellar radiation field
  double gamma_isrf;
  // for arbitrary grain size distribution
  double gamma_isrf2;

  /* Gas/grain energy transfer. */
  double *gas_grain;
  // for arbitrary grain size distribution
  double *gas_grain2;

  /* CIE cooling rate (Yoshida et al. 2006) */
  double *cieY06;

  /* H2 and HD cooling rates (collision with HI; Hollenbach & McKee 1979) */
  int    *LH2_N, LH2_Size;
  double *LH2_D, *LH2_T, *LH2_H, LH2_dD, LH2_dT, LH2_dH, *LH2_L;
  int    *LHD_N, LHD_Size;
  double *LHD_D, *LHD_T, *LHD_H, LHD_dD, LHD_dT, LHD_dH, *LHD_L;

  /* Fine-structure cooling rates (collision with HI; Maio et al. 2007) */
  int    *LCI_N, LCI_Size;
  double *LCI_D, *LCI_T, *LCI_H, LCI_dD, LCI_dT, LCI_dH, *LCI_L;
  int    *LCII_N, LCII_Size;
  double *LCII_D, *LCII_T, *LCII_H, LCII_dD, LCII_dT, LCII_dH, *LCII_L;
  int    *LOI_N, LOI_Size;
  double *LOI_D, *LOI_T, *LOI_H, LOI_dD, LOI_dT, LOI_dH, *LOI_L;

  /* metal molecular cooling rates (collision with H2I; UMIST table) */
  int    *LCO_N, LCO_Size;
  double *LCO_D, *LCO_T, *LCO_H, LCO_dD, LCO_dT, LCO_dH, *LCO_L;
  int    *LOH_N, LOH_Size;
  double *LOH_D, *LOH_T, *LOH_H, LOH_dD, LOH_dT, LOH_dH, *LOH_L;
  int    *LH2O_N, LH2O_Size;
  double *LH2O_D, *LH2O_T, *LH2O_H, LH2O_dD, LH2O_dT, LH2O_dH, *LH2O_L;

  /* primordial opacity */
  int    *alphap_N, alphap_Size;
  double *alphap_D, *alphap_T, alphap_dD, alphap_dT;
  double *alphap_Data;

  /* metal/dust abundance */
  int    *gr_N, gr_Size;
  double gr_dT, *gr_Td;
  int     SN0_N;
  double *SN0_XC , *SN0_XO , *SN0_XMg, *SN0_XAl, *SN0_XSi, *SN0_XS , *SN0_XFe;
  double *SN0_fC , *SN0_fO , *SN0_fMg, *SN0_fAl, *SN0_fSi, *SN0_fS , *SN0_fFe;
  double *SN0_fSiM, *SN0_fFeM, *SN0_fMg2SiO4, *SN0_fMgSiO3, *SN0_fFe3O4
       , *SN0_fAC, *SN0_fSiO2D, *SN0_fMgO, *SN0_fFeS, *SN0_fAl2O3
       , *SN0_freforg , *SN0_fvolorg , *SN0_fH2Oice;
  double *SN0_r0SiM, *SN0_r0FeM, *SN0_r0Mg2SiO4, *SN0_r0MgSiO3, *SN0_r0Fe3O4
       , *SN0_r0AC, *SN0_r0SiO2D, *SN0_r0MgO, *SN0_r0FeS, *SN0_r0Al2O3
       , *SN0_r0reforg , *SN0_r0volorg , *SN0_r0H2Oice;
  double *SN0_kpSiM, *SN0_kpFeM, *SN0_kpMg2SiO4, *SN0_kpMgSiO3, *SN0_kpFe3O4
       , *SN0_kpAC, *SN0_kpSiO2D, *SN0_kpMgO, *SN0_kpFeS, *SN0_kpAl2O3
       , *SN0_kpreforg , *SN0_kpvolorg , *SN0_kpH2Oice;

  /* UV background data */
  UVBtable UVbackground_table;

  /* Primordial cooling table */
  cloudy_data cloudy_primordial;

  /* Metal cooling table */
  cloudy_data cloudy_metal;

  /* New/old cloudy data flag */
  int cloudy_data_new;

  /* tracks the initial value of the code-units */
  code_units initial_units;
} chemistry_data_storage;

/**************************
 *** photo-rate storage ***
 **************************/

typedef struct
{

    /* Radiative rates for 6-species. */
    double k24;
    double k25;
    double k26;
    /* Radiative rates for 9-species. */
    double k27;
    double k28;
    double k29;
    double k30;
    double k31;

    double piHI;
    double piHeI;
    double piHeII;

    // spectrum averaged absorption cross sections
    double crsHI;
    double crsHeI;
    double crsHeII;

    // X-ray compton heating coefficient
    double comp_xray;
    // X-ray compton heating temperature (K)
    double temp_xray;

} photo_rate_storage;

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __CHEMISTRY_DATA_H__ */
