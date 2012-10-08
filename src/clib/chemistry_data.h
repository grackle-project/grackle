#define CLOUDY_COOLING_MAX_DIMENSION 5

struct chemistry_data
{

  // adiabatic index
  float Gamma;

  /****************************************
   *** chemistry and cooling parameters ***
   ****************************************/

  int use_chemistry;
  int primordial_chemistry; // 1) HI, HII, HeI, HeII, HeIII, e
                            // 2) + H2, H2I+, H-
                            // 3) + D, D+, HD
  int metal_cooling;        // 0) off, 1) on using Cloudy tables
  int h2_on_dust;  // should be left off for now

  // Use a CMB temperature floor.
  int cmb_temperature_floor;

  // Flag to control whether or not to include heating from Cloudy.
  int include_metal_heating;

  // Cooling grid file.
  char *cloudy_table_file;

  /* additional H2 chemistry parameters
     best left unchanged. */
  
  int three_body_rate;
  int cie_cooling;
  int h2_optical_depth_approximation;

  /* photo-electric heating from irradiated dust */

  int photoelectric_heating;
  float photoelectric_heating_rate; // in CGS

  /***************************************
   *** radiation background parameters ***
   ***************************************/

  int RadiationFieldType;
  int AdjustUVBackground; 
  int AdjustUVBackgroundHighRedshift; 
  float SetUVBAmplitude;
  float SetHeIIHeatingScale;
  int RadiationFieldLevelRecompute;
  int RadiationXRaySecondaryIon;
  int RadiationXRayComptonHeating;
  int TabulatedLWBackground;
  float RadiationFieldRedshift;

  FLOAT TimeFieldLastUpdated;
  int  RadiationShield;

  /* Frequency dependant information (constant). */

  int NumberOfFrequencyBins;
  float FrequencyBinWidth;    // in log10(eV)

  float *HICrossSection;
  float *HeICrossSection;
  float *HeIICrossSection;

  float *StellarSpectrum;
  float *QuasarSpectrum;

  float *Spectrum[4];
  float *Emissivity[4];

  /**************************************
   *** primordial chemistry rate data ***
   **************************************/

  int NumberOfTemperatureBins;   
  int CaseBRecombination;
  float TemperatureStart;        // range of temperature in K
  float TemperatureEnd;

  /* 6 species rates */

  float *k1;
  float *k2;
  float *k3;
  float *k4;
  float *k5;
  float *k6;

  /* 9 species rates (including H2) */

  float *k7;
  float *k8;
  float *k9;
  float *k10;
  float *k11;
  float *k12;
  float *k13;
  float *k14;
  float *k15;
  float *k16;
  float *k17;
  float *k18;
  float *k19;
  float *k20;  /* currently not used */
  float *k21;  /* currently not used */
  float *k22;  /* 3-body H2 formation */
  float *k23;  /* H2-H2 dissociation */

  float *k13dd;  /* density dependent version of k13 (collisional H2
                    dissociation); actually 7 functions instead of 1. */

  /* Radiative rates for 6-species (for external field). */

  float k24;
  float k25;
  float k26;

  /* Radiative rates for 9-species (for external field). */

  float k27;
  float k28;
  float k29;
  float k30;
  float k31;

  /* 12 species rates (with Deuterium). */

  float *k50;
  float *k51;
  float *k52;
  float *k53;
  float *k54;
  float *k55;
  float *k56;

  /* H2 formation on dust. */

  int NumberOfDustTemperatureBins;   
  float DustTemperatureStart;        // range of temperature in K
  float DustTemperatureEnd;
  float *h2dust;                     // function of Tgas and Tdust

  /* Chemical heating from H2 formation. */
  /* numerator and denominator of Eq 23 of Omukai ea. 2000. */

  float *n_cr_n;
  float *n_cr_d1;
  float *n_cr_d2;

  /********************
   *** cooling data ***
   ********************/

  int ih2co;                     // flag for H2 cooling (0-off/1-on)
  int ipiht;                     // flag for photoionization cooling

  /* Radiation parameters (should be put in own spot). */

  float alpha0;
  float f3;
  float f0to3;
  float RadiationRedshiftOn;
  float RadiationRedshiftOff;
  float RadiationRedshiftFullOn;
  float RadiationRedshiftDropOff;
  float HydrogenFractionByMass;
  float DeuteriumToHydrogenRatio;
  float SolarMetalFractionByMass;

  /* 6 species rates */

  float *ceHI;                   // collisional excitation rates
  float *ceHeI;
  float *ceHeII;
  float *ciHI;                   // collisional ionization
  float *ciHeI;
  float *ciHeIS;
  float *ciHeII;
  float *reHII;                  // recombination
  float *reHeII1;
  float *reHeII2;
  float *reHeIII;
  float *brem;                   // free-free (Bremsstrahlung)
  float comp;                    // Compton cooling
  float comp_xray;               // X-ray compton heating coefficient
  float temp_xray;               // X-ray compton heating temperature (K)
  float gammah;                  // Photoelectric heating (code units)

  /* radiative rates (external field). */

  float piHI;                    // photo-ionization cooling
  float piHeI;                   //    (no temperature dependance)
  float piHeII;

  /* 9 species rates (including H2) 
       The first five are for the Lepp & Shull rates.
       The next two are for the (better) Galli & Palla 1999 rates. 
       The selection is controlled by a flag in cool1d_multi.src. */

  float *hyd01k;
  float *h2k01;
  float *vibh;
  float *roth;
  float *rotl;

  float *GP99LowDensityLimit;
  float *GP99HighDensityLimit;

  /* Revised H2 cooling rates from Glover & Abel 2008 */
  float *GAHI;
  float *GAH2;
  float *GAHe;
  float *GAHp;
  float *GAel;

  /* 12 species rates (including HD) */

  float *HDlte;
  float *HDlow;
  float *HDcool;

  /* CIE cooling */
  float *cieco;

  /* Metal line cooling (tabulated by temperature and electron
     fraction) */

  int NumberOfElectronFracBins;
  float ElectronFracStart;
  float ElectronFracEnd;
  float *metals;

  /* Gas/grain energy transfer. */
  float *gas_grain;

  /* For analysis, ratios of metal fine structure line emission is
     desired sometimes. */

  float *metal_ratios;
  char *MetalCoolingTable;
  int MR_NumberOfTemperatureBins;
  int MR_NumberOfElectronFracBins;
  float MR_TemperatureStart;
  float MR_TemperatureEnd;
  float MR_ElectronFracStart;
  float MR_ElectronFracEnd;

  /***************************
   *** cloudy cooling data ***
   ***************************/

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  float CloudyElectronFractionFactor;

  // Rank of Cloudy dataset.
  int CloudyCoolingGridRank;

  // Dimension of Cloudy dataset.
  int *CloudyCoolingGridDimension;

  // Dataset parameter values.
  //  float *CloudyCoolingGridParameters[CLOUDY_COOLING_MAX_DIMENSION];
  float **CloudyCoolingGridParameters;

  // Heating values
  float *CloudyHeating;

  // Cooling values
  float *CloudyCooling;

  // Length of 1D flattened Cloudy data
  int CloudyDataSize;
};
