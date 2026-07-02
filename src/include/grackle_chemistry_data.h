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

// this should go before the header-guard
#if !defined(GRIMPL_PUBLIC_INCLUDE) && !defined(GRIMPL_COMPILING_CORE_LIB)
  #include "grackle_misc.h"
  GRIMPL_COMPTIME_WARNING(
    "You are using a deprecated header file; include the public \"grackle.h\" "
    "header file instead! In a future Grackle version, "
    "\"grackle_chemistry_data.h\" may cease to exist (or contents may change "
    "in an incompatible manner)."
  );
#endif

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

  /* minimum temperature to include tabulated metal cooling.
     This is used to supplement the non-equilibrium metal cooling
     from setting metal_chemistry=1, which is only valid up to T~1e4 K.
     - A value of -1.0 is used to specify no minimum temperature, i.e.,
       tabulated cooling is always used.
     - A value of -2.0 indicates the parameter is unset. If unset, the
       following behavior is applied:
       - if metal_chemistry = 0, this is set to -1.0 (i.e., disabled)
       - if metal_chemistry = 1, this is set to 1e4 K. */
  double tabulated_cooling_minimum_temperature;

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

  /* Alternative formulations of interstellar radiation heating
     rate of grains. Both are based on Goldsmith (2001) and
     Krumholz (2014) and are in fact very similar. See
     rate_functions.c for more details.
   * 0: 3.9e-24 / mh / fgr, where fgr is local dust-to-gas ratio
   * 1: 8.60892e-24 / (2.0 * mh) / fgr
   */
  int uniform_grain_isrf_heating_rate;

  /* maximum number of subcycle iterations for solve_chemistry */
  int max_iterations;

  /* flag to exit if max iterations exceeded */
  int exit_after_iterations_exceeded;

  /* number of OpenMP threads, if supported */
  int omp_nthreads;

   /* solver to be used 
   * 1: Gauss-Seidel-Newton-Raphson
   * 2: Gauss-Seidel-only
   * 3: Newton-Raphson-only
   */
  int solver_method;

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
 ******* Generic Interpolation Table ******
 ******************************************/
// as with the other components of chemistry_data_storage, this struct and its
// contents should be treated as an implementation detail
// -> in the future, it would be nice to unify this with cloudy_data
// -> to help facillitate this goal, we track the grid properties in a data
//    structure (that we can probably reuse) that doesn't hold the values that
//    are actually interpolated.
// -> this may be useful since there may be a variable number of grids that
//    need to be interpolated (e.g. cloudy-primordial tables have 3 grids,
//    cloudy-metal tables have 2 grids, generic interpolation tables have 1
//    grid). But maybe we should revisit this in the future?

typedef struct gr_interp_grid_props
{
  /// Rank of dataset
  ///
  /// TODO: do we need this attribute? In most cases, we know the rank
  ///       of a table ahead of time
  long long rank;

  /// Dimension of dataset.
  long long dimension[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Dataset parameter values (in the common case where there there is
  /// constant spacing, we could probably track less data).
  double *parameters[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Value of the constant paramter spacing
  double parameter_spacing[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Length of 1D flattened data grid
  long long data_size;

} gr_interp_grid_props;

typedef struct gr_interp_grid
{
  /// properties of the interpolation grid
  gr_interp_grid_props props;
  /// the actual data that gets interpolated
  double* data;
} gr_interp_grid;


/******************************************
 ********** Forward declarations **********
 ******************************************/

// the following is an opaque type. The definition is intentionally opaque to
// consumers of the Grackle library. Ideally, the entire chemistry_data_storage
// struct will become opaque in the future
struct gr_opaque_storage;

/******************************************
 *** Chemistry and cooling data storage ***
 ******************************************/
typedef struct
{

  /**********************************
   * primordial chemistry rate data *
   **********************************/

  // the (non-radiative) collisional rates, that are handled in a standard
  // manner, are now tracked within opaque_storage

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

  /* The next 8 attributes hold interpolation grids representing collision
   * rates of a species with HI or H2I. The parameters along each axis are:
   *   grid.parameter[0]: log10( ndens_{species} / (dv/dr) )
   *   grid.parameter[1]: log10( T )
   *   grid.parameter[2]: log10( ndens_HI ) OR log10( ndens_H2I )
   */

  /* H2 and HD cooling rates (collision with HI; Hollenbach & McKee 1979) */
  gr_interp_grid LH2;
  gr_interp_grid LHD;

  /* Fine-structure cooling rates (collision with HI; Maio et al. 2007) */
  gr_interp_grid LCI;
  gr_interp_grid LCII;
  gr_interp_grid LOI;

  /* metal molecular cooling rates (collision with H2I; UMIST table) */
  gr_interp_grid LCO;
  gr_interp_grid LOH;
  gr_interp_grid LH2O;

  /* primordial opacity table
   * -> alphap.parameters[0] is log10(mass density)
   * -> alphap.parameters[1] is log10(temperature) */
  gr_interp_grid alphap;

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

  /// holds data that in a manner opaque to consumers of grackle
  struct gr_opaque_storage* opaque_storage;
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
