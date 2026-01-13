/******************************************************************************/
/*
 * ----------Function to Calculate Multispecies Rate Lookup Table---------------
 * 
 * Written by: Ewan Jones
 * Date of Creation: 18/2/2021
 * 
 * Purpose:
 *  Calculates lookup tables for the rate coefficients (and some cooling functions) as a function 
 *  of log(temperature).
 * 
 * Inputs:
 *  *my_chemistry : Pointer to the structure of type chemistry_data which 
 *                   holds flags for the calculations
 *  *my_rates : Pointer to the structure of type chemistry_data_storage
 *               which holds pointers to the rates.
 *  *my_units : Pointer to the structure which defines the units for 
 *               the variables.
 *  *co_length_unit : Pointer to the comoving length unit.
 *  *co_density_unit : Pointer to the comoving density unit.
 *  *iOutput : Integer value which is modified when function finishes.
 * 
 * Outputs:* 
 *  
 * Units:
 *  Most rate coefficients are calculated using cgs units and then converted
 *  for convenience. Units and constants are included explicitly where the 
 *  code units are not used. Also be aware that comoving densities are NOT
 *  used in the rate equations.
 * 
 * Rate Coefficients:
 *   We only list the rates explicitly mentioned by name in this file
 *   @  k13 @     H2I + HI --> 3HI
 *   @  h2dust @  2H + grain --> H2 + grain
 *
 * The logic in this file was originally transcribed from Fortran
 */
/******************************************************************************/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "grackle.h"
#include "grackle_macros.h"
#include "grackle_rate_functions.h"
#include "collisional_rate_props.hpp"  // init_extra_collisional_rates
#include "dust/grain_species_info.hpp"
#include "init_misc_species_cool_rates.hpp"  // init_misc_species_cool_rates
#include "initialize_dust_yields.hpp"  // initialize_dust_yields
#include "initialize_rates.hpp"
#include "internal_types.hpp" // new_CollisionalRxnRateCollection
#include "LUT.hpp" // CollisionalRxnLUT
#include "opaque_storage.hpp" // gr_opaque_storage
#include "phys_constants.h"
#include "status_reporting.h"

// this function pointer type is defined inside an extern "C" block because all
// described functions have C linkage
extern "C" {

//Define the type of a scalar rate function.
typedef double (*scalar_rate_function)(double, chemistry_data*);

}  // extern "C"

namespace { // stuff inside an anonymous namespace are only locally visible

//Define a function which calculates start temperature and increments in log space.
void logT_spacing(double *logT_start, double *d_logT, chemistry_data *my_chemistry)
{
    *logT_start = log(my_chemistry->TemperatureStart);
    *d_logT = (log(my_chemistry->TemperatureEnd) - *logT_start)
                / (my_chemistry->NumberOfTemperatureBins - 1);
}

//Define a function which calculates start temperature and incremenets in log space for dust.
void logT_spacing_dust(double *logT_start_dust, double *d_logT_dust, chemistry_data *my_chemistry)
{
    *logT_start_dust = log(my_chemistry->DustTemperatureStart);
    *d_logT_dust = (log(my_chemistry->DustTemperatureEnd) - *logT_start_dust)
                                / (my_chemistry->NumberOfDustTemperatureBins - 1);
}

//Define a function which is able to add a scalar rate to the calculation.
int add_scalar_reaction_rate(double *rate_ptr, scalar_rate_function my_function, double units,
                        chemistry_data *my_chemistry)
{
    //As these are scalar they have no temperature-dependence.
    (*rate_ptr) = my_function(units, my_chemistry);

    return SUCCESS;
}

/// Initialize values in a pre-allocated rate_ptr
static int init_preallocated_rate(double *rate_ptr, rate_function my_function,
                                  double units, chemistry_data *my_chemistry)
{
  // Calculate temperature spacing.
  double logT_start, d_logT;
  logT_spacing(&logT_start, &d_logT, my_chemistry);

  //Calculate rate at each temperature.
  for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++){
    // Set rate to tiny for safety.
    rate_ptr[i] = tiny;

    // Calculate bin temperature.
    double logT = logT_start + i*d_logT;
    double T = exp(logT);

    // Calculate rate and store.
    rate_ptr[i] = my_function(T, units, my_chemistry);
  }

  return GR_SUCCESS;
}

/// Allocate and initialize a reaction rate
static int add_reaction_rate(double **rate_ptr, rate_function my_function,
                             double units,
                             chemistry_data *my_chemistry)
{
  *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins
                              * sizeof(double));
  return init_preallocated_rate(*rate_ptr, my_function, units, my_chemistry);
}

//Define a function which will calculate the k13dd rates.
int add_k13dd_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate array memory to the pointer.
    *rate_ptr = (double*)malloc(14 * my_chemistry->NumberOfTemperatureBins * sizeof(double));

    //Create a temporary array to retrieve the 14 coefficients calculated by the k13dd function.
    double temp_array[14];

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT;
    logT_spacing(&logT_start, &d_logT, my_chemistry);

    //Calculate k13dd for all temperatures. Store in results array.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++){
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        //Obtain 14 coefficients at this temperature and place them in results array.
        k13dd_rate(T, units, temp_array, my_chemistry);
        for (int j = 0; j < 14; j++){
            (*rate_ptr)[i + (my_chemistry->NumberOfTemperatureBins * j)] = temp_array[j];
        }
    }

    return SUCCESS;
}
    
//Define a function which will calculate h2dust rates.
int add_h2dust_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

// Define a function which will calculate h2dust_C rates.
int add_h2dust_C_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_C_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

// Define a function which will calculate h2dust_S rates.
int add_h2dust_S_reaction_rate(double **rate_ptr, double units, chemistry_data *my_chemistry)
{
    //Allocate memory for h2dust.
    *rate_ptr = (double*)malloc(my_chemistry->NumberOfTemperatureBins * my_chemistry->NumberOfDustTemperatureBins
                        * sizeof(double));

    //Calculate temperature spacing.
    double T, logT, logT_start, d_logT, T_dust, logT_dust, logT_start_dust, d_logT_dust;
    logT_spacing(&logT_start, &d_logT, my_chemistry);
    logT_spacing_dust(&logT_start_dust, &d_logT_dust, my_chemistry);
    
    //Calculate h2dust.
    for (int i = 0; i < my_chemistry->NumberOfTemperatureBins; i++)
    {   
        //Calculate bin temperature.
        logT = logT_start + i*d_logT;
        T = exp(logT);

        for (int j = 0; j < my_chemistry->NumberOfDustTemperatureBins; j++)
        {
            //Calculate dust bin temperature.
            logT_dust = logT_start_dust + j*d_logT_dust;
            T_dust = exp(logT_dust);

            //Calculate rate and store.
            (*rate_ptr)[i + my_chemistry->NumberOfTemperatureBins*j] = h2dust_S_rate(T, T_dust, units, my_chemistry);
        }
    }
    return GR_SUCCESS;
}

/// sets up the species-specific h2dust grain coefficient grids
int setup_h2dust_grain_rates(chemistry_data* my_chemistry,
                             chemistry_data_storage *my_rates,
                             double kUnit) {

  //H2 formation on dust grains with C and S compositions
  if (
    (add_h2dust_C_reaction_rate(&my_rates->h2dustC, kUnit, my_chemistry)
      != GR_SUCCESS) ||
    (add_h2dust_S_reaction_rate(&my_rates->h2dustS, kUnit, my_chemistry)
      != GR_SUCCESS)
  ) {
    return GR_FAIL;
  }

  // initialize my_rates->opaque_storage->h2dust_grain_interp_props
  long long n_Tdust = (long long)(my_chemistry->NumberOfDustTemperatureBins);
  long long n_Tgas = (long long)(my_chemistry->NumberOfTemperatureBins);
  double* d_Td = (double*)malloc(n_Tdust * sizeof(double));
  double* d_Tg = (double*)malloc(n_Tgas * sizeof(double));

  const double logtem_start = std::log(my_chemistry->TemperatureStart);
  const double dlogtem = (std::log(my_chemistry->TemperatureEnd) -
                          std::log(my_chemistry->TemperatureStart)) /
                         (double)(my_chemistry->NumberOfTemperatureBins - 1);

  const double logTdust_start = std::log(my_chemistry->DustTemperatureStart);
  const double dlogTdust = (std::log(my_chemistry->DustTemperatureEnd) -
                            std::log(my_chemistry->DustTemperatureStart)) /
                           (double)(my_chemistry->NumberOfDustTemperatureBins - 1);

  for (long long idx = 0; idx < n_Tdust; idx++) {
    d_Td[idx] = logTdust_start + (double)idx * dlogTdust;
  }
  for (long long idx = 0; idx < n_Tgas; idx++) {
    d_Tg[idx] = logtem_start + (double)idx * dlogtem;
  }

  gr_interp_grid_props* interp_props
    = &(my_rates->opaque_storage->h2dust_grain_interp_props);
  interp_props->rank = 2ll;
  interp_props->dimension[0] = n_Tdust;
  interp_props->dimension[1] = n_Tgas;
  interp_props->parameters[0] = d_Td;
  interp_props->parameters[1] = d_Tg;
  interp_props->parameter_spacing[0] = dlogTdust;
  interp_props->parameter_spacing[1] = dlogtem;
  interp_props->data_size = n_Tdust*n_Tgas;

  return GR_SUCCESS;
}

// Down below we define functionality to initialize the table of ordinary
// collisional rates. If we more fully embraced C++ (and used templates rather
// than C-style function pointers), this could all be a lot more concise.

/// this is a callback function that simply counts up the number of rates that
/// will be used from a grackle::impl::CollisionalRxnRateCollection instance
/// during a given calculation
///
/// @param ctx A pointer to an unsigned int
void tables_counter_callback(grackle::impl::KColProp /* ignored */, void* ctx) {
  *((unsigned int*)ctx) = 1u + *((unsigned int*)ctx);
}


/// tracks state between calls to @ref table_init_callback
struct TablesInitCallbackCtx{
  grackle::impl::CollisionalRxnRateCollection* tables;
  int* used_kcol_rate_indices;
  chemistry_data* my_chemistry;
  double kunit_2bdy;
  double kunit_3bdy;
  unsigned int counter;
};

/// this function is repeatedly called once for each rate table that will be
/// held by a grackle::impl::CollisionalRxnRateCollection instance
///
/// @param rate_prop Specifies information about the current rate
/// @param ctx A pointer to an instance of @ref TablesInitCallbackCtx
void tables_init_callback(struct grackle::impl::KColProp rate_prop,
                          void* ctx) {
  TablesInitCallbackCtx& my_ctx = *((TablesInitCallbackCtx*)ctx);
  init_preallocated_rate(
    my_ctx.tables->data[rate_prop.kcol_lut_index],
    rate_prop.fn_ptr,
    (rate_prop.is_2body) ? my_ctx.kunit_2bdy : my_ctx.kunit_3bdy,
    my_ctx.my_chemistry
  );
  my_ctx.used_kcol_rate_indices[my_ctx.counter] = rate_prop.kcol_lut_index;
  my_ctx.counter++;
}

/// allocate and initialize the table of standard collisional rates that is
/// stored within 
int init_kcol_rate_tables(
  gr_opaque_storage* opaque_storage, chemistry_data* my_chemistry,
  double kunit_2bdy, double kunit_3bdy
) {
  // allocate storage for kcol_rate_tables
  grackle::impl::CollisionalRxnRateCollection* tables =
    new grackle::impl::CollisionalRxnRateCollection;

  // allocate storage within kcol_rate_tables
  *tables = grackle::impl::new_CollisionalRxnRateCollection(
    my_chemistry->NumberOfTemperatureBins);

  // count up the number of rates & allocate used_kcol_rate_indices
  unsigned int n_kcol_rate_indices = 0;
  grackle::impl::visit_rate_props(my_chemistry, &tables_counter_callback,
                                  (void*)(&n_kcol_rate_indices));
  if (n_kcol_rate_indices > (unsigned int)CollisionalRxnLUT::NUM_ENTRIES) {
    GRIMPL_ERROR("something went very wrong when counting up rates");
  } else if (n_kcol_rate_indices == 0) {
    return GrPrintAndReturnErr("this logic shouldn't be executed in a config "
                               "with 0 \"ordinary\" collisional rxn rates");
  }
  int* used_kcol_rate_indices = new int[n_kcol_rate_indices];

  // now its time to initialize the storage
  TablesInitCallbackCtx ctx{
    tables, used_kcol_rate_indices, my_chemistry, kunit_2bdy, kunit_3bdy, 0
  };
  grackle::impl::visit_rate_props(my_chemistry, &tables_init_callback,
                                  (void*)(&ctx));

  // wrap things up!
  opaque_storage->kcol_rate_tables = tables;
  opaque_storage->used_kcol_rate_indices = used_kcol_rate_indices;
  opaque_storage->n_kcol_rate_indices = (int)n_kcol_rate_indices;
  return GR_SUCCESS;
}


/// a function that encodes the algorithm for looking up the `i`th pointer
/// collisional reaction rate from an instance of `chemistry_data_storage`.
///
/// This is intended to be registered with RegBuilder in order to provide
/// access to the various collisional reaction rates.
///
/// @param my_rates The object from which the rate Entry is loaded
/// @param i the index of the queried rate
grackle::impl::ratequery::Entry get_CollisionalRxn_Entry
    (chemistry_data_storage* my_rates, int i) {
  // sanity check! (this shouldn't actually happen)
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr) ||
      (my_rates->opaque_storage->kcol_rate_tables == nullptr)) {
    return grackle::impl::ratequery::mk_invalid_Entry();
  }

  // import new_Entry into the current scope (so we don't need the full name)
  using ::grackle::impl::ratequery::new_Entry;

  double** data = my_rates->opaque_storage->kcol_rate_tables->data;

  // this implementation leverages the following properties of the
  // CollisionalRxnLUT enum:
  // - each entry of collisional_rxn_rate_members.def has a corresponding
  //   enumeration-constant
  // - the very first enumeration constant has a value of 0 (since a value
  //   wasn't explicitly specified)
  // - the value of each other enumeration constants is 1 larger than the value
  //   of the previous value (if a value isn't explicitly specified)
  // - CollisionalRxnLUT::NUM_ENTRIES specifies the number of other enumeration
  //   constants (excluding CollisionalRxnLUT::NUM_ENTRIES) in the enum
  switch (i) {
#define TO_STR(s) #s
#define ENTRY(NAME)                                                            \
  case CollisionalRxnLUT::NAME: { return new_Entry(data[i], TO_STR(NAME)); }
#include "collisional_rxn_rate_members.def"

#undef ENTRY
#undef TO_STR
    default: {
      return grackle::impl::ratequery::mk_invalid_Entry();
    }
  }
}


/// a function that encodes the algorithm for looking up the k13dd
/// collisional reaction rate from an instance of `chemistry_data_storage`.
///
/// This is intended to be registered with RegBuilder in order to provide
/// access to the various collisional reaction rates.
///
/// @param my_rates The object from which the rate Entry is loaded
/// @param i the index of the queried rate
grackle::impl::ratequery::Entry get_k13dd_Entry(
    chemistry_data_storage* my_rates, int i) {
  if (i == 0) {
    double* ptr = (my_rates == nullptr) ? nullptr : my_rates->k13dd;
    return grackle::impl::ratequery::new_Entry(ptr, "k13dd");
  } else {
    return grackle::impl::ratequery::mk_invalid_Entry();
  }
}


} // anonymous namespace

//Definition of the initialise_rates function.
int grackle::impl::initialize_rates(
  chemistry_data *my_chemistry, chemistry_data_storage *my_rates,
  code_units *my_units, double co_length_unit, double co_density_unit,
  ratequery::RegBuilder* reg_builder)
{ 
    // TODO: we REALLY need to do an error check that
    //    my_chemistry->NumberOfTemperatureBins >= 2
    //  is satisfied. We should give some thought about whether this should
    //  produce errors when primordial_chemistry == 0

    // TODO: also add error-checks for:
    //   TemperatureStart, TemperatureEnd, NumberOfDustTemperatureBins,
    //   DustTemperatureStart, DustTemperatureEnd

    int anyDust;
    if ( my_chemistry->h2_on_dust > 0 || my_chemistry->dust_chemistry > 0 || my_chemistry->dust_recombination_cooling > 0) {
        anyDust = TRUE;
    } else {
        anyDust = FALSE;
    }

    //* Obtain the conversion factors which convert between code and physical units. We define a_value = 1 at z = zInit such that a = a_value * [a].
    double timeBase1 = my_units->time_units;
    double lengthBase1 = co_length_unit / (my_units->a_value * my_units->a_units);
    double densityBase1 = co_density_unit * pow(my_units->a_value * my_units->a_units, 3);

    /*
    * 1) Set the dimensions of the non-radiative rate coefficients.
    * 
    *   --Following explanation taken from original fortran code--
    * 
    * 
    *   Note that we have included the units that convert density to 
    *   number density, so the rate equations should look like:
    * 
    *       d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
    * 
    *   where k~ is the dimenionless rate coefficients and d0-2~ are three
    *   dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
    *   dimensionless expansion coefficient.
    * 
    *   rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
    *   rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
    *   rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
    * 
    *   So, k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] ) [~]
    *   
    *   Reminder: the number densities here are normalized with [dens] which is not
    *   a constant (it has a factor a^3), so the number densities must be converted
    *   from comoving to proper.
    */ 
    double kUnit = (pow(my_units->a_units, 3) * mh) / (densityBase1 * timeBase1);
    double kUnit_3Bdy = kUnit * (pow(my_units->a_units, 3) * mh) / densityBase1;

    /*
    * 2) Set the dimensions of the cooling coefficients.
    * 
    *   --Following explanation taken from original fortran code--
    * 
    * 
    *   This equation has a rho because e is the specific energy, not energy/unit volume).
    * 
    *   delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
    *   [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
    *   delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
    * 
    *   So, L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
    *   but [e] = ([a]*[x])**2 / [time]**2, and [a] = 1 / (1 + zri)
    *   
    *   Therefore, [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
    * 
    * 
    *   Note: some of the coffiecients have only one power of n.  These do not have the /a^3
    *   factor, also they have units:
    *     
    *     [L1] = ([a]**2 * [x]**2 * mh) / [time]**3  =  [L] * [dens] * [a]**3 / mh
    *       
    *   This is done through the dom variable in cool.src (some have three powers of n and they
    *   are different by the reciprocal of the above factor multiplying [L]).
    * 
    */
    double coolingUnits = (pow(my_units->a_units, 5) * pow(lengthBase1, 2) * pow(mh, 2))
                          / (densityBase1 * pow(timeBase1, 3));

    // These always need to be allocated since we define other variables by them.
    my_rates->gr_N = (int*)calloc(2, sizeof(int));

    if (my_chemistry->use_primordial_continuum_opacity == 1) {
      initialize_primordial_opacity(my_chemistry, my_rates);
    }

    //* 3) Units for radiative transfer coefficients are 1/[time].

    //* Compute rates for primordial chemistry.

    //* 1) Rate Coefficients (excluding the external radiation field)
    if (my_chemistry->primordial_chemistry > 0){

      // handle all "standard" collisional rates
      if (init_kcol_rate_tables(
            my_rates->opaque_storage, my_chemistry, kUnit, kUnit_3Bdy
          ) != GR_SUCCESS) {
        fprintf(stderr, "Error in init_kcol_rate_tables.\n");
        return GR_FAIL;
      }

      // register the recipe for looking up the "standard" collisional rates
      if (grackle::impl::ratequery::RegBuilder_recipe_1d(
              reg_builder, CollisionalRxnLUT::NUM_ENTRIES,
              &get_CollisionalRxn_Entry, my_chemistry->NumberOfTemperatureBins
          ) != GR_SUCCESS) {
        return GrPrintAndReturnErr("error registering standard collisional "
                                   "reaction rates");
      }

        //--------Calculate coefficients for density-dependent collisional H2 dissociation rate--------
        //
        // idt = 0 calculates coefficients for direct collisional dissociation idt = 1 
        // calculates coefficients for dissociative tunneling as prescribed in 
        // (Martin et al, 1996, ApJ, 461, 265). Code is ran with both values of idt
        // so that coefficients for both cases are calculated.
        // 
        // Results are stored within the k13dd array in chemistry_data_storage. The coefficients
        // corresponding to idt=0 come first, followed by those corresponding to idt=1. There are
        // seven coefficients for each case, and each coefficient has a value at each temperature bin.
        // Therefore the array is structured as follows:
        // k13dd = {coeff1(idt=0, Tbin1), coeff1(idt=0, Tbin2), ..., coeff1(idt=0, TbinFinal), coeff2(idt=0, Tbin1), ..., 
        //          coeff7(idt=0, TbinFinal), coeff1(idt=1, Tbin1), ..., coeff7(idt=1, TbinFinal)}
        add_k13dd_reaction_rate(&my_rates->k13dd, kUnit, my_chemistry);

        // maybe k13dd should be considered multi-dimensional?
        if (grackle::impl::ratequery::RegBuilder_recipe_1d(
                reg_builder, 1, &get_k13dd_Entry,
                my_chemistry->NumberOfTemperatureBins * 14) != GR_SUCCESS) {
          return GrPrintAndReturnErr("error registering k13dd rate");
        }


        //H2 formation on dust grains requires loop over the dust temperature.
        add_h2dust_reaction_rate(&my_rates->h2dust, kUnit, my_chemistry);

        //H2 formation heating terms from Equation 23, Omuaki (2000).
        add_reaction_rate(&my_rates->n_cr_n, n_cr_n_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d1, n_cr_d1_rate, kUnit, my_chemistry);
        add_reaction_rate(&my_rates->n_cr_d2, n_cr_d2_rate, kUnit, my_chemistry);

        //* 2) Cooling and heating Rates (excluding external radiation field)

        //* a) Collisional excitations (Black 1981; Cen 1992).

        add_reaction_rate(&my_rates->ceHI, ceHI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeI, ceHeI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ceHeII, ceHeII_rate, coolingUnits, my_chemistry);

        //* b) Collisional ionizations (Cen 1992, Abel 1996).

        add_reaction_rate(&my_rates->ciHeIS, ciHeIS_rate, coolingUnits, my_chemistry);
        //Collisional ionizations. Polynomial fits from Tom Abel.
        add_reaction_rate(&my_rates->ciHI, ciHI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeI, ciHeI_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->ciHeII, ciHeII_rate, coolingUnits, my_chemistry);

        //* c) Recombinations (Hui & Gnedin 1997, except where noted).

        add_reaction_rate(&my_rates->reHII, reHII_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeII1, reHeII1_rate, coolingUnits, my_chemistry);
        //Dielectronic recombination (Cen, 1992).
        add_reaction_rate(&my_rates->reHeII2, reHeII2_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->reHeIII, reHeIII_rate, coolingUnits, my_chemistry);

        //* d) Bremsstrahlung (Black, 1981 and Spitzer & Hart, 1979)
        add_reaction_rate(&my_rates->brem, brem_rate, coolingUnits, my_chemistry);

        //* e) Molecular hydrogen cooling

        //------Calculate Lepp & Shull rates------

        add_reaction_rate(&my_rates->vibh, vibh_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->hyd01k, hyd01k_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->h2k01, h2k01_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->rotl, rotl_rate, coolingUnits, my_chemistry);
        add_reaction_rate(&my_rates->roth, roth_rate, coolingUnits, my_chemistry);
         
        //------Calculate Galli & Palla (1999) Rates------ 
        //             (as fit by Tom Abel)
        
        //Low density limit (Galli & Palla).
        add_reaction_rate(&my_rates->GP99LowDensityLimit, GP99LowDensityLimit_rate, coolingUnits, my_chemistry);
        //High density limit from HM79.
        add_reaction_rate(&my_rates->GP99HighDensityLimit, GP99HighDensityLimit_rate, coolingUnits, my_chemistry);

        //------Low density rates (Glover & Abel, 2008)------

        //Excitation by HI.
        add_reaction_rate(&my_rates->GAHI, GAHI_rate, coolingUnits, my_chemistry);
        //Excitation by H2.
        add_reaction_rate(&my_rates->GAH2, GAH2_rate, coolingUnits, my_chemistry);
        //Excitation by He.
        add_reaction_rate(&my_rates->GAHe, GAHe_rate, coolingUnits, my_chemistry);
        //Excitation by HII.
        //Collisional excitation rates from Honvault et al (2011, Phys. Rev. Lett., 107, 023201) and 
        //Honvault et al (2012, Phys. Rev. Lett., 108, 109903).
        add_reaction_rate(&my_rates->GAHp, GAHp_rate, coolingUnits, my_chemistry);
        //Excitation by electrons.
        //Based on data from  Yoon et al (2008, J. Phys. Chem. Ref. Data, 37, 913).
        add_reaction_rate(&my_rates->GAel, GAel_rate, coolingUnits, my_chemistry);
        //------New fit to LTE rate------ 
        //from Glover (2015, MNRAS, 451, 2082)
        add_reaction_rate(&my_rates->H2LTE, H2LTE_rate, coolingUnits, my_chemistry);

        // Chiaki & Wise 2019 rates
        // Note: these are still defined in initialize_metal_chemistry_rates.c.
        // They should be moved someday.
        initialize_cooling_rate_H2(my_chemistry, my_rates, coolingUnits);
        initialize_cooling_rate_HD(my_chemistry, my_rates, coolingUnits);

        //* f) HD cooling.

        //HD cooling function has units of ergs cm^3 / s
        //Fit from Coppola et al 2011. LTE (ergs/s) -> hdlte (ergs cm3/s)
        add_reaction_rate(&my_rates->HDlte, HDlte_rate, coolingUnits, my_chemistry);
        //Rate based on (Wrathmall, Gusdorf & Flower, 2007) HD-H collisional excitation rates.
        add_reaction_rate(&my_rates->HDlow, HDlow_rate, coolingUnits, my_chemistry);

        //* g) CIE cooling (Ripamonti & Abel, 2003).

        //Below explanation originally written by Matt Turk:
        //  The above function returns the rate with units of ergs * s^-1 * cm^3 * gram^-1 *
        //  (gram in H2 molecules)^-1. To reproduce equation 5 in RA04, divide c_t_c_r by mh,
        //  so to get erg/s cm^3 we multiply by mh. Divide by 2 for the H2 mass --> number.
        add_reaction_rate(&my_rates->cieco, cieco_rate, coolingUnits, my_chemistry); 

    }//End of primordial chemistry if-statement.

    //* Dust-related processes.
    if (anyDust == TRUE) {

        //Energy transfer from gas to dust grains.
        //(Equation 2.15, Hollenbach & McKee, 1989)
        add_reaction_rate(&my_rates->gas_grain, gasGrain_rate, coolingUnits, my_chemistry); 

        //Electron recombination onto dust grains.
        //(Equation 9, Wolfire et al., 1995)
        add_reaction_rate(&my_rates->regr, regr_rate, coolingUnits, my_chemistry);

        // set up the species-specific h2dust grain coefficient grids
        setup_h2dust_grain_rates(my_chemistry, my_rates, kUnit);

        //Heating of dust by interstellar radiation field, with an arbitrary grain size distribution
        add_scalar_reaction_rate(&my_rates->gamma_isrf2, gamma_isrf2_rate, coolingUnits, my_chemistry);

        //Gas-grain energy transfer, with an arbitrary grain size distribution
        add_reaction_rate(&my_rates->gas_grain2, gasGrain2_rate, coolingUnits, my_chemistry);

        //Grain growth rate
        add_reaction_rate(&my_rates->grain_growth_rate, grain_growth_rate, kUnit, my_chemistry);

    }//End of anyDust if-statement.

    //Below rates are scalar -- temperature independent.

    //* i) Compton cooling (Peebles 1971).
    add_scalar_reaction_rate(&my_rates->comp, comp_rate, coolingUnits, my_chemistry); 

    //* j) Photoelectric heating by UV-irradiated dust (Wolfire 1995).

    add_scalar_reaction_rate(&my_rates->gammah, gammah_rate, coolingUnits, my_chemistry); 
    //Heating of dust by interstellar radiation field.
    //(Equation B15, Krumholz, 2014)
    add_scalar_reaction_rate(&my_rates->gamma_isrf, gamma_isrf_rate, coolingUnits, my_chemistry); 

    // Miscellaneous species-based cooling rates
    if (grackle::impl::init_misc_species_cool_rates(my_chemistry, my_rates, my_units) != GR_SUCCESS) {
      fprintf(stderr, "Error in initialize_metal_chemistry_rates.\n");
      return GR_FAIL;
    }

    // Dust Grain Species Information
    // (it may make sense want to handle more of the dust separately)
    if (my_chemistry->dust_species > 0) {
      my_rates->opaque_storage->grain_species_info = 
        new grackle::impl::GrainSpeciesInfo;
      *(my_rates->opaque_storage->grain_species_info) =
        grackle::impl::new_GrainSpeciesInfo(my_chemistry->dust_species);
    }

    // Dust rates
    if (grackle::impl::initialize_dust_yields(my_chemistry, my_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in initialize_dust_yields.\n");
      return FAIL;
    }

    return SUCCESS;
}

