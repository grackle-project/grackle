/***********************************************************************
/
/ Grackle function prototypes
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_H__
#define __GRACKLE_H__

// this logic must occur before including Grackle's headers
#ifndef GRIMPL_PUBLIC_INCLUDE
  #define GRIMPL_PUBLIC_INCLUDE 0
#elif GRIMPL_PUBLIC_INCLUDE == 1
  // we temporarily allow this case for internal use within the Grackle lib to
  // avoid merge conflicts (external Grackle users shouldn't depend on this!)
#else
  #error "it is an error for GRIMPL_PUBLIC_INCLUDE to be defined"
#endif /* GRIMPL_PUBLIC_INCLUDE */


// include the private headers
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// standard error codes returned by the grackle functions
#define GR_SUCCESS 1
#define GR_FAIL 0

#define GR_SPECIFY_INITIAL_A_VALUE -1

extern int grackle_verbose;

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

double get_velocity_units(const code_units *my_units);

void set_velocity_units(code_units *my_units);

double get_temperature_units(const code_units *my_units);

double gr_query_units(const chemistry_data_storage * my_rates,
                      const char* units_name, double current_a_value);

int set_default_chemistry_parameters(chemistry_data *my_grackle);

int local_initialize_chemistry_parameters(chemistry_data *my_chemistry);

int initialize_chemistry_data(code_units *my_units);

int local_initialize_chemistry_data(chemistry_data *my_chemistry, 
                                    chemistry_data_storage *my_rates,
                                    code_units *my_units);

int* local_chemistry_data_access_int(chemistry_data* my_chemistry,
                                     const char* param_name);
double* local_chemistry_data_access_double(chemistry_data* my_chemistry,
                                           const char* param_name);
char** local_chemistry_data_access_string(chemistry_data* my_chemistry,
                                          const char* param_name);

const char* param_name_int(unsigned int i);
const char* param_name_double(unsigned int i);
const char* param_name_string(unsigned int i);

unsigned int grackle_num_params(const char* type_name);

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value);

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value);

int calculate_cooling_time(code_units *my_units,
                           grackle_field_data *my_fields,
                           gr_float *cooling_time);

int local_calculate_cooling_time(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates,
                                 code_units *my_units,
                                 grackle_field_data *my_fields,
                                 gr_float *cooling_time);

int calculate_dust_temperature(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *dust_temperature);

int local_calculate_dust_temperature(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *dust_temperature);

int calculate_gamma(code_units *my_units,
                    grackle_field_data *my_fields,
                    gr_float *my_gamma);

int local_calculate_gamma(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *my_gamma);

int calculate_pressure(code_units *my_units,
                       grackle_field_data *my_fields,
                       gr_float *pressure);

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure);

int calculate_temperature(code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *temperature);

int local_calculate_temperature(chemistry_data *my_chemistry,
                                chemistry_data_storage *my_rates,
                                code_units *my_units,
                                grackle_field_data *my_fields,
                                gr_float *temperature);

int free_chemistry_data(void);

int local_free_chemistry_data(chemistry_data *my_chemistry, chemistry_data_storage *my_rates);

grackle_version get_grackle_version(void);

// the following function is NOT part of the public API (and it NEVER will be)
// - it does a lot of heavy lifting for `gr_check_consistency`
// - the signature may change in the future
int grimpl_check_consistency_(int gr_float_hdrsize);

static inline int gr_check_consistency(void) {
  return grimpl_check_consistency_((int)sizeof(gr_float));
}

int gr_initialize_field_data(grackle_field_data *my_fields);

/// unstable/experimental machinery
/// ===============================
/// entities defined down below (prefixed with `grunstable_`) are experimental
/// features that are not yet part of the public API. We reserve the right to
/// change, modify, remove these features at ANY time (we may provide more
/// details on a function by function basis)

/// @defgroup Experimental API for dynamically Accessing Rates
///
/// This group of entities defines the experimental API for accessing rates in
/// a manner that is API/ABI stable without needing to access implementation
/// details (i.e. the internals if `chemistry_data_storage`).
///
/// In the longer term it would be REALLY nice to support a truly dynamic set
/// of rates (i.e. the collection of allowed rates need not be known at
/// compile-time). If this is a serious goal, it may be nice if this API is
/// general enough so that it doesn't need to change in that transition.
///
/// > [!note]
/// > The API currently uses a handle called rateid. You lookup the handle
/// > from a string and then use the handle to lookup various properties.
/// > Currently, we just support looking up the array, but we could add
/// > functions to query other properties like shape or rate-type (i.e. is it
/// > a cooling rate or rxn rate).
/// > - while this may seem like overkill, this is a lesson learned from a
/// >   implementing dynamic api for parameters. Mapping a string to a value
/// >   (or entity) without any mutable dynamic memory allocation in a
/// >   simple, portable-threadsafe manner (without undefined behavior),
/// >   requires a `O(N)` linear search. (Using the API to iterate
/// >   over all rates has `O(N^2)` complexity)
/// > - but, if we're willing to allocate a dict/map when initializing Grackle,
/// >   then I would be happy to get rid of rateid (we'd need to do this
/// >   anyways to support any variant of this API for dynamic rates)
///
/// > [!note]
/// > There are a number of considerations before stabilization. Some are
/// > discussed in docstrings. Below we summarize more generic considerations:
/// > 1. Do we want to be able to use this API to eventually support a dynamic
/// >    set of rates? If so, then:
/// >    - we are free to get rid of the rateid handle (there are still some
/// >      reasons to retain it, but we need to decide whether simplicity is
/// >      better)
/// >    - we should probably update all of the function signatures so that a
/// >      pointer to chemistry_data_storage* is passed as an argument to each
/// >      of the functions.
/// > 2. We should generally consider whether the API is consistent enough with
/// >    APIs for accessing other data structures.
/// >    - currently the dynamic parameter API is grackle's only other
/// >      data-access API, but we could imagine creating an API for the fields
/// >      (in order to achieve ABI stability)
/// > 3. I think it may be beneficial to follow the example of yt and make keys
/// >    2 element strings (rather than just a single string).
/// >    - this is attractive because it would let us group things together
/// >    - for example, we could group k1, k2, k3, ... under a group of
/// >      collisional reaction rates. It might be useful to list yields for
/// >      different injection pathways (e.g. used in the multi-dust grain
/// >      species model) under a separate group.
/// > 4. We should consider improving "ergonomics" of this API
/// >    - accessing arrays of strings are currently very clunky.
/// >      - Maybe we should avoid requiring a deep-copy? And just require the
/// >        user to allocate a buffer to hold pointers to internal strings? I
/// >        don't love this.
/// >      - There's a piece of me that wonders if we could entirely eliminate
/// >        arrays of strings (while still supporting single strings) by
/// >        adopting keys composed of 2 strings and refactoring
/// >        grackle_field_data to only provide access to strings through keys
/// >        (as in PR 271)
/// >    - property-querying also feels a little clunky... See the note
/// >      attached to that docstring for more details
/// > 5. We need to come up with a consistent strategy when rates are not
/// >    defined.
/// >    - In a lot of cases, I think it's ok to simply "not" define the rate.
/// >      For example, when primordial_chemistry=0, I don't think we should
/// >      define rates like k1, k2, k3, ... If we ever want to support a
/// >      dynamic API, I think this approach makes a lot of sense in a lot of
/// >      cases.
/// >    - There's a separate, related question of whether we should support the
/// >      idea of an empty entry (i.e. a 1D array with 0 entries). This may
/// >      come up in the context of arrays of strings
/// >      - for example, if we normally list the names of all dust grain
/// >        species, we may want an empty string when there are no dust grain
/// >        grain species
/// >      - but does it make sense to support an empty scalar?
/// >      - as we discuss earlier in this list, it's plausible certain changes
/// >        could make string-datasets unnecessary. In that case, this point
/// >        may become moot.
/** @{ */

/// @typedef grunstable_rateid_type
/// @brief Type of rate-ids returned to the user
///
/// > [!note]
/// > Before stabilizing, we should consider:
/// > - if we should use `int64_t` rather than `long long`
/// > - if we want to make this a more generic name like `gr_idtype` (i.e. if
/// >   we plan to use rates in other parts of the code)
typedef long long grunstable_rateid_type;

/// Lookup the rate_id associated with a rate of a given name.
///
/// If the rate is unknown, the same value will be returned as when you pass a
/// `NULL` pointer (the precise returned value may change at any given time)
///
/// > [!important]
/// > The mapping of variable names to rate_id values is not stable. You should
/// > only rely upon the mapping for a given version of Grackle.
///
/// > [!note]
/// > While this is unstable, not all rates may be known yet (but, the number
/// > rates are all accessible).
grunstable_rateid_type grunstable_ratequery_id(
    const chemistry_data_storage* my_rates, const char* name);

/// Copy data associated with the @p rate_id in @p my_rates into buf
///
/// The behavior is currently undefined if:
/// - the @p my_rates argument is a NULL pointer
/// - the @p my_rates argument isn't configured to use the specified rate
/// - the @p buf argument is NULL
/// - the @p buf isn't large enough to store the copied rates
///
/// @param[in] my_rates The object from which data is retrieved
/// @param[in] rate_id The id of the rate for which the data is retrieved
/// @param[out] buf The buffer that the function writes to
///
/// @return Returns GR_SUCCESS if successful. If an invalid @p rate_id is
///    specified, this returns a different value
///
/// > [!note]
/// > Before stablizing, we should consider whether we want to add an argument
/// > that specifies the supplied size of `buf` (and provide well-defined
/// > behavior when `buf` is too small
int grunstable_ratequery_get_f64(
  chemistry_data_storage* my_rates, grunstable_rateid_type rate_id,
  double* buf
);


/// Overwrite the data associated with @p rate_id in @p my_rates with the
/// values provided by @p buf
///
/// The behavior is currently undefined if:
/// - the @p my_rates argument is a NULL pointer
/// - the @p my_rates argument isn't configured to use the specified rate
/// - the @p buf argument is NULL
/// - the @p buf isn't large enough to store the copied rates
///
/// @param[out] my_rates The object from which data is retrieved
/// @param[in] rate_id The id of the rate for which the data is retrieved
/// @param[in] buf The buffer that the function writes to
///
/// @return Returns GR_SUCCESS if successful. If an invalid @p rate_id is
///    specified, this returns a different value
///
/// > [!note]
/// > Before stabilizing, we should consider whether we want to add an argument
/// > that specifies the supplied size of `buf` (and provide well-defined
/// > behavior when `buf` is too small
int grunstable_ratequery_set_f64(
  chemistry_data_storage* my_rates, grunstable_rateid_type rate_id,
  const double* buf
);

int grunstable_ratequery_get_str(
  chemistry_data_storage* my_rates, grunstable_rateid_type rate_id,
  char* const * buf
);


/// Describe types known to grackle
///
/// > [!note]
/// > Before stabilizing, we should consider:
/// > 1. whether this is general enough to be useful throughout grackle.
/// >    - For example, we could imagine adding a function in the future to the
/// >      chemistry_data dynamic-access API that generically tries to query
/// >      datatype.
/// >    - In that case, we would want to reuse these macros
/// > 2. if these should actually be defined as macros (that may make it easier
/// >    to support them from Fortran
enum grunstable_types {
  GRUNSTABLE_TYPE_F64 = 1,
  GRUNSTABLE_TYPE_STR = 2,
};

/// Describe Rate-Query Property Types
///
/// > [!note]
/// > It may make more sense to use macros if we want to support these from
/// > Fortran
///
/// > [!important]
/// > Users should obviously avoid hardcoding values in their codebase.
enum grunstable_ratequery_prop_kind {
  GRUNSTABLE_QPROP_NDIM = 1,
  GRUNSTABLE_QPROP_SHAPE = 2,
  GRUNSTABLE_QPROP_MAXITEMSIZE = 3,
  GRUNSTABLE_QPROP_WRITABLE = 4,
  GRUNSTABLE_QPROP_DTYPE = 5,
};

/// Query a property of the specified rate
///
/// @param[in]  my_rates The object being queried
/// @param[in]  rate_id The id of the rate for which the property is queried
/// @param[in]  prop_kind The property to query
/// @param[out] ptr The pointer where the property is recorded
///
/// @returns GR_SUCCESS if successful. Otherwise, a different value is returned.
///
/// The behavior is undefined when @p my_rates is a `nullptr`, @p ptr is a
/// nullptr or @p ptr doesn't have enough space to store the queried property
///
/// @note
/// It turns out that using this interface is a little clunky...
/// - before stabilization, it may be better to change this function so that
///   - it accepts an additional argument specifying the total size of the
///     provided buffer (and report an error if the buffer isn't long enough)
///   - it would also be great (but not essential) to change the return-value
///     to specify:
///     - the number of entries that were filled written to ptr by this
///       function
///     - aside: we need to return 0 for GRUNSTABLE_QPROP_SHAPE when
///       querying a scalar quantity
///     - a negative value would denote an error
///   - optionally, when ptr is a nullptr, we could have the function return
///     the number of required ptr needs (similar to the interface for
///     snprintf), but this isn't necessary
/// - Doing this would allow users to write extra code to handle queries of
///   GRUNSTABLE_QPROP_SHAPE in a special way (i.e. currently you **NEED** to
///   query GRUNSTABLE_QPROP_NDIM, first to make sure you allocate enough space)
int grunstable_ratequery_prop(const chemistry_data_storage* my_rates,
                              grunstable_rateid_type rate_id,
                              enum grunstable_ratequery_prop_kind prop_kind,
                              long long* ptr);

/// Query the name (and optionally the rate_id) of the ith registered rate
///
/// > [!warning]
/// > The order of parameters may change between different versions of Grackle
///
/// @param[in]  my_rates The object being queried
/// @param[in]  i the index of the accessed rate
/// @param[out] out_rate_id A pointer to store the rate of the queried rate_id.
///    The behavior is **NOT** currently well defined when there are `i` or
///    fewer registered rates.
/// @result Pointer to the string-literal specifying the rate's name. This is
///    `NULL`, if there are `i` or fewer registered rates.
const char* grunstable_ith_rate(
  const chemistry_data_storage* my_rates, unsigned long long i,
  grunstable_rateid_type* out_rate_id
);

/// Query the number of rates accessible through the ratequery API
unsigned long long grunstable_ratequery_nrates(
    const chemistry_data_storage* my_rates);

/** @}*/ // end of group

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#if !defined(GRIMPL_COMPILING_CORE_LIB)
// avoid leaking these macros used for defining implementation details out of
// the core library
// -> while this not essential, it's considered a good practice to do this
//    (people can't use the macros if they aren't defined)
// -> note: these definitions will leak if external codes include deprecated
//    (but that's probably ok)

#undef GRACKLE_CLOUDY_TABLE_MAX_DIMENSION
#undef GRIMPL_MAX_INJ_PATHWAYS
#endif

// this is important for letting us diagnose improper use of private headrs
#if GRIMPL_PUBLIC_INCLUDE == 0
  #undef GRIMPL_PUBLIC_INCLUDE
#endif

#endif /* __GRACKLE_H__ */
