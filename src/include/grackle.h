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
/// > There are a number of considerations before stablization. Some are
/// > discussed in docstrings. Below we summarize more generic considerations:
/// > 1. Do we want to be able to use this API to eventually support a dynamic
/// >    set of rates? If so, then:
/// >    - we are free to get rid of the rateid handle (there are still some
/// >      reasons to retain it, but we need to decide whether simplicity is
/// >      better)
/// >    - we should probably update all of the function signatures so that a
/// >      pointer to chemistry_data_storage* is passed as an argument to each
/// >      of the functions.
/// > 2. Should we prevent direct access to the underlying data pointers?
/// >    Instead we would provide getter and setters. We discuss this further
/// >    down below (in the relevant function's API). This is important for
/// >    different backends.
/// > 3. We should generally consider whether the API is consistent enough with
/// >    APIs for accessing other data structures.
/// >    - currently the dynamic parameter API is grackle's only other
/// >      data-access API, but we could imagine creating an API for the fields
/// >      (in order to achieve ABI stability)
/** @{ */

/// @typedef grunstable_rateid_type
/// @brief Type of rate-ids returned to the user
///
/// > [!note]
/// > Before stablizing, we should consider:
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

/// Access the pointer associated with the rateid from myrates
///
/// The behavior is **NOT** currently defined if the `my_rates` struct isn't
/// configured to use the specified rate.
///
/// @return The pointer to the accessed data (whether the rate data corresponds
///    to an array or scalar). If an invalid `rate_id` is specified, this
///    returns NULL.
///
/// > [!note]
/// > Before stablizing, we should consider whether we actually want this
/// > function. It may be better to replace it with a getter (that copies the
/// > rate data to the supplied buffer) and a setter (that copies the provided
/// > data out of the provided buffer). This may be appealing for a number of
/// > reasons:
/// > 1. When we add GPU support, the rate data may live permanently on the GPU
/// >    (or some data may live on the GPU while other data lives on the CPU).
/// >    With dedicated setters/getters, we could always require that the
/// >    provided buffer data is on the CPU.
/// > 2. We have the option to conditionally disable the setter. For example,
/// >    if we choose to support custom rates, I assume we will want to specify
/// >    all custom choices as part of initialization (for simplicity) and then
/// >    we may want to deprecate/remove the setter. (One could also imagine,
/// >    disabling the setter when using GPUs).
/// > 3. It eliminates a certain class of memory-lifetime bugs (People could
/// >    try use the pointer returned by the incarnation of this function after
/// >    destroying Grackle. This problem can't happen with the proposed
/// >    getter/setter)
double* grunstable_ratequery_get_ptr(
  chemistry_data_storage* my_rates, grunstable_rateid_type rate_id
);

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
  GRUNSTABLE_QPROP_TYPE = 3,
  GRUNSTABLE_QPROP_MAXITEMSIZE = 4,
  // I don't like the next one
  GRUNSTABLE_QPROP_WRITABLE = 5,
};

/// Query a property of the specified rate
///
/// @param[in]  my_rates The object being queried
/// @param[in]  rate_id The id of the rate for which the property is queried
/// @param[in]  prop_kind The proprty to query
/// @param[out] ptr The pointer where the property is recorded
///
/// @returns GR_SUCCESS if successful. Otherwise, a different value is returned.
///
/// The behavior is undefined when @p my_rates is a `nullptr`, @p ptr is a
/// nullptr or @p ptr doesn't have enough space to store the queried property
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

/** @}*/ // end of group

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */


// this is important for letting us diagnose improper use of private headrs
#if GRIMPL_PUBLIC_INCLUDE == 0
  #undef GRIMPL_PUBLIC_INCLUDE
#endif

#endif /* __GRACKLE_H__ */
