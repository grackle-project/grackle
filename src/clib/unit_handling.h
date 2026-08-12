/***********************************************************************
/
/ Implement and declare automatic unit-handling
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef UNIT_HANDLING_H
#define UNIT_HANDLING_H

#include <stdio.h>  // provides fprintf
#include "grackle.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



/// Automatically compute a code_units instance when Grackle is configured
/// with automatic code-units.
///
/// @param specified_units Pointer to the specified units
/// @param my_rates Pointer to the chemistry_data_storage object that contains
///     the initial rates that were specified when the grackle_solver was
///     initialized
/// @param current_a_value The specified current scale-factor
/// @param unit_handling The value of Grackle's unit handling flag
/// @param func_name Name of the function that is calling the error check
///
/// @returns a fully initialized code_units struct. If there was any kind of
///     issue, then the a_units member will be negative.
///
/// @note
/// This is declared inline to reduce as much overhead as possible
///
/// @note
/// This is currently a proof-of-concept, some improvements are definitely
/// possible. For example, we could:
///   * directly pass in my_rates->initial_units rather than my_rates
///   * refactor gr_query_units so that we could directly call into its
///     internals (the error-checking and string comparisons introduce a bunch
///     of unnecessary overhead)
static inline code_units determine_code_units (
    const code_units* specified_units,
    const chemistry_data_storage* my_rates,
    double current_a_value,
    int unit_handling,
    const char* func_name)
{
  code_units out;
  out.a_units = -1.0;

  if (unit_handling == GR_UNIT_HANDLING_LEGACY) {
    // handle error-cases
    if (current_a_value != GR_SPECIFY_INITIAL_A_VALUE) {
      fprintf(stderr,
          "ERROR in %s: current_a_value shouldn't be specified when using "
          "Grackle's legacy unit handling\n", func_name);
      out.a_units = -1.0;
      return out;
    } else if (specified_units == NULL) {
      fprintf(stderr,
          "ERROR in %s: code_units argument can't be NULL when using "
          "Grackle's legacy unit handling\n", func_name);
      return out;
    }

    // this is the happy path
    out = *specified_units;
    return out;

  } else if (unit_handling == GR_UNIT_HANDLING_AUTOMATIC) {
    // handle error-case
    if (specified_units != NULL) {
      fprintf(stderr,
          "ERROR in %s: code_units argument must be NULL when using Grackle's "
          "automatic unit handling\n", func_name);
      return out;
    }

    // this is the happy path
    out = my_rates->initial_units;
    out.a_value = current_a_value;
    out.density_units = gr_query_units(my_rates, "density_units",
                                       current_a_value);
    out.length_units = gr_query_units(my_rates, "length_units",
                                      current_a_value);
    out.time_units = gr_query_units(my_rates, "time_units",
                                    current_a_value);
    set_velocity_units(&out);
    return out;

  } else {
    fprintf(stderr,
        "ERROR in %s: Grackle's unit_handling parameter has an unrecognized "
        "value\n", func_name);
    return out;
  }

}

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* UNIT_HANDLING_H */

