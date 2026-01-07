/***********************************************************************
/
/ Declare utility functions used internally by Grackle (across routines)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "grackle_chemistry_data.h"
#include "grackle_types.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/// Perform an error check related to self-shielding. If the check fails, the
/// function returns FAIL and prints an error message to stderr
///
/// @param my_chemistry Holds configuration of chemistry solver
/// @param fields Specify the field names
/// @param func_name Name of the function that is calling the error check
int self_shielding_err_check(const chemistry_data *my_chemistry,
                             const grackle_field_data *fields,
                             const char* func_name);


#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus

#endif  // UTILS_H
