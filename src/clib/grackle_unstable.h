/***********************************************************************
/
/ Unstable Grackle function prototypes
/
/ This header exposes a subset of functionality to external programs.
/ As a rule, no guarantees are made about the stability/long-term
/ compatability of functions exposed in this header (more details may be
/ provided on a function-by-function basis). Functionality in this header
/ might include:
/   - experimental functionality
/   - debugging functionality
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

// this header should NOT directly be included by other public headers.
// Downstream users should directly include this file

#include "grackle_types.h"


// debugging functionality
// =======================


/// function used to dump state of various grackle objects for debugging
/// purposes. The information will be written to a newly created HDF5 file (at
/// the path specified by `fname`) OR to a location specified within an already
/// open hdf5 file (at the location specified by `dest_hid`)
///
/// @param[in] fname When not NULL, specifies the path to the file that will be
///     created by this function. When NULL, the data will be written to the
///     location specified by dest_hid
/// @param[in] dest_hid Used to optionally specify an hdf5 handle that
///     corresponds to the location where the data will be written. This should
///     be set to `-1` when `fname` is not `NULL`.
/// @param[in] chemistry_data Pointer to chemistry_data. Intended to be the
///     `chemistry_data` instance used during initialization of the
///     `chemistry_data_storage` struct (i.e. within `initialize_chemistry_data`
///     or `local_initialize_chemistry_data`)
/// @param[in] initial_code_units Pointer to the `code_units` instance used
///     during initialization of the `chemistry_data_storage` struct. If any
///     member (namely `a_value`) has been modified since then, you should pass
///     this `NULL`.
/// @param[in] current_code_units Pointer to the `code_units` instance currently
///     being used (The main purpose of this argument is to track changes in
///     `a_value`).
/// @param[in] current_code_units Pointer to the `code_units` instance currently
///     being used (The main purpose of this argument is to track changes in
///     `a_value`).
///
/// @note
/// We explicitly elect to use `long long` as the datatype of `dest_hid` rather
/// than `hid_t` to ensure that the hdf5-library is not a public dependency of
/// this library (i.e. to avoid including the hdf5 header in public headers)
/// 
int grunstable_h5dump_state(const char* fname, long long dest_hid,
                            const chemistry_data* chemistry_data,
                            const code_units* initial_code_units,
                            const code_units* current_code_units,
                            const grackle_field_data* my_fields);
