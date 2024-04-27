/***********************************************************************
/
/ Private internal header that declare "visitor functionality" designed to 
/ help dump Grackle's internal state (for debugging purposes)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/


// Originally, this functionality and the functions for dumping state were all
// defined in one place, but that produced an overwhelming amount of code

// General Design Philosophy
// =========================
//
// For some of Grackle's public structs, we want to dump the state in multiple
// ways (e.g. to a human readable format and a binary format).
//
// For a given format and a given struct, we need to apply formatting rules to
// each member of that struct. This can be implemented by:
// - defining a set of formatting-functions (specific to the dump-format)
// - iterating over each struct-member: for each struct-member, pass the name
//   and value of the struct-member to one of these formatting-functions (we
//   choose the struct-member based on the type of that struct-member)
//
// This is the strategy we choose to adopt for dumping Grackle's simplest
// public structs. We choose to generalize this process:
// - we adopt a struct used to store the various formatting-functions pointers
// - we define a function that applies the function pointers to the each member
//   of the struct
//
// This control-flow resembles the Visitor design pattern.
//
// Our use of function pointers slows things down a little bit, but that's
// reasonable given that this is purely for debugging/informational purposes
// (i.e. it is not invoked inside of any calculation loop). If speed is a
// concern we could make use of macros.
//
// NOTE: I am very hesitant to expose this functionality (at least in its
//       current form) as part of the public API right now
// - Some context: at this time of writing, we are preparing to transcribe a
//   lot of Grackle's internals from Fortran to C/C++. I think it might be
//   useful to internally reorganize some things (e.g. reorganize some members
//   of large structs so that they are composed of smaller structs)
//
// - It's worth mentioning that we can adapt this functionality for the sake of
//   deserializing data (it currently follows a control-flow very similar to
//   the ones used by various C++ serialization frameworks, like Charm++'s pup
//   framework (of course we rely on function pointers with context objects and
//   they use classes with function-overloading). However, if we support such a
//   case, users could start relying on this functionality even if it isn't
//   part of the public API
//
// - With all of that said, if users express an interest in having this
//   functionality, we could definitely add this functionality into the public
//   API (I'm primarily concerned with locking into the current implementation,
//   which isn't very polished)

#ifndef VISITOR_UTILS
#define VISITOR_UTILS

#include <stdio.h>
#include <hdf5.h>

#include "grackle_types.h" // gr_float

// when rank is less than 3, don't fill in unused dimensions
typedef struct array_props_{ int rank; int dimensions[3]; } array_props_;

/// A struct used to temporarily hold the function pointers for the purpose of
/// dumping state (usually: it is used to "visit" members of a struct)
///
/// The common control-flow related to this function is similar to a class --
/// there are usually funcs for setup/constructor and teardown/destruction.
/// We could implement things a lot more concisely in C++ (with inheritance
/// or templates)
///
/// @note
/// If we ultimately decide to expose this as part of the public API, we should
/// be very careful. In particular, I think we should expose an opaque type and
/// not provide direct access to this struct (this will help us with adding
/// new members to this struct - such as to support new types)
///
/// @note
/// Rather than holding 5 function pointer, this could have held just 1
/// function pointer with an extra argument to denote type. We opt not to do
/// that for extra type-checking (e.g. if we change the type of a struct and
/// we forget to update the visit-function, the compiler should complain)
struct member_visitor_ {
  void* context;
  int (*fn_INT)(void*, const char*, int);
  int (*fn_DOUBLE)(void*, const char*, double);
  int (*fn_STRING)(void*, const char*, const char*);
  int (*fn_INTARR)(void*, const char*, const int*, array_props_);
  int (*fn_GRFLOATARR)(void*, const char*, const gr_float*, array_props_);
};


/// Create a visitor object with a JSON plugin.
///
/// The visitor is configured to write the members of an object out to a file,
/// `fp`, in the JSON format
struct member_visitor_ create_json_visitor_(FILE *fp);

/// delete/cleanup the visitor object (with hdf5 plugin)
///
/// In addition to freeing memory, this writes some extra information to the
/// underlying file.
void free_json_visitor_(struct member_visitor_* visitor);

/// The idea is to track additional useful information inside this struct and
/// pass it around throughout our function stack...
///
/// To start out, we will just use it to track custom commonly used dtypes (so
/// that we don't need to constantly create new versions of this dtype)
typedef struct contextH5_ {
  hid_t gr_floattype;
  hid_t var_strtype;
} contextH5_;

int initialize_contextH5_(contextH5_* h_ctx);
int cleanup_contextH5_(contextH5_* p);

/// Create a visitor object with an hdf5 plugin.
///
/// The visitor is configured to write the members of an object out to an hdf5
/// file at a location specified by `loc_id`.
struct member_visitor_ create_h5_visitor_(hid_t loc_id, contextH5_* h_ctx,
                                          const char* group_name);

/// delete/cleanup the visitor object (with hdf5 plugin)
///
/// In addition to freeing memory, this closes the group that all attributes
/// were written into.
///
/// @param[in] visitor This visitor to cleanup
/// @param[in] objname_errfmt Optionally specify the name of the object that
///     we were trying to dump (only used for error-formatting)
///
/// @return denotes whether the visitor ever encountered any issues
int free_h5_visitor_(struct member_visitor_* visitor,
                     const char* objname_errfmt);

/// Creates an "annotated" HDF5 group called `name`
///
/// In more detail:
/// - this is just an ordinary hdf5 group with some 1 or more standardized
///   attributes that are used to annotate metadata (the names of these
///   attributes are prefixed with "@:?GRACKLEDUMP_")
/// - The group should be closed by `h5dump_create_annotated_grp_`
/// - The presence of the "@:?GRACKLEDUMP_VERSION" attribute denotes that the
///   group is an annotated group
/// - The "@:?GRACKLEDUMP_DONE" attribute is not initially written when we open
///   the group; it's ONLY written when we properly close the group. 
///   - when the group is properly closed without any perceived issues, the
///     associated value will be explicitly overwritten with a value of 1.
///   - when the group is properly closed, but issues have been identified, the
///     associated value will be overwritten with some arbitrary negative value
///   - the idea is this provides at lease some assurance about the validity of
///     the group's contents. If there is no attribute OR the value is not 1,
///     then there was a logical error or something went wrong during execution
///     (a segmentation fault, the program was interupted, etc.).
///
/// @note
/// This may seem like overkill, but in my experience this stuff usually comes
/// in handy
hid_t h5dump_create_annotated_grp_(hid_t loc_id, contextH5_* h_ctx,
                                   const char* name);

/// Closes an "annotated" HDF5 group
///
/// @note
/// See the docstring for h5dump_create_annotated_grp_ for more details
int h5dump_close_annotated_grp_(hid_t grp_id, contextH5_* h_ctx,
                                int indicate_errs);

#endif /* VISITOR_UTILS */
