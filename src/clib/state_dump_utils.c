/***********************************************************************
/
/ Utilities related to dumping state (for debugging purposes)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // strlen
#include <hdf5.h>

#include "grackle.h"
#include "grackle_macros.h"

// when rank is less than 3, don't fill in unused dimensions
typedef struct array_props_{ int rank; int dimensions[3]; } array_props_;

// General Design Philosophy
// =========================

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
// - It's worth mentioning that we can easily adapt this functionality for the
//   sake of deserializing data (it currently follows a control-flow very
//   similar to the ones used by various C++ serialization frameworks, like
//   Charm++'s pup framework (of course we rely on function pointers with
//   context objects and they use classes with function-overloading). However,
//   if we support such a case, users could start relying on this functionality
//   even if it isn't part of the public API
//
// - With all of that said, if users express an interest in having this
//   functionality, we could definitely add this functionality into the public
//   API (I'm primarily concerned with locking into the current implementation,
//   which isn't very polished)


/// A struct used to temporarily hold the function pointers for the purpose of
/// dumping state.
///
/// @note
/// If we ultimately decide to expose this as part of the public API, we should
/// be very careful. In particular, I think we should expose an opaque type and
/// not provide direct access to this struct (this will help us with adding
/// new members to this struct - such as to support new types)
struct member_visitor_ {
  void* context;
  int (*fn_INT)(void*, const char*, int);
  int (*fn_DOUBLE)(void*, const char*, double);
  int (*fn_STRING)(void*, const char*, const char*);
  int (*fn_INTARR)(void*, const char*, const int*, array_props_);
  int (*fn_GRFLOATARR)(void*, const char*, const gr_float*, array_props_);
};

static void visit_parameters_(const chemistry_data * my_chemistry,
                              struct member_visitor_ * visitor){
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) \
    visitor->fn_ ## TYPE (visitor->context, #FIELD, my_chemistry->FIELD);
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
}

static void visit_version_(const grackle_version * gversion,
                           struct member_visitor_ * visitor) {
  visitor->fn_STRING(visitor->context, "version", gversion->version);
  visitor->fn_STRING(visitor->context, "branch", gversion->branch);
  visitor->fn_STRING(visitor->context, "revision", gversion->revision);
}

static void visit_code_units_(const code_units * units,
                              struct member_visitor_ * visitor)
{
  visitor->fn_INT(visitor->context, "comoving_coordinates",
                  units->comoving_coordinates);
  visitor->fn_DOUBLE(visitor->context, "density_units", units->density_units);
  visitor->fn_DOUBLE(visitor->context, "length_units", units->length_units);
  visitor->fn_DOUBLE(visitor->context, "time_units", units->time_units);
  visitor->fn_DOUBLE(visitor->context, "velocity_units", units->velocity_units);
  visitor->fn_DOUBLE(visitor->context, "a_units", units->a_units);
  visitor->fn_DOUBLE(visitor->context, "a_value", units->a_value);
}

static void visit_field_data_(const grackle_field_data * my_fields,
                              struct member_visitor_ * v)
{
  // part 0: copy some data into fixed size arrays
  array_props_ field_prop = {my_fields->grid_rank, {0,0,0}};
  int grid_start[3] = {0, 0, 0};
  int grid_end[3] = {0, 0, 0};
  for (int i = 0; i < my_fields->grid_rank; i++){
    field_prop.dimensions[i] = my_fields->grid_dimension[i];
    grid_start[i] = my_fields->grid_start[i];
    grid_end[i] = my_fields->grid_end[i];
  }

  // part 1: dump the generic description about the field and grid
  {
    const array_props_ array_prop_3elem = {1, {3,0,0}};

    v->fn_INTARR(v->context, "grid_dimension", field_prop.dimensions,
                 array_prop_3elem);
    v->fn_INTARR(v->context, "grid_start", grid_start, array_prop_3elem);
    v->fn_INTARR(v->context, "grid_end", grid_end, array_prop_3elem);
    v->fn_INT(v->context, "grid_rank", my_fields->grid_rank);
    v->fn_DOUBLE(v->context, "grid_dx", my_fields->grid_dx);
  }

  // part 2: dump the field-data (we use X-Macros)
  #define ENTRY(MEMBER, _1)                                                \
    v->fn_GRFLOATARR(v->context, #MEMBER, my_fields->MEMBER, field_prop);
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY
}

// define some functionality to help write json objects
// ====================================================

struct json_obj_writer {
  FILE* fp;
  const char* member_spacing; // holds the whitespace between members
  int num_separations;
};

static void json_write_separation_(struct json_obj_writer* writer){
  if (writer->num_separations > 0) fputc(',', writer->fp);
  fprintf(writer->fp, "%s", writer->member_spacing);
  writer->num_separations++;
}

// Define helpers for the show_parameters function
static int json_field_INT(void* json_ptr, const char* field, int val)
{
  struct json_obj_writer* writer = (struct json_obj_writer*)json_ptr;
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %d", field, val);
  return SUCCESS;
}
static int json_field_DOUBLE(void* json_ptr, const char* field, double val)
{
  struct json_obj_writer* writer = (struct json_obj_writer*)json_ptr;
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %.17g", field, val);
  return SUCCESS;
}
static int json_field_STRING(void* json_ptr, const char* field, const char* val)
{
  struct json_obj_writer* writer = (struct json_obj_writer*)json_ptr;
  json_write_separation_(writer);
  if (val == NULL){
    fprintf(writer->fp, "\"%s\" : null", field);
  } else {
    fprintf(writer->fp, "\"%s\" : \"%s\"", field, val);
  }
  return SUCCESS;
}

static int json_field_INTARR(void* json_ptr, const char* field,
                             const int* val, array_props_ array_prop)
{
  fprintf(stderr, "Currently no support for json arrays. ABORTING NOW!\n");
  abort();
}

static int json_field_GRFLOATARR(void* json_ptr, const char* field,
                                 const gr_float* val, array_props_ array_prop)
{
  fprintf(stderr, "Currently no support for json arrays. ABORTING NOW!\n");
  abort();
}

struct member_visitor_ create_json_visitor_(FILE *fp) {
  // setup the json_obj_writer instance and do initial work
  struct json_obj_writer temporary = {fp, "\n  ", 0};
  fputc('{', fp); // intentionally omit '\n'

  // allocate a pointer to store in the member_visitor_ struct
  struct json_obj_writer* ptr = malloc(sizeof(struct json_obj_writer));
  *ptr = temporary;
  struct member_visitor_ out = {ptr, json_field_INT, json_field_DOUBLE,
                                json_field_STRING, json_field_INTARR,
                                json_field_GRFLOATARR};
  return out;
}

static void free_json_visitor(struct member_visitor_* visitor) {
  struct json_obj_writer* writer = (struct json_obj_writer*)(visitor->context);
  if (writer->num_separations == 0) {
    fprintf(writer->fp, "}\n");
  } else {
    fprintf(writer->fp, "\n}\n");
  }
  free(writer);
}

// functions to write json representations of various structs to disk
// ==================================================================

void show_parameters_(FILE *fp, const chemistry_data *my_chemistry){
  struct member_visitor_ visitor = create_json_visitor_(fp);
  visit_parameters_(my_chemistry, &visitor);
  free_json_visitor(&visitor); // end the json object
}

void show_version_(FILE *fp, const grackle_version* gversion)
{
  struct member_visitor_ visitor = create_json_visitor_(fp);
  visit_version_(gversion, &visitor);
  free_json_visitor(&visitor); // end the json object
}

// if we make the following function part of the stable API:
// - we should probably move the function to a different file (currently it's
//   here because it's implicitly required to properly dump the field_data)
// - we could theoretically start explicitly testing whether users specify the
//   correct fields

int grunstable_initialize_field_data(grackle_field_data *my_fields,
                                     int only_init_datafields)
{
  if (my_fields == NULL) {
    fprintf(stderr, "gr_initial_field_data was passed a NULL pointer\n");
    return FAIL;
  }

  // branching logic based on information encoded by only_init_datafields
  if (only_init_datafields == 0) { // initialize all members
    my_fields->grid_rank = -1;
    my_fields->grid_dimension = NULL;
    my_fields->grid_start = NULL;
    my_fields->grid_end = NULL;
    my_fields->grid_dx = -1.0;
  } else if (only_init_datafields != 1) {
    fprintf(stderr, "gr_initial_field_data received an invalid arg\n");
    return FAIL;
  }

  // now, modify all members holding datafields to have values of NULL
  // (we use X-Macros to do this)
  #define ENTRY(MEMBER_NAME, _1) my_fields->MEMBER_NAME = NULL;
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY

  return SUCCESS;
}

// implement generic hdf5-dumping machinery:
// =========================================

/// The idea is to track additional useful information inside this struct and
/// pass it around throughout our function stack...
///
/// To start out, we will just use it to track custom commonly used dtypes (so
/// that we don't need to constantly create new versions of this dtype)
typedef struct contextH5_ {
  hid_t gr_floattype;
  hid_t var_strtype;
} contextH5_;

static int cleanup_contextH5_(contextH5_* p) {
  if (p == NULL)  return FAIL;

  // ignore h_ctx->gr_floattype (we didn't make a new type here)

  // close h_ctx->var_strtype (if its in a valid state)
  if ((p->var_strtype != H5I_INVALID_HID) && (H5Tclose(p->var_strtype) < 0)) {
    fprintf(stderr, "error closing contextH5_->var_strtype\n");
    return FAIL;
  }
  p->var_strtype = H5I_INVALID_HID;

  return SUCCESS;
}

static int initialize_contextH5_(contextH5_* h_ctx) {
  if (h_ctx == NULL) {
    fprintf(stderr, "initialize_contextH5_ can't initialize a NULL ptr\n");
    return FAIL;
  }

  // Step 1: determine datatype equivalent to gr_float
  h_ctx->gr_floattype = (sizeof(gr_float) == 8) ?
    H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

  // Step 2: create a datatype to represent a variable length string
  h_ctx->var_strtype = H5Tcopy(H5T_C_S1);
  if (H5Tset_size(h_ctx->var_strtype, H5T_VARIABLE) < 0) {
    fprintf(stderr, "Error while initializing contextH5_ instance");
    return FAIL;
  }

  return SUCCESS;
}


/// write out a simple individual attribute
///
/// @param[in] loc_id location where the attribute will be attached to
/// @param[in] name the name of the attribute
/// @param[in] ptr Specifies the data to be written
/// @param[in] type_id the type of the attribute
/// @param[in] length The number of elements in the attribute (if it's a 1D
///     array). If writing a scalar, this should be 0
static int h5_write_attr_(hid_t loc_id, const char* name, const void* ptr,
                          hid_t type_id, int length) {
  if ((ptr == NULL) || (length < 0)) return FAIL;

  int ret_val = FAIL;
  hsize_t casted_length = (hsize_t)length;

  // create the dataspace
  hid_t space_id = (length == 0)
    ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &casted_length, NULL);
  if (space_id == H5I_INVALID_HID) goto fail_mkspace;

  // create the attribute
  hid_t attr_id = H5Acreate1(loc_id, name, type_id, space_id, H5P_DEFAULT);
  if (attr_id == H5I_INVALID_HID) goto fail_mkattr;

  // write the attribute
  if (H5Awrite(attr_id, type_id, ptr) >= 0) {
    ret_val = SUCCESS;
  }

  // cleanup time
  H5Aclose(attr_id);
fail_mkattr:
  H5Sclose(space_id);
fail_mkspace:
  // nothing to cleanup

  if (ret_val != SUCCESS){
    fprintf(stderr, "problem writing attribute: %s\n", name);
  }

  return ret_val;
}

/// Write a contiguous array as a simple dataset
///
/// @param[in] loc_id location where the attribute will be attached to
/// @param[in] name the name of the dataset
/// @param[in] ptr Specifies the data to be written
/// @param[in] type_id the elment-type of the array
/// @param[in] array_prop Specifies the expected dimensions of the data (you
///     don't need to do anything special when ptr is NULL)
static int h5dump_write_dset_(hid_t loc_id, const char* name, const void* ptr,
                              hid_t type_id, array_props_ array_prop)
{
  int ret_val = FAIL;

  // create dataspace
  hsize_t dims[3] = {0, 0, 0};
  for (int i = 0; i < array_prop.rank; i++) {
    dims[i] = array_prop.dimensions[i];
  }
  hid_t space_id = (ptr == NULL)
    ? H5Screate(H5S_NULL) // (in this case, it's an empty member)
    : H5Screate_simple(array_prop.rank, dims, NULL);
  if (space_id == H5I_INVALID_HID) goto fail_mkspace;

  // create dataset
  hid_t dset_id = H5Dcreate1(loc_id, name, type_id, space_id, H5P_DEFAULT);
  if (space_id == H5I_INVALID_HID) goto fail_mkdset;

  // write to the dataset
  if (H5Dwrite(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr) >= 0) {
    ret_val = SUCCESS;
  }

  // cleanup:
  H5Dclose(dset_id);
fail_mkdset:
  H5Sclose(space_id);
fail_mkspace:
  // nothing to cleanup

  return ret_val;
}

static const char * const H5_PROTOCOL_VERS_ATTR_ = "@:?GRACKLEDUMP_VERSION";
static const char * const H5_COMPLETION_ATTR_ = "@:?GRACKLEDUMP_DONE";

/// Closes the "annotated" HDF5 group
///
/// @note
/// See the docstring for h5dump_create_annotated_grp_ for more details
static int h5dump_close_annotated_grp_(hid_t grp_id, contextH5_* h_ctx,
                                       int indicate_errs)
{
  int ret_val = SUCCESS;

  if (H5Aexists(grp_id, H5_PROTOCOL_VERS_ATTR_) <= 0) {
    // includes cases where attribute does not exist or there was a failure
    fprintf(stderr, "ERROR: the \"%s\" attr does not seem to exist\n",
            H5_PROTOCOL_VERS_ATTR_);
    ret_val = FAIL;
  } else if (H5Aexists(grp_id, H5_COMPLETION_ATTR_) != 0) {
    // includes cases where attribute already exist or there was a failure
    fprintf(stderr,
            "ERROR: \"%s\" attr check failed or indicates it already exists\n",
            H5_PROTOCOL_VERS_ATTR_);
    ret_val = FAIL;
  }

  // when giving an indication that everything succeeded correctly, we are
  // going to try to flush the buffers before writing this attribute. Nothing
  // is guaranteed (we are at the mercy of the OS), but this decreases the
  // chance that an interuption to the program could leave us in a state, where
  // we indicate success and the group is only partially written to disk
  if ((ret_val == SUCCESS) && (indicate_errs == 0) && (H5Gflush(grp_id) < 0)) {
    fprintf(stderr, "ERROR: something went wrong during a group-flush\n");
    ret_val = FAIL;
  }

  // now, write out the H5_COMPLETION_ATTR_ attribute
  const int attr_val = ((ret_val == SUCCESS) && (indicate_errs == 0)) ? 1 : -1;

  if (h5_write_attr_(grp_id, H5_COMPLETION_ATTR_, &attr_val,
                     H5T_NATIVE_INT, 0) != SUCCESS) {
    fprintf(stderr, "ERROR while writing the \"%s\" attr\n",
            H5_COMPLETION_ATTR_);
    ret_val = FAIL;
  }

  // finally, let's actually close the group!
  H5Gclose (grp_id);
  return ret_val;
}

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
static hid_t h5dump_create_annotated_grp_(hid_t loc_id, contextH5_* h_ctx,
                                          const char* name)
{
  if (name == NULL) {
    fprintf(stderr, "error: the name of the desired group is a NULL ptr\n");
    return H5I_INVALID_HID;
  } else if (H5Lexists(loc_id, name, H5P_DEFAULT) != 0) {
    fprintf(stderr,
            "error: the check for whether the \"%s\" group exists "
            "failed or it indicates that the group already exists\n", name);
    return H5I_INVALID_HID;
  }

  // create the group
  hid_t grp_id = H5Gcreate2(loc_id, name, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
  if (grp_id == H5I_INVALID_HID){
    fprintf(stderr, "problem creating the \"%s\" group\n", name);
    return H5I_INVALID_HID;
  }

  // write the version number
  int version = 1;
  if (h5_write_attr_(grp_id, H5_PROTOCOL_VERS_ATTR_, &version,
                     H5T_NATIVE_INT, 0) != SUCCESS){
    fprintf(stderr, "problem creating the \"%s\" group\n", name);
    // close the group and denote an error
    h5dump_close_annotated_grp_(grp_id, h_ctx, 0);
    return H5I_INVALID_HID;
  }

  return grp_id;
}

// implement functions specific to dumping other structs to hdf5 file
// ==================================================================

struct visitor_h5plugin_ {hid_t loc_id; contextH5_* h_ctx; int err; };

static int visitorh5_INT(void* h5plugin, const char* field, int val)
{
  struct visitor_h5plugin_* p = (struct visitor_h5plugin_*)h5plugin;
  if (p->err == SUCCESS) {
    p->err = h5_write_attr_(p->loc_id, field, &val, H5T_NATIVE_INT, 0);
  }
  return p->err;
}

static int visitorh5_DOUBLE(void* h5plugin, const char* field, double val)
{
  struct visitor_h5plugin_* p = (struct visitor_h5plugin_*)h5plugin;
  if (p->err == SUCCESS) {
    p->err = h5_write_attr_(p->loc_id, field, &val, H5T_NATIVE_DOUBLE, 0);
  }
  return p->err;
}

static int visitorh5_STRING(void* h5plugin, const char* field, const char* val)
{
  struct visitor_h5plugin_* p = (struct visitor_h5plugin_*)h5plugin;
  if (p->err == SUCCESS) {
    p->err = h5_write_attr_(p->loc_id, field, &val, p->h_ctx->var_strtype, 0);
  }
  return p->err;
}

static int visitorh5_INTARR(void* h5plugin, const char* field,
                            const int* val, array_props_ array_prop)
{
  struct visitor_h5plugin_* p = (struct visitor_h5plugin_*)h5plugin;
  if (p->err == SUCCESS) {
    p->err = h5dump_write_dset_(p->loc_id, field, val, H5T_NATIVE_INT,
                                array_prop);
  }
  return p->err;
}

static int visitorh5_GRFLOATARR(void* h5plugin, const char* field,
                                const gr_float* val, array_props_ array_prop)
{
  struct visitor_h5plugin_* p = (struct visitor_h5plugin_*)h5plugin;
  if (p->err == SUCCESS) {
    p->err = h5dump_write_dset_(p->loc_id, field, val, p->h_ctx->gr_floattype,
                                array_prop);
  }
  return p->err;
}

struct member_visitor_ create_h5_visitor_(hid_t loc_id, contextH5_* h_ctx,
                                          const char* group_name) {
  // allocate and initialize the pluggin-pointer
  struct visitor_h5plugin_* ptr = malloc(sizeof(struct visitor_h5plugin_));
  {
    struct visitor_h5plugin_ tmp = {H5I_INVALID_HID, h_ctx, SUCCESS};
    *ptr = tmp;
  }

  // create the new group (that the visitor will write to)
  ptr->loc_id = h5dump_create_annotated_grp_(loc_id, h_ctx, group_name);
  if (ptr->loc_id == H5I_INVALID_HID) ptr->err = FAIL;

  // now package up member_visitor_ and return it!
  struct member_visitor_ out = {ptr, visitorh5_INT, visitorh5_DOUBLE,
                                visitorh5_STRING, visitorh5_INTARR,
                                visitorh5_GRFLOATARR};
  return out;
}

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
static int free_h5_visitor(struct member_visitor_* visitor,
                           const char* objname_errfmt) {
  // argument preprocessing
  struct visitor_h5plugin_* ptr = (struct visitor_h5plugin_*)visitor->context;
  const char* objname = (objname_errfmt == NULL) ? "<missing>" : objname_errfmt;

  // close the group (that the visitor previously wrote data to)
  int ret_val = ptr->err;
  int close_status = SUCCESS;
  if (ptr->loc_id != H5I_INVALID_HID) {
    close_status = h5dump_close_annotated_grp_(ptr->loc_id, ptr->h_ctx,
                                               ret_val != SUCCESS);
  }

  // error handling:
  if (ret_val != SUCCESS) {
    fprintf(stderr, "problem dumping contents of the %s object\n", objname);
  } else if (close_status != SUCCESS) {
    fprintf(stderr, "problem closing group for the %s object\n", objname);
    ret_val = FAIL;
  }

  free(ptr); // free memory

  return ret_val;
}

static int h5dump_chemistry_data_(hid_t loc_id, contextH5_* h_ctx,
                                  const char* name,
                                  const chemistry_data *chemistry_data) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (chemistry_data != NULL)  visit_parameters_(chemistry_data, &visitor);
  free_h5_visitor(&visitor, name);
  return SUCCESS;
}

static int h5dump_code_units_(hid_t loc_id, contextH5_* h_ctx,
                              const char* name, const code_units *units) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (units != NULL)  visit_code_units_(units, &visitor);
  free_h5_visitor(&visitor, name);
  return SUCCESS;
}

static int h5dump_version_(hid_t loc_id, contextH5_* h_ctx,
                           const char* name, const grackle_version *version) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (version != NULL)  visit_version_(version, &visitor);
  free_h5_visitor(&visitor, name);
  return SUCCESS;
}

static int h5dump_field_data_(hid_t loc_id, contextH5_* h_ctx,
                              const char* name,
                              const grackle_field_data *my_fields) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (my_fields != NULL)  visit_field_data_(my_fields, &visitor);
  free_h5_visitor(&visitor, name);
  return SUCCESS;
}

// implement the publicly exposed hdf5 dumper function
// ===================================================

/// performs some basic argument handling. Handles 3 cases:
///   1. Identify invalid combinations of fname and dest_hid
///   2. Open a new file (and specify that it must be closed)
///   3. Don't open anything (and specify that it should not be closed)
static int handle_fname_hid_args_(const char* fname, long long* dest_hid,
                                  int* require_file_close_hid) {

  // there's only 1 case where require_file_close_hid should have a value of 1
  // so we set it to 0 by default
  *require_file_close_hid = 0;
  
  if ((fname == NULL) && (*dest_hid != -1)) {
    // totally invalid state
    fprintf(stderr,
            "Can't handle case where fname and dest_hid are both provided\n");
    return FAIL;

  } else if (fname == NULL) {
    // in this case, the caller wants to specify a location within an existing
    // hdf5 file where data should be written

    H5I_type_t h5type = H5Iget_type(*dest_hid);
    if ((h5type != H5I_FILE) && (h5type != H5I_GROUP)){
      fprintf(stderr, "no fname specified & dest_hid doesn't specify a valid "
              "file or group\n");
      return FAIL;
    }
    return SUCCESS;

  } else {
    // try to open up a new hdf5 file at specified path 
    *dest_hid = (long long)(H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT,
                                      H5P_DEFAULT));
    if (*dest_hid == H5I_INVALID_HID) {
      fprintf(stderr, "unable to open an hdf5 file at %s\n", fname);
      return FAIL;
    }
    *require_file_close_hid = 1;
    return SUCCESS;
  }
}

int grunstable_h5dump_state(const char* fname, long long dest_hid,
                            const chemistry_data* my_chemistry,
                            const code_units* initial_code_units,
                            const code_units* current_code_units,
                            const grackle_field_data* my_fields)
{
  int ret_val = FAIL;
  int require_file_close_hid = 0;

  if (handle_fname_hid_args_(fname, &dest_hid, &require_file_close_hid)
      != SUCCESS) {
    goto fail_fname_hid_args;
  }

  // setup the contextH5_ instance
  contextH5_ h_ctx;
  if (initialize_contextH5_(&h_ctx) != SUCCESS)  goto fail_init_contextH5_;

  // create a group to contain everything
  const char* group_name = "grackle_statedump";
  hid_t grp_id = h5dump_create_annotated_grp_(dest_hid, &h_ctx, group_name);
  if (grp_id == H5I_INVALID_HID){
    fprintf(stderr, "problem creating the %s group\n", group_name);
    goto fail_create_grp;
  }

  // let's start dumping stuff!
  grackle_version vers = get_grackle_version();
  if (h5dump_version_(grp_id, &h_ctx, "grackle_version", &vers) != SUCCESS) {
    goto general_cleanup;
  } else if (h5dump_chemistry_data_(grp_id, &h_ctx, "chemistry_data",
                                    my_chemistry) != SUCCESS) {
    goto general_cleanup;
  } else if (h5dump_code_units_(grp_id, &h_ctx, "initial_code_units",
                                initial_code_units) != SUCCESS) {
    goto general_cleanup;
  } else if (h5dump_code_units_(grp_id, &h_ctx, "current_code_units",
                                current_code_units) != SUCCESS) {
    goto general_cleanup;
  } else if (h5dump_field_data_(grp_id, &h_ctx, "grackle_field_data",
                                my_fields) != SUCCESS) {
    goto general_cleanup;
  }

  ret_val = SUCCESS;

general_cleanup:
  if (h5dump_close_annotated_grp_(grp_id, &h_ctx,
                                  ret_val != SUCCESS) != SUCCESS) {
    fprintf(stderr, "problem closing the \"%s\" group\n", group_name);
    ret_val = FAIL;
  }
fail_create_grp:
  cleanup_contextH5_(&h_ctx);
fail_init_contextH5_:
fail_fname_hid_args:
  if (require_file_close_hid) H5Fclose(dest_hid);

  return ret_val;
}
