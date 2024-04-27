/***********************************************************************
/
/ Implement vistor utilities
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <stdlib.h> // malloc

#include "grackle_macros.h" // SUCCESS, FAIL
#include "visitor_utils.h"

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

void free_json_visitor_(struct member_visitor_* visitor) {
  struct json_obj_writer* writer = (struct json_obj_writer*)(visitor->context);
  if (writer->num_separations == 0) {
    fprintf(writer->fp, "}\n");
  } else {
    fprintf(writer->fp, "\n}\n");
  }
  free(writer);
}

// implement generic hdf5-dumping machinery:
// =========================================

int initialize_contextH5_(contextH5_* h_ctx) {
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

int cleanup_contextH5_(contextH5_* p) {
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

hid_t h5dump_create_annotated_grp_(hid_t loc_id, contextH5_* h_ctx,
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

int h5dump_close_annotated_grp_(hid_t grp_id, contextH5_* h_ctx,
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
int free_h5_visitor_(struct member_visitor_* visitor,
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
