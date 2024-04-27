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
#include <hdf5.h>

#include "grackle.h"
#include "grackle_macros.h"

#include "grackle_unstable.h" // forward declare grunstable_h5dump_state
#include "state_dump_utils.h" // more forward-declarations
#include "visitor_utils.h"

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

// functions to write json representations of various structs to disk
// ==================================================================

void show_parameters_(FILE *fp, const chemistry_data *my_chemistry){
  struct member_visitor_ visitor = create_json_visitor_(fp);
  visit_parameters_(my_chemistry, &visitor);
  free_json_visitor_(&visitor); // end the json object
}

void show_version_(FILE *fp, const grackle_version* gversion)
{
  struct member_visitor_ visitor = create_json_visitor_(fp);
  visit_version_(gversion, &visitor);
  free_json_visitor_(&visitor); // end the json object
}

// functions to write hdf5 representations of various structs to disk
// ==================================================================

static int h5dump_chemistry_data_(hid_t loc_id, contextH5_* h_ctx,
                                  const char* name,
                                  const chemistry_data *chemistry_data) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (chemistry_data != NULL)  visit_parameters_(chemistry_data, &visitor);
  free_h5_visitor_(&visitor, name);
  return SUCCESS;
}

static int h5dump_code_units_(hid_t loc_id, contextH5_* h_ctx,
                              const char* name, const code_units *units) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (units != NULL)  visit_code_units_(units, &visitor);
  free_h5_visitor_(&visitor, name);
  return SUCCESS;
}

static int h5dump_version_(hid_t loc_id, contextH5_* h_ctx,
                           const char* name, const grackle_version *version) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (version != NULL)  visit_version_(version, &visitor);
  free_h5_visitor_(&visitor, name);
  return SUCCESS;
}

static int h5dump_field_data_(hid_t loc_id, contextH5_* h_ctx,
                              const char* name,
                              const grackle_field_data *my_fields) {
  struct member_visitor_ visitor = create_h5_visitor_(loc_id, h_ctx, name);
  if (my_fields != NULL)  visit_field_data_(my_fields, &visitor);
  free_h5_visitor_(&visitor, name);
  return SUCCESS;
}

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
