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


// define some functionality to help write json objects
// ====================================================

struct json_obj_writer {
  FILE* fp;
  const char* member_spacing; // holds the whitespace between members
  int num_separations;
};

static struct json_obj_writer json_create_writer_(FILE *fp){
  fputc('{', fp); // intentionally omit '\n'
  struct json_obj_writer out = {fp, "\n  ", 0};
  return out;
}

void json_write_separation_(struct json_obj_writer* writer){
  if (writer->num_separations > 0) fputc(',', writer->fp);
  fprintf(writer->fp, "%s", writer->member_spacing);
  writer->num_separations++;
}

void json_finish_(struct json_obj_writer* writer){
  if (writer->num_separations == 0) {
    fprintf(writer->fp, "}\n");
  } else {
    fprintf(writer->fp, "\n}\n");
  }
}

// Define helpers for the show_parameters function
static void json_field_INT(struct json_obj_writer* writer, const char* field,
                           int val)
{
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %d", field, val);
}
static void json_field_DOUBLE(struct json_obj_writer* writer,
                              const char* field, double val)
{
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %.17g", field, val);
}
static void json_field_STRING(struct json_obj_writer* writer,
                              const char* field, const char* val)
{
  json_write_separation_(writer);
  if (val == NULL){
    fprintf(writer->fp, "\"%s\" : null", field);
  } else {
    fprintf(writer->fp, "\"%s\" : \"%s\"", field, val);
  }
}

// functions to write json representations of various structs to disk
// ==================================================================

void show_parameters_(FILE *fp, const chemistry_data *my_chemistry){
  struct json_obj_writer writer = json_create_writer_(fp);

  // write all of the fields
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) \
    json_field_ ## TYPE (&writer, #FIELD, my_chemistry->FIELD);
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY

  // end the json object
  json_finish_(&writer);
}

void show_version_(FILE *fp, const grackle_version* gversion)
{
  struct json_obj_writer writer = json_create_writer_(fp);
  json_field_STRING(&writer, "version", gversion->version);
  json_field_STRING(&writer, "branch", gversion->branch);
  json_field_STRING(&writer, "revision", gversion->revision);
  json_finish_(&writer);
}

void show_code_units_(FILE *fp, const code_units* units)
{
  struct json_obj_writer writer = json_create_writer_(fp);
  json_field_INT(&writer, "comoving_coordinates", units->comoving_coordinates);
  json_field_DOUBLE(&writer, "density_units", units->density_units);
  json_field_DOUBLE(&writer, "length_units", units->length_units);
  json_field_DOUBLE(&writer, "time_units", units->time_units);
  json_field_DOUBLE(&writer, "velocity_units", units->velocity_units);
  json_field_DOUBLE(&writer, "a_units", units->a_units);
  json_field_DOUBLE(&writer, "a_value", units->a_value);
  json_finish_(&writer);
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

/// write out a simple individual attribute
///
/// @param[in] loc_id location where the attribute will be attached to
/// @param[in] name the name of the attribute
/// @param[in] type_id the type of the attribute
/// @param[in] length The number of elements in the attribute (if it's a 1D
///     array). If writing a scalar, this should be 0
/// @param[in] ptr Specifies the data to be written
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

// implement functions specific to dumping grackle_field_data to hdf5 file
// =======================================================================

/// Dumps subset of grackle_field_data's members into hdf5 attributes
static int h5dump_field_data_attributes_(hid_t loc_id,
                                         const grackle_field_data* my_fields)
{

  // for simplicity: always make grid_start and grid_end into 3D arrays
  int grid_start[3] = {0, 0, 0};
  int grid_end[3] = {0, 0, 0};
  for (int i = 0; i < my_fields->grid_rank; i++){
    grid_start[i] = my_fields->grid_start[i];
    grid_end[i] = my_fields->grid_end[i];
  }

  // here we actually write out the attributes
  if (h5_write_attr_(loc_id, "grid_start", grid_start,
                     H5T_NATIVE_INT, 3) != SUCCESS) {
    return FAIL;
  } else if (h5_write_attr_(loc_id, "grid_end", grid_end, H5T_NATIVE_INT, 3)
             != SUCCESS) {
    return FAIL;
  } else if (h5_write_attr_(loc_id, "grid_rank", &(my_fields->grid_rank),
                            H5T_NATIVE_INT, 0) != SUCCESS) {
    return FAIL;
  } else if (h5_write_attr_(loc_id, "grid_dx", &(my_fields->grid_dx),
                            H5T_NATIVE_DOUBLE, 0) != SUCCESS) {
    return FAIL;
  }

  return SUCCESS;
}

/// Write the field_data member (that holds field_data) to the container
static int h5dump_field_data_single_dset_(hid_t container,
                                          const grackle_field_data* my_fields,
                                          const char* field_name,
                                          const gr_float* values)
{
  int ret_val = FAIL;
  hid_t dtype = (sizeof(gr_float) == 4) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;

  // create dataspace
  hsize_t dims[3] = {0, 0, 0};
  for (int i = 0; i < my_fields->grid_rank; i++) {
    dims[i] = my_fields->grid_dimension[i];
  }
  hid_t space_id = (values == NULL)
    ? H5Screate(H5S_NULL) // (in this case, it's an empty member)
    : H5Screate_simple(my_fields->grid_rank, dims, NULL);
  if (space_id == H5I_INVALID_HID) goto fail_mkspace;

  // create dataset
  hid_t dset_id = H5Dcreate1(container, field_name, dtype, space_id,
                             H5P_DEFAULT);
  if (space_id == H5I_INVALID_HID) goto fail_mkdset;

  // write to the dataset
  if (H5Dwrite(dset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) >= 0) {
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

/// dump a grackle_field_data instance to loc_id
static int h5dump_field_data_(hid_t loc_id, const void* ptr)
{
  const grackle_field_data *my_fields = ptr;
  // arg-checking:
  if (my_fields == NULL) {
    fprintf(stderr, "gr_initial_field_data was passed a NULL pointer\n");
    return FAIL;
  }

  // write out some attributes:
  if (h5dump_field_data_attributes_(loc_id, my_fields) != SUCCESS) {
    return FAIL;
  }

  // now dump members of my_fields holding field-data (we use X-Macros)
  #define ENTRY(MEMBER_NAME, _1)                                              \
    if (h5dump_field_data_single_dset_(loc_id, my_fields, #MEMBER_NAME,       \
                                       my_fields->MEMBER_NAME) != SUCCESS) {  \
      return FAIL;                                                            \
    }
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY

  return SUCCESS;
}

// implement functions specific to dumping other structs to hdf5 file
// ==================================================================

// for now, these are just placeholders.
// - We are going to start with a simple strategy: writing the json
//   representation to a single attribute.
// - In the future, we could generalize the json writing machinery and
//   write individual key-value pairs as separate attributes... (but that may
//   not be worth the effort)

// to get a json-string representation, we will temporarily dump the data to a
// file and then load the data back into memory. This is a dumb strategy, but
// it's good enough for debugging purposes. Faster approaches might rely upon
// platform-specific c extensions (like open_memstream or fopencookie)

/// returns a newly allocated null-terminated string that contains a copy of
/// all of the characters in fp
///
/// @note
/// requires that fp was openned in "wb+" mode (this is consistent with using
/// an object produced by tmpfile(). To my knowledge, there is no way to
/// explicitly check this
static char* copy_f_contents_to_str_(FILE* fp) {
  // TODO: it would be nice to have some additional error handling

  fseek(fp, 0, SEEK_END); // advance file-position to end of the file
  const long length = ftell(fp);
  if (length == -1L)  return NULL;
  fseek(fp, 0, SEEK_SET); // advance file-position to start of the file

  const size_t str_len_with_nul = (size_t)length + 1;
  char* buf = malloc(sizeof(char) * str_len_with_nul);
  if ((buf != NULL) && (str_len_with_nul > 1)) {
    size_t num_chars = fread(buf, sizeof(char), str_len_with_nul, fp);
    if ((num_chars+1) != str_len_with_nul) {
      free(buf);
      buf = NULL;
    }
  }

  if (buf != NULL)  buf[str_len_with_nul - 1] = '\0';
  return buf;
}

/// this is just a helper function
/// - essentially the idea is that you call this after you've dumped a long
///   string to fp
/// - this function then reads that data into memory and then immediately
///   writes it to an hdf5 attribute
static int copy_fcontents_to_attr_(hid_t loc_id, const char* attr_name,
                                   FILE* fp) {
  char* str = copy_f_contents_to_str_(fp);
  if (str == NULL) {
    fprintf(stderr,
            "something went wrong while loading file buffer into memory (for "
            "the sake of writing the \"%s\" hdf5 attr)", attr_name);
    return FAIL;
  }

  // we are "cheating" here. We are writing the string as 1D array of
  // characters (it would probably be more robust to write it as
  // a variable-length string)

  int out = h5_write_attr_(loc_id, attr_name, str, H5T_NATIVE_CHAR,
                           strlen(str)+1);
  free(str);
  return out;
}


static int h5dump_chemistry_data_(hid_t loc_id, const void* ptr){
  const chemistry_data *chemistry_data = ptr;
  // dump json representation to tmpfile
  FILE* fp = tmpfile();
  if (fp == NULL) return FAIL;
  show_parameters_(fp, chemistry_data);

  // copy json representation to hdf5 attribute
  int out = copy_fcontents_to_attr_(loc_id, "json_str", fp);
  fclose(fp);
  return SUCCESS;
}

static int h5dump_code_units_(hid_t loc_id, const void* ptr){
  const code_units *units = ptr;
  FILE* fp = tmpfile();
  if (fp == NULL) return FAIL;
  show_code_units_(fp, units);

  // copy json representation to hdf5 attribute
  int out = copy_fcontents_to_attr_(loc_id, "json_str", fp);
  fclose(fp);
  return SUCCESS;
}

static int h5dump_version_(hid_t loc_id, const void* ptr){
  const grackle_version *version = ptr;

  // dump json representation to tmpfile
  FILE* fp = tmpfile();
  if (fp == NULL) return FAIL;
  show_version_(fp, version);

  // copy json representation to hdf5 attribute
  int out = copy_fcontents_to_attr_(loc_id, "json_str", fp);
  fclose(fp);
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

struct state_dump_entry_{
  const char* group_name;             ///< name of the group to create
  const void* obj;                    ///< object to actually dump
  int (*fn_ptr)(hid_t, const void*);  ///< function for dumping state
};

int grunstable_h5dump_state(const char* fname, long long dest_hid,
                            const chemistry_data* my_chemistry,
                            const code_units* initial_code_units,
                            const code_units* current_code_units,
                            const grackle_field_data* my_fields)
{
  int ret_val = SUCCESS;
  int require_file_close_hid = 0;

  if (handle_fname_hid_args_(fname, &dest_hid, &require_file_close_hid)
      != SUCCESS) {
    ret_val = FAIL;
    goto general_cleanup;
  }

  // TODO: create an attribute called "grackle_dump_finished" with a value of 0

  grackle_version version = get_grackle_version();

  // Do not make entry_l as a static variable (otherwise, that requires that
  // the expressions in the initializer list are constant-expressions)
  const struct state_dump_entry_ entry_l[] =
    {
     {"grackle_version", &version, h5dump_version_},
     {"chemistry_data", my_chemistry, h5dump_chemistry_data_},
     {"initial_code_units", initial_code_units, h5dump_code_units_},
     {"current_code_units", current_code_units, h5dump_code_units_},
     {"grackle_field_data", my_fields, h5dump_field_data_}
    };

  int num_entries = (int)(sizeof(entry_l) / sizeof(struct state_dump_entry_));
  for (int i = 0; i < num_entries; i++) {
    if (entry_l[i].obj == NULL)  continue;

    // create a group
    hid_t group_id = H5Gcreate2(dest_hid, entry_l[i].group_name,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id == H5I_INVALID_HID){
      fprintf(stderr, "problem creating the %s group\n", entry_l[i].group_name);
      ret_val = FAIL;
      break;
    }

    int cur_dump_result = entry_l[i].fn_ptr(group_id, entry_l[i].obj);
    H5Gclose (group_id);
    if (cur_dump_result != SUCCESS) {
      fprintf(stderr, "problem dumping contents in the %s group\n",
              entry_l[i].group_name);
      ret_val = FAIL;
      break;
    }
  }

  if (ret_val == SUCCESS) {
    // overwrite the "grackle_dump_finished" attribute with a value of 1
  }

general_cleanup:
  if (require_file_close_hid) H5Fclose(dest_hid);

  return ret_val;
}
