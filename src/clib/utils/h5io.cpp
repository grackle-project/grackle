//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements functions to help read from hdf5 files
///
//===----------------------------------------------------------------------===//
#include <climits>
#include <cstdio>
#include <cstring>

#include "hdf5.h"
#include "grackle.h"
#include "grackle_macros.h"
#include "h5io.hpp"
#include "../status_reporting.h"

int grackle::impl::h5io::read_str_attribute(hid_t attr_id, int bufsz,
                                            char* buffer) {
  if (bufsz < 0 || attr_id == H5I_INVALID_HID) {
    return -1;
  }

  // I'm a little suspicious of this line of code...
  // -> if we ever see Parameter names get clipped, my guess is that this line
  //    doesn't handle variable length strings properly
  // -> with that said, I *think* this is correct
  hsize_t storage_size = H5Aget_storage_size(attr_id);
  if (storage_size > static_cast<hsize_t>(INT_MAX - 1)) {
    std::fprintf(
        stderr,
        "the string attribute requires a bigger buffer than expected\n");
    return -1;
  }
  int required_bufsz = static_cast<int>(storage_size) + 1;

  if (bufsz == 0) {
    return required_bufsz;
  } else if (required_bufsz > bufsz) {
    return -1;
  }

  // we inspect a few properties about the datatype on disk
  // - its unfortunate that we need to worry about utf8 encoding... This is
  //   relevant for the _high_density table (I suspect that it may be a default
  //   within h5py). We'll validate that the characters are compatible with
  //   ASCII before we return
  bool is_variable, uses_utf8_encoding;
  {
    hid_t disk_typeid = H5Aget_type(attr_id);
    if (disk_typeid == H5I_INVALID_HID) {
      std::fprintf(stderr, "Failed to get type of attribute.\n");
      return -1;
    } else if (H5Tdetect_class(disk_typeid, H5T_STRING) <= 0) {
      std::fprintf(stderr, "Expected to be a string.\n");
      H5Tclose(disk_typeid);
      return -1;
    }

    htri_t is_variable_rslt = H5Tis_variable_str(disk_typeid);
    if (is_variable_rslt < 0) {
      std::fprintf(stderr, "Error in H5Tis_variable_str.\n");
      H5Tclose(disk_typeid);
      return -1;
    }
    is_variable = (is_variable_rslt > 0);

    H5T_cset_t cset = H5Tget_cset(disk_typeid);
    if (cset == H5T_CSET_ERROR) {
      std::fprintf(stderr, "Error in H5Tget_cset.\n");
      H5Tclose(disk_typeid);
      return -1;
    }
    uses_utf8_encoding = (cset == H5T_CSET_UTF8);

    H5Tclose(disk_typeid);
  }

  // construct the disk_typeid that we use for describing `buffer` to the hdf5
  // library
  // - it might be nice if we actually cached this type and carried it around
  //   so that we can avoid remaking the type every time we load a string
  hid_t memory_typeid = H5Tcopy(H5T_C_S1);
  if (H5Tset_size(memory_typeid, is_variable ? H5T_VARIABLE : bufsz) < 0) {
    fprintf(stderr, "error in H5Tset_size (for memory datatype)\n");
    H5Tclose(memory_typeid);
    return -1;
  }
  if (uses_utf8_encoding && (H5Tset_cset(memory_typeid, H5T_CSET_UTF8) < 0)) {
    fprintf(stderr, "error in H5Tset_size (for memory datatype)\n");
    H5Tclose(memory_typeid);
    return -1;
  }
  if (H5Tset_strpad(memory_typeid, H5T_STR_NULLTERM) < 0) {
    fprintf(stderr, "error in H5Tset_strpad (for memory datatype)\n");
    H5Tclose(memory_typeid);
    return -1;
  }
  if (is_variable) {
    char* tmp_str = nullptr;
    if (H5Aread(attr_id, memory_typeid, static_cast<void*>(&tmp_str)) < 0) {
      // is there a memory leak?
      fprintf(stderr, "error in H5Aread for a variable length string\n");
      H5Tclose(memory_typeid);
      return -1;
    }
    std::strncpy(buffer, tmp_str, static_cast<std::size_t>(required_bufsz));
    H5free_memory(tmp_str);
  } else {
    if (H5Aread(attr_id, memory_typeid, buffer) < 0) {
      H5Tclose(memory_typeid);
      fprintf(stderr, "error in H5Aread\n");
      return -1;
    }
  }
  H5Tclose(memory_typeid);

  if (uses_utf8_encoding) {
    // if the string used utf8 encoding on disk, we need to confirm that the
    // characters are all compatible with ASCII characterset.
    //
    // For the uninitiated
    // - for the uninitiated, all 128 standard ASCII characters are encoded in
    //   7bits and have an identical representation in utf8-encoding.
    // - all other utf8 codepoints are represented by 2 or more bytes &
    //   the leading byte has a value doesn't correspond to an ASCII character
    const char max_val = static_cast<char>(127);
    for (int i = 0; i < bufsz; i++) {
      if (buffer[i] == '\0') {
        break;
      } else if (buffer[i] < '\0' || buffer[i] > max_val) {
        // I'm pretty sure we must verfiy that buffer[i] BOTH isn't smaller
        // than '\0' AND doesn't exceed max_val since the standard doesn't
        // specify whether char is signed or unsigned
        std::fprintf(stderr,
                     "Error: read a utf8-encoded string from an attribute that "
                     "contains non-ASCII characters\n");
        return -1;
      }
    }
  }

  return required_bufsz;
}

namespace {  // stuff inside an anonymous namespace is local to this file

grackle::impl::h5io::ArrayShape shape_from_space(hid_t space_id) {
  grackle::impl::h5io::ArrayShape out;
  out.ndim = -2;  // set up an output value that denotes an error
  if (space_id == H5I_INVALID_HID) {
    return out;
  }

  H5S_class_t space_type = H5Sget_simple_extent_type(space_id);
  switch (space_type) {
    case H5S_NULL: {
      out.ndim = -1;
      return out;
    }
    case H5S_SCALAR: {
      out.ndim = 0;
      return out;
    }
    case H5S_SIMPLE: {
      int ndim = H5Sget_simple_extent_ndims(space_id);
      hsize_t tmp[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
      if (ndim <= 0) {  // unclear what ndim of 0 would mean
        std::fprintf(stderr, "error in call to H5Sget_simple_extent_ndims\n");
        return out;
      } else if (ndim > GRACKLE_CLOUDY_TABLE_MAX_DIMENSION) {
        std::fprintf(
            stderr,
            "hdf5 array has %d dimensions, exceeding the hardcoded maximum\n",
            ndim);
        return out;
      } else if (ndim != H5Sget_simple_extent_dims(space_id, tmp, nullptr)) {
        std::fprintf(stderr, "error in call to H5Sget_simple_extent_dims\n");
        return out;
      }

      for (int i = 0; i < ndim; i++) {
        if (tmp[i] > INT64_MAX) {
          std::fprintf(stderr, "shape[%d] is too big\n", i);
          return out;
        }
        out.shape[i] = static_cast<std::int64_t>(tmp[i]);
      }
      out.ndim = ndim;  // make sure to do this last!
      return out;
    }
    default: {
      std::fprintf(stderr, "unexpected dataspace issue\n");
      return out;
    }
  }
}

grackle::impl::h5io::ArrayShape mk_invalid_array_shape() {
  return shape_from_space(H5I_INVALID_HID);
}

}  // anonymous namespace

grackle::impl::h5io::ArrayShape grackle::impl::h5io::read_dataset_shape(
    hid_t file_id, const char* name) {
  if (name == nullptr) {
    std::fprintf(stderr, "name is a nullptr");
    return mk_invalid_array_shape();
  }

  hid_t dset_id = H5Dopen(file_id, name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Failed to open dataset \"%s\".\n", name);
    return mk_invalid_array_shape();
  }

  hid_t space_id = H5Dget_space(dset_id);
  H5Dclose(dset_id);

  ArrayShape out = shape_from_space(space_id);
  H5Sclose(space_id);

  if (ArrayShape_is_null(out)) {
    std::fprintf(stderr, "Error reading the dataspace of \"%s\".\n", name);
  }
  return out;
}

/// read the dataset named dset_name from file_id into buffer
int grackle::impl::h5io::read_dataset(
    hid_t file_id, const char* dset_name, double* buffer,
    const grackle::impl::h5io::ArrayShape* expected_shape) {
  if (dset_name == nullptr) {
    std::fprintf(stderr, "dset_name is a nullptr");
    return GR_FAIL;
  }

  hid_t dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Failed to open dataset \"%s\".\n", dset_name);
    return GR_FAIL;
  }

  // we may perform an extra quick check about the dataset's shape
  if (expected_shape != nullptr) {
    hid_t space_id = H5Dget_space(dset_id);
    ArrayShape actual_shape = shape_from_space(space_id);
    H5Sclose(space_id);

    if (!ArrayShape_is_same(*expected_shape, actual_shape)) {
      H5Dclose(dset_id);
      std::fprintf(stderr, "The \"%s\" dataset has an unexpected shape.\n",
                   dset_name);
      return GR_FAIL;
    }
  }

  herr_t status =
      H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  H5Dclose(dset_id);
  if (status < 0) {
    std::fprintf(stderr, "Failed to read dataset \"%s\".\n", dset_name);
    return GR_FAIL;
  }

  return GR_SUCCESS;
}

namespace {  // stuff inside an anonymous namespace is local to this file

static constexpr int attr_name_set_capacity = 30;

/// Acts as a crude set of attribute names
///
/// @note
/// if we ever choose to embrace C++, it would be **MUCH** nicer to use
/// std::set<std::string>
struct AttrNameSet {
  char* vals[attr_name_set_capacity];
  int capacity;
  int length;
};

/// helper function that calculates ArrayShape from the "Rank" and "Dimension"
/// attributes
///
/// @note
/// we only carry around dset_name argument for nicer error messages
grackle::impl::h5io::ArrayShape shape_from_grid_attrs(hid_t dset_id,
                                                      const char* dset_name) {
  hid_t attr_id = H5Aopen_name(dset_id, "Rank");
  if (attr_id == H5I_INVALID_HID) {
    std::fprintf(stderr,
                 "Failed to open \"Rank\" attribute of \"%s\" dataset.\n",
                 dset_name);
    return mk_invalid_array_shape();
  }
  long long rank;
  if (H5Aread(attr_id, HDF5_I8, &rank) < 0) {
    H5Aclose(attr_id);
    std::fprintf(stderr,
                 "Failed to read \"Rank\" attribute of \"%s\" dataset.\n",
                 dset_name);
    return mk_invalid_array_shape();
  }
  H5Aclose(attr_id);
  if ((rank < 0) || (GRACKLE_CLOUDY_TABLE_MAX_DIMENSION < rank)) {
    H5Dclose(dset_id);
    std::fprintf(
        stderr,
        "\"Rank\" attribute of \"%s\" dataset, %lld, is negative or exceeds "
        "the hardcoded maximum, %d.\n",
        dset_name, rank, GRACKLE_CLOUDY_TABLE_MAX_DIMENSION);
    return mk_invalid_array_shape();
  }

  // read in the Dimension attribute
  std::int64_t grid_dimensions[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    grid_dimensions[i] = -1;
  }
  attr_id = H5Aopen_name(dset_id, "Dimension");
  // todo: it would be nice to validate that the shape of dimension matches
  //       our expectation
  if (attr_id == H5I_INVALID_HID) {
    std::fprintf(stderr,
                 "Failed to open \"Dimension\" attribute in \"%s\" dataset.\n",
                 dset_name);
    return mk_invalid_array_shape();
  }
  if (H5Aread(attr_id, HDF5_I8, grid_dimensions) < 0) {
    std::fprintf(stderr,
                 "Failed to read \"Dimension\" attribute in \"%s\" dataset.\n",
                 dset_name);
    return mk_invalid_array_shape();
  }
  H5Aclose(attr_id);
  // a quick check
  for (int i = 0; i < rank; i++) {
    if (grid_dimensions[i] < 1) {
      std::fprintf(
          stderr,
          "\"Dimension\" attribute of \"%s\" dataset has non-positive value\n",
          dset_name);
      return mk_invalid_array_shape();
    }
  }

  // finally, let's format the output
  grackle::impl::h5io::ArrayShape out;
  out.ndim = static_cast<int>(rank);
  for (int i = 0; i < out.ndim; i++) {
    out.shape[i] = grid_dimensions[i];
  }
  return out;
}

static constexpr int max_attr_name_length = 64;

/// updates axes with parsed parameters
///
/// returns the number of uniqued attributes accessed by this function
int set_grid_axes_props(hid_t dset_id, const char* dset_name,
                        grackle::impl::h5io::ArrayShape grid_shape,
                        grackle::impl::h5io::GridTableAxis* axes) {
  using grackle::impl::h5io::read_str_attribute;
  int accessed_attr_count = 0;

  for (int i = 0; i < grid_shape.ndim; i++) {
    bool is_last_axis = (i + 1) == grid_shape.ndim;

    // Step 1: get the name of the attribute holding the values along axis i
    //         AND fill axes[i].name with the name of the varying quantity
    char val_attr_name[max_attr_name_length];
    std::snprintf(val_attr_name, max_attr_name_length, "Parameter%d", i + 1);

    if (H5Aexists(dset_id, val_attr_name) > 0) {
      char tmp_attr_name[max_attr_name_length];
      std::snprintf(tmp_attr_name, max_attr_name_length, "Parameter%d_Name",
                    i + 1);
      hid_t attr_id = H5Aopen_name(dset_id, tmp_attr_name);

      if (attr_id == H5I_INVALID_HID) {
        std::fprintf(stderr,
                     "Failed to open \"%s\" attribute of \"%s\" dataset.\n",
                     tmp_attr_name, dset_name);
        return -1;
      }

      int min_buf_length = read_str_attribute(attr_id, 0, nullptr);
      if (min_buf_length < 0) {
        std::fprintf(
            stderr,
            "can't determine string size held by \"%s\" attr of \"%s\" dset.\n",
            tmp_attr_name, dset_name);
        return -1;
      }

      axes[i].name = new char[min_buf_length];
      if (read_str_attribute(attr_id, min_buf_length, axes[i].name) < 0) {
        std::fprintf(stderr,
                     "error loading the \"%s\" attr from \"%s\" dataset\n",
                     tmp_attr_name, dset_name);
        return -1;
      }
      H5Aclose(attr_id);
      accessed_attr_count++;

    } else if (is_last_axis && H5Aexists(dset_id, "Temperature") > 0) {
      axes[i].name = new char[12];
      std::snprintf(val_attr_name, max_attr_name_length, "Temperature");
      std::snprintf(axes[i].name, 12, "Temperature");

    } else {
      const char* extra_detail =
          (is_last_axis) ? "" : " Neither is \"Temperature\".";
      std::fprintf(
          stderr,
          "Failed to determine attribute associated with axis %d of the "
          "\"%s\" dataset. \"%s\" is not a known attribute.%s\n",
          i, val_attr_name, dset_name, extra_detail);
      return -1;
    }

    // step 2: load the values associated with axes[i].name and store them in
    //         axes[i].values

    hid_t attr_id = H5Aopen_name(dset_id, val_attr_name);
    if (attr_id == H5I_INVALID_HID) {
      std::fprintf(stderr, "Failed to open \"%s\" attr of \"%s\" dataset.\n",
                   val_attr_name, dset_name);
      return -1;
    }

    hid_t space_id = H5Aget_space(attr_id);
    grackle::impl::h5io::ArrayShape axis_shape = shape_from_space(space_id);
    H5Sclose(space_id);
    if (axis_shape.ndim != 1 || axis_shape.shape[0] != grid_shape.shape[i]) {
      std::fprintf(
          stderr,
          "The \"%s\" attr of the \"%s\" dataset isn't an array with a shape "
          "that is consistent with the dataset's \"Dimension\" attr.\n",
          val_attr_name, dset_name);
      return -1;
    }

    axes[i].values = new double[grid_shape.shape[i]];
    if (H5Aread(attr_id, HDF5_R8, axes[i].values) < 0) {
      H5Aclose(attr_id);
      std::fprintf(stderr, "failed to read \"%s\" attr of \"%s\" dataset.\n",
                   val_attr_name, dset_name);
      return -1;
    }
    H5Aclose(attr_id);
    accessed_attr_count++;
  }
  return accessed_attr_count;
}

int get_num_attrs(hid_t dset_id, const char* dset_name) {
#if H5_VERSION_LE(1, 10, 2)
  return -2;
#else
#if H5_VERSION_GE(1, 12, 0)
  H5O_info2_t info;
  herr_t status = H5Oget_info3(dset_id, &info, H5O_INFO_NUM_ATTRS);
#else
  H5O_info_t info;
  herr_t status = H5Oget_info2(dset_id, &info, H5O_INFO_NUM_ATTRS);
#endif
  if (status < 0) {
    std::fprintf(stderr, "Can't get num_attrs for \"%s\" dataset.\n",
                 dset_name);
    return -1;
  }
  if (info.num_attrs > INT_MAX) {
    std::fprintf(stderr, "\"%s\" dataset has too many attrs.\n", dset_name);
    return -1;
  } else {
    return static_cast<int>(info.num_attrs);
  }
#endif
}

}  // anonymous namespace

grackle::impl::h5io::GridTableProps grackle::impl::h5io::parse_GridTableProps(
    hid_t file_id, const char* dset_name) {
  // setup the output object so that we're always prepared to return an object
  // - the default object is constructed to denote a failure
  grackle::impl::h5io::GridTableProps out;
  out.table_shape = shape_from_space(H5I_INVALID_HID);
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    out.axes[i].name = nullptr;
    out.axes[i].values = nullptr;
  }

  if (dset_name == nullptr) {  // sanity check!
    std::fprintf(stderr, "dset_name is a nullptr");
    return out;
  }

  hid_t dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Can't open \"%s\" dataset.\n", dset_name);
    return out;
  }

  // infer the shape of the table from the attributes
  grackle::impl::h5io::ArrayShape inferred_shape =
      shape_from_grid_attrs(dset_id, dset_name);
  if (!ArrayShape_is_valid(inferred_shape)) {
    H5Dclose(dset_id);
    // shape_from_grid_attrs already printed error messages in this case
    return out;
  }
  int shape_attr_count = 2;

  // parse the quantities along each axis
  int axes_prop_attr_count =
      set_grid_axes_props(dset_id, dset_name, inferred_shape, out.axes);
  if (axes_prop_attr_count < 0) {
    H5Dclose(dset_id);
    drop_GridTableProps(&out);
    // shape_from_grid_attrs already printed error messages in this case
    return out;
  }

  int total_accessed_attrs_count = shape_attr_count + axes_prop_attr_count;

  // perform a check to confirm that we have accessed every available attribute
  // -> we may want to disable this check...
  int num_attrs = get_num_attrs(dset_id, dset_name);
  if (num_attrs == -2) {
    num_attrs = total_accessed_attrs_count;
  } else if (num_attrs == -1) {
    H5Dclose(dset_id);
    drop_GridTableProps(&out);
    // get_num_attrs already printed error messages in this case
    return out;
  }

  H5Dclose(dset_id);

  if (num_attrs != total_accessed_attrs_count) {
    drop_GridTableProps(&out);
    std::fprintf(
        stderr,
        "the number of accessed attributes, %d, for the \"%s\" dataset isn't "
        "equal to the number of available attributes, %d. This may be a sign "
        "of a problem\n",
        total_accessed_attrs_count, dset_name, num_attrs);
    return out;
  }

  out.table_shape = inferred_shape;  // we intentionally do this as the very
                                     // last step!
  return out;
}
