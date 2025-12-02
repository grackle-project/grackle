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
/// The error-handling in these functions leaves a lot to be desired...
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

// HDF5 definitions
// -> these were moved out of the grackle_macros.h header
// -> most of these are unnecessary!

#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE

#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_I8 H5T_NATIVE_LLONG
#define HDF5_R4 H5T_NATIVE_FLOAT
#define HDF5_R8 H5T_NATIVE_DOUBLE

namespace {  // stuff inside an anonymous namespace is local to this file

/// Determines whether the specified buffer contains a null-terminated byte
/// string of ascii characters.
///
/// @returns true if buffer only contains ascii characters and contains a null
///   character at buffer[i] for `i < bufsz`. Otherwise it returns true.
///
/// @par Purpose
/// The function is primarily intended for checking whether a UTF-8 encoded
/// string is also a valid ASCII string (i.e. it only includes UTF-8 code
/// points that are valid ASCII characters). For the uninitiated:
/// - all 128 standard ASCII characters are encoded in 7 bits. UTF-8 was
///   designed for backwards compatibility with these characters; the first
///   128 UTF-8 code points are all encoded using a single byte and have an
///   **EXACT** one-to-one correspondence with the ascii characters.
/// - all other utf8 codepoints are represented by 2 or more bytes &
///   the leading byte has a value doesn't correspond to an ASCII character
bool is_ascii_string(const char* buffer, int bufsz) {
  const char max_val = static_cast<char>(127);
  for (int i = 0; i < bufsz; i++) {
    if (buffer[i] == '\0') {
      return true;
    } else if (buffer[i] < '\0' || buffer[i] > max_val) {
      // I'm pretty sure we must verfiy that buffer[i] BOTH isn't smaller
      // than '\0' AND doesn't exceed max_val since the standard doesn't
      // specify whether char is signed or unsigned
      return false;
    }
  }

  return false;  // (buffer doesn't contain a null character)
}

/// does the heavy lifting for read_str_attribute and read_str_dataset
int read_str_data_helper_(hid_t id, bool is_attr, int bufsz, char* buffer) {
  if (bufsz < 0 || id == H5I_INVALID_HID) {
    return -1;
  }

  // I'm a little suspicious of this line of code...
  // -> if we ever see Parameter names get clipped, my guess is that this line
  //    doesn't handle variable length strings properly
  // -> with that said, I *think* this is correct
  hsize_t storage_size =
      (is_attr) ? H5Aget_storage_size(id) : H5Dget_storage_size(id);
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
    hid_t disk_typeid = (is_attr) ? H5Aget_type(id) : H5Dget_type(id);
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

  // actually read the data
  if (is_variable) {
    char* tmp_str = nullptr;
    herr_t status =
        (is_attr) ? H5Aread(id, memory_typeid, static_cast<void*>(&tmp_str))
                  : H5Dread(id, memory_typeid, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            static_cast<void*>(&tmp_str));
    if (status < 0) {
      H5Tclose(memory_typeid);
      if (tmp_str != nullptr) {
        H5free_memory(tmp_str);
      }
      fprintf(stderr, "error in %s for a variable length string\n",
              (is_attr) ? "H5Aread" : "H5Dread");
      return -1;
    }
    std::strncpy(buffer, tmp_str, static_cast<std::size_t>(required_bufsz));
  } else {
    herr_t status = (is_attr) ? H5Aread(id, memory_typeid, buffer)
                              : H5Dread(id, memory_typeid, H5S_ALL, H5S_ALL,
                                        H5P_DEFAULT, buffer);
    if (status < 0) {
      H5Tclose(memory_typeid);
      fprintf(stderr, "error in %s for a fixed length string\n",
              (is_attr) ? "H5Aread" : "H5Dread");
      return -1;
    }
  }
  H5Tclose(memory_typeid);

  if (uses_utf8_encoding) {
    // if the string was stored on the disk with a utf-8 encoding, we need to
    // confirm that the string *only* includes the subset of utf-8 code-points
    // that are also valid characters are all compatible with ASCII.
    if (!is_ascii_string(buffer, bufsz)) {
      std::fprintf(stderr,
                   "Error: read a UTF-8 encoded string from an attribute that "
                   "contains non-ASCII code points\n");
      return -1;
    }
  }

  return required_bufsz;
}

}  // anonymous namespace

int grackle::impl::h5io::read_str_attribute(hid_t attr_id, int bufsz,
                                            char* buffer) {
  return read_str_data_helper_(attr_id, true, bufsz, buffer);
}

int grackle::impl::h5io::read_str_dataset(hid_t file_id, const char* dset_name,
                                          int bufsz, char* buffer) {
  if (dset_name == nullptr) {
    std::fprintf(stderr, "dset_name is a nullptr");
    return GR_FAIL;
  }

  hid_t dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Failed to open dataset \"%s\".\n", dset_name);
    return GR_FAIL;
  }

  int out = read_str_data_helper_(dset_id, false, bufsz, buffer);
  H5Dclose(dset_id);
  return out;
}

namespace {  // stuff inside an anonymous namespace is local to this file

/// Construct an ArrayShape instance from a HDF5 dataspace
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

/// Construct an ArrayShape object that represents an invalid instance
///
/// @note
/// This is useful for encoding that a function produced an error
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

/// @defgroup AttrNameRecorderGroup Attribute Recording Group
/// This is a collection of crude machinery used when we parse hdf5 dataset
/// attributes that specify Grid Table properties. This machinery is designed
/// to operate in 2 modes:
/// 1. A lightweight counter of attribute names
/// 2. A heavier weight mode that actually tracks alls attribute names. The
///    premise is to only use this mode (when we already know that there is an
///    error) in order to provide a descriptive error message.
///@{
static constexpr int attr_name_recorder_capacity =
    (GRACKLE_CLOUDY_TABLE_MAX_DIMENSION * 3 + 4);

/// Acts as a crude tool for checking the set of attribute names
///
/// @note
/// if we ever choose to embrace C++, it would be **MUCH** nicer to use
/// std::set<std::string>
struct AttrNameRecorder {
  bool only_count;
  int length;
  int capacity;
  char* names[attr_name_recorder_capacity];
};

AttrNameRecorder new_AttrNameRecorder(bool only_count) {
  AttrNameRecorder out;
  out.only_count = only_count;
  out.capacity = attr_name_recorder_capacity;
  out.length = 0;
  for (int i = 0; i < attr_name_recorder_capacity; i++) {
    out.names[i] = nullptr;
  }
  return out;
}

/// acts as a destructor for all data contained within ptr
void drop_AttrNameRecorder(AttrNameRecorder* ptr) {
  if (!ptr->only_count) {
    for (int i = 0; i < ptr->length; i++) {
      delete[] ptr->names[i];
    }
  }
  ptr->length = 0;
}

/// records the attribute name
///
/// @returns GR_SUCCESS if successful. Any other value denotes an issue
///
/// @note
/// There are implicit assumptions that ptr and name are not nullptrs. It also
/// makes some sense to assume that name wasn't already recorded.
int AttrNameRecorder_record_name(AttrNameRecorder* ptr, const char* name) {
  if (ptr->length == ptr->capacity) {
    std::fprintf(stderr,
                 "encountered more attribute names than hardcoded maximum\n");
    return GR_FAIL;
  }

  if (!ptr->only_count) {
    std::size_t len_without_nullchr = std::strlen(name);
    std::size_t buf_len = len_without_nullchr + 1;
    char* buf = new char[buf_len];
    std::memcpy(buf, name, buf_len);
    ptr->names[ptr->length] = buf;
  }

  ptr->length++;
  return GR_SUCCESS;
}

int AttrNameRecorder_length(AttrNameRecorder* ptr) { return ptr->length; }

/// writes a representation of all recorded attribute names to a character
/// string @p buffer.
///
/// Just like ``std::snprintf``, this function writes at most ``buf_size - 1``
/// characters and the resulting character string is terminated with a null
/// character. As with ``std::snprintf``, when @p buf_size is zero, nothing is
/// written, but the return value is still calculated.
///
/// @param[in] buffer Pointer to the buffer where characters are written. This
///   can **only** be a nullptr if `buf_size` is 0.
/// @param[in] buf_size Specifies that `buf_size - 1` characters can be written
///   to `buffer`, in addition to the null terminator character.
/// @param[in] obj The object whose contents are being printed
///
/// \returns If successful, returns number of characters that would be written
/// for a sufficiently large buffer. This value always excludes the terminating
/// null character. Otherwise, a negative value denotes an error.
///
/// @note
/// Obviously, this could be simplified substantially if we fully embrace C++
/// (e.g. we could simply return a `std::string` rather than mimic the
/// interface of `std::snprintf`)
int AttrNameRecorder_stringify_attr_names(char* buffer, std::size_t buf_size,
                                          const AttrNameRecorder* obj) {
  int length = obj->length;

  if (obj->only_count) {
    return std::snprintf(buffer, buf_size, "{<%d unspecified attrs>}", length);
  }

  // this is inefficient, but probably okay since its only for error reporting

  const int total_bufsz = (buf_size > static_cast<std::size_t>(INT_MAX))
                              ? INT_MAX
                              : static_cast<int>(buf_size);
  int total_offset = 0;  // accumulated offset from start of buffer

  auto helper_fn = [&total_offset, total_bufsz, buffer](const char* s) -> bool {
    bool write = total_bufsz > total_offset;
    int rslt = std::snprintf((write) ? buffer + total_offset : nullptr,
                             (write) ? total_bufsz - total_offset : 0, "%s", s);
    bool is_err = rslt < 0;
    total_offset = rslt + ((is_err) ? 0 : total_offset);
    fprintf(stderr, "component, total_offset: `%s`, %d\n", s, total_offset);
    return is_err;
  };

  if (helper_fn("{")) {
    return total_offset;
  }

  for (int i = 0; i < length; i++) {
    if (i > 0 && helper_fn(", ")) {
      return total_offset;
    }
    if (helper_fn("\"") || helper_fn(obj->names[i]) || helper_fn("\"")) {
      return total_offset;
    }
  }

  helper_fn("}");
  return total_offset;
}

///@}

/// helper function that calculates ArrayShape from the "Rank" and "Dimension"
/// attributes
///
/// @param[in] dset_id The dataset identifier
/// @param[in] dset_name Name associated with dset_id (only used for producing
///     nicer error messages)
/// @param[inout] name_recorder Updated to track each attribute involved
///     in parsing a dataset's grid table properties.
///
/// @returns the parsed array shape for the grid shape properties. The caller
///     should ensure that this function was successful by calling
///     ArrayShape_is_valid on the returned value.
grackle::impl::h5io::ArrayShape shape_from_grid_attrs(
    hid_t dset_id, const char* dset_name, AttrNameRecorder* name_recorder) {
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
  if (AttrNameRecorder_record_name(name_recorder, "Rank") != GR_SUCCESS) {
    std::fprintf(
        stderr,
        "issue recording access of \"Rank\" attribute of \"%s\" dataset.",
        dset_name);
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
  if (AttrNameRecorder_record_name(name_recorder, "Dimension") != GR_SUCCESS) {
    std::fprintf(
        stderr,
        "issue recording access of \"Dimension\" attribute of \"%s\" dataset.",
        dset_name);
    return mk_invalid_array_shape();
  }

  // finally, let's format the output
  // -> we use mk_invalid_array_shape to suppress compiler warnings that out
  //    may be uninitialized
  grackle::impl::h5io::ArrayShape out = mk_invalid_array_shape();
  out.ndim = static_cast<int>(rank);
  for (int i = 0; i < out.ndim; i++) {
    out.shape[i] = grid_dimensions[i];
  }
  return out;
}

static constexpr int max_attr_name_length = 64;

/// updates axes with parsed parameters
///
/// @param[in] dset_id The dataset identifier
/// @param[in] dset_name Name associated with dset_id (only used for producing
///     nicer error messages)
/// @param[in] grid_shape Specifies the pre-parsed grid shape
/// @param[out] axes Buffers theat will be updated with the grid axes
/// @param[inout] name_recorder Updated to track each attribute involved
///     in parsing a dataset's grid table properties.
///
/// @returns GR_SUCCESS if successful. Any other value denotes an issue
int set_grid_axes_props(hid_t dset_id, const char* dset_name,
                        grackle::impl::h5io::ArrayShape grid_shape,
                        grackle::impl::h5io::GridTableAxis* axes,
                        AttrNameRecorder* name_recorder) {
  using grackle::impl::h5io::read_str_attribute;

  bool legacy_mode = H5Aexists(dset_id, "Temperature") > 0;

  for (int i = 0; i < grid_shape.ndim; i++) {
    // Step 1: fill axes[i].name with the name of the quantity that varies
    //         along axis i and fill val_attr_name with the name of the
    //         attribute containing the values of the quantitiy
    char val_attr_name[max_attr_name_length];

    if (legacy_mode && ((i + 1) == grid_shape.ndim)) {
      axes[i].name = new char[12];
      std::snprintf(axes[i].name, 12, "Temperature");

      std::snprintf(val_attr_name, max_attr_name_length, "Temperature");
    } else {
      char tmp_attr_name[max_attr_name_length];
      std::snprintf(tmp_attr_name, max_attr_name_length, "Parameter%d_Name",
                    i + 1);
      hid_t attr_id = H5Aopen_name(dset_id, tmp_attr_name);

      if (attr_id == H5I_INVALID_HID) {
        std::fprintf(stderr,
                     "Failed to open \"%s\" attribute of \"%s\" dataset.\n",
                     tmp_attr_name, dset_name);
        return GR_FAIL;
      }

      int min_buf_length = read_str_attribute(attr_id, 0, nullptr);
      if (min_buf_length < 0) {
        std::fprintf(
            stderr,
            "can't determine string size held by \"%s\" attr of \"%s\" dset.\n",
            tmp_attr_name, dset_name);
        return GR_FAIL;
      }
      axes[i].name = new char[min_buf_length];
      if (read_str_attribute(attr_id, min_buf_length, axes[i].name) < 0) {
        std::fprintf(stderr,
                     "error loading the \"%s\" attr from \"%s\" dataset\n",
                     tmp_attr_name, dset_name);
        return GR_FAIL;
      }
      H5Aclose(attr_id);

      if (AttrNameRecorder_record_name(name_recorder, tmp_attr_name) !=
          GR_SUCCESS) {
        std::fprintf(stderr,
                     "issue recording access of \"%s\" attr of \"%s\" dataset",
                     tmp_attr_name, dset_name);
        return GR_FAIL;
      }

      std::snprintf(val_attr_name, max_attr_name_length, "Parameter%d", i + 1);
    }

    // step 2: load the values associated with axis i and store them in
    //         axes[i].values
    hid_t attr_id = H5Aopen_name(dset_id, val_attr_name);
    if (attr_id == H5I_INVALID_HID) {
      std::fprintf(stderr, "Failed to open \"%s\" attr of \"%s\" dataset.\n",
                   val_attr_name, dset_name);
      return GR_FAIL;
    }

    hid_t space_id = H5Aget_space(attr_id);
    grackle::impl::h5io::ArrayShape axis_shape = shape_from_space(space_id);
    H5Sclose(space_id);
    if (axis_shape.ndim != 1) {
      std::fprintf(stderr, "The \"%s\" dataset's \"%s\" attr isn't 1D\n",
                   val_attr_name, dset_name);
      return GR_FAIL;
    } else if (axis_shape.shape[0] != grid_shape.shape[i]) {
      std::fprintf(
          stderr,
          "The \"%s\" dataset's \"Dimension\" attr suggests that axis %d has "
          "%lld values. The \"%s\" attribute actually %lld values\n",
          dset_name, i, static_cast<long long>(grid_shape.shape[i]),
          val_attr_name, static_cast<long long>(axis_shape.shape[0]));
      return GR_FAIL;
    }

    axes[i].values = new double[grid_shape.shape[i]];
    if (H5Aread(attr_id, HDF5_R8, axes[i].values) < 0) {
      H5Aclose(attr_id);
      std::fprintf(stderr, "failed to read \"%s\" attr of \"%s\" dataset.\n",
                   val_attr_name, dset_name);
      return GR_FAIL;
    }
    H5Aclose(attr_id);
    if (AttrNameRecorder_record_name(name_recorder, val_attr_name) !=
        GR_SUCCESS) {
      std::fprintf(stderr,
                   "issue recording access of \"%s\" attr of \"%s\" dataset",
                   val_attr_name, dset_name);
      return GR_FAIL;
    }
  }
  return GR_SUCCESS;
}

grackle::impl::h5io::GridTableProps mk_invalid_GridTableProps() {
  grackle::impl::h5io::GridTableProps out;
  out.table_shape = mk_invalid_array_shape();
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    out.axes[i].name = nullptr;
    out.axes[i].values = nullptr;
  }
  return out;
}

/// helper function that parses the GridTableProps from dataset attributes
///
/// @param[in] dset_id The dataset identifier
/// @param[in] dset_name Name associated with dset_id (only used for producing
///     nicer error messages)
/// @param[inout] name_recorder Updated to track each attribute involved
///     in parsing a dataset's grid table properties.
///
/// @returns Returns the appropriate GridTableProps object. The caller should
///     use the GridTableProps_is_valid function to confirm that the function
///     was successful.
grackle::impl::h5io::GridTableProps parse_GridTableProps_helper(
    hid_t dset_id, const char* dset_name, AttrNameRecorder* name_recorder) {
  // setup the output object so that we're always prepared to return an object
  // - the default object is constructed to denote a failure
  grackle::impl::h5io::GridTableProps out = mk_invalid_GridTableProps();

  // if optional Description attribute is present, record that we accessed it
  // (this is done purely for error-handling purposes).
  if ((H5Aexists(dset_id, "Description") > 0) &&
      (AttrNameRecorder_record_name(name_recorder, "Description") !=
       GR_SUCCESS)) {
    return out;
  }

  // infer the shape of the table from the attributes
  grackle::impl::h5io::ArrayShape inferred_shape =
      shape_from_grid_attrs(dset_id, dset_name, name_recorder);
  if (!ArrayShape_is_valid(inferred_shape)) {
    return out;
  }

  // parse the quantities along each axis
  if (set_grid_axes_props(dset_id, dset_name, inferred_shape, out.axes,
                          name_recorder) != GR_SUCCESS) {
    drop_GridTableProps(&out);
    return out;
  }

  out.table_shape = inferred_shape;  // we intentionally do this as the very
                                     // last step!
  return out;
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
  if (dset_name == nullptr) {  // sanity check!
    std::fprintf(stderr, "dset_name is a nullptr");
    return mk_invalid_GridTableProps();
  }

  hid_t dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Can't open \"%s\" dataset.\n", dset_name);
    return mk_invalid_GridTableProps();
  }

  // construct object to count accessed attributes
  AttrNameRecorder attr_counter = new_AttrNameRecorder(true);
  // actually parse GridTableProps
  GridTableProps out =
      parse_GridTableProps_helper(dset_id, dset_name, &attr_counter);
  // get the attribute count and cleanup the counter
  int total_accessed_attrs_count = AttrNameRecorder_length(&attr_counter);
  drop_AttrNameRecorder(&attr_counter);

  if (!GridTableProps_is_valid(out)) {
    H5Dclose(dset_id);
    // parse_GridTableProps_helper already printed appropriate errors messages
    return out;
  }

  // now, we will perform a check validating that dataset doesn't have
  // unexpected attributes (i.e. we want to avoid the situation where the
  // creator expects new attributes to introduce some behavior that Grackle
  // doesn't know about)
  int num_attrs = get_num_attrs(dset_id, dset_name);
  if (num_attrs == -2) {
    num_attrs = total_accessed_attrs_count;
  } else if (num_attrs == -1) {
    H5Dclose(dset_id);
    drop_GridTableProps(&out);
    // get_num_attrs already printed error messages in this case
    return out;
  } else if (num_attrs != total_accessed_attrs_count) {
    drop_GridTableProps(&out);
    // to provide a detailed error message that will make it straight-forward
    // to debug the underlying problem, we are going to call
    // parse_GridTableProps_helper, but this time, we are going to actually
    // record every accessed name.
    AttrNameRecorder attr_name_recorder = new_AttrNameRecorder(false);
    // let's parse GridTableProps (again)
    GridTableProps tmp =
        parse_GridTableProps_helper(dset_id, dset_name, &attr_name_recorder);

    int bufsz_without_nullchr =
        AttrNameRecorder_stringify_attr_names(nullptr, 0, &attr_name_recorder);
    int bufsz = bufsz_without_nullchr + 1;
    char* stringified_attr_list = nullptr;
    if (bufsz > 0) {
      stringified_attr_list = new char[bufsz];
      if (AttrNameRecorder_stringify_attr_names(stringified_attr_list,
                                                static_cast<std::size_t>(bufsz),
                                                &attr_name_recorder) < 0) {
        delete[] stringified_attr_list;
        stringified_attr_list = nullptr;
      }
    }

    if (stringified_attr_list == nullptr) {
      std::fprintf(
          stderr,
          "encountered error while stringifying list of recognized attrs for "
          "the \"%s\" dataset\n",
          dset_name);
    } else {
      std::fprintf(
          stderr,
          "while parsing the %dD grid table properties from attributes of the "
          "\"%s\" dataset, encountered (one or more) unrecognized attributes. "
          "Recognized attributes include: %s\n",
          tmp.table_shape.ndim, dset_name, stringified_attr_list);
      delete[] stringified_attr_list;
    }
    drop_GridTableProps(&tmp);
    drop_AttrNameRecorder(&attr_name_recorder);
    H5Dclose(dset_id);
    return mk_invalid_GridTableProps();
  }

  // check for consistency between the GridTableProps and the dataset shape
  hid_t space_id = H5Dget_space(dset_id);
  ArrayShape actual_shape = shape_from_space(space_id);
  H5Sclose(space_id);
  if (!ArrayShape_is_null(actual_shape) &&
      !ArrayShape_is_same(out.table_shape, actual_shape)) {
    drop_GridTableProps(&out);
    H5Dclose(dset_id);
    std::fprintf(
        stderr,
        "The grid table properties parsed from the attributes of the \"%s\" "
        "dataset are inconsistent with the dataset's shape.\n",
        dset_name);
    return mk_invalid_GridTableProps();
  }

  H5Dclose(dset_id);
  return out;
}
