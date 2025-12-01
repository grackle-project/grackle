//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares functions to help read from hdf5 files
///
//===----------------------------------------------------------------------===//

#ifndef UTILS_H5IO_HPP
#define UTILS_H5IO_HPP

#include <cstdint>

#include "hdf5.h"
#include "grackle.h"

namespace grackle::impl::h5io {

/// copies the string encoded in the specified hdf5 dataset into ``buffer`` as
/// a null-terminated string, and returns ``min_req_bufsz`` (if successful).
///
/// @param[in] attr_id Attribute identifier
/// @param[in] bufsz Number of ascii characters (including the null terminating
///   character) that can be written to @p buffer
/// @param[out] buffer Pointer to the buffer where characters are written. This
///   can **only** be a nullptr if @p bufsz is 0.
///
/// @returns If successful, returns ``min_req_bufsz`` (see below). Otherwise,
///   returns a negative value.
///
/// ``min_req_bufsz`` is the minimum required @p bufsz that this function must
/// receive for it to attempt to load the string.
/// - this is the maximum length of the string (including the null character).
///   Thus, after succesfully calling this function
///   ``std::strlen(buffer) + 1 <= min_req_bufsz``.
/// - this function's behavoir is described in terms of ``min_req_bufsz``,
///   rather than the exact required buffer length because the exact length
///   can't be determined without loading the buffer.
///
/// This function fails if @p bufsz is smaller than ``min_req_bufsz``, unless
/// @p bufsz is zero. In that case, nothing is written to @p buffer and the
/// returns ``min_req_bufsz``. The function reports an error if the user tries
/// to reads a utf8-encoded string that contains non-ASCII characters.
///
/// @note
/// If we are more willing to embrace C++, we could return a std::string or
/// std::vector rather than requiring a pre-allocated buffer
int read_str_attribute(hid_t attr_id, int bufsz, char* buffer);

/// copies the string encoded in the specified hdf5 dataset into ``buffer`` as
/// a null-terminated string, and returns ``min_req_bufsz`` (if successful).
///
/// @param[in] file_id File identifier
/// @param[in] dset_name The name of the dataset to read attributes from.
/// @param[in] bufsz Number of ascii characters (including the null terminating
///   character) that can be written to @p buffer
/// @param[out] buffer Pointer to the buffer where characters are written. This
///   can **only** be a nullptr if @p bufsz is 0.
///
/// @returns If successful, returns ``min_req_bufsz`` (see below). Otherwise,
///   returns a negative value.
///
/// ``min_req_bufsz`` is the minimum required @p bufsz that this function must
/// receive for it to attempt to load the string.
/// - this is the maximum length of the string (including the null character).
///   Thus, after succesfully calling this function
///   ``std::strlen(buffer) + 1 <= min_req_bufsz``.
/// - this function's behavoir is described in terms of ``min_req_bufsz``,
///   rather than the exact required buffer length because the exact length
///   can't be determined without loading the buffer.
///
/// This function fails if @p bufsz is smaller than ``min_req_bufsz``, unless
/// @p bufsz is zero. In that case, nothing is written to @p buffer and the
/// returns ``min_req_bufsz``. The function reports an error if the user tries
/// to reads a utf8-encoded string that contains non-ASCII characters.
///
/// @note
/// If we are more willing to embrace C++, we could return a std::string or
/// std::vector rather than requiring a pre-allocated buffer
///
/// @note
/// The choice to accept @p file_id and @p dset_name, rather than an already
/// open dataset identifier, was made for consistency with the interfaces of
/// read_dataset and read_dataset_shape. However, the choice is a little
/// "clunky" because this function is usually called twice (once to query
/// ``min_req_bufsz`` and once to load the string).
int read_str_dataset(hid_t file_id, const char* dset_name, int bufsz,
                     char* buffer);

/// represents a contiguous array shape
///
/// @note
/// An ndim of -1 corresponds to a null dataset (i.e. ``H5S_NULL``). Any other
/// negative value denotes an invalid shape. An ndim of 0 denotes a scalar
struct ArrayShape {
  int ndim;
  std::int64_t shape[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
};

/// checks whether shape is valid
inline bool ArrayShape_is_valid(ArrayShape shape) { return shape.ndim >= -1; }

/// checks whether shape refers to a scalar
inline bool ArrayShape_is_scalar(ArrayShape shape) { return shape.ndim == 0; }

/// checks whether shape is null
inline bool ArrayShape_is_null(ArrayShape shape) { return shape.ndim == -1; }

/// calculates the total number of elements in the array
inline std::int64_t ArrayShape_elem_count(ArrayShape shape) {
  if (shape.ndim < 0) {
    return -1;
  } else {  // this works even if shape.ndim is 0
    std::int64_t product = 1;
    for (int i = 0; i < shape.ndim; i++) {
      product *= shape.shape[i];
    }
    return product;
  }
}

/// checks whether shape_a and shape_b are the same
///
/// returns false if the either shape is invalid
inline bool ArrayShape_is_same(ArrayShape shape_a, ArrayShape shape_b) {
  if ((!ArrayShape_is_valid(shape_a)) || (shape_a.ndim != shape_b.ndim)) {
    return false;
  }
  for (int i = 0; i < shape_a.ndim; i++) {
    if (shape_a.shape[i] != shape_b.shape[i]) {
      return false;
    }
  }
  return true;
}

/// load the shape of the dataset
ArrayShape read_dataset_shape(hid_t file_id, const char* dset_name);

/// read the dataset named dset_name from file_id into buffer
///
/// When expected_shape is provided, the dataset's shape is checked before
/// the data is read
int read_dataset(hid_t file_id, const char* dset_name, double* buffer,
                 const ArrayShape* expected_shape = nullptr);

struct GridTableAxis {
  char* name;
  double* values;
};

/// Used to represent properties of an interpolation table
struct GridTableProps {
  ArrayShape table_shape;
  GridTableAxis axes[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
};

inline bool GridTableProps_is_valid(GridTableProps grid_props) {
  return ArrayShape_is_valid(grid_props.table_shape);
}

/// acts as a destructor for the contents within ptr
inline void drop_GridTableProps(GridTableProps* ptr) {
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    if (ptr->axes[i].name != nullptr) {
      delete[] ptr->axes[i].name;
    }
    if (ptr->axes[i].values != nullptr) {
      delete[] ptr->axes[i].values;
    }
    ptr->axes[i].name = nullptr;
    ptr->axes[i].values = nullptr;
  }
}

/// parses the GridTableProps from dataset attributes
///
/// @param[in] file_id File identifier
/// @param[in] dset_name The name of the dataset to read attributes from.
///
/// @returns Returns the appropriate GridTableProps object. The caller should
///     use the GridTableProps_is_valid function to confirm that the function
///     was successful.
GridTableProps parse_GridTableProps(hid_t file_id, const char* dset_name);

}  // namespace grackle::impl::h5io

#endif /* UTILS_H5IO_HPP */
