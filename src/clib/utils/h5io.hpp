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

/// copies the string encoded in the specified hdf5 attribute and into
/// ``buffer`` as a null-terminated string, and returns ``min_req_bufsz``.
///
/// ``min_req_bufsz`` is the minimum required ``bufsz`` (i.e. ``bufsz`` is the
/// the length of ``buffer``) where this function will try to load the
/// attribute.
/// - this is the maximum length of the string (including the null character).
///   To put it another way, after this function succeeds, std::strlen(buffer),
///   will return a ``x`` that is bounded by ``0 <= x <= min_req_bufsz``.
/// - this function's behavoir is described in terms of ``min_req_bufsz``,
///   rather than the exact required buffer length because the exact length
///   can't be determined without loading the buffer.
///
/// This function fails if ``bufsz`` is smaller than ``min_req_bufsz``, unless
/// unless ``bufsz`` is zero. In the event that ``bufsz`` is zero, nothing is
/// written and ``buffer`` may be a ``nullptr``. In this case, the function
/// returns ``min_req_bufsz``. If this function fails, a negative value is
/// returned.
///
/// @note
/// The function reports an error if it reads a utf8-encoded string that
/// contains non-ASCII characters
int read_str_attribute(hid_t attr_id, int bufsz, char* buffer);

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
/// @param[in] dset_name The name of the dataset that the parameters are
///     attached to.
///
/// @returns Returns the appropriate GridTableProps object. The caller should
///     use the GridTableProps_is_valid function to confirm that the function
///     was successful.
GridTableProps parse_GridTableProps(hid_t file_id, const char* dset_name);

}  // namespace grackle::impl::h5io

#endif /* UTILS_H5IO_HPP */
