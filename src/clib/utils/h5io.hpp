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

/// checks whether shape is valid
inline bool ArrayShape_is_scalar(ArrayShape shape) { return shape.ndim == 0; }

/// checks whether shape is null
inline bool ArrayShape_is_null(ArrayShape shape) { return shape.ndim == -1; }

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

}  // namespace grackle::impl::h5io

#endif /* UTILS_H5IO_HPP */
