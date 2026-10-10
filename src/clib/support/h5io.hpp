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

#ifndef SUPPORT_H5IO_HPP
#define SUPPORT_H5IO_HPP

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "hdf5.h"
#include "grackle.h"
#include "config.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace h5io {

/// @brief retrieve a copy the string encoded in the specified hdf5 attribute
///
/// @param attr_id Attribute identifier
///
/// @note
/// The function reports an error if the user tries to reads a utf8-encoded
/// string that contains non-ASCII characters.
std::optional<std::string> read_str_attribute(hid_t attr_id);

/// @brief retrieve a copy of the string encoded in the specified hdf5 dataset
///
/// @param file_id File identifier
/// @param dset_name The name of the dataset to read attributes from.
///
/// @note
/// The function reports an error if the user tries to reads a utf8-encoded
/// string that contains non-ASCII characters.
std::optional<std::string> read_str_dataset(hid_t file_id,
                                            const char* dset_name);

/// @brief represents a contiguous array shape
///
/// @note
/// An ndim of -1 corresponds to a null dataset (i.e. ``H5S_NULL``). An ndim of
/// 0 denotes a scalar
struct ArrayShape {
  int ndim;
  std::int64_t shape[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  // this only exists to avoid warnings about potentially uninitialized members
  ArrayShape() : ndim{0} {
    for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
      shape[i] = 0;
    }
  }

  /// @brief checks whether shape refers to a scalar
  bool is_scalar() const { return ndim == 0; }

  /// @brief checks whether the shape is null
  bool is_null() const { return ndim == -1; }

  /// @brief calculates the total number of elements in the array
  int64_t elem_count() const {
    if (ndim < 0) {
      return -1;
    } else {  // this works even if shape.ndim is 0
      int64_t product = 1;
      for (int i = 0; i < ndim; i++) {
        product *= shape[i];
      }
      return product;
    }
  }

  /// @brief overloads the equality comparison (``==``) operation
  bool operator==(const ArrayShape& other) const {
    for (int i = 0; i < ndim; i++) {
      if (shape[i] != other.shape[i]) {
        return false;
      }
    }
    return true;
  }
};

/// load the shape of the dataset
std::optional<ArrayShape> read_dataset_shape(hid_t file_id,
                                             const char* dset_name);

/// read the dataset named dset_name from file_id into buffer
///
/// When expected_shape is provided, the dataset's shape is checked before
/// the data is read
int read_dataset(hid_t file_id, const char* dset_name, double* buffer,
                 const ArrayShape* expected_shape = nullptr);

struct GridTableAxis {
  std::string name;
  std::vector<double> values;
};

/// Used to represent properties of an interpolation table
struct GridTableProps {
  ArrayShape table_shape;
  GridTableAxis axes[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// @brief overloads the equality comparison (``==``) operation
  ///
  /// @note returns false if the either object is invalid
  bool operator==(const GridTableProps& o) const;

  /// @brief checks whether the specified dataset has consistent grid properties
  ///
  /// @returns `true` indicates that the dataset has consistent properties
  bool assert_is_consistent(hid_t file_id, const char* dset_name) const;
};

/// parses the GridTableProps from dataset attributes
///
/// @param[in] file_id File identifier
/// @param[in] dset_name The name of the dataset to read attributes from.
///
/// @returns Returns the appropriate GridTableProps object.
std::optional<GridTableProps> parse_GridTableProps(hid_t file_id,
                                                   const char* dset_name);

}  // namespace h5io
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_H5IO_HPP
