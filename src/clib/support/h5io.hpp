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

/// @brief represents a contiguous array shape
///
/// @note
/// An ndim of -1 corresponds to a null dataset (i.e. ``H5S_NULL``). Any other
/// negative value denotes an invalid shape. An ndim of 0 denotes a scalar
struct ArrayShape {
  int ndim;
  std::int64_t shape[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// @brief checks whether shape is valid
  bool is_valid() const { return ndim >= -1; }

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
  ///
  /// @note returns false if the either shape is invalid
  bool operator==(const ArrayShape& other) const {
    if ((!is_valid()) || (ndim != other.ndim)) {
      return false;
    }
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

  bool is_valid() const { return table_shape.is_valid(); }

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
