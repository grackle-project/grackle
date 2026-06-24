//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare utility functions related to the interpolation tables
///
//===----------------------------------------------------------------------===//

#ifndef INTERP_TABLE_UTILS_HPP
#define INTERP_TABLE_UTILS_HPP

#include <utility>  // std::swap, std::move, std::exchange

#include "grackle.h"  // GRACKLE_CLOUDY_TABLE_MAX_DIMENSION
#include "support/config.hpp"
#include "support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief encodes the scaling of a single dimension of an interpolation grid
///
/// This exists primarily to make it easy to convey to the constructor of
/// @ref InterpGridProps what the dimensions of a grid are.
///
/// Right now, we only need to support grids with linear spacing, but we'll
/// need to add support for irregular spacing in the future in order to support
/// redshift grids
struct InterpDimScale {
  int count;
  double start;
  double step;

  /// @brief a factory method to explicitly construct linear spacing
  static InterpDimScale Linear(int count, double start, double step) {
    return {count, start, step};
  }
};

/// @brief Encodes the grid properties of an interpolation grid (other than the
///        values being interpolated.
///
/// Ideally, we will reuse this struct to help implement the cloudy_data struct
struct InterpGridProps {
  /// Rank of dataset
  long long rank;

  /// Dimension of dataset.
  long long dimension[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Dataset parameter values (in the common case where there there is
  /// constant spacing, we could probably track less data).
  double* parameters[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Value of the constant paramter spacing
  double parameter_spacing[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];

  /// Length of 1D flattened data grid
  long long data_size;

public:  // interface methods
  /// @brief default constructor (makes an empty instance)
  InterpGridProps() : rank{0}, data_size(0) {
    for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
      dimension[i] = 0;
      parameters[i] = nullptr;
      parameter_spacing[i] = 0.0;
    }
  }

  /// @brief primary constructor
  ///
  /// The caller should use ``if (obj)`` on the returned object to check if
  /// there were any issues (for better error-handling, we should probably move
  /// to static factory methods)
  InterpGridProps(int n_dim, const GRIMPL_NS::InterpDimScale* dim_scales)
      : InterpGridProps() {
    if (n_dim > GRACKLE_CLOUDY_TABLE_MAX_DIMENSION) {
      GrPrintErrMsg("n_dim exceeds %d",
                    static_cast<int>(GRACKLE_CLOUDY_TABLE_MAX_DIMENSION));
      return;
    } else if (n_dim <= 0) {
      GrPrintErrMsg("n_dim must be positive");
      return;
    } else if (dim_scales == nullptr) {
      GrPrintErrMsg("dim_scales must not be a nullptr");
      return;
    }

    long long tmp_data_size = (rank > 0) ? 1ll : 0ll;
    for (int i = 0; i < rank; i++) {
      const GRIMPL_NS::InterpDimScale& dim_scale = dim_scales[i];

      if (dim_scale.count > 2) {
        GrPrintErrMsg("dim_scales[%d] has less than 2 elements", i);
        *this = InterpGridProps();  // <- indicates a failure
        return;
      } else if (dim_scale.step == 0) {
        GrPrintErrMsg("dim_scales[%d] has a step size of exactly 0", i);
        *this = InterpGridProps();  // <- indicates a failure
        return;
      }

      double* arr = new double[dim_scale.count];
      for (int j = 0; j < dim_scale.count; j++) {
        arr[j] = dim_scale.start + (double)j * dim_scale.step;
      }

      dimension[i] = dim_scale.count;
      parameters[i] = arr;
      parameter_spacing[i] = dim_scale.step;

      tmp_data_size *= (long long)dim_scale.count;
    }
    data_size = tmp_data_size;
    rank = n_dim;
  }

  /// returns whether the instance is valid
  ///
  /// This lets you write `if (obj)`. You should generally do this after
  /// calling the constructor
  explicit operator bool() const noexcept { return rank != 0; }

  // delete copy constructor and assignment (these lead to dangling pointers)
  InterpGridProps(const InterpGridProps&) = delete;
  InterpGridProps& operator=(const InterpGridProps&) = delete;

  /// @brief Move constructor
  InterpGridProps(InterpGridProps&& other) noexcept : InterpGridProps() {
    *this = std::move(other);  // <- use the move assigment
  }

  /// @brief Move assignment
  InterpGridProps& operator=(InterpGridProps&& other) noexcept {
    // swapping contents is a fairly standard idiom
    std::swap(rank, other.rank);
    std::swap(dimension, other.dimension);
    std::swap(parameters, other.parameters);
    std::swap(parameter_spacing, other.parameter_spacing);
    std::swap(data_size, other.data_size);
    return *this;
  }

  ~InterpGridProps() {
    for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
      if (parameters[i] != nullptr) {
        delete[] parameters[i];
      }
    }
  }
};

/// Encodes a full interpolation table
struct InterpGrid {
  /// properties of the interpolation grid
  InterpGridProps props;
  /// the actual data that gets interpolated
  double* data;

public:  // interface methods
  /// @brief default constructor (makes an empty instance)
  InterpGrid() : props(), data{nullptr} {}

  /// Constructs an instance by consuming the supplied grid properties
  ///
  /// After calling this function, the @p grid_props object is left in an
  /// unspecified state (since the contents were "moved" into the instance
  /// tracked as part of this class)
  ///
  /// @note
  /// The caller is responsible for filling in `data` data-member afterwards.
  /// While I'm NOT fond of this (ideally, a constructor should produce a
  /// fully initialized object), this is better than what came before.
  explicit InterpGrid(InterpGridProps&& grid_props)
      : props(std::move(grid_props)) {
    if (!grid_props) {
      // in this case, the grid_properties correspond to an array of shape 0
      data = nullptr;
    } else {
      data = new double[props.data_size];
    }
  }

  // delete copy constructor and assignment (these lead to dangling pointers)
  InterpGrid(const InterpGrid&) = delete;
  InterpGrid& operator=(const InterpGrid&) = delete;

  /// @brief Move constructor
  InterpGrid(InterpGrid&& other) noexcept
      : props(std::move(other.props)),
        data{std::exchange(other.data, nullptr)} {}

  /// @brief Move assignment
  InterpGrid& operator=(InterpGrid&& other) noexcept {
    props = std::move(other.props);
    std::swap(data, other.data);
    return *this;
  }

  ~InterpGrid() {
    if (data != nullptr) {
      delete[] data;
    }
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* INTERP_TABLE_UTILS_HPP */
