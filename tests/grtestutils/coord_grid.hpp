//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare and implement machinery related to iterating over grid coordinates
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_COORD_GRID_HPP
#define GRTESTUTILS_COORD_GRID_HPP

#include <iostream>  // needed to teach googletest how to print Coordinate
#include <vector>

#include "support/status_reporting.hpp"
#include "./iterator_adaptor.hpp"

namespace grtest {

/// @brief Represents a single coordinate on a grid modelled by a @ref CoordGrid
struct Coordinate {
  double components[5];
  int rank;

  /// overloads array access operation
  ///
  /// Give an instance `c` and an integer `i`, this method lets us write c[i],
  /// which is equivalent to `c.components[i]`
  double operator[](int index) const noexcept { return components[index]; }
};

// teach GoogleTest how to print grtest::Coordinate
inline void PrintTo(const Coordinate& c, std::ostream* os) {
  *os << '{' << c.components[0];
  for (int i = 1; i < c.rank; i++) {
    *os << ", " << c.components[i];
  }
  *os << '}';
}

namespace gridcoord_detail {

/// @brief The max rank of a @ref CoordGrid
static constexpr int MAX_RANK = 5;

/// The plugin type used to implement to implement the iterator over a
/// @ref Coordinate (using @ref IteratorAdaptor)
///
/// @note
/// We only chose to parameterize this in terms of a template out of
/// convenience. Alternatively, we could have simply forward-declared
/// @ref CoordGrid and implemented the member functions of this type in
/// a source file.
template <typename CoordGridType>
class CoordGridItrPlugin {
  const CoordGridType& grid_props_;

public:
  using ValueType = Coordinate;

  explicit CoordGridItrPlugin(const CoordGridType& grid_props)
      : grid_props_(grid_props) {}

  Coordinate operator()(unsigned long long i) const {
    return grid_props_.from_flatidx(i);
  }

  bool operator==(const CoordGridItrPlugin<CoordGridType>& other) const {
    return grid_props_ == other.grid_props_;
  }
};

struct IdxCycle {
  int extent;
  int element_reps;
  int cycle_len;
};

}  // namespace gridcoord_detail

/// Represents a coordinate grid
///
/// @note
/// this class assumes that the very first coordinate component is the "slow
/// axis" and the very last component is the "fast axis". When iterating over
/// the coordinates, the component along the "slow axis" changes with the least
/// frequency. This convention is based on the convention currently used for
/// Grackle's interpolation tables
class CoordGrid {
  using iterator =
      IteratorAdaptor<gridcoord_detail::CoordGridItrPlugin<CoordGrid>>;

  std::vector<std::vector<double>> coords_;
  gridcoord_detail::IdxCycle cycle_info_[gridcoord_detail::MAX_RANK];
  int rank_;

  iterator setup_iterator_(bool make_start) const {
    int size = this->size();
    gridcoord_detail::CoordGridItrPlugin<CoordGrid> plugin(*this);
    return iterator(make_start ? 0 : size, size, plugin);
  }

public:
  explicit CoordGrid(std::vector<std::vector<double>> coords)
      : coords_(coords), rank_(coords.size()) {
    GRIMPL_REQUIRE(rank_ >= 1 && rank_ <= gridcoord_detail::MAX_RANK,
                   "passed an invalid rank");

    for (int i = rank_ - 1; i >= 0; i--) {
      // setup cycle_info_[i]
      int extent = coords[i].size();
      int element_reps = (i + 1 == rank_) ? 1 : cycle_info_[i + 1].cycle_len;
      int cycle_len = extent * element_reps;
      cycle_info_[i] = {extent, element_reps, cycle_len};
    }
  }

  /// @brief overloads the equality check operation
  bool operator==(const CoordGrid& other) const {
    // no need to compare cycle_info_ since it is derived from the other data
    // members
    return rank_ == other.rank_ && coords_ == other.coords_;
  }

  /// @brief get the number of dimensions
  int rank() const noexcept { return rank_; }

  int axis_len(int axis) const {
    GRIMPL_REQUIRE(axis >= 0 && axis < rank_, "passed an invalid arg");
    return coords_[axis].size();
  }

  /// The number of locations in a grid
  int size() const {
    int out = 1;
    for (int i = 0; i < rank_; i++) {
      out *= coords_[i].size();
    }
    return out;
  }

  /// @note Returns the appropriate
  Coordinate from_flatidx(int flat_idx) const {
    Coordinate out;
    out.rank = rank_;
    for (int axis = 0; axis < rank_; axis++) {
      const gridcoord_detail::IdxCycle& tmp = cycle_info_[axis];
      out.components[axis] =
          coords_[axis][(flat_idx % tmp.cycle_len) / tmp.element_reps];
    }
    return out;
  }

  /// Returns an iterator to the first element of the CoordGrid
  ///
  /// @warning
  /// Mutating the described grid coordinates will invalidate the iterator
  iterator begin() const { return setup_iterator_(true); }

  /// Returns iterator to the element after the last element in CoordGrid
  ///
  /// @warning
  /// Mutating the described grid coordinates will invalidate the iterator
  iterator end() const { return setup_iterator_(false); }
};

}  // namespace grtest

#endif  // GRTESTUTILS_COORD_GRID_HPP