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

#include <vector>

#include "support/status_reporting.hpp"
#include "./iterator_adaptor.hpp"

namespace grtest {

/// represents numbers evenly spaced on a geometric scale
///
/// @note
/// I feel like I've written this exact functionality more times than I
/// can count.
class LinSpaced {
  double first_;
  double last_;
  int n_entries_;
  double step_;

public:
  LinSpaced() = delete;
  LinSpaced(const LinSpaced&) = default;
  LinSpaced(double first, double last, int n_entries)
      : first_{first}, last_{last}, n_entries_{n_entries} {
    GRIMPL_REQUIRE(n_entries >= 2, "must hold at least 2 entries");
    step_ = (last_ - first_) / n_entries_;
  }

  bool is_increasing() const { return last_ > first_; }
  int size() const { return n_entries_; }

  // overloads the elementwise array access operation
  // -> if an instance is called `obj`, this method is invoked by `obj[3]`
  [[gnu::always_inline]] inline double operator[](int idx) const noexcept {
    bool is_last = idx + 1 == n_entries_;
    double tmp = first_ + step_ * idx;
    return (!is_last) * tmp + is_last * last_;
  }

  /// writes values in the sequence to the range specified by [`start`, `stop`)
  void fill_buf(double* start, double* stop) const {
    GRIMPL_REQUIRE(start + n_entries_ == stop, "buffer has incorrect size");
    for (int i = 0; i < n_entries_; i++) {
      start[i] = first_ + step_ * i;
    }
    start[n_entries_ - 1] = last_;  // avoid round-off errors
  }

  // this is a convenience
  std::vector<double> to_vec() const {
    std::vector<double> out(n_entries_);
    fill_buf(out.data(), out.data() + n_entries_);
    return out;
  }

  bool operator==(const LinSpaced& other) const {
    return (
        (first_ == other.first_) && (last_ == other.last_) &&
        (n_entries_ == other.n_entries_) &&
        (step_ == other.step_)  // <- the last comparison shouldn't be needed...
    );
  }
};

/// @brief Represents a single coordinate on a grid modelled by a @ref CoordGrid
struct Coordinate {
  double components[5];
};

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
/// this class assumes that the very first coordinate component is the "fast
/// axis" and the very last component is the "slow axis". When iterating over
/// the coordinates, the component along the "slow axis" changes with the least
/// frequency. This convention is based on the convention currently used for
/// Grackle's interpolation tables
class CoordGrid {
  using iterator =
      IteratorAdaptor<gridcoord_detail::CoordGridItrPlugin<CoordGrid>>;

  std::vector<LinSpaced> coords_;
  gridcoord_detail::IdxCycle cycle_info_[gridcoord_detail::MAX_RANK];
  int rank_;

  iterator setup_iterator_(bool make_start) const {
    int size = this->size();
    gridcoord_detail::CoordGridItrPlugin<CoordGrid> plugin(*this);
    return iterator(make_start ? 0 : size, size, plugin);
  }

public:
  CoordGrid(int rank, const LinSpaced* coords) : rank_{rank} {
    GRIMPL_REQUIRE(rank >= 1 && rank <= gridcoord_detail::MAX_RANK,
                   "passed an invalid rank");
    GRIMPL_REQUIRE(coords != nullptr, "passed a nullptr");

    int running_product = 1;
    for (int i = 0; i < rank; i++) {
      coords_.push_back(coords[i]);

      // setup cycle_info_[i]
      if (i > 0) {
        running_product *= cycle_info_[i - 1].cycle_len;
      }
      int extent = coords[i].size();
      int element_reps = running_product;
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
  Coordinate from_flatidx(int flat_idx) {
    Coordinate out;
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