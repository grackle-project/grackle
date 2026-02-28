//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares functionality for working with/testing multidimensional "views"
///
/// For additional context, we use the word "view" in the same sense as
/// std::string_view or a std::span (it is a view onto externally owned data)
///
//===----------------------------------------------------------------------===//

#ifndef GRTESTUTILS_VIEW_HPP
#define GRTESTUTILS_VIEW_HPP

#include <cassert>
#include <cstddef>  // std::size_t
#include <optional>
#include <ostream>
#include <string>
#include <utility>  // std::pair

#include "./harness/status.hpp"

#include "status_reporting.h"

namespace grtest {

/// To be used with @ref IdxMapping
enum struct DataLayout {
  LEFT,  ///< the leftmost dimension has a stride 1
  RIGHT  ///< the rightmost dimension has a stride 1
};

/// Maps multi-dimensional indices to a 1D pointer offset
///
/// Broader Context
/// ===============
/// To best describe this type, it's insightful to draw comparisons with C++
/// conventions.
///
/// Background
/// ----------
/// For some background, C++23 introduced `std::mdspan` to describe
/// multi-dimensional views. A `std::mdspan` is parameterized by
/// - the data's extents (aka the shape)
/// - the data's layout, which dictates how a multidimensional index is mapped
///   to a 1D pointer offset
///
/// For views of contiguous data there are 2 obvious layouts:
/// 1. layout-right: where the stride is `1` along the rightmost extent.
///    - for extents `{a,b,c}`, an optimal nested for-loop will iterates from
///      `0` up to `a` in the outermost loop and from `0` up to `c`
///      in the innermost loop
///    - this is the "natural layout" for a multidimensional c-style array
///      `arr[a][b][c]`
/// 2. layout-left: where the stride is `1` along the leftmost extent
///    - for extents `{a,b,c}`, an optimal nested for-loop will iterates from
///      `0` up to `c` in the outermost loop and from `0` up to `a`
///      in the innermost loop
///    - this is the "natural layout" for a multidimensional fortran array
///      `arr[a][b][c]`
///
/// > Aside: More sophisticated layouts are possible when strides along each
/// > axis aren't directly tied to extents (this comes up when making subviews).
///
/// About this type
/// ---------------
/// This type specifies the data layout, extents, and provides the mapping. You
/// draw analogies with C++23 types:
/// - `IdxMapping<DataLayout::LEFT>`   <--->  `std::layout_left::mapping`
/// - `IdxMapping<DataLayout::RIGHT>`  <--->  `std::layout_right::mapping`
///
/// At the time of writing, the type **only** represents contiguous data layouts
///
/// @note
/// The template specialization using @ref DataLayout::RIGHT mostly exists for
/// exposition purposes
/// - it's useful to talk about this scenario since it is the "natural" layout
///   for C and C++
/// - in practice, Grackle was written assuming @ref DataLayout::LEFT. Thus,
///   most loops will initially get written assuming that layout (but we always
///   have the option to conditionally provide the better kind of loop)
template <DataLayout Layout>
struct IdxMapping {
  static_assert(Layout == DataLayout::LEFT || Layout == DataLayout::RIGHT,
                "A new layout type was introduced that we don't yet support");
  static constexpr DataLayout layout = Layout;

  /// max rank value
  static constexpr int MAX_RANK = 3;

private:
  // attributes:
  int rank_;               ///< the number of dimensions
  int extents_[MAX_RANK];  ///< dimensions of the multi-dimensional index-space

  // private functions:

  // default constructor is only invoked by factory methods
  // -> this explicitly set each extent_ to a value of 0
  IdxMapping() : rank_(0), extents_{} {}

  /// factory method that aborts with an error message upon failure
  static IdxMapping<Layout> create_or_abort_(int rank, const int* extents) {
    std::pair<IdxMapping<Layout>, Status> pair =
        IdxMapping<Layout>::try_create(rank, extents);
    if (pair.second.is_err()) {
      std::string tmp = pair.second.to_string();
      GR_INTERNAL_ERROR("%s", tmp.c_str());
    }
    return pair.first;
  }

public:
  /// factory method that tries to create an instance
  static std::pair<IdxMapping<Layout>, Status> try_create(
      int rank, const int* extents) noexcept {
    // arg checking
    if ((rank < 1) || (rank > MAX_RANK)) {
      return {IdxMapping<Layout>(), error::Adhoc("rank is invalid")};
    } else if (extents == nullptr) {
      return {IdxMapping<Layout>(), error::Adhoc("extents is a nullptr")};
    }
    for (int i = 0; i < rank; i++) {
      if (extents[i] < 1) {
        return {IdxMapping<Layout>(),
                error::Adhoc("extents must hold positive vals")};
      }
    }

    // build and return the mapping
    IdxMapping<Layout> mapping;
    mapping.rank_ = rank;
    for (int i = 0; i < rank; i++) {
      mapping.extents_[i] = extents[i];
    }
    return {mapping, OkStatus()};
  }

  explicit IdxMapping(int extent0) {
    // for less experienced C++ devs: the `explicit` kwarg prevents the use of
    // this constructor for implicit casts
    int extents[1] = {extent0};
    *this = create_or_abort_(1, extents);
  }

  IdxMapping(int extent0, int extent1) {
    int extents[2] = {extent0, extent1};
    *this = create_or_abort_(2, extents);
  }

  IdxMapping(int extent0, int extent1, int extent2) {
    int extents[3] = {extent0, extent1, extent2};
    *this = create_or_abort_(3, extents);
  }

  IdxMapping(const IdxMapping<Layout>&) = default;
  IdxMapping(IdxMapping<Layout>&&) = default;
  IdxMapping<Layout>& operator=(const IdxMapping<Layout>&) = default;
  IdxMapping<Layout>& operator=(IdxMapping<Layout>&&) = default;
  ~IdxMapping() = default;

  /// access the rank
  int rank() const noexcept { return rank_; }

  /// access the extents pointer (the length is given by @ref rank)
  const int* extents() const noexcept { return extents_; }

  /// construct an equivalent 3d IdxMapping
  IdxMapping<Layout> to_3d_mapping() const noexcept {
    int out_rank = 3;
    int rank_diff = out_rank - this->rank_;
    GR_INTERNAL_REQUIRE(rank_diff >= 0, "current rank exceeds new rank");

    IdxMapping<Layout> out;
    out.rank_ = out_rank;
    for (int i = 0; i < out_rank; i++) {
      out.extents_[i] = 1;
    }
    int offset = (Layout == DataLayout::RIGHT) ? rank_diff : 0;
    for (int i = out; i < this->rank_; i++) {
      out.extents_[i + offset] = this->extents_[i];
    }
    return out;
  }

  static constexpr bool is_contiguous() { return true; }

  int n_elements() const noexcept {
    int product = 1;
    for (int i = 0; i < this->rank_; i++) {
      product *= this->extents_[i];
    }
    return product;
  }

  /** @{ */  // <- open the group of member functions with a shared docstring
  /// compute the 1D pointer offset associated with the multidimensional index
  ///
  /// @note
  /// For less experienced C++ devs, these methods overloads the "function call
  /// operator". In python they would be named `__call__(self, ...)`
  ///
  /// Behavior is undefined if the number of arguments doesn't match the value
  /// returned by `this->rank()`
  [[gnu::always_inline]] int operator()(int i) const noexcept {
    assert(this->rank_ == 1);
    return i;
  }

  [[gnu::always_inline]] int operator()(int i, int j) const noexcept {
    assert(this->rank_ == 2);
    if constexpr (Layout == DataLayout::LEFT) {
      return i + this->extents_[0] * j;
    } else {  // Layout == DataLayout::RIGHT
      return j + this->extents_[1] * i;
    }
  }

  [[gnu::always_inline]] int operator()(int i, int j, int k) const noexcept {
    assert(this->rank_ == 3);
    if constexpr (Layout == DataLayout::LEFT) {
      return i + this->extents_[0] * (j + k * this->extents_[1]);
    } else {  // Layout == DataLayout::RIGHT
      return k + this->extents_[2] * (j + i * this->extents_[1]);
    }
  }
  /** @} */  // <- close the group of member functions with a shared docstring

  /// Convert a 1D pointer offset to a multidimensional index
  ///
  /// The number of components for the output index is given by calling the
  /// @ref rank method. Behavior is undefined if @p out doesn't have enough
  /// space.
  ///
  /// @param[in] offset The 1D pointer offset (must be non-negative)
  /// @param[out] out Buffer where result is stored
  void offset_to_md_idx(int offset, int* out) const noexcept {
    assert(out != nullptr);

    const int contig_ax = (Layout == DataLayout::LEFT) ? 0 : (this->rank_ - 1);
    const int slowest_ax = (Layout == DataLayout::LEFT) ? (this->rank_ - 1) : 0;

    switch (this->rank_) {
      case 1:
        out[0] = offset;
        return;
      case 2:
        out[slowest_ax] = offset / this->extents_[contig_ax];
        out[contig_ax] = offset % this->extents_[contig_ax];
        return;
      case 3: {
        int largest_stride = this->extents_[contig_ax] * this->extents_[1];
        out[slowest_ax] = offset / largest_stride;
        int remainder = offset % largest_stride;
        out[1] = remainder / this->extents_[contig_ax];
        out[contig_ax] = remainder % this->extents_[contig_ax];
        return;
      }
      default:
        GR_INTERNAL_ERROR("should be unreachable");
    }
  }

  // teach googletest how to print this type
  friend void PrintTo(const IdxMapping<Layout>& mapping, std::ostream* os) {
    const char* layout =
        (Layout == DataLayout::LEFT) ? "DataLayout::LEFT" : "DataLayout::RIGHT";
    *os << "IdxMapping<" << layout << ">(";
    for (int i = 0; i < mapping.rank_; i++) {
      *os << mapping.extents_[i];
      *os << ((i + 1) < mapping.rank_ ? ',' : ')');
    }
  }
};

/// equivalent of converting output of printf("%g", val) to std::string
std::string to_pretty_string(float val);
std::string to_pretty_string(double val);

/// formats a pointer as a string
std::string ptr_to_string(const float* ptr, std::size_t len);
std::string ptr_to_string(const double* ptr, std::size_t len);

/// formats a pointer as a string (current implementation prints a 1D array)
template <DataLayout Layout>
std::string ptr_to_string(const float* ptr, IdxMapping<Layout> idx_mapping) {
  static_assert(IdxMapping<Layout>::is_contiguous());
  return ptr_to_string(ptr, idx_mapping.n_elements());
}

/// formats a pointer as a string (current implementation prints a 1D array)
template <DataLayout Layout>
std::string ptr_to_string(const double* ptr, IdxMapping<Layout> idx_mapping) {
  static_assert(IdxMapping<Layout>::is_contiguous());
  return ptr_to_string(ptr, idx_mapping.n_elements());
}

}  // namespace grtest

#endif  // GRTESTUTILS_VIEW_HPP