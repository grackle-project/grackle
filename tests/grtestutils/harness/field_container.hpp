//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare the FieldContainer class
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_FIELD_CONTAINER_HPP
#define GRTESTUTILS_FIELD_CONTAINER_HPP

#include "grackle.h"

#include "./grackle_ctx_pack.hpp"
#include "./status.hpp"
#include "status_reporting.h"

#include <functional>  // std::less
#include <iosfwd>
#include <map>
#include <memory>  // std::map
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <utility>  // std::pair

namespace grtest {

namespace field_detail {
struct CorePack;
}  // namespace field_detail

/// Represents Grid Properties
///
/// This type is probably a little over-engineered
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
//       `0` up to `a` in the outermost loop and from `0` up to `c`
///      in the innermost loop
///    - this is the "natural layout" for a multidimensional c-style array
///      `arr[a][b][c]`
/// 2. layout-left: where the stride is `1` along the leftmost extent
///    - for extents `{a,b,c}`, an optimal nested for-loop will iterates from
//       `0` up to `c` in the outermost loop and from `0` up to `a`
///      in the innermost loop
///    - this is the "natural layout" for a multidimensional fortran array ()
///      `arr[a][b][c]`
///
/// > Aside: More sophisticated layouts are possible when strides along each
/// > axis aren't directly tied to extents (this comes up when making subviews).
///
/// About this type
/// ---------------
/// This type specifies array extents for a layout-left mapping. It also
/// specifies the region of valid values. For context, arrays of
/// fluid-quantities in mesh-based hydro codes have an outer layer "ghost zones"
/// (Computer Scientists sometimes call this a "halo") that may not contain
/// valid values when calling Grackle.
///
/// > Aside: We only describe the concept of "ghost zones" because we want to
/// > to support tests for validating that ghost zones aren't modified.
/// > Otherwise, we could reframe this description in terms of strided layouts.
class GridLayout {
  friend field_detail::CorePack;

  /// the number of dimensions
  int rank_;
  /// number of elements along an axis
  int dim_[3];

  /// the first active-zone index along an axis
  int start_[3];

  /// the last active-zone index along an axis
  ///
  /// @note
  /// We track this because Grackle natively understands this
  int end_[3];

  /// elements are one larger than the corresponding value in @ref end
  ///
  /// @note
  /// This is the more natural than @ref end for 0-based indexing
  int stop_[3];

  // only invoked by factory methods
  GridLayout() = default;

  static std::pair<GridLayout, Status> create_(int rank, const int* dim,
                                               const int* start,
                                               const int* stop) {
    GridLayout out;
    out.rank_ = rank;
    if ((rank < 1) || (rank > 3)) {
      return {out, error::Adhoc("rank is invalid")};
    } else if (dim == nullptr) {
      return {out, error::Adhoc("dim is a nullptr")};
    }

    for (int i = 0; i < 3; i++) {
      out.dim_[i] = (i < rank) ? dim[i] : 1;
      out.start_[i] = (start != nullptr && i < rank) ? start[i] : 0;
      out.stop_[i] = (stop != nullptr && i < rank) ? stop[i] : out.dim_[i];
      out.end_[i] = out.stop_[i] - 1;

      if (out.dim_[i] < 1) {
        return {out, error::Adhoc("dim must hold positive vals")};
      } else if (out.start_[i] < 0) {
        return {out, error::Adhoc("start must hold non-negative vals")};
      } else if (out.stop_[i] <= out.start_[i]) {
        return {out, error::Adhoc("stop must exceed start")};
      } else if (out.stop_[i] > out.dim_[i]) {
        return {out, error::Adhoc("stop must not exceed dim")};
      }
    }

    return {out, OkStatus()};
  }

public:
  /// factory method
  ///
  /// This is for the common case where we want a 1D layout with a hardcoded
  /// number of entries (and we know at compile-time that it can't fail)
  template <int N>
  static GridLayout create_1d() noexcept {
    static_assert(N >= 1, "N must be positive");
    int dim[3] = {N, 0, 0};
    return GridLayout::create_(1, dim, nullptr, nullptr).first;
  }

  /// factory method
  static std::pair<GridLayout, Status> try_from_dim(int rank, const int* dim) {
    return GridLayout::create_(rank, dim, nullptr, nullptr);
  }

  /// factory method (aborts for invalid arguments)
  static GridLayout from_dim(int rank, const int* dim) {
    std::pair<GridLayout, Status> tmp = GridLayout::try_from_dim(rank, dim);
    if (tmp.second.is_err()) {
      std::string msg = tmp.second.to_string();
      GR_INTERNAL_ERROR("%s", msg.c_str());
    }
    return tmp.first;
  }

  /// factory method (aborts for invalid arguments)
  static GridLayout from_ghostdepth_and_dims(int rank, const int* ghostdepth,
                                             const int* dim) {
    GR_INTERNAL_REQUIRE(rank == 1 || rank == 2 || rank == 3, "rank is invalid");
    GR_INTERNAL_REQUIRE(ghostdepth != nullptr, "ghostdepth is a nullptr");
    GR_INTERNAL_REQUIRE(dim != nullptr, "dim is a nullptr");

    int start[3] = {0, 0, 0};
    int stop[3] = {1, 1, 1};
    for (int i = 0; i < rank; i++) {
      GR_INTERNAL_REQUIRE(ghostdepth[i] >= 0, "ghostdepth must be nonnegative");
      GR_INTERNAL_REQUIRE(2 * ghostdepth[i] < dim[i],
                          "dim[i] must exceed 2*ghostdepth[i]")
      start[i] = ghostdepth[i];
      stop[i] = dim[i] - ghostdepth[i];
    }

    std::pair<GridLayout, Status> tmp =
        GridLayout::create_(rank, dim, start, stop);
    if (tmp.second.is_err()) {
      std::string msg = tmp.second.to_string();
      GR_INTERNAL_ERROR("%s", msg.c_str());
    }
    return tmp.first;
  }

  int rank() const noexcept { return rank_; }
  const int* dim() const noexcept { return dim_; }
  const int* start() const noexcept { return start_; }
  const int* stop() const noexcept { return stop_; }
  const int* end() const noexcept { return end_; }

  friend bool operator==(const GridLayout& lhs, const GridLayout& rhs) noexcept;
  friend bool operator!=(const GridLayout& lhs,
                         const GridLayout& rhs) noexcept {
    return !(lhs == rhs);
  }

  int n_elements(bool exclude_inactive = false) const noexcept {
    int total = 1;
    for (int i = 0; i < rank_; i++) {
      total *= (exclude_inactive) ? stop_[i] - start_[i] : dim_[i];
    }
    return total;
  }

  // teach googletest how to print this type
  friend void PrintTo(const GridLayout& layout, std::ostream* os);

  /// create a string representation
  std::string to_string() const;
};

namespace field_detail {

inline int elements_per_field_(const grackle_field_data& my_fields) {
  int total = 1;
  for (int i = 0; i < my_fields.grid_rank; i++) {
    total *= my_fields.grid_dimension[i];
  }
  return total;
}

// by including std::less<> in the configuration, we ensure that we can perform
// a lookup with a string
using MapType = std::map<std::string, gr_float*, std::less<>>;

/// Holds the core data of a FieldContainer
///
/// @par Implementation Notes
/// This holds a @ref grackle_field_data instance, associated buffers, and a
/// string-buffer mapping.
///
/// Consider a case where @ref my_fields is fully initialized:
/// - each pointer data member for specifying grid properties refers to a
///   memory tracked by @ref grid_layout
/// - each pointer data member for specifying an unused grackle-field is holds
///   a nullptr
/// - each pointer data member for specifying an actively used grackle-field:
///   - holds an associated pointer to a non-overlapping segment of the buffer
///     tracked by @ref data_buf. The numer of elements per segment is given by
///     `elements_per_field_(*my_fields)`.
///   - has a corresponding key-value pair tracked by @ref map, where the key
///     is the name of the grackle-field, and the value is associated pointer
struct CorePack {
  /// the buffer used for holding field shape and ghost zones
  std::unique_ptr<GridLayout> layout;
  /// the buffer used to hold all field data
  std::unique_ptr<gr_float[]> data_buf;
  /// the struct understood by Grackle
  std::unique_ptr<grackle_field_data> my_fields;
  /// maps field names to pointers
  MapType map;

  /// prefer this method over directly modifying layout
  ///
  /// @note This obviously requires that this->layout and this->my_fields are
  /// not nullptr.
  void override_layout(const GridLayout& new_layout) noexcept {
    // overriding this->layout implicitly updates the members of `my_fields`,
    // `my_fields->grid_(dimension|start|end)`, because these members all point
    // to statically sized C array members of GridLayout
    *this->layout = new_layout;
    // make sure grid_rank remains up to date
    this->my_fields->grid_rank = this->layout->rank();
  }

  /// factory method that consumes a @p premade_map
  ///
  /// @param premade_map A string to pointer mapping. The keys of this argument
  ///     should specify the names for each desired Grackle field. We assume
  ///     that pointers associated with each key hold meaningless garbage values
  /// @param buf_size The number of elements to allocate per field
  static std::pair<CorePack, Status> setup_1d(MapType&& premade_map,
                                              int buf_size);
};

}  // namespace field_detail

/// A container wrapping grackle_field_data that owns the underlying data
///
/// @par Implementation Notes
/// The current implementation attempts to be a general-purpose test-harness
/// that can be used for testing and for benchmarking. At the moment,
/// construction is a little slow. That's ok for benchmarking since @ref
/// copy_into is fast (it's effectively just a single memcpy)
///
/// @note
/// At the moment, I have explicitly avoided addressing grid_dx
class FieldContainer {
  field_detail::CorePack data_;

  /// Can only be invoked by a factory method
  FieldContainer() = default;

  void copy_into_helper_(FieldContainer& other) const;

public:
  using MapType = field_detail::MapType;

  FieldContainer(const FieldContainer& other) = delete;
  FieldContainer& operator=(const FieldContainer& other) = delete;
  FieldContainer(FieldContainer&& other) = default;
  FieldContainer& operator=(FieldContainer&& other) = default;
  ~FieldContainer() = default;

  /// A factory method to make a simple 1d container
  ///
  /// @param ctx_pack The Grackle Configuration used for initialization
  /// @param buf_size The positive number of elements in the container
  /// @param exclude_fields Names of fields that should be excluded
  ///
  /// @note This is provided as a convenience. Do we really need it?
  static std::pair<FieldContainer, Status> create_1d(
      const GrackleCtxPack& ctx_pack, int buf_size,
      const std::set<std::string>& exclude_fields = {});

  /// A factory method to make a container
  ///
  /// @param ctx_pack The Grackle Configuration used for initialization
  /// @param layout The Grid layout to use
  /// @param exclude_fields Names of fields that should be excluded
  static std::pair<FieldContainer, Status> create(
      const GrackleCtxPack& ctx_pack, const GridLayout& layout,
      const std::set<std::string>& exclude_fields = {});

  /// Create a clone of FieldContainer
  FieldContainer clone() const;

  /// Overwrite the field data in @p dest with the field data from `this`
  ///
  /// This doesn't perform any allocations.
  ///
  /// @warning
  /// Setting bypass_check to `true` is risky. It primarily exists for the
  /// case where you call this method in a loop and we already know that 2
  /// containers are compatible
  Status copy_into(FieldContainer& dest, bool bypass_check = false) const {
    if (!(bypass_check || this->same_grid_props(dest))) {
      return error::Adhoc("grid properties are incompatible");
    } else if (!(bypass_check || this->same_fields(dest))) {
      return error::Adhoc("field sets are incompatible");
    }
    this->copy_into_helper_(dest);
    return OkStatus();
  }

  /// returns whether `this` and @p other contains the same grid properties
  bool same_grid_props(const FieldContainer& other) const {
    return (*this->data_.layout) == (*other.data_.layout);
  }

  /// returns whether `this` and @p other contains the same set of fields
  bool same_fields(const FieldContainer& other) const;

  /**@{*/
  /// get the pointer to the wrapped @ref grackle_field_data instance
  ///
  /// @note
  /// This primarily exists to support Grackle API function. Avoid mutating
  /// the returned pointer
  const grackle_field_data* get_ptr() const { return data_.my_fields.get(); }

  // NOLINTNEXTLINE(readability-make-member-function-const)
  grackle_field_data* get_ptr() { return data_.my_fields.get(); }
  /**@}*/

  /**@{*/
  /// finds the field data pointer for the field with the specified name
  std::optional<gr_float*> find(std::string_view key) {
    auto s = data_.map.find(key);
    return (s != data_.map.end()) ? std::optional{s->second} : std::nullopt;
  }

  std::optional<const gr_float*> find(std::string_view key) const {
    auto s = data_.map.find(key);
    return (s != data_.map.end()) ? std::optional{s->second} : std::nullopt;
  }
  /**@}*/

  /**@{*/
  /// finds the field data pointer or aborts the program
  gr_float* get_or_abort(std::string_view key) {
    std::optional<gr_float*> tmp = find(key);
    if (!tmp.has_value()) {
      std::string key_copy{key};  // required b/c key isn't '\0' terminated
      GR_INTERNAL_ERROR("\"%s\" field wasn't found", key_copy.c_str());
    }
    return *tmp;
  }

  const gr_float* get_or_abort(std::string_view key) const {
    std::optional<const gr_float*> tmp = find(key);
    if (!tmp.has_value()) {
      std::string key_copy{key};  // required b/c key isn't '\0' terminated
      GR_INTERNAL_ERROR("\"%s\" field wasn't found", key_copy.c_str());
    }
    return *tmp;
  }
  /**@}*/

  const GridLayout& grid_layout() const { return *data_.layout; }

  int n_fields() const { return static_cast<int>(data_.map.size()); }

  // we define both of the following methods to support the writing of
  // range-based for-loops to iterate over key-buffer pairs
  MapType::const_iterator begin() const { return data_.map.begin(); }
  MapType::const_iterator end() const { return data_.map.end(); }

  // teach googletest how to print this type
  friend void PrintTo(const FieldContainer& fc, std::ostream* os);
};
}  // namespace grtest

#endif  // GRTESTUTILS_FIELD_CONTAINER_HPP
