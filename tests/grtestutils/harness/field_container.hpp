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

#include <functional>  // std::less
#include <map>
#include <memory>  // std::map
#include <optional>
#include <string>
#include <string_view>
#include <utility>  // std::pair

namespace grtest {

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
///   non-overlapping segment of the buffer tracked by @ref prop_buf
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
  std::unique_ptr<int[]> prop_buf;
  /// the buffer used to hold all field data
  std::unique_ptr<gr_float[]> data_buf;
  /// the struct understood by Grackle
  std::unique_ptr<grackle_field_data> my_fields;
  /// maps field names to pointers
  MapType map;
};

}  // namespace field_detail

/// Encapsulates an arbitrary GridLayout
struct GridLayout {
  int rank;
  int dimension[3];
  int ghost_depth[3];
};

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
  static std::pair<FieldContainer, Status> create_1d(
      const GrackleCtxPack& ctx_pack, int buf_size, bool disable_metal = false);

  /// A factory method to make a container
  static std::pair<FieldContainer, Status> create(
      const GrackleCtxPack& ctx_pack, const GridLayout& layout,
      bool disable_metal = false);

  /// Create a clone of FieldContainer
  FieldContainer clone() const;

  /// Overwrite the field data in @p dest with the field data from `this`
  ///
  /// This doesn't perform any allocations.
  ///
  /// @warning
  /// Setting bypass_check to `true` is risky. It primarily exists for the
  /// case where you call this method in a loop
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
  bool same_grid_props(const FieldContainer& other) const;

  /// returns whether `this` and @p other contains the same set of fields
  bool same_fields(const FieldContainer& other) const;

  const grackle_field_data* get_ptr() const { return data_.my_fields.get(); }
  grackle_field_data* get_ptr() { return data_.my_fields.get(); }

  std::optional<gr_float*> find(std::string_view key) {
    auto s = data_.map.find(key);
    return (s != data_.map.end()) ? std::optional{s->second} : std::nullopt;
  }

  std::optional<const gr_float*> find(std::string_view key) const {
    auto s = data_.map.find(key);
    return (s != data_.map.end()) ? std::optional{s->second} : std::nullopt;
  }

  int n_fields() const { return static_cast<int>(data_.map.size()); }

  int rank() const { return data_.my_fields->grid_rank; }

  /// the number of elements per field (unaffected by ghost zones)
  int elements_per_field() const {
    return field_detail::elements_per_field_(*data_.my_fields);
  }

  // we define both of the following methods to support the writing of
  // range-based for-loops to iterate over key-buffer pairs
  MapType::const_iterator begin() const { return data_.map.begin(); }
  MapType::const_iterator end() const { return data_.map.end(); }
};
}  // namespace grtest

#endif  // GRTESTUTILS_FIELD_CONTAINER_HPP
