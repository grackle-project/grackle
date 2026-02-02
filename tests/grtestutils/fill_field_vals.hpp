//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// define fill_field_vals
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_FILL_FIELD_VALS_HPP
#define GRTESTUTILS_FILL_FIELD_VALS_HPP

#include "GrackleCtxPack.hpp"
#include "grackle.h"
#include "status_reporting.h"
#include "./field_container.hpp"
#include "./field_info_detail.hpp"
#include "./status.hpp"

#include <optional>
#include <string_view>
#include <vector>

namespace grtest {

/// @defgroup fillfieldsgrp Logic for filling field values
///
/// This group of entities defines the logic for initializing values in a
/// @ref FluidContainer.
///
/// Since our testing tools may be used to test (or benchmark) scenarios with
/// an arbitrarily large number of elements, we adopt a general-purpose scheme
/// where we essentially repeat a "tile" of field values 1 or more times.
///
/// In general, a "tile" refers to an instance of a generic type that implements
/// the methods illustrated in the following snippet.
///
/// @code{C++}
/// struct MySampleTileType {
///   /// number of elements in the tile
///   int get_n_vals() const noexcept;
///
///   // iF field_name is known, fill buf with associated values associated
///   // and return true. Otherwise, return false.
///   //
///   // buf should have a length given by get_n_vals()
///   bool fill_buf(gr_float* buf, std::string_view field_name) const noexcept;
/// };
/// @endcode
///
/// This approach was picked to maximize flexibility.
/// - it can work for initializing a single value or an arbitrary number of
///   values. This reduces the maintenance burden (at the cost of some
///   performance)
/// - it should be easy enough to implement a new tile type
/// - it is conceivable that we could reuse this machinery for initializing
///   a field from serialized values (maybe from a json or toml or hdf5 file)
/**@{*/  // open the doxygen group

/// Fill in the values of @p fc by repeating the values from @p field_tile
///
/// @param[in] field_tile Specifies the values to use
/// @param[in, out] fc The field container that will be filled
template <class Tile>
Status fill_field_vals(const Tile& field_tile, FieldContainer& fc) {
  // access field grid index properties
  const GridLayout& layout = fc.grid_layout();
  int rank = layout.rank();
  int ix_start = layout.start()[0];
  int ix_stop = layout.stop()[0];
  int iy_start = (rank >= 2) ? layout.start()[1] : 0;
  int iy_stop = (rank >= 2) ? layout.stop()[1] : 1;
  int iz_start = (rank == 3) ? layout.start()[2] : 0;
  int iz_stop = (rank == 3) ? layout.stop()[2] : 1;

  int mx = fc.get_ptr()->grid_dimension[0];
  int my = (rank >= 2) ? fc.get_ptr()->grid_dimension[1] : 1;

  int n_x_vals = ix_stop - ix_start;

  // access properties about the tile
  int tile_len = field_tile.get_n_vals();

  // a few error consistency checks
  if (tile_len <= 0) {
    return error::Adhoc("tile_len must be positive");
  } else if (tile_len > n_x_vals) {
    // we may want to revisit this in the future
    return error::Adhoc("tile_len can't exceed length of contiguous axis");
  }

  // allocate the buffer
  std::vector<gr_float> tile_buf(tile_len);

  for (auto [field_name, field_ptr] : fc) {
    if (!field_tile.fill_buf(tile_buf.data(), field_name)) {
      return error::Adhoc("field_tile doesn't recognize field name: " +
                          field_name);
    }

    for (int iz = iz_start; iz < iz_stop; iz++) {
      for (int iy = iy_start; iy < iy_stop; iy++) {
        for (int ix = ix_start; ix < ix_stop; ix++) {
          int i = ix + mx * (iy + my * iz);
          field_ptr[i] = tile_buf[ix % tile_len];
        }
      }
    }
  }

  return OkStatus();
}

/// Create a @ref FieldContainer and fill in values by repeating the values from
/// @p field_tile
///
/// @param field_tile Specifies the values to use
/// @param ctx_pack The Grackle Configuration used for initialization
/// @param layout The Grid layout to use
/// @param exclude_fields Field names that should be excluded during creation
///
/// This function is simple: it just calls @ref FieldContainer::create and
/// @ref fill_field_vals. But because this sequence of events is relatively
/// common, the aggregation of @ref Status instance is quite convenient.
template <class Tile>
std::pair<FieldContainer, Status> create_and_fill_FieldContainer(
    const Tile& field_tile, const GrackleCtxPack& ctx_pack,
    const GridLayout& layout,
    const std::set<std::string>& exclude_fields = {}) {
  std::pair<grtest::FieldContainer, grtest::Status> out =
      grtest::FieldContainer::create(ctx_pack, layout, exclude_fields);
  if (out.second.is_ok()) {
    out.second = fill_field_vals(field_tile, out.first);
  }
  return out;
}

/// Represents a very simplistic set of conditions
struct SimpleFieldTile {
  double mfrac_metal;
  double mfrac_H;
  double mfrac_He;
  double mfrac_D;

  // both of these are in code-units
  double common_density;
  double common_eint;

  /// number of values specified by the tile
  int get_n_vals() const noexcept { return 1; }

  /// fill @p buf with values associated with the specified @p field_name
  ///
  /// @param[out] buf This is filled. The length is given by @ref get_n_vals()
  /// @param[in] field_name Name of the field
  /// @returns `true` if the field is known. Otherwise, returns `false`
  bool fill_buf(gr_float* buf, std::string_view field_name) const noexcept {
    if (field_name == "density") {
      buf[0] = common_density;
    } else if (field_name == "internal_energy") {
      buf[0] = common_eint;
    } else if (field_name == "metal_density") {
      buf[0] = mfrac_metal * common_density;
    } else if (field_name == "HI_density") {
      buf[0] = mfrac_H * common_density;
    } else if (field_name == "HeI_density") {
      buf[0] = mfrac_He * common_density;
    } else if (field_name == "DI_density") {
      buf[0] = mfrac_D * common_density;
    } else {
      field_detail::Kind kind = field_detail::Kind::UNSET;
      std::optional<field_detail::FieldInfo> maybe_field_info =
          field_detail::query_field_info(field_name);
      if (maybe_field_info.has_value()) {
        kind = maybe_field_info->kind;
      }

      if (kind == field_detail::Kind::ELEMENTWISE_SOLVER_ARG) {
        buf[0] = 0.0;
      } else if (kind == field_detail::Kind::PRIMORDIAL_SPECIES) {
        gr_float tiny_number = 1.e-20;
        buf[0] = tiny_number * common_density;
      } else {
        return false;
      }
    }
    return true;
  }
};

inline SimpleFieldTile make_simple_tile(const GrackleCtxPack& ctx_pack,
                                        code_units units, double metallicity) {
  GR_INTERNAL_REQUIRE(ctx_pack.is_initialized(), "ctx_pack was initialized");

  const chemistry_data* my_chem = ctx_pack.my_chemistry();

  // the internal energy roughly corresponds to 1000 K
  double common_eint = 1000.0 / get_temperature_units(&units);

  return SimpleFieldTile{
      /* mfrac_metal = */ metallicity * my_chem->SolarMetalFractionByMass,
      /* mfrac_H = */ my_chem->HydrogenFractionByMass,
      /* mfrac_He = */ (1.0 - my_chem->HydrogenFractionByMass),
      /* mfrac_D = */ 2.0 * 3.4e-5,

      // both of these are in code-units
      /* common_density = */ 1.0,
      /* common_eint = */ common_eint};
}

/**@}*/  // close the doxygen group

}  // namespace grtest

#endif  // GRTESTUTILS_SET_FIELD_CONDITIONS_HPP
