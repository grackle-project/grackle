//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement logic pertaining to @ref FieldContainer
///
//===----------------------------------------------------------------------===//

#include "./GrackleCtxPack.hpp"
#include "./field_container.hpp"
#include "./field_info_detail.hpp"
#include "./status.hpp"

#include "grackle.h"
#include "status_reporting.h"

#include <cstring>  // std::memcmp, std::memcpy
#include <utility>  // std::pair, std::move

namespace grtest {
namespace field_detail {

/// Construct a name-pointer mapping based on @p ctx_pack
///
/// The returned where the keys are the names of every active Grackle field
/// and the associated values are all nullptr
static MapType make_nullptr_map_(const GrackleCtxPack& ctx_pack,
                                 bool exclude_metal) {
  MapType m;
  // fill up m with (field, nullptr) pairs for each field enabled by ctx_pack
  auto fn = [&m](const char* name, const FieldInfo& info) {
    if (!m.emplace(name, nullptr).second) {
      GR_INTERNAL_ERROR("%s was inserted multiple times", name);
    }
  };
  for_each_named_field(ctx_pack, fn);

  if (exclude_metal) {
    auto search = m.find("metal_density");
    if (search != m.end()) {
      m.erase(search);
    }
  }
  return m;
}

/// Construct a `CorePack` instance by consuming @p premade_map
///
/// @param premade_map A string to pointer mapping. The keys of this argument
///     should specify the names for each desired Grackle field. We assume that
///     the pointers associated with each key hold meaningless garbage values
/// @param buf_size The number of elements to allocate per field
static std::pair<CorePack, Status> setup_1d_CorePack(MapType&& premade_map,
                                                     int buf_size) {
  if (buf_size <= 0) {
    return {CorePack{}, error::Adhoc("buf_size must be positive")};
  }
  std::size_t data_sz = premade_map.size() * buf_size;

  CorePack pack{
      // the trailing parentheses when allocating an array of integers or
      // floating point values sets the initial values to 0
      /* prop_buf */ std::unique_ptr<int[]>(new int[9]()),
      /* data_buf */ std::unique_ptr<gr_float[]>(new gr_float[data_sz]()),
      /* my_fields */
      std::unique_ptr<grackle_field_data>(new grackle_field_data),
      /* map */ std::move(premade_map)};
  if (gr_initialize_field_data(pack.my_fields.get()) != GR_SUCCESS) {
    // this really should never fail
    return {CorePack{}, error::Adhoc("gr_initialize_field_data failed")};
  }

  // setup grid properties
  pack.my_fields->grid_dimension = pack.prop_buf.get();
  pack.my_fields->grid_start = pack.prop_buf.get() + 3;
  pack.my_fields->grid_end = pack.prop_buf.get() + 6;

  pack.my_fields->grid_rank = 1;
  pack.my_fields->grid_dimension[0] = buf_size;
  pack.my_fields->grid_end[0] = buf_size - 1;

  // distribute allocated field-memory among all the appropriate slots of
  // pack.my_fields and pack.map
  gr_float* ptr = pack.data_buf.get();
  int counter = 0;
  auto fn = [ptr, buf_size, &counter, &pack](const char* name,
                                             const FieldInfo& finfo) {
    // return immediately if finfo.name isn't a key of `pack.map`
    MapType::iterator search = pack.map.find(name);
    if (search == pack.map.end()) {
      return;  // finfo.name isn't a key of `map`
    }
    // get the pointer to the memory reserved for the current field
    gr_float* cur_ptr = ptr + (counter * buf_size);
    counter++;
    // update `pack.my_fields` to associate the current field with `cur_ptr`
    pack.my_fields.get()->*finfo.relative_addr = cur_ptr;
    // update `map` to associate the current field with `cur_ptr`
    search->second = cur_ptr;
  };
  // this acts like a for-loop that passes a FieldInfo struct for every known
  // grackle-field into `fn`
  for_each_named_field(fn);
  return {std::move(pack), OkStatus()};
}

}  // namespace field_detail

std::pair<FieldContainer, Status> FieldContainer::create_1d(
    const GrackleCtxPack& ctx_pack, int buf_size, bool disable_metal) {
  if (!ctx_pack.is_initialized()) {
    return {FieldContainer{}, error::Adhoc("ctx_pack isn't initialized")};
  }

  // construct a map for each relevant field (the values are all nullptr)
  field_detail::MapType m =
      field_detail::make_nullptr_map_(ctx_pack, disable_metal);
  std::pair<field_detail::CorePack, Status> tmp =
      field_detail::setup_1d_CorePack(std::move(m), buf_size);
  if (tmp.second.is_err()) {
    return {FieldContainer(), std::move(tmp.second)};
  } else {
    FieldContainer fc;
    fc.data_ = std::move(tmp.first);
    return {std::move(fc), OkStatus()};
  }
}

std::pair<FieldContainer, Status> FieldContainer::create(
    const GrackleCtxPack& ctx_pack, const GridLayout& layout,
    bool disable_metal) {
  if (layout.rank < 1 || layout.rank > 3) {
    return {FieldContainer(), error::Adhoc("layout.rank must be 1, 2, or 3")};
  }
  int total_count = 1;
  for (int i = 0; i < layout.rank; i++) {
    if (layout.dimension[i] <= 0) {
      return {FieldContainer(),
              error::Adhoc("layout.dimension[i] must be positive")};
    } else if (layout.ghost_depth[i] < 0) {
      return {FieldContainer(),
              error::Adhoc("layout.ghost_depth[i] is negative")};
    } else if (layout.ghost_depth[i] * 2 >= layout.dimension[i]) {
      return {FieldContainer(),
              error::Adhoc(
                  "layout.dimension[i] must exceed 2*layout.ghost_depth[i]")};
    }
    total_count *= layout.dimension[i];
  }

  std::pair<FieldContainer, Status> tmp =
      FieldContainer::create_1d(ctx_pack, total_count, disable_metal);

  if (tmp.second.is_ok()) {
    grackle_field_data* my_fields = tmp.first.data_.my_fields.get();
    tmp.first.data_.my_fields->grid_rank = layout.rank;
    for (int i = 0; i < layout.rank; i++) {
      int dim = layout.dimension[i];
      int ghost_depth = layout.ghost_depth[i];
      my_fields->grid_dimension[i] = dim;
      my_fields->grid_start[i] = ghost_depth;
      my_fields->grid_end[i] = dim - ghost_depth - 1;
    }
  }
  return tmp;
}

FieldContainer FieldContainer::clone() const {
  // make a copy of m, and use it construct a new CorePack
  // -> the keys of `m` tell `setup_1d_CorePack` which Grackle fields get used
  // -> right after we create
  field_detail::MapType m = this->data_.map;
  std::pair<field_detail::CorePack, Status> tmp =
      field_detail::setup_1d_CorePack(std::move(m), this->elements_per_field());
  // it shouldn't be possible for tmp.second.is_err() to return true
  FieldContainer out;
  out.data_ = std::move(tmp.first);

  // copy over layout properties
  int rank = this->rank();
  for (int i = 0; i < rank; i++) {
    out.data_.my_fields->grid_start[i] = this->data_.my_fields->grid_start[i];
    out.data_.my_fields->grid_end[i] = this->data_.my_fields->grid_end[i];
    out.data_.my_fields->grid_dimension[i] =
        this->data_.my_fields->grid_dimension[i];
  }

  copy_into_helper_(out);
  return out;
}

void FieldContainer::copy_into_helper_(FieldContainer& other) const {
  // I'm still not entirely sure how to best handle grid_dx, but this is here
  // to make sure we don't forget to propagate the logic when we ultimately
  // make a decision
  if (this->get_ptr()->grid_dx != other.get_ptr()->grid_dx) {
    GR_INTERNAL_ERROR("did you forget to update grid_dx handling?");
  }

  const gr_float* src = this->data_.data_buf.get();
  gr_float* dst = other.data_.data_buf.get();
  int length = this->elements_per_field() * this->n_fields();
  std::memcpy(dst, src, length * sizeof(gr_float));
}

bool FieldContainer::same_grid_props(const FieldContainer& other) const {
  const grackle_field_data* a = this->data_.my_fields.get();
  const grackle_field_data* b = other.data_.my_fields.get();

  if (a->grid_rank != b->grid_rank) {
    return false;
  }
  int rank = a->grid_rank;
  std::size_t sz = rank * sizeof(int);
  bool s = std::memcmp(a->grid_start, b->grid_start, sz) == 0;
  bool e = std::memcmp(a->grid_end, b->grid_end, sz) == 0;
  bool d = std::memcmp(a->grid_dimension, b->grid_dimension, sz) == 0;
  return s && d && e;
}

bool FieldContainer::same_fields(const FieldContainer& other) const {
  // this takes advantage of the fact that MapType is sorted
  MapType::const_iterator it_a = this->data_.map.begin();
  MapType::const_iterator stop_a = this->data_.map.end();
  MapType::const_iterator it_b = other.data_.map.begin();

  for (; it_a != stop_a; ++it_a, ++it_b) {
    if (it_a->first != it_b->first) {
      return false;
    }
  }
  return true;
}

}  // namespace grtest