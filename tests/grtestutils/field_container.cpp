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
#include "./utils.hpp"

#include "grackle.h"
#include "status_reporting.h"

#include <cstring>  // std::memcmp, std::memcpy
#include <ostream>
#include <sstream>
#include <utility>  // std::pair, std::move

namespace grtest {

bool operator==(const GridLayout& lhs, const GridLayout& rhs) noexcept {
  bool rank = lhs.rank_ == rhs.rank_;
  bool start = std::memcmp(lhs.start_, rhs.start_, 3 * sizeof(int)) == 0;
  bool stop = std::memcmp(lhs.stop_, rhs.stop_, 3 * sizeof(int)) == 0;
  bool dim = std::memcmp(lhs.dim_, rhs.dim_, 3 * sizeof(int)) == 0;
  return rank && start && stop && dim;
}

void PrintTo(const GridLayout& layout, std::ostream* os) {
  auto write_arr_ = [&](const int* vals) -> void {
    for (int i = 0; i < layout.rank_; i++) {
      *os << ((i == 0) ? "{" : ", ");
      *os << vals[i];
    }
    *os << '}';
  };
  *os << "GridLayout{ dim: ";
  write_arr_(layout.dim_);
  *os << ", start: ";
  write_arr_(layout.start_);
  *os << ", stop: ";
  write_arr_(layout.stop_);
  *os << '}';
}

std::string GridLayout::to_string() const {
  std::stringstream s;
  PrintTo(*this, &s);
  return s.str();
}

namespace field_detail {

/// Construct a name-pointer mapping based on @p ctx_pack
///
/// The returned where the keys are the names of every active Grackle field
/// and the associated values are all nullptr
static MapType make_nullptr_map_(const GrackleCtxPack& ctx_pack,
                                 const std::set<std::string>& exclude_fields) {
  MapType m;
  // fill up m with (field, nullptr) pairs for each field enabled by ctx_pack
  auto fn = [&m](const char* name, const FieldInfo& info) {
    if (!m.emplace(name, nullptr).second) {
      GR_INTERNAL_ERROR("%s was inserted multiple times", name);
    }
  };
  for_each_named_field(ctx_pack, fn);

  for (const std::string& name : exclude_fields) {
    auto search = m.find(name);
    if (search != m.end()) {
      m.erase(search);
    }
  }
  return m;
}

std::pair<CorePack, Status> CorePack::setup_1d(MapType&& premade_map,
                                               int buf_size) {
  if (buf_size <= 0) {
    return {CorePack{}, error::Adhoc("buf_size must be positive")};
  }
  std::size_t data_sz = premade_map.size() * buf_size;

  GridLayout layout = GridLayout::from_dim(1, &buf_size);

  CorePack pack{
      // the trailing parentheses when allocating an array of integers or
      // floating point values sets the initial values to 0
      /* layout */ std::unique_ptr<GridLayout>(new GridLayout(layout)),
      /* data_buf */ std::unique_ptr<gr_float[]>(new gr_float[data_sz]()),
      /* my_fields */
      std::unique_ptr<grackle_field_data>(new grackle_field_data),
      /* map */ std::move(premade_map)};
  if (gr_initialize_field_data(pack.my_fields.get()) != GR_SUCCESS) {
    // this really should never fail
    return {CorePack{}, error::Adhoc("gr_initialize_field_data failed")};
  }

  // setup grid properties
  pack.my_fields->grid_dimension = pack.layout->dim_;
  pack.my_fields->grid_start = pack.layout->start_;
  pack.my_fields->grid_end = pack.layout->end_;

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
    const GrackleCtxPack& ctx_pack, int buf_size,
    const std::set<std::string>& exclude_fields) {
  if (!ctx_pack.is_initialized()) {
    return {FieldContainer{}, error::Adhoc("ctx_pack isn't initialized")};
  }

  // construct a map for each relevant field (the values are all nullptr)
  field_detail::MapType m =
      field_detail::make_nullptr_map_(ctx_pack, exclude_fields);
  std::pair<field_detail::CorePack, Status> tmp =
      field_detail::CorePack::setup_1d(std::move(m), buf_size);
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
    const std::set<std::string>& exclude_fields) {
  int total_count = layout.n_elements();

  std::pair<FieldContainer, Status> tmp =
      FieldContainer::create_1d(ctx_pack, total_count, exclude_fields);

  if (tmp.second.is_ok()) {
    tmp.first.data_.override_layout(layout);
  }
  return tmp;
}

FieldContainer FieldContainer::clone() const {
  // make a copy of m, and use it construct a new CorePack
  // -> the keys of `m` tell `setup_1d_CorePack` which Grackle fields get used
  // -> right after we create
  field_detail::MapType m = this->data_.map;
  std::pair<field_detail::CorePack, Status> tmp =
      field_detail::CorePack::setup_1d(std::move(m),
                                       this->grid_layout().n_elements());
  // it shouldn't be possible for tmp.second.is_err() to return true
  FieldContainer out;
  out.data_ = std::move(tmp.first);

  // copy over layout properties
  out.data_.override_layout(*this->data_.layout);

  // now copy over field values
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
  int length = this->grid_layout().n_elements() * this->n_fields();
  std::memcpy(dst, src, length * sizeof(gr_float));
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

void PrintTo(const FieldContainer& fc, std::ostream* os) {
  const GridLayout& layout = fc.grid_layout();
  *os << "FieldContainer{\n"
      << "  " << layout.to_string() << ",\n"
      << "  grid_dx = " << fc.data_.my_fields->grid_dx << ",\n"
      << "  fields = {\n";
  std::size_t n_elements = layout.n_elements();
  for (const auto& it : fc) {
    const std::string& field_name = it.first;
    const gr_float* ptr = it.second;
    *os << "    " << field_name << " =\n";
    *os << "      " << ptr_to_string(ptr, n_elements) << '\n';
  }
  *os << "  }\n}\n";
}

}  // namespace grtest