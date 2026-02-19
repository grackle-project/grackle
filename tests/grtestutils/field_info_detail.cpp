//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// help implement logic pertaining to query_field_info
///
//===----------------------------------------------------------------------===//

#include "./field_info_detail.hpp"

#include <array>
#include <string_view>
#include <unordered_map>
#include <utility>  // std::pair

namespace grtest::field_detail {

namespace {  // stuff inside an anonymous namespace is local to this file

constexpr int get_known_field_count_() {
  int i = 0;
  for_each_named_field([&i](const char* name, const FieldInfo& finfo) { i++; });
  return i;
}

/// Holds the number of known grackle fields (known at compile-time)
constexpr int N_KNOWN_FIELDS = get_known_field_count_();

using FieldNameInfoArray =
    std::array<std::pair<std::string_view, FieldInfo>, N_KNOWN_FIELDS>;

constexpr FieldNameInfoArray make_name_info_array() {
  FieldNameInfoArray out;
  std::size_t i = 0;
  for_each_named_field([&out, &i](const char* name, const FieldInfo& finfo) {
    std::pair<std::string_view, FieldInfo>& cur_pair = out[i++];
    cur_pair.first = std::string_view(name);
    cur_pair.second = finfo;
  });
  return out;
}

/// Holds the number of known grackle fields
///
/// @note Computed at compile-time
constexpr FieldNameInfoArray name_info_array = make_name_info_array();

/// A mapping between known all known field names and FieldInfo structs
///
/// @par Startup Overhead
/// Each time an executable linked against this file starts up, this mapping is
/// constructed before the main function gets executed. It's unfortunate that
/// we have this overhead:
/// - To reduce the overhead, we compute the key-value pairs at compile-time
/// - The best alternative would involve constructing this mapping the very
///   first time we invoke query_field_info. We probably would want to use
///   something like std::call_once (but we need to be careful about colliding
///   threading runtimes)
///
/// @par Implementation Note
/// Ordinarily, C++ containers should avoid holding std::string_view instances
/// since a std::string_view isn't tied to an underlying allocation, which can
/// lead to unexpected lifetime issues. In this case
/// - the use of std::string_view is ok because
///   1. each value is a string-literal (i.e. there is no way for the underlying
///      memory to be de-allocated)
///   2. this variable is constant (i.e. there is no way for keys to be added
///      that don't wrap string-literals)
/// - moreover, the only alternative is to use std::string and that would
///   involve heap allocations (adding to startup overhead)
///
/// @note
const std::unordered_map<std::string_view, FieldInfo> name_info_map_(
    name_info_array.begin(), name_info_array.end());

}  // anonymous namespace

std::optional<FieldInfo> query_field_info(std::string_view name) noexcept {
  auto search = name_info_map_.find(name);
  return (search == name_info_map_.end())
             ? std::optional<FieldInfo>{}
             : std::optional<FieldInfo>{search->second};
}

}  // namespace grtest::field_detail