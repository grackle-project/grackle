//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Logic pertaining to ParamVal
///
//===----------------------------------------------------------------------===//
#include "status_reporting.h"

#include "./param.hpp"
#include "./grackle_ctx_pack.hpp"  // param_detail::StrAllocTracker

#include <cstring>             // std::strlen
#include <ostream>             // std::ostream
#include <type_traits>         // std::false_type, std::is_same_v, std::decay_t

namespace {  // stuff inside an anonymous namespace is local to this file
// this is used to report a compile-time error in the lambda function that
// we use in std::visit (this is standard practice)
template <class...>
constexpr std::false_type always_false_{};
}  // anonymous namespace

namespace grtest {

std::string ParamVal::to_string(bool unwrap) const {
  std::string tmp;
  // this lambda function is passed a reference to the value within this->val_
  auto get_string_repr_ = [&tmp](const auto& v) -> void {
    using T = std::decay_t<decltype(v)>;
    if constexpr (std::is_same_v<T, int> || std::is_same_v<T, double>) {
      tmp = std::to_string(v);
    } else if constexpr (std::is_same_v<T, std::string>) {
      tmp.reserve(v.size() + 2);
      tmp += '"';
      tmp += v;
      tmp += '"';
    } else if constexpr (std::is_same_v<T, std::nullptr_t>) {
      tmp = "nullptr";
    } else {
      static_assert(always_false_<T>, "encountered unhandled type");
    }
  };
  std::visit(get_string_repr_, this->val_);

  if (unwrap) {
    return tmp;
  }
  std::string out;
  out.reserve(tmp.size() + 10);
  out = "ParamVal(";
  out += tmp;
  out += ")";
  return out;
}

void PrintTo(const ParamVal& p, std::ostream* os) { *os << p.to_string(); }

bool ParamVal::is_equal(int val) const {
  return std::holds_alternative<int>(val_) && val == std::get<int>(val_);
}

bool ParamVal::is_equal(double val) const {
  return std::holds_alternative<double>(val_) && val == std::get<double>(val_);
}

bool ParamVal::is_equal(std::string_view val) const {
  return std::holds_alternative<std::string>(val_) &&
         val == std::get<std::string>(val_);
}

bool ParamVal::is_equal(const char* val) const {
  if (val == nullptr) {
    return std::holds_alternative<std::nullptr_t>(val_);
  }
  return is_equal(std::string_view(val));
}

namespace {  // stuff inside an anonymous namespace is local to this file

bool set_str_(chemistry_data& my_chem, const std::string& name, const char* val,
              std::size_t sz_with_nul,
              param_detail::StrAllocTracker* str_allocs) {
  // NOTE: we should NOT directly modify characters held by field_ptr
  const char** dest = const_cast<const char**>(
      local_chemistry_data_access_string(&my_chem, name.c_str()));
  if (dest == nullptr) {
    return false;  // field is either not known or not a string
  }

  // the following acts as a compiler hint (to suppress a warning)
  GR_INTERNAL_REQUIRE((val == nullptr) != (sz_with_nul > 0), "compiler-hint");

  if (str_allocs != nullptr) {
    // (if applicable) allocate a new buffer and deallocate the old buffer
    char* new_alloc = str_allocs->alloc_buf_and_free_old(sz_with_nul, *dest);
    if (sz_with_nul > 0) {
      std::memcpy(new_alloc, val, sz_with_nul - 1);
      new_alloc[sz_with_nul - 1] = '\0';
    }
    (*dest) = new_alloc;
  } else {
    (*dest) = val;
  }
  return true;
}

/// Tries to set a string parameter tracked by @p val
///
/// @param[in,out] my_chem Tracks various Grackle parameters
/// @param[in] name The name of the parameter getting updated
/// @param[in] val The value of the parameter
/// @param[in] str_allocs Tracks allocations of the strings held by @p my_chem
/// @returns true if successful and `false` if there was an error (e.g. @p name
///     isn't a known parameter or is a parameter that doesn't expect a string)
///
/// @note
/// Ordinarily, this function will update @p str_allocs to hold a copy of
/// @p val and @p my_chem will be updated to store a pointer to that copy. When
/// @p str_allocs is a nullptr, @p my_chem is directly updated to track the data
/// referenced by @p val.
bool set_str(chemistry_data& my_chem, const std::string& name, const char* val,
             param_detail::StrAllocTracker* str_allocs) {
  std::size_t sz_with_nul = (val == nullptr) ? 0 : std::strlen(val) + 1;
  return set_str_(my_chem, name, val, sz_with_nul, str_allocs);
}

bool set_str(chemistry_data& my_chem, const std::string& name,
             const std::string& val,
             param_detail::StrAllocTracker* str_allocs) {
  return set_str_(my_chem, name, val.data(), val.size() + 1, str_allocs);
}

bool set_int(chemistry_data& my_chem, const std::string& name, int val) {
  int* dest = local_chemistry_data_access_int(&my_chem, name.c_str());
  if (dest == nullptr) {
    return false;
  }
  (*dest) = val;
  return true;
}

bool set_double(chemistry_data& my_chem, const std::string& name, double val) {
  double* dest = local_chemistry_data_access_double(&my_chem, name.c_str());
  if (dest == nullptr) {
    return false;
  }
  (*dest) = val;
  return true;
}

}  // anonymous namespace

bool set_param(chemistry_data& my_chem, const std::string& name,
               const ParamVal& par_val,
               param_detail::StrAllocTracker* str_allocs) {
  bool success = false;

  // this lambda function is passed a reference to the value within par_val.val_
  auto set_val = [&success, &my_chem, &name,
                  str_allocs](const auto& v) -> void {
    using T = std::decay_t<decltype(v)>;
    if constexpr (std::is_same_v<T, int>) {
      success = set_int(my_chem, name, v);
    } else if constexpr (std::is_same_v<T, double>) {
      success = set_double(my_chem, name, v);
    } else if constexpr (std::is_same_v<T, std::string> ||
                         std::is_same_v<T, std::nullptr_t>) {
      success = set_str(my_chem, name, v, str_allocs);
    } else {
      static_assert(always_false_<T>, "encountered unhandled type");
    }
  };

  std::visit(set_val, par_val.val_);
  return success;
}

}  // namespace grtest