//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implement logic pertaining to GrackleCtxPack
///
//===----------------------------------------------------------------------===//

#include "GrackleCtxPack.hpp"
#include "status_reporting.h"

#include <cstring>

namespace grtest::param_detail {

static bool set_str_(chemistry_data& my_chem, const std::string& name,
                     const char* val, std::size_t sz_with_nul,
                     StrAllocTracker* str_allocs) {
  // NOTE: we should NOT directly modify characters held by field_ptr
  const char** dest = const_cast<const char**>(
      local_chemistry_data_access_string(&my_chem, name.c_str()));
  if (dest == nullptr) {
    return false;  // field is either not known or not a string
  }

  // the following acts as a compiler hint (to suppress a warning)
  GR_INTERNAL_REQUIRE((val == nullptr) != (sz_with_nul == 0), "compiler-hint");

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

bool set_str(chemistry_data& my_chem, const std::string& name, const char* val,
             StrAllocTracker* str_allocs) {
  std::size_t sz_with_nul = (val == nullptr) ? 0 : std::strlen(val) + 1;
  return set_str_(my_chem, name, val, sz_with_nul, str_allocs);
}

bool set_str(chemistry_data& my_chem, const std::string& name,
             const std::string val, StrAllocTracker* str_allocs) {
  return set_str_(my_chem, name, val.data(), val.size() + 1, str_allocs);
}

}  // namespace grtest::param_detail