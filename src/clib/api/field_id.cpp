//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements machinery for querying field_id and field_names
///
//===----------------------------------------------------------------------===//

#include "grackle.h"
#include "support/error.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

static Expected<int64_t, Error> get_field_id_(
    const chemistry_data_storage& my_rates, int64_t category,
    const char* name) {
  return Unexpected(Error::msg_literal("not implemented yet!"));
}

extern "C" int gr_get_field_id(const chemistry_data_storage* my_rates,
                               int64_t category, const char* name,
                               int64_t* out) {
  Expected<int64_t, Error> rslt = get_field_id_(*my_rates, category, name);
  if (rslt.has_value()) {
    *out = rslt.value();
    return GR_SUCCESS;
  }
  return GR_FAIL;
}

static Expected<const char*, Error> get_field_name_(
    const chemistry_data_storage& my_rates, int64_t category, int64_t id) {
  return Unexpected(Error::msg_literal("not implemented yet!"));
}

extern "C" int gr_get_field_name(const chemistry_data_storage* my_rates,
                                 int64_t category, int64_t id,
                                 const char** out) {
  Expected<const char*, Error> rslt = get_field_name_(*my_rates, category, id);
  if (rslt.has_value()) {
    *out = rslt.value();
    return GR_SUCCESS;
  }
  return GR_FAIL;
}

}  // namespace GRIMPL_NAMESPACE_DECL