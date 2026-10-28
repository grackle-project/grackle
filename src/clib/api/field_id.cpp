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

#include <cstdint>
#include "grackle.h"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief determines the appropriate bimap
static Expected<const FrozenKeyIdxBiMap*, Error> find_bimap_(
    const chemistry_data_storage& my_rates, int64_t category) {
  // in the future, we'll fetch the bimap from my_rates->opaque_storage
  return Unexpected(Error::msg_literal("recieved unknown category"));
}

/// does the heavy lifting
static Expected<int64_t, Error> get_field_id_(
    const chemistry_data_storage& my_rates, int64_t category,
    const char* name) {
  Expected<const FrozenKeyIdxBiMap*, Error> tmp =
      find_bimap_(my_rates, category);
  if (!tmp.has_value()) {
    return Unexpected(tmp.error());
  }
  const FrozenKeyIdxBiMap& bimap = *tmp.value();
  if (name == nullptr) {
    return bimap.size();
  }
  std::optional<uint16_t> search_rslt = bimap.find(name);
  if (search_rslt.has_value()) {
    return static_cast<int64_t>(search_rslt.value());
  } else {
    return -1;  // <- a negative value indicates that name couldn't be found
  }
}

extern "C" int gr_get_field_id(const chemistry_data_storage* my_rates,
                               int64_t category, const char* name,
                               int64_t* out) {
  Expected<int64_t, Error> rslt = get_field_id_(*my_rates, category, name);
  if (rslt.has_value()) {
    *out = rslt.value();
    return GR_SUCCESS;
  }
  rslt.error().write(stderr);
  return GR_FAIL;
}

static Expected<const char*, Error> get_field_name_(
    const chemistry_data_storage& my_rates, int64_t category, int64_t id) {
  Expected<const FrozenKeyIdxBiMap*, Error> tmp =
      find_bimap_(my_rates, category);
  if (!tmp.has_value()) {
    return Unexpected(tmp.error());
  }
  const FrozenKeyIdxBiMap& bimap = *tmp.value();
  // it's worth considering whether the following 2 failure-modes should be
  // undefined behavior.
  if (id < 0) {
    return Unexpected(Error::msg_literal("id is negative"));
  } else if (id >= bimap.size()) {
    return Unexpected(Error::msg_literal("id is too big"));
  } else {
    return bimap.inverse_find(static_cast<int>(id));
  }
}

extern "C" int gr_get_field_name(const chemistry_data_storage* my_rates,
                                 int64_t category, int64_t id,
                                 const char** out) {
  Expected<const char*, Error> rslt = get_field_name_(*my_rates, category, id);
  if (rslt.has_value()) {
    *out = rslt.value();
    return GR_SUCCESS;
  }
  rslt.error().write(stderr);
  return GR_FAIL;
}

}  // namespace GRIMPL_NAMESPACE_DECL