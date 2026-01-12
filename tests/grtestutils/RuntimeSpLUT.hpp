//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines RuntimeSpLUT, which is a nice C++ wrapper around some internal
/// constructs. This is very useful fore testing-purposes
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_RUNTIMESPLUT_HPP
#define GRTESTUTILS_RUNTIMESPLUT_HPP

#include <optional>
#include <string>

#include "runtime_splut.hpp"
#include "dust/grain_species_info.hpp"
#include "support/config.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"

namespace grtest {

class RuntimeSpLUT {
  GRIMPL_NS::PartMap partition_map_;
  GRIMPL_NS::FrozenKeyIdxBiMap bimap_;

public:
  void swap(RuntimeSpLUT& other) noexcept {
    GRIMPL_NS::PartMap tmp_pm = this->partition_map_;
    GRIMPL_NS::FrozenKeyIdxBiMap tmp_bimap = this->bimap_;
    this->partition_map_ = other.partition_map_;
    this->bimap_ = other.bimap_;
    other.partition_map_ = tmp_pm;
    other.bimap_ = tmp_bimap;
  }

  RuntimeSpLUT()
      : partition_map_(GRIMPL_NS::new_PartMap(nullptr, nullptr, 0)),
        bimap_(GRIMPL_NS::new_FrozenKeyIdxBiMap(
            nullptr, 0, GRIMPL_NS::BiMapMode::COPIES_KEYDATA)) {}

  // copy constructor
  RuntimeSpLUT(const RuntimeSpLUT& other)
      : partition_map_(other.partition_map_),
        bimap_(GRIMPL_NS::FrozenKeyIdxBiMap_clone(&other.bimap_)) {}

  // copy-assignment
  RuntimeSpLUT& operator=(const RuntimeSpLUT& other) noexcept {
    RuntimeSpLUT temp = other;
    swap(temp);
    return *this;
  }

  // move-constructor
  RuntimeSpLUT(RuntimeSpLUT&& o) noexcept : RuntimeSpLUT() { swap(o); }

  // move-assignment
  RuntimeSpLUT& operator=(RuntimeSpLUT&& other) noexcept {
    RuntimeSpLUT::swap(other);
    return *this;
  }

  ~RuntimeSpLUT() {
    GRIMPL_NS::drop_PartMap(&partition_map_);
    GRIMPL_NS::drop_FrozenKeyIdxBiMap(&bimap_);
  }

  static std::optional<RuntimeSpLUT> create(
      int primordial_chemistry, int metal_chemistry,
      const GRIMPL_NS::GrainSpeciesInfo* grain_info);

  static std::optional<RuntimeSpLUT> create(int primordial_chemistry,
                                            int metal_chemistry,
                                            int dust_species);

  // this is a little silly (but it's useful)
  static std::optional<RuntimeSpLUT> create(
      const chemistry_data* my_chemistry) {
    return RuntimeSpLUT::create(my_chemistry->primordial_chemistry,
                                my_chemistry->metal_chemistry,
                                my_chemistry->dust_species);
  }

  int size() const { return GRIMPL_NS::FrozenKeyIdxBiMap_size(&bimap_); }

  std::optional<int> find(const char* key) const {
    GRIMPL_NS::bimap::AccessRslt tmp =
        GRIMPL_NS::FrozenKeyIdxBiMap_find(&bimap_, key);
    return (tmp.has_value) ? std::optional{tmp.value} : std::nullopt;
  }

  std::optional<int> find(const std::string& key) const {
    return find(key.c_str());
  }

  std::optional<std::string> inverse_find(int index) const {
    const char* tmp = GRIMPL_NS::FrozenKeyIdxBiMap_inverse_find(&bimap_, index);
    return (tmp == nullptr) ? std::nullopt : std::optional{std::string{tmp}};
  }

  std::optional<int> sp_kind(int index) const {
    GRIMPL_NS::partmap::IdxSearch tmp =
        PartMap_search_idx(&partition_map_, index);
    return (tmp.has_val) ? std::optional{tmp.pd} : std::nullopt;
  }
  std::optional<int> sp_kind(const char* key) const {
    GRIMPL_NS::partmap::IdxSearch tmp =
        GRIMPL_NS::key_partition_search(&partition_map_, &bimap_, key);
    return (tmp.has_val) ? std::optional{tmp.pd} : std::nullopt;
  }
  std::optional<int> sp_kind(const std::string& key) const {
    return sp_kind(key.c_str());
  }

  std::optional<int> kind_index(int index) const {
    GRIMPL_NS::partmap::IdxSearch tmp =
        PartMap_search_idx(&partition_map_, index);
    return (tmp.has_val) ? std::optional{tmp.start_offset} : std::nullopt;
  }
  std::optional<int> kind_index(const char* key) const {
    GRIMPL_NS::partmap::IdxSearch tmp =
        GRIMPL_NS::key_partition_search(&partition_map_, &bimap_, key);
    return (tmp.has_val) ? std::optional{tmp.start_offset} : std::nullopt;
  }
  std::optional<int> kind_index(const std::string& key) const {
    return kind_index(key.c_str());
  }

  /// expects GRIMPL_NS::SpKind::CHEMICAL or GRIMPL_NS::SpKind::DUST
  GRIMPL_NS::IdxInterval kind_interval(int kind) const {
    return grackle::impl::PartMap_part_bounds(&partition_map_, kind);
  }
};

}  // namespace grtest

#endif  // GRTESTUTILS_RUNTIMESPLUT_HPP