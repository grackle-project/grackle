//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare the GrackleCtxType
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_GRACKLECTXPACK_HPP
#define GRTESTUTILS_GRACKLECTXPACK_HPP

#include "./preset.hpp"
#include "./status.hpp"
#include <algorithm>  // std::remove_if
#include <memory>     // std::unique_ptr

namespace grtest {

namespace param_detail {

/// Tracks allocations of string parameters stored by a @ref chemistry_data
///
/// The basic premise is that when you construct a @ref chemistry_data instance,
/// you would also construct a StrAllocTracker instance that has the same
/// lifetime as @ref chemistry_data. And then you would pass both instances to
/// @ref set_str when you want to store a string parameter
///
/// @par Implementation Notes
/// * We explicitly avoid tracking the allocations within this type using
///   std::string because SSO (small string optimization) could cause weird
///   bugs in the future during refactoring if we aren't super careful
///
/// * We assume max number of strings is small (i.e. <=5). If the number grows,
///   use std::unordered_set instead of std::vector
class StrAllocTracker {
  /// manages the lifetime of string allocations
  std::vector<std::unique_ptr<char[]>> bufs_;

public:
  /// allocate a new buffer and try to free the old buffer
  ///
  /// @param tot_len Total length of new buffer (including the nul terminator)
  /// @param old The old buffer that we are trying to replace
  char* alloc_buf_and_free_old(std::size_t tot_len, const char* old) {
    // erase the entry in bufs_ holding old (this also frees memory). If old
    // isn't found, it's a nullptr or we assume that it's a string literal
    for (std::size_t i = 0; i < bufs_.size(); i++) {
      if (bufs_[i].get() == old) {
        bufs_.erase(bufs_.begin() + i);
        break;
      }
    }

    if (tot_len == 0) {
      return nullptr;
    }
    const std::unique_ptr<char[]>& p = bufs_.emplace_back(new char[tot_len]);
    return p.get();
  }
};

}  // namespace param_detail

/// Tracks the group of Grackle objects needed for executing API functions
///
/// The primary motivation for this object's existence is making sure that
/// the allocations get cleaned up when a test fails
///
/// @note
/// Ideally, we would only make it possible to create a fully initialized
/// instance (in that case, we would delete the default constructor), but that
/// involves a bunch more work
///
/// We choose to implement this in terms of std::unique_ptr (rather than raw
/// pointers) since it implements proper move semantics for us
class GrackleCtxPack {
  /// units used for initializing chemistry_data
  code_units initial_units_;
  /// tracks string allocations for @ref my_chemistry_
  param_detail::StrAllocTracker str_allocs_;
  /// the fully initialized chemistry_data instance
  std::unique_ptr<chemistry_data> my_chemistry_;
  /// the fully initialized chemistry_data_storage_instance
  std::unique_ptr<chemistry_data_storage> my_rates_;

public:
  /// Construct an uninitialized instance
  GrackleCtxPack()
      : initial_units_{}, my_chemistry_(nullptr), my_rates_(nullptr) {}

  GrackleCtxPack(GrackleCtxPack&&) = default;
  GrackleCtxPack& operator=(GrackleCtxPack&&) = default;

  // forbid copy and assignment operations...
  // -> we could re-enable them if we wanted to be able to add to clone the
  //    types (or if we wanted to internally use std::shared_ptr
  GrackleCtxPack(const GrackleCtxPack&) = delete;
  GrackleCtxPack& operator=(const GrackleCtxPack&) = delete;

  ~GrackleCtxPack() {
    if (!this->is_initialized()) {
      return;
    }
    local_free_chemistry_data(this->my_chemistry_.get(), this->my_rates_.get());
    // unique_ptr destructor will handle calls to delete
  }

  bool is_initialized() const { return this->my_chemistry_ != nullptr; }

  // getter functions
  const code_units& initial_units() const { return this->initial_units_; }
  chemistry_data* my_chemistry() { return this->my_chemistry_.get(); }
  const chemistry_data* my_chemistry() const {
    return this->my_chemistry_.get();
  }
  chemistry_data_storage* my_rates() { return this->my_rates_.get(); }

  /// create an initialized instance
  static std::pair<GrackleCtxPack, Status> create(const ParamConf& param_conf);
};

}  // namespace grtest

#endif  // GRTESTUTILS_GRACKLECTXPACK_HPP
