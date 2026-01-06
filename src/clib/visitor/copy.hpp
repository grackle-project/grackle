// See LICENSE file for license and copyright information

/// @file copy.hpp
/// @brief declare/implement the CopyMembers visitor

#ifndef VISITOR_COPY_HPP
#define VISITOR_COPY_HPP

#include "common.hpp"
#include <cstring>  // std::memcpy

namespace grackle::impl::visitor {

/// A visitor that copies the contents of a type
///
/// @note
/// This is for demonstration purposes. We will need to create a type that is
/// a lot like this in order to support GPUs
class CopyMembers {
  VisitorCtx ctx;

public:
  CopyMembers() = delete;
  explicit CopyMembers(VisitorCtx ctx) : ctx(ctx) {}

  void operator()(const char* name, int*& src, int*& dst,
                  const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    std::memcpy(dst, src, sizeof(int) * n_elem);
  }

  void operator()(const char* name, long long*& src, long long*& dst,
                  const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    std::memcpy(dst, src, sizeof(long long) * n_elem);
  }

  void operator()(const char* name, double*& src, double*& dst,
                  const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    std::memcpy(dst, src, sizeof(double) * n_elem);
  }
};

}  // namespace grackle::impl::visitor

#endif  // VISITOR_COPY_HPP
