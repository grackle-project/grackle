/// @file memory.hpp
/// @brief declare/implement the Printer visitor
///
/// We can implement multiple different kinds of strategies. Examples include:
/// 1. the current strategy where we allocate/free on a member-by-member basis
/// 2. a strategy where we (i) use a "Sizer" visitor to count up the total
///    number of elements with a given type, across all members, (ii) allocate
///    the memory for each type all at once, (iii) use a visitor to distribute
///    the memory of each type to individual members.
///
/// Strategy #2 This should be measurably faster. Most immediately it will be
/// faster when we have relatively few elements per buffer. It will also be
/// faster because we have fewer calls to free.

#ifndef VISITOR_MEMORY_HPP
#define VISITOR_MEMORY_HPP

#include "common.hpp"
#include "grackle_macros.h"  // GRACKLE_FREE

#include <cstdlib>  // std::malloc

namespace grackle::impl::visitor {

/// Visitor that separately allocates memory for each visited member
class AllocateMembers {
  VisitorCtx ctx;

public:
  AllocateMembers() = delete;
  explicit AllocateMembers(VisitorCtx ctx) : ctx{ctx} {}

  void operator()(const char* name, int*& mem, const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    mem = static_cast<int*>(malloc(sizeof(int) * n_elem));
  }

  void operator()(const char* name, long long*& mem,
                  const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    mem = static_cast<long long*>(malloc(sizeof(long long) * n_elem));
  }

  void operator()(const char* name, double*& mem,
                  const BufLenSpec& spec) const {
    std::size_t n_elem = get_buf_len(spec, ctx);
    mem = static_cast<double*>(malloc(sizeof(double) * n_elem));
  }
};

// no need to specialize begin_visit or end_visit for AllocateMembers

/// Visitor that separately deallocates memory for each visited member
struct FreeMembers {
  FreeMembers() = default;

  void operator()(const char* name, int*& mem, const BufLenSpec& spec) {
    GRACKLE_FREE(mem);
  }

  void operator()(const char* name, long long*& mem,
                  const BufLenSpec& spec) const {
    GRACKLE_FREE(mem);
  }

  void operator()(const char* name, double*& mem,
                  const BufLenSpec& spec) const {
    GRACKLE_FREE(mem);
  }
};

// no need to specialize begin_visit or end_visit for FreeMembers

}  // namespace grackle::impl::visitor

#endif  // VISITOR_MEMORY_HPP
