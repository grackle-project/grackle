//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define the FullRxnRateBuf type
///
//===----------------------------------------------------------------------===//
#ifndef GRACKLE_FULL_RXN_RATE_BUF_HPP
#define GRACKLE_FULL_RXN_RATE_BUF_HPP

#include "LUT.hpp"
#include "internal_types.hpp"
#include "support/config.hpp"
#include "visitor/common.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// encloses some constants used for defining @ref FullFxnRateBuf
///
/// For some context, @ref FullRxnRateBuf holds various 1D rate buffers for
/// various kinds of rates. As explained in the @ref FullRxnRateBuf docstring,
/// we aggregate all of these buffers into single array, where contiguous
/// sections are dedicated to different kinds of rates. This namespace defines
/// the max known number of each kind of rate (and other derived quantities).
///
/// The known reaction rate kinds are currently:
/// - COLLISIONAL (collisional reaction rates)
/// - PHOTO (radiative reaction rates)
/// - GRAIN_GROWTH (tracks grain species growth rates)
/// - MISC (namely just tracks the h2dust buffer)
namespace rxn_rate_buf_detail {

/// holds constants specifying the maximum number of rates of a given kind
/// that may be stored the shared buffer
struct MaxSize {
  enum {
    COLLISIONAL = CollisionalRxnLUT::NUM_ENTRIES,
    PHOTO = PhotoRxnLUT::NUM_ENTRIES,
    GRAIN_GROWTH = OnlyGrainSpLUT::NUM_ENTRIES,
    MISC = 1
  };
};

/// holds constants specifying the offset from the start of the buffer to the
/// section holding a given kind of rate
struct SectionOffset {
  enum {
    COLLISIONAL = 0,
    PHOTO = COLLISIONAL + MaxSize::COLLISIONAL,
    GRAIN_GROWTH = PHOTO + MaxSize::PHOTO,
    MISC = GRAIN_GROWTH + MaxSize::GRAIN_GROWTH
  };
};

/// the max number of rates that FullRxnRateBuf::data may hold
inline constexpr int MAX_LEN = SectionOffset::MISC + MaxSize::MISC;

}  // namespace rxn_rate_buf_detail

/// Collection of buffers that store all computed reaction rates
///
/// For some context, Grackle evolves species densities by computing batches of
/// values (i.e. for 1 or more spatial location) for a variety of reaction
/// rates, and subsequently use the rates to determine time derivatives.
/// Instances of this type hold a buffer to hold the batch of values for every
/// kind of rate.
///
/// To access a buffer, use the associated functions
/// (@ref FullRxnRateBuf_kcol_bufs, @ref FullRxnRateBuf_kph_bufs,
/// @ref FullRxnRateBuf_grain_growth_bufs, or @ref FullRxnRateBuf_h2dust).
///
/// @par Implementation Strategy
/// Our current strategy involves aggregating as many of the rate buffers in
/// the @ref data member as possible. Different contiguous slices of @ref data
/// are dedicated to storing rates of a given kind.
///
/// @par Allocation Strategy
/// Currently, we reserve memory for the maximum number of rates known to
/// Grackle. This has been the historical convention adopted by Grackle.
/// Obviously this should be addressed in the future.
///
/// To move away from this over-allocation strategy, we'll need to move away
/// from using compile-time constants to index/size @ref data. Instead, we
/// would probably use the data structure proposed in #493
///
/// @par Longer Term Plan
/// The longer term goal is to adopt a dynamic strategy like the one outlined
/// in the docstring of @ref chemistry::species_density_updates_gauss_seidel.
/// This will be a lot easier if the various rate buffers are aggregated behind
/// accessible by indexing into a single array
struct FullRxnRateBuf {
  double* data[rxn_rate_buf_detail::MAX_LEN];
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template <class BinaryFn>
void visit_member_pair(FullRxnRateBuf& obj0, FullRxnRateBuf& obj1, BinaryFn f) {
  visitor::begin_visit("FullRxnRateBuf", f);
  for (int i = 0; i < rxn_rate_buf_detail::MAX_LEN; i++) {
    f(VIS_MEMBER_NAME("data[...]"), obj0.data[i], obj1.data[i],
      visitor::idx_range_len_multiple(1));
  }

  visitor::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(FullRxnRateBuf* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, FullRxnRateBuf, ptr, fn)
}

/// allocates the contents of a new FullRxnRateBuffer
///
/// @param nelem The number of elements in each buffer
inline FullRxnRateBuf new_FullRxnRateBuf(int nelem) {
  FullRxnRateBuf out;
  double* ptr = new double[nelem * rxn_rate_buf_detail::MAX_LEN];
  for (int i = 0; i < rxn_rate_buf_detail::MAX_LEN; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

/// performs cleanup of the contents of FullRxnRateBuf
///
/// This effectively invokes a destructor
inline void drop_FullRxnRateBuffer(FullRxnRateBuf* ptr) {
  // since we only allocate a single pointer for the data member, we only need
  // to call delete on the very first entry of data
  delete[] ptr->data[0];
}

// helper functions to retrieve the relevant sections of the buffers

/// fetch array of buffers reserved for holding collisional reaction rates
///
/// The returned array is expected to be indexed with @ref CollisionalRxnLUT.
/// The buffer at index `CollisionalRxnLUT::<rate>`, is reserved for a rate
/// called `<rate>`.
inline double* const* FullRxnRateBuf_kcol_bufs(const FullRxnRateBuf* ptr) {
  return ptr->data + rxn_rate_buf_detail::SectionOffset::COLLISIONAL;
}

/// fetch array of buffers reserved for holding radiative reaction rates
///
/// The returned array is expected to be indexed with @ref PhotoRxnLUT.
/// The buffer at index `PhotoRxnLUT::<rate>`, is reserved for a rate
/// called `<rate>`.
inline double* const* FullRxnRateBuf_kph_bufs(const FullRxnRateBuf* ptr) {
  return ptr->data + rxn_rate_buf_detail::SectionOffset::PHOTO;
}

/// fetch array of buffers reserved for holding grain growth rates
///
/// The returned array is expected to be indexed with @ref OnlyGrainSpLUT. The
/// buffer at index `OnlyGrainSpLUT::<name>`, is reserved for the growth rate
/// of the grain species called `<name>`.
inline double* const* FullRxnRateBuf_grain_growth_bufs(
    const FullRxnRateBuf* ptr) {
  return ptr->data + rxn_rate_buf_detail::SectionOffset::GRAIN_GROWTH;
}

/// fetch the buffer reserved for holding h2dust (the rate of molecular
/// Hydrogen formation on dust grains).
inline double* FullRxnRateBuf_h2dust(const FullRxnRateBuf* ptr) {
  return ptr->data[rxn_rate_buf_detail::SectionOffset::MISC];
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // GRACKLE_FULL_RXN_RATE_BUF_HPP
