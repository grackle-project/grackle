
//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define some logic pertaining to ln(T) buffers
///
//===----------------------------------------------------------------------===//

#include "lnT_prep.hpp"
#include "visitor/memory.hpp"

namespace GRIMPL_NAMESPACE_DECL {

LnTLinInterpBuf new_LnTLinInterpBuf(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  LnTLinInterpBuf out;
  visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  visit_member(&out, visitor::AllocateMembers{ctx});
  return out;
}

void drop_LnTLinInterpBuf(LnTLinInterpBuf* ptr) {
  visit_member(ptr, visitor::FreeMembers{});
}

}  // namespace GRIMPL_NAMESPACE_DECL