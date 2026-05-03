
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

LogTLinInterpScratchBuf new_LogTLinInterpScratchBuf(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  LogTLinInterpScratchBuf out;
  visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  visit_member(&out, visitor::AllocateMembers{ctx});
  return out;
}

void drop_LogTLinInterpScratchBuf(LogTLinInterpScratchBuf* ptr) {
  visit_member(ptr, visitor::FreeMembers{});
}

}  // namespace GRIMPL_NAMESPACE_DECL