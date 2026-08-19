//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/implements @ref DustSolver
///
//===----------------------------------------------------------------------===//

#ifndef DUST_SOLVER_HPP
#define DUST_SOLVER_HPP

#include "../support/config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Provides a standard interface for computing the contributions of
///     dust chemistry to the overall chemistry solve
///
/// The basic premise is that:
/// - while performing the full chemistry calculation, the rest of grackle
///   doesn't need to know anything at all about the dust model beyond this
///   object (this is somewhat aspirational, right now make_consistent and
///   the step_rate_ functions still require some knowledge, but we're moving
///   away from that)
/// - objects of this type are immutable. After they are created, they have no
///   mutable state (any scratch buffers will need to be passed in)
class DustSolver {
  // in the future, we'll add some configuration data
  // - example: when using the Chiaki multi-grain-species-growth model, we may
  //   track injection pathway information as part of this type (after we
  //   decouple metal-chemistry from injection pathway information)

public:
  /// @brief default constructor
  DustSolver() = default;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // DUST_SOLVER_HPP