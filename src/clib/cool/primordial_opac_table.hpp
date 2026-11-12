//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares @ref init_primordial_opacity
///
//===----------------------------------------------------------------------===//

#ifndef COOL_PRIMORDIAL_OPAC_TABLE_HPP
#define COOL_PRIMORDIAL_OPAC_TABLE_HPP

#include "../interp_grid.hpp"
#include "../support/config.hpp"
#include "../support/error.hpp"
#include "../support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<InterpGrid, Error> init_primordial_opacity();

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // COOL_PRIMORDIAL_OPAC_TABLE_HPP