//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the function for loading injection pathway data
///
//===----------------------------------------------------------------------===//

#ifndef INJECT_MODEL_LOAD_DATA_HPP
#define INJECT_MODEL_LOAD_DATA_HPP

#include "grackle.h"
#include "../ratequery.hpp"
#include "../support/config.hpp"
#include "../support/error.hpp"
#include "../support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// loads the model data for the various injection pathways and update
/// @p my_rates, accordingly
///
/// @todo
/// It might be nice if we returned the constructed injection_pathway object
Expected<void, Error> load_inject_path_data(const chemistry_data* my_chemistry,
                                            chemistry_data_storage* my_rates,
                                            ratequery::RegBuilder* reg_builder);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif /* INJECT_MODEL_LOAD_DATA_HPP */
