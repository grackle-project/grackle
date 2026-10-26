//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares (and in several cases also defines) field-scaling functions to
/// account for cosmology
///
//===----------------------------------------------------------------------===//

#ifndef SCALE_FIELDS_HPP
#define SCALE_FIELDS_HPP

#include "grackle.h"
#include "chem_model/infer_species.hpp"
#include "field_adaptor.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "support/View.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Scales fields related to computing dust temperature
///
/// @param[in] factor The factor that is multiplied by the fields
/// @param[in] idx_range The index range that will be modified
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[inout] my_fields holds the fields that will be updated in-place
/// @param[inout] sp_densities Specifies the densities of the various species
///     that Grackle evolves (if any) in a format that allows the values to be
///     accessed with the index lookup table. Wherever possible, data should be
///     be accessed through this argument, rather than with @p my_fields
/// @param[in] n_inj_path_ptrs The number of pointers tracked by
///     `my_fields->inject_pathway_metal_density`
/// @param[in] sp_kind_map Specifies species index ranges that correspond
///     to particular categories of species
void scale_fields_dust(gr_float factor, IndexRange idx_range,
                       chemistry_data* my_chemistry,
                       grackle_field_data* my_fields,
                       SpeciesMultiView<gr_float> sp_densities,
                       int n_inj_path_ptrs, const PartMap<SpKind>& sp_kind_map);

/// @brief Scales density-related fields for the specified index range
///
/// @param[in] factor The factor that is multiplied by the fields
/// @param[in] idx_range The index range that will be modified
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[inout] my_fields holds the fields that will be updated in-place
/// @param[inout] sp_densities Specifies the densities of the various species
///     that Grackle evolves (if any) in a format that allows the values to be
///     accessed with the index lookup table. Wherever possible, data should be
///     be accessed through this argument, rather than with @p my_fields
/// @param[in] n_inj_path_ptrs The number of pointers tracked by
///     `my_fields->inject_pathway_metal_density`
/// @param[in] sp_kind_map Specifies species index ranges that correspond
///     to particular categories of species
void scale_fields(gr_float factor, IndexRange idx_range,
                  chemistry_data* my_chemistry, grackle_field_data* my_fields,
                  SpeciesMultiView<gr_float> sp_densities, int n_inj_path_ptrs,
                  const PartMap<SpKind>& sp_kind_map);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SCALE_FIELDS_HPP
