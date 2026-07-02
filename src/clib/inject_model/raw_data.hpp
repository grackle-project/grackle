//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares some data structures used for organizing the raw input injection
/// model data, which is currently embedded as part of the Grackle library) and
/// declares a function for accessing this data.
///
/// @note
/// It's worth emphasizing that all machinery declared in this file is only
/// intended to be temporary. The long-term plan is to drop all of the
/// machinery once we start encoding the data within HDF5 files. The main
/// barrier to adopting HDF5 for this purpose is that I would like to first
/// remove all usage of the data currently encoded in
/// MetalNuclideYieldProps::total_yield (so that we can avoid encoding it in
/// the HDF5 file).
///
//===----------------------------------------------------------------------===//

#ifndef INJECT_MODEL_RAW_DATA_HPP
#define INJECT_MODEL_RAW_DATA_HPP

#include "grackle.h"  // GR_SUCCESS

namespace grackle::impl::inj_model_input {

/// the number of opacity-related coefficients that are grouped together
/// in the opacity table
static constexpr int N_Opacity_Coef = 4;

/// the number of Tdust values in an opacity table for a given grain species
/// and injection pathway. For each Tdust, there are N_Opacity_Coef coefficients
static constexpr int N_Tdust_Opacity_Table = 35;

/// the number of hard-coded injection pathways
static constexpr int N_Injection_Pathways = GRIMPL_MAX_INJ_PATHWAYS;

struct MetalNuclideYieldProps;
struct GrainSpeciesYieldProps;

/// This is used to organize the storage of the injection pathway data
///
/// Essentially, this dictates how the data is organized on disk (as part of
/// the Grackle library). This data will be somewhat reordered to actually
/// perform calculations
///
/// There are a few reasons why it makes some sense to store this data in a
/// a slightly different format than the way is organized during the core
/// calculations:
/// - the way we organize it here (which has some similarities to a struct of
///   arrays) is a lot less error-prone.
/// - the way we organize it here allows us to use a more sparse representation
///   (i.e. we omit yields from pathways when the yields are 0)
///
/// @note
/// The plan is to eventually shift the data into a file loaded at runtime
struct InjectionPathwayInputData {
  const MetalNuclideYieldProps* metal_nuclide_yields;
  int n_metal_nuclide_yields;

  const GrainSpeciesYieldProps* initial_grain_props;
  int n_injected_grain_species;
};

struct MetalNuclideYieldProps {
  /// name of the nuclide
  const char* name;

  /// fraction of non-primordial injection corresponding to the nuclide in
  /// the gas phase
  double gas_yield;

  /// total fraction of non-primordial injection corresponding to the nuclide
  /// - I'm pretty sure we can get rid of this information in the future
  double total_yield;
};

/// Each injection array will hold an array of these structs
struct GrainSpeciesYieldProps {
  /// name of the Grain Species
  const char* name;

  /// Initial yield
  double nonprimoridal_yield_frac;

  /// the 1st, 2nd, and 3rd order moments of the initial size distribution for
  /// the grain species.
  ///
  /// - Element 0 specifies \f$\langle r^1 \rangle_j\f$. This is the grain
  ///   species's average **INITIAL** radius. Has units of cm.
  /// - Element 1 specifies \f$\langle r^2 \rangle_j\f$. The product of this
  ///   value and π gives the grain species's average **INITIAL**
  ///   cross-section. Has units of centimeters squared.
  /// - Element 2 specifies \f$\langle r^3 \rangle_j\f$. The product of this
  ///   value and (4π/3) is the grain species's average **INITIAL** volume.
  ///   Has units of centimeters cubed.
  double size_moments[3];

  double opacity_coef_table[N_Tdust_Opacity_Table][N_Opacity_Coef];
};

extern "C" {
/// protoype for callback function used with input_inject_model_iterate
typedef int (*inj_iterate_t)(const char* name,
                             const InjectionPathwayInputData* input,
                             void* op_data);
}

/// iterates over injection model input data packs with user callback routine
///
/// @param[in] op Callback function
/// @param[inout] op_data User-defined callback function context
///
/// @returns GR_SUCCESS if all calls to @p op returns GR_SUCCESS. Any other
///     value denotes an error.
///
/// If any call to @p op is not `GR_SUCCESS`, this function exits immediately
///
/// @note
/// we choose to provide access to each injection model's raw input data via
/// callback routines since the equivalent HDF5 (when we eventually shift) will
/// also involve callbacks (i.e. via `H5Literate`). By initially writing the
/// code using callbacks, the transition to HDF5 should be easier.
int input_inject_model_iterate(inj_iterate_t op, void* op_data);

}  // namespace grackle::impl::inj_model_input

#endif  // INJECT_MODEL_RAW_DATA_HPP
