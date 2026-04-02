//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implments the function for loading injection pathway data
///
//===----------------------------------------------------------------------===//

#include <cstring>  // std::strcmp
#include "../dust/grain_species_info.hpp"
#include "grackle_chemistry_data.h"
#include "load_data.hpp"  // forward declarations
#include "raw_data.hpp"
#include "../LUT.hpp"
#include "../opaque_storage.hpp"
#include "../ratequery.hpp"
#include "../status_reporting.h"  // GrPrintAndReturnErr
#include "../support/FrozenKeyIdxBiMap.hpp"

namespace {  // stuff inside an anonymous namespace is local to this file

/// a function that encodes the algorithm for looking up the pointer to the
/// gas yield fractions (as in fractions of the total injected non-primordial
/// material) from each injection pathway for the `i`th metal nuclide
///
/// This is intended to be registered with RegBuilder in order to provide
/// access to these quantities
///
/// @param my_rates The object from which the rate Entry is loaded
/// @param i the index of the queried rate
grackle::impl::ratequery::Entry nuclide_gas_yield_recipe(
    chemistry_data_storage* my_rates, int i) {
  namespace rateq = grackle::impl::ratequery;
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr) ||
      (my_rates->opaque_storage->inject_pathway_props == nullptr)) {
    return rateq::mk_invalid_Entry();
  }
  const grackle::impl::yields::MetalTables& tab =
      my_rates->opaque_storage->inject_pathway_props->gas_metal_nuclide_yields;

  switch (i) {
    case 0:
      return rateq::new_Entry(tab.C, "inject_path_gas_yield_frac.C");
    case 1:
      return rateq::new_Entry(tab.O, "inject_path_gas_yield_frac.O");
    case 2:
      return rateq::new_Entry(tab.Mg, "inject_path_gas_yield_frac.Mg");
    case 3:
      return rateq::new_Entry(tab.Al, "inject_path_gas_yield_frac.Al");
    case 4:
      return rateq::new_Entry(tab.Si, "inject_path_gas_yield_frac.Si");
    case 5:
      return rateq::new_Entry(tab.S, "inject_path_gas_yield_frac.S");
    case 6:
      return rateq::new_Entry(tab.Fe, "inject_path_gas_yield_frac.Fe");
    default:
      return rateq::mk_invalid_Entry();
  }
}

/// a function that encodes the algorithm for looking up the pointer to the
/// yield fractions (as in fractions of the total injected non-primordial
/// material) from each injection pathway for the `i`th dust grain species
///
/// This is intended to be registered with RegBuilder in order to provide
/// access to these quantities
///
/// @param my_rates The object from which the rate Entry is loaded
/// @param i the index of the queried rate
///
/// @note
/// This function makes 2 things clear:
/// 1. We *may* want to configure recipes to accept an arbitrary callback
///    argument so that we don't have to hardcode the names of every dust
///    grain species
/// 2. It would be advantageous to adopt key-names that are composed of 2
///    strings (that way we could reuse string-literals holding the grain
///    species names)
grackle::impl::ratequery::Entry grain_yield_recipe(
    chemistry_data_storage* my_rates, int i) {
  namespace rateq = grackle::impl::ratequery;
  if ((my_rates == nullptr) || (my_rates->opaque_storage == nullptr) ||
      (my_rates->opaque_storage->inject_pathway_props == nullptr)) {
    return rateq::mk_invalid_Entry();
  }
  const grackle::impl::GrainSpeciesCollection& tab =
      my_rates->opaque_storage->inject_pathway_props->grain_yields;
  double* const* data = tab.data;

  switch (i) {
    case 0:
      return rateq::new_Entry(data[OnlyGrainSpLUT::MgSiO3_dust],
                              "inject_path_grain_yield_frac.MgSiO3_dust");
    case 1:
      return rateq::new_Entry(data[OnlyGrainSpLUT::AC_dust],
                              "inject_path_grain_yield_frac.AC_dust");
    case 2:
      return rateq::new_Entry(data[OnlyGrainSpLUT::SiM_dust],
                              "inject_path_grain_yield_frac.SiM_dust");
    case 3:
      return rateq::new_Entry(data[OnlyGrainSpLUT::FeM_dust],
                              "inject_path_grain_yield_frac.FeM_dust");
    case 4:
      return rateq::new_Entry(data[OnlyGrainSpLUT::Mg2SiO4_dust],
                              "inject_path_grain_yield_frac.Mg2SiO4_dust");
    case 5:
      return rateq::new_Entry(data[OnlyGrainSpLUT::Fe3O4_dust],
                              "inject_path_grain_yield_frac.Fe3O4_dust");
    case 6:
      return rateq::new_Entry(data[OnlyGrainSpLUT::SiO2_dust],
                              "inject_path_grain_yield_frac.SiO2_dust");
    case 7:
      return rateq::new_Entry(data[OnlyGrainSpLUT::MgO_dust],
                              "inject_path_grain_yield_frac.MgO_dust");
    case 8:
      return rateq::new_Entry(data[OnlyGrainSpLUT::FeS_dust],
                              "inject_path_grain_yield_frac.FeS_dust");
    case 9:
      return rateq::new_Entry(data[OnlyGrainSpLUT::Al2O3_dust],
                              "inject_path_grain_yield_frac.Al2O3_dust");
    case 10:
      return rateq::new_Entry(data[OnlyGrainSpLUT::ref_org_dust],
                              "inject_path_grain_yield_frac.ref_org_dust");
    case 11:
      return rateq::new_Entry(data[OnlyGrainSpLUT::vol_org_dust],
                              "inject_path_grain_yield_frac.vol_org_dust");
    case 12:
      return rateq::new_Entry(data[OnlyGrainSpLUT::H2O_ice_dust],
                              "inject_path_grain_yield_frac.H2O_ice_dust");
    default:
      return rateq::mk_invalid_Entry();
  }
}

int configure_RegBuilder(const chemistry_data_storage* my_rates,
                         grackle::impl::ratequery::RegBuilder* reg_builder,
                         const char* const* inj_path_name_l, int n_pathways) {
  namespace rateq = grackle::impl::ratequery;

  // make list of pathway names available to users through the ratequery API
  if (rateq::RegBuilder_copied_str_arr1d(reg_builder, "inject_model_names",
                                         inj_path_name_l,
                                         n_pathways) != GR_SUCCESS) {
    return GrPrintAndReturnErr(
        "There was an issue making names of inject pathways queryable");
  }

  // the length of each gas nuclide yield array is equal to the number of
  // injection pathways
  if (rateq::RegBuilder_recipe_1d(reg_builder, 7, &nuclide_gas_yield_recipe,
                                  n_pathways) != GR_SUCCESS) {
    return GrPrintAndReturnErr(
        "There was an issue making nuclide gas yield fractions (for each "
        "injection pathway) queryable");
  }

  if ((my_rates->opaque_storage != nullptr) &&
      (my_rates->opaque_storage->grain_species_info != nullptr)) {
    int n_grain_species =
        my_rates->opaque_storage->grain_species_info->n_species;

    // the length of each grain species yield array is equal to the number of
    // injection pathways
    if (rateq::RegBuilder_recipe_1d(reg_builder, n_grain_species,
                                    &grain_yield_recipe,
                                    n_pathways) != GR_SUCCESS) {
      return GrPrintAndReturnErr(
          "There was an issue making nuclide gas yield fractions (for each "
          "injection pathway) queryable");
    }
  }

  return GR_SUCCESS;
}

/// names of all injection pathways known to grackle listed in the order
/// consistent with the logic for implementing InjectPathFieldPack
///
/// @see SetupCallbackCtx::inj_path_names - The entity that is initialized with
///   the values in this variable. Its docstring provides more details.
constexpr const char* const known_inj_path_names[] = {
    "local_ISM", "ccsn13", "ccsn20", "ccsn25",  "ccsn30",  "fsn13",
    "fsn15",     "fsn50",  "fsn80",  "pisn170", "pisn200", "y19",
};

static_assert(sizeof(known_inj_path_names) ==
                  (sizeof(char*) *
                   grackle::impl::inj_model_input::N_Injection_Pathways),
              "inconsistency b/t number of entries in known_inj_path_names and "
              "grackle::impl::inj_model_input::N_Injection_Pathways");

/// a crude map-like function
bool lookup_metal_yield_ptrs(
    grackle::impl::GrainMetalInjectPathways* inject_pathway_props,
    const char* metal_nuclide_name, double** total_yield_ptr,
    double** gasonly_yield_ptr) {
  if (std::strcmp("C", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.C;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.C;
  } else if (std::strcmp("O", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.O;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.O;
  } else if (std::strcmp("Mg", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.Mg;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.Mg;
  } else if (std::strcmp("Al", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.Al;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.Al;
  } else if (std::strcmp("Si", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.Si;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.Si;
  } else if (std::strcmp("S", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.S;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.S;
  } else if (std::strcmp("Fe", metal_nuclide_name) == 0) {
    (*total_yield_ptr) = inject_pathway_props->total_metal_nuclide_yields.Fe;
    (*gasonly_yield_ptr) = inject_pathway_props->gas_metal_nuclide_yields.Fe;
  } else {
    return false;
  }
  return true;
}

/// the context object for setup_yield_table_callback
struct SetupCallbackCtx {
  /// The object that gets updated by the values that the callback loads
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props;

  /// a counter that is incremented every time setup_yield_table_callback loads
  /// data from an injection pathway
  int counter;

  /// maps the names of known injection pathways to a unique index
  ///
  /// @par Why This Is Needed
  /// For some background context:
  /// - to make use of the injection pathway information, Grackle requires
  ///   users to provide mass densities fields specifying the total amount of
  ///   injected metals by each pathway
  /// - currently, the names of the user accessed fields have hard-coded names
  /// - currently, Grackle loops over the data specified by these fields in a
  ///   standardized order and assumes data within an instance of
  ///   grackle::impl::GrainMetalInjectPathways has the same order
  ///
  /// @par How this should change
  /// Before any functionality involving injection pathways is included in a
  /// public release of Grackle, the plan is to stop hard-coding the names of
  /// injection pathway density fields, and to load in injection pathway data
  /// from HDF5 files.
  /// - At that point, we'll need to tweak the callback function to load in
  ///   data for models with arbitrary names.
  /// - I suspect that we'll want to adopt a policy for ensuring that the order
  ///   of models is well-defined (to make results bitwise reproducible). In
  ///   that scenario, we might use this to enforce an alphanumberic ordering
  const grackle::impl::FrozenKeyIdxBiMap* inj_path_names;

  /// maps the names of the grain species for which data will be loaded to the
  /// appropriate grain species index
  const grackle::impl::FrozenKeyIdxBiMap* grain_species_names;
};

/// a callback function that sets up the appropriate parts of
/// GrainMetalInjectPathways given data for a particular injection pathway
///
/// @param[in] name Name of the current injection pathway
/// @param[in] input Holds the data for the current injection pathway.
/// @param[inout] ctx A pointer to a SetupCallbackCtx instance
///
/// @note
/// This function has C linkage because that's the expectation of the function
/// receiving this callback
extern "C" int setup_yield_table_callback(
    const char* name,
    const grackle::impl::inj_model_input::InjectionPathwayInputData* input,
    void* ctx) {
  namespace inj_input = ::grackle::impl::inj_model_input;
  namespace bimap = ::grackle::impl::bimap;

  SetupCallbackCtx* my_ctx = static_cast<SetupCallbackCtx*>(ctx);

  // lookup the pathway associated with the current injection pathway name
  // and report an error if there is one
  // -> see the docstring for SetupCallbackCtx::inj_path_names for how this
  //    behavior will change when we start loading data from HDF5 files
  bimap::AccessRslt maybe_pathway_idx =
      FrozenKeyIdxBiMap_find(my_ctx->inj_path_names, name);
  if (!maybe_pathway_idx.has_value) {
    return GR_SUCCESS;
  }
  int pathway_idx = static_cast<int>(maybe_pathway_idx.value);

  // load the object that we update with the data we read
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props =
      my_ctx->inject_pathway_props;

  // record the yields for each metal nuclide
  for (int i = 0; i < input->n_metal_nuclide_yields; i++) {
    const inj_input::MetalNuclideYieldProps& yield_info =
        input->metal_nuclide_yields[i];

    double *total_yield, *gas_yield;
    if (!lookup_metal_yield_ptrs(inject_pathway_props, yield_info.name,
                                 &total_yield, &gas_yield)) {
      return GrPrintAndReturnErr("`%s` not a known metal nuclide",
                                 yield_info.name);
    }

    total_yield[pathway_idx] = yield_info.total_yield;
    gas_yield[pathway_idx] = yield_info.gas_yield;
  }

  if (my_ctx->grain_species_names != nullptr) {
    // loop over chunks of yield data. Each chunk holds data for the current
    // pathway that corresponds to a distinct grain species
    for (int yield_idx = 0; yield_idx < input->n_injected_grain_species;
         yield_idx++) {
      const inj_input::GrainSpeciesYieldProps& yield_info =
          input->initial_grain_props[yield_idx];

      bimap::AccessRslt maybe_grain_idx =
          FrozenKeyIdxBiMap_find(my_ctx->grain_species_names, yield_info.name);
      if (!maybe_grain_idx.has_value) {
        continue;
      }
      int grain_species_idx = static_cast<int>(maybe_grain_idx.value);

      // copy the nonprimordial yield fraction
      inject_pathway_props->grain_yields.data[grain_species_idx][pathway_idx] =
          yield_info.nonprimoridal_yield_frac;

      // copy the 1st, 2nd, and 3rd moments of the size distribution
      // (the 0th moment isn't recorded anywhere
      double* size_mom_table =
          inject_pathway_props->size_moments.data[grain_species_idx];
      for (int i = 0; i < 3; i++) {
        size_mom_table[pathway_idx * 3 + i] = yield_info.size_moments[i];
      }

      // copy over the opacity coefficients table
      double* opac_coef_table =
          inject_pathway_props->opacity_coef_table.data[grain_species_idx];
      int n_Td = grackle::impl::inj_model_input::N_Tdust_Opacity_Table;
      int n_coef = grackle::impl::inj_model_input::N_Opacity_Coef;

      for (int i_Td = 0; i_Td < n_Td; i_Td++) {
        for (int i_coef = 0; i_coef < n_coef; i_coef++) {
          int i = (pathway_idx * n_Td * n_coef) + (i_Td * n_coef) + i_coef;
          opac_coef_table[i] = yield_info.opacity_coef_table[i_Td][i_coef];
        }
      }
    }
  }

  (my_ctx->counter)++;

  return GR_SUCCESS;
}

/// a helper function to override the values of all dust inject properties
///
/// @note
/// to start, this only handles the grain yields
void zero_out_dust_inject_props(
    grackle::impl::GrainMetalInjectPathways* inject_pathway_props) {
  const int n_grain_species = OnlyGrainSpLUT::NUM_ENTRIES;

  int n_pathways = inject_pathway_props->n_pathways;

  // handle the grain yields
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
       grain_species_idx++) {
    double* arr = inject_pathway_props->grain_yields.data[grain_species_idx];
    for (int iSN = 0; iSN < n_pathways; iSN++) {
      arr[iSN] = 0.0;
    }
  }

  // handle the size moments
  int moment_table_len = n_pathways * 3;  // there are 3 size moments
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
       grain_species_idx++) {
    double* arr = inject_pathway_props->size_moments.data[grain_species_idx];
    for (int i = 0; i < moment_table_len; i++) {
      arr[i] = 0.0;
    }
  }

  // handle the opacity coefficient table
  long long opac_table_len =
      (static_cast<long long>(n_pathways) *
       static_cast<long long>(inject_pathway_props->n_opac_poly_coef) *
       inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
       grain_species_idx++) {
    double* arr =
        inject_pathway_props->opacity_coef_table.data[grain_species_idx];
    for (long long i = 0LL; i < opac_table_len; i++) {
      arr[i] = 0.0;
    }
  }
}

}  // anonymous namespace

int grackle::impl::load_inject_path_data(const chemistry_data* my_chemistry,
                                         chemistry_data_storage* my_rates,
                                         ratequery::RegBuilder* reg_builder) {
  // Currently, this function "loads" injection pathway data from data that is
  // directly embedded as part of the Grackle library in a manner controlled
  // by raw_data.hpp and raw_data.cpp
  //
  // Longer term, the goal is to "load" the data from an external HDF5 file

  if (my_chemistry->metal_chemistry == 0) {
    return GR_SUCCESS;
  }

  // an upper bound on the number of allowed injection pathways
  int max_n_pathways = grackle::impl::inj_model_input::N_Injection_Pathways;

  // get the list of injection pathways
  // -> currently this requires us to look at the my_chemistry->multi_metals
  //    and my_chemstry->metal_abundances variables.
  // -> when we move away from hardcoded names and start loading from HDF5,
  //    we'll remove these variables
  const char* const* inj_path_name_l = nullptr;
  int n_pathways = 0;
  bool valid_metal_abundances =
      ((0 <= my_chemistry->metal_abundances) &&
       (my_chemistry->metal_abundances < max_n_pathways));

  if ((my_chemistry->multi_metals == 0) && !valid_metal_abundances) {
    return GrPrintAndReturnErr(
        "the metal_abundances parameter must not be negative or exceed %d",
        max_n_pathways - 1);
  } else if (my_chemistry->multi_metals == 0) {
    inj_path_name_l = known_inj_path_names + my_chemistry->metal_abundances;
    n_pathways = 1;
  } else if (my_chemistry->multi_metals == 1) {
    inj_path_name_l = known_inj_path_names;
    n_pathways = max_n_pathways;
  } else {
    return GrPrintAndReturnErr("the multi_metals parameter isn't 0 or 1");
  }

  // construct a mapping of the injection pathway names
  // -> right now, since the strings are statically allocated, we use
  //    BiMapMode::REFS_KEYDATA to instruct the map to avoid making copies.
  // -> In the future, when model names are dynamically specified by an HDF5
  //    file, we'll need to use BiMapMode::COPIES_KEYDATA.
  FrozenKeyIdxBiMap inj_path_names = new_FrozenKeyIdxBiMap(
      inj_path_name_l, n_pathways, BiMapMode::REFS_KEYDATA);
  if (!FrozenKeyIdxBiMap_is_ok(&inj_path_names)) {
    return GrPrintAndReturnErr(
        "there was a problem building the map of model names");
  }

  // initialize the object that will hold the loaded data
  int n_log10Tdust_vals = grackle::impl::inj_model_input::N_Tdust_Opacity_Table;
  int n_opac_poly_coef = grackle::impl::inj_model_input::N_Opacity_Coef;
  my_rates->opaque_storage->inject_pathway_props =
      new grackle::impl::GrainMetalInjectPathways;
  *(my_rates->opaque_storage->inject_pathway_props) =
      new_GrainMetalInjectPathways(n_pathways, n_log10Tdust_vals,
                                   n_opac_poly_coef);

  grackle::impl::GrainMetalInjectPathways* inject_pathway_props =
      my_rates->opaque_storage->inject_pathway_props;

  if (!GrainMetalInjectPathways_is_valid(inject_pathway_props)) {
    drop_FrozenKeyIdxBiMap(&inj_path_names);
    return GR_FAIL;
  }

  // initialize the grid of dust temperatures associated with the opacity table
  double log10Tdust_lo = 0.0;
  double log10Tdust_step = 0.1;

  inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0] =
      log10Tdust_step;
  double* log10Tdust_vals =
      inject_pathway_props->log10Tdust_interp_props.parameters[0];
  for (int iTd = 0; iTd < n_log10Tdust_vals; iTd++) {
    log10Tdust_vals[iTd] = log10Tdust_lo + (double)iTd * log10Tdust_step;
  }

  // zero-out all metal injection yield fractions and dust grain properties
  yields::MetalTables_zero_out(
      &inject_pathway_props->total_metal_nuclide_yields, n_pathways);
  yields::MetalTables_zero_out(&inject_pathway_props->gas_metal_nuclide_yields,
                               n_pathways);
  zero_out_dust_inject_props(inject_pathway_props);

  // actually load in the data for each injection pathway
  GrainSpeciesInfo* grain_species_info =
      my_rates->opaque_storage->grain_species_info;
  const FrozenKeyIdxBiMap* grain_species_names =
      (grain_species_info == nullptr) ? nullptr : &grain_species_info->name_map;

  SetupCallbackCtx ctx = {
      /* inject_pathway_props = */ my_rates->opaque_storage
          ->inject_pathway_props,
      /* counter = */ 0,
      /* inj_path_names = */ &inj_path_names,
      /* grain_species_names = */ grain_species_names,
  };

  int ret = grackle::impl::inj_model_input::input_inject_model_iterate(
      &setup_yield_table_callback, static_cast<void*>(&ctx));

  drop_FrozenKeyIdxBiMap(&inj_path_names);
  if (ret != GR_SUCCESS) {
    return GrPrintAndReturnErr(
        "some kind of unspecified error occured when loading data from each "
        "injection pathway");
  } else if (ctx.counter != n_pathways) {
    return GrPrintAndReturnErr(
        "Only loaded data for %d of the %d available pathways", ctx.counter,
        n_pathways);
  }

  // Before we finish, let's register some of these it can be queried by
  // Grackle-users
  if (configure_RegBuilder(my_rates, reg_builder, inj_path_name_l,
                           n_pathways) != GR_SUCCESS) {
    return GrPrintAndReturnErr("error making injection pathway info queryable");
  }

  return GR_SUCCESS;
}
