//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the functions to initialize the dust yields
///
//===----------------------------------------------------------------------===//

#include <cstring>
#include <stdlib.h>
#include <hdf5.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_chemistry_data.h"
#include "initialize_dust_yields.hpp" // forward declarations
#include "inject_model/raw_data.hpp"
#include "LUT.hpp"
#include "opaque_storage.hpp"
#include "status_reporting.h" // GrPrintAndReturnErr

namespace {  // stuff inside an anonymous namespace is local to this file

int lookup_pathway_idx(const char* name) {
  const char * known_names[12] = {
    "local_ISM", "ccsn13", "ccsn20", "ccsn25", "ccsn30", "fsn13",
    "fsn15", "fsn50", "fsn80", "pisn170", "pisn200", "y19",
  };

  for (int i = 0; i < 12; i++){
    if (std::strcmp(known_names[i], name) == 0){
      return i;
    }
  }
  return -1;
}

/// the context object for setup_yield_table_callback
struct SetupCallbackCtx{
  chemistry_data_storage *my_rates;
  int counter;
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
    const grackle::impl::inj_model_input::InjectionPathwayInputData *input,
    void* ctx)
{
  namespace inj_input = ::grackle::impl::inj_model_input;


  int pathway_idx = lookup_pathway_idx(name);

  if (pathway_idx < 0) {
    return GrPrintAndReturnErr(
      "`%s` is an unexpected injection pathway name", name);
  }

  chemistry_data_storage *my_rates =
    static_cast<SetupCallbackCtx*>(ctx)->my_rates;

  grackle::impl::GrainMetalInjectPathways* inject_pathway_props
    = my_rates->opaque_storage->inject_pathway_props;

  // record each metal nuclide yield
  // -> there is less value to using string keys in this case, but it makes
  //    some sense to be semi-consistent with the handling of the dust species
  //    yields
  for (int i = 0; i < input->n_metal_nuclide_yields; i++) {
    const inj_input::MetalNuclideYieldProps& yield_info =
      input->metal_nuclide_yields[i];

    double* total_yield = nullptr;
    double* gas_yield = nullptr;

    // todo: refactor to use a map (I have a rough plan)
    if (std::strcmp("C", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.C;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.C;
    } else if (std::strcmp("O", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.O;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.O;
    } else if (std::strcmp("Mg", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.Mg;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.Mg;
    } else if (std::strcmp("Al", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.Al;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.Al;
    } else if (std::strcmp("Si", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.Si;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.Si;
    } else if (std::strcmp("S", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.S;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.S;
    } else if (std::strcmp("Fe", yield_info.name) == 0) {
      total_yield = inject_pathway_props->total_metal_nuclide_yields.Fe;
      gas_yield = inject_pathway_props->gas_metal_nuclide_yields.Fe;
    } else {
      return GrPrintAndReturnErr(
        "`%s` not a known metal nuclide", yield_info.name);
    }

    total_yield[pathway_idx] = yield_info.total_yield;
    gas_yield[pathway_idx] = yield_info.gas_yield;
  }

  // record each grain species yield
  for (int yield_idx = 0; yield_idx < input->n_injected_grain_species;
       yield_idx++) {
    const inj_input::GrainSpeciesYieldProps& yield_info =
      input->initial_grain_props[yield_idx];

    int grain_species_idx = -1;

    // with a little refactoring, this will get a lot more concise
    if (std::strcmp("MgSiO3_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::MgSiO3_dust;
    } else if (std::strcmp("AC_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::AC_dust;
    } else if (std::strcmp("SiM_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::SiM_dust;
    } else if (std::strcmp("FeM_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::FeM_dust;
    } else if (std::strcmp("Mg2SiO4_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Mg2SiO4_dust;
    } else if (std::strcmp("Fe3O4_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Fe3O4_dust;
    } else if (std::strcmp("SiO2_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::SiO2_dust;
    } else if (std::strcmp("MgO_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::MgO_dust;
    } else if (std::strcmp("FeS_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::FeS_dust;
    } else if (std::strcmp("Al2O3_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Al2O3_dust;
    } else if (std::strcmp("ref_org_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::ref_org_dust;
    } else if (std::strcmp("vol_org_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::vol_org_dust;
    } else if (std::strcmp("H2O_ice_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::H2O_ice_dust;
    } else {
      return GrPrintAndReturnErr(
        "`%s` not a known grain species", yield_info.name);
    }

    // copy the nonprimordial yield fraction
    inject_pathway_props->grain_yields.data[grain_species_idx][pathway_idx]
      = yield_info.nonprimoridal_yield_frac;

    // copy the 1st, 2nd, and 3rd moments of the size distribution
    // (the 0th moment isn't recorded anywhere
    double* size_mom_table =
      inject_pathway_props->size_moments.data[grain_species_idx];
    for (int i = 0; i < 3; i++) {
      size_mom_table[pathway_idx*3+i] = yield_info.size_moments[i];
    }

    // copy over the opacity coefficients table
    {
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

  (static_cast<SetupCallbackCtx*>(ctx)->counter)++;

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
      grain_species_idx++){
    double* arr = inject_pathway_props->grain_yields.data[grain_species_idx];
    for(int iSN = 0; iSN < n_pathways; iSN++) {
      arr[iSN] = 0.0;
    }
  }

  // handle the size moments
  int moment_table_len = n_pathways * 3; // there are 3 size moments
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
      grain_species_idx++){
    double* arr = inject_pathway_props->size_moments.data[grain_species_idx];
    for(int i = 0; i < moment_table_len; i++) {
      arr[i] = 0.0;
    }
  }

  // handle the opacity coefficient table
  long long opac_table_len = (
      static_cast<long long>(n_pathways) *
      static_cast<long long>(inject_pathway_props->n_opac_poly_coef) *
      inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
      grain_species_idx++){
    double* arr =
      inject_pathway_props->opacity_coef_table.data[grain_species_idx];
    for(long long i = 0LL; i < opac_table_len; i++) {
      arr[i] = 0.0;
    }
  }
}

}  // anonymous namespace

int grackle::impl::initialize_dust_yields(chemistry_data *my_chemistry,
                                          chemistry_data_storage *my_rates,
                                          code_units *my_units)
{

  if (my_chemistry->metal_chemistry == 0) {
    return SUCCESS;
  }

  int n_pathways = 12;
  int n_log10Tdust_vals = grackle::impl::inj_model_input::N_Tdust_Opacity_Table;
  int n_opac_poly_coef = grackle::impl::inj_model_input::N_Opacity_Coef;
  my_rates->opaque_storage->inject_pathway_props =
    new grackle::impl::GrainMetalInjectPathways;
  *(my_rates->opaque_storage->inject_pathway_props) =
    new_GrainMetalInjectPathways(n_pathways, n_log10Tdust_vals,
                                 n_opac_poly_coef);

  grackle::impl::GrainMetalInjectPathways* inject_pathway_props
    = my_rates->opaque_storage->inject_pathway_props;

  if (!GrainMetalInjectPathways_is_valid(inject_pathway_props)) {
    return GR_FAIL;
  }

  double log10Tdust_lo = 0.0;
  double log10Tdust_step = 0.1;

  inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0]
    = log10Tdust_step;
  double* log10Tdust_vals =
    inject_pathway_props->log10Tdust_interp_props.parameters[0];
  for(int iTd = 0; iTd < n_log10Tdust_vals; iTd++) {
    log10Tdust_vals[iTd] = log10Tdust_lo + (double)iTd * log10Tdust_step;
  }


  // zero-out all metal injection yield fractions and dust grain properties
  yields::MetalTables_zero_out(
      &inject_pathway_props->total_metal_nuclide_yields, n_pathways);
  yields::MetalTables_zero_out(
      &inject_pathway_props->gas_metal_nuclide_yields, n_pathways);
  zero_out_dust_inject_props(inject_pathway_props);

  SetupCallbackCtx ctx = {my_rates, 0};

  int ret = grackle::impl::inj_model_input::input_inject_model_iterate(
      &setup_yield_table_callback, static_cast<void*>(&ctx));
  if (ret != GR_SUCCESS) {
    return GR_FAIL;
  }

  GR_INTERNAL_REQUIRE(ctx.counter == 12, "sanity-checking");

  return GR_SUCCESS;
}
