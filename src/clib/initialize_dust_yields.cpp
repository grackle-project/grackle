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
    double* opac_coef_table = nullptr;

    // with a little refactoring, this will get a lot more concise
    if (std::strcmp("MgSiO3_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::MgSiO3_dust;
      opac_coef_table = my_rates->SN0_kpMgSiO3;
    } else if (std::strcmp("AC_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::AC_dust;
      opac_coef_table = my_rates->SN0_kpAC;
    } else if (std::strcmp("SiM_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::SiM_dust;
      opac_coef_table = my_rates->SN0_kpSiM;
    } else if (std::strcmp("FeM_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::FeM_dust;
      opac_coef_table = my_rates->SN0_kpFeM;
    } else if (std::strcmp("Mg2SiO4_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Mg2SiO4_dust;
      opac_coef_table = my_rates->SN0_kpMg2SiO4;
    } else if (std::strcmp("Fe3O4_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Fe3O4_dust;
      opac_coef_table = my_rates->SN0_kpFe3O4;
    } else if (std::strcmp("SiO2_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::SiO2_dust;
      opac_coef_table = my_rates->SN0_kpSiO2D;
    } else if (std::strcmp("MgO_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::MgO_dust;
      opac_coef_table = my_rates->SN0_kpMgO;
    } else if (std::strcmp("FeS_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::FeS_dust;
      opac_coef_table = my_rates->SN0_kpFeS;
    } else if (std::strcmp("Al2O3_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::Al2O3_dust;
      opac_coef_table = my_rates->SN0_kpAl2O3;
    } else if (std::strcmp("ref_org_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::ref_org_dust;
      opac_coef_table = my_rates->SN0_kpreforg;
    } else if (std::strcmp("vol_org_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::vol_org_dust;
      opac_coef_table = my_rates->SN0_kpvolorg;
    } else if (std::strcmp("H2O_ice_dust", yield_info.name) == 0) {
      grain_species_idx = OnlyGrainSpLUT::H2O_ice_dust;
      opac_coef_table = my_rates->SN0_kpH2Oice;
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

/// a helper function override the values of all metal injection properties
void override_metal_inject_props(
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props,
  double value, int n_pathways)
{
  for(int iSN = 0; iSN < n_pathways; iSN++) {
    inject_pathway_props->total_metal_nuclide_yields.C [iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.O [iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.Mg[iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.Al[iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.Si[iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.S [iSN] = value;
    inject_pathway_props->total_metal_nuclide_yields.Fe[iSN] = value;

    inject_pathway_props->gas_metal_nuclide_yields.C [iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.O [iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.Mg[iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.Al[iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.Si[iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.S [iSN] = value;
    inject_pathway_props->gas_metal_nuclide_yields.Fe[iSN] = value;
  }
}

/// a helper function to override the values of all dust inject properties
///
/// @note
/// to start, this only handles the grain yields
void override_dust_inject_props(chemistry_data_storage *my_rates,
                                double value, int n_pathways) {

  GRIMPL_REQUIRE(my_rates != nullptr && my_rates->opaque_storage != nullptr,
                 "sanity check -- should never ever fail!");
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props
    = my_rates->opaque_storage->inject_pathway_props;

  GRIMPL_REQUIRE(inject_pathway_props != nullptr,
                 "sanity check -- should never ever fail!");

  const int n_grain_species = OnlyGrainSpLUT::NUM_ENTRIES;

  // handle the grain yields
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
      grain_species_idx++){
    double* arr = inject_pathway_props->grain_yields.data[grain_species_idx];
    for(int iSN = 0; iSN < n_pathways; iSN++) {
      arr[iSN] = value;
    }
  }

  // handle the size moments
  int moment_table_len = n_pathways * 3; // there are 3 size moments
  for (int grain_species_idx = 0; grain_species_idx < n_grain_species;
      grain_species_idx++){
    double* arr = inject_pathway_props->size_moments.data[grain_species_idx];
    for(int i = 0; i < moment_table_len; i++) {
      arr[i] = value;
    }
  }

  if (true) { // NOLINT(readability-simplify-boolean-expr)

    // handle the opacity coefficient table
    int opac_table_len = n_pathways * my_rates->gr_Size;
    for(int i = 0; i < opac_table_len; i++) {
      my_rates->SN0_kpSiM     [i] = value;
      my_rates->SN0_kpFeM     [i] = value;
      my_rates->SN0_kpMg2SiO4 [i] = value;
      my_rates->SN0_kpMgSiO3  [i] = value;
      my_rates->SN0_kpFe3O4   [i] = value;
      my_rates->SN0_kpAC      [i] = value;
      my_rates->SN0_kpSiO2D   [i] = value;
      my_rates->SN0_kpMgO     [i] = value;
      my_rates->SN0_kpFeS     [i] = value;
      my_rates->SN0_kpAl2O3   [i] = value;
      my_rates->SN0_kpreforg  [i] = value;
      my_rates->SN0_kpvolorg  [i] = value;
      my_rates->SN0_kpH2Oice  [i] = value;
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


  int NSN = n_pathways;  // todo: delete me!
  my_rates->SN0_N = n_pathways;

  // crude hack (holding onto this momentarily)
  my_rates->SN0_fC  = inject_pathway_props->gas_metal_nuclide_yields.C;
  my_rates->SN0_fO  = inject_pathway_props->gas_metal_nuclide_yields.O;
  my_rates->SN0_fMg = inject_pathway_props->gas_metal_nuclide_yields.Mg;
  my_rates->SN0_fAl = inject_pathway_props->gas_metal_nuclide_yields.Al;
  my_rates->SN0_fSi = inject_pathway_props->gas_metal_nuclide_yields.Si;
  my_rates->SN0_fS  = inject_pathway_props->gas_metal_nuclide_yields.S;
  my_rates->SN0_fFe = inject_pathway_props->gas_metal_nuclide_yields.Fe;

  my_rates->SN0_r0SiM     = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::SiM_dust];
  my_rates->SN0_r0FeM     = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::FeM_dust];
  my_rates->SN0_r0Mg2SiO4 = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Mg2SiO4_dust];
  my_rates->SN0_r0MgSiO3  = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::MgSiO3_dust];
  my_rates->SN0_r0Fe3O4   = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Fe3O4_dust];
  my_rates->SN0_r0AC      = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::AC_dust];
  my_rates->SN0_r0SiO2D   = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::SiO2_dust];
  my_rates->SN0_r0MgO     = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::MgO_dust];
  my_rates->SN0_r0FeS     = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::FeS_dust];
  my_rates->SN0_r0Al2O3   = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::Al2O3_dust];
  my_rates->SN0_r0reforg  = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::ref_org_dust];
  my_rates->SN0_r0volorg  = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::vol_org_dust];
  my_rates->SN0_r0H2Oice  = inject_pathway_props->size_moments.data[OnlyGrainSpLUT::H2O_ice_dust];

  // write out the opacity related quantities

  // todo: consider renaming Nmom -> Ncoef
  int Nmom = n_opac_poly_coef; // todo: remove me!
  // todo: more this into GrainMetalInjectPathways
  // - essentially, each SN0_kpGRSP array is a 3D array of shape
  //   (n_pathways, NTd, Nmom). This shape uses numpy conventions (i.e. Nmom is
  //   the fast-axis).
  // - In more detail:
  //   - n_pathways is the number of injection pathways
  //   - NTd corresponds to the number of dust temperatures we consider
  //   - Nmom corresponds to the number of coefficients that we track for
  //     computing an opacity related quantity. This quantity is given by
  //       4πζ/3 * ( SN0_kpGRSP[path_j, iTd, 3] +
  //                 SN0_kpGRSP[path_j, iTd, 2] * δr(t) +
  //                 SN0_kpGRSP[path_j, iTd, 1] * δr²(t) +
  //                 SN0_kpGRSP[path_j, iTd, 0] * δr³(t))
  //     where
  //       - ζ is the mass density of a single grain (in g/cm³)
  //       - δr(t) refers to the derived "size increment" (it is a central
  //         quantity in the model)
  //       - I **think** the resulting quantity is the optical cross-section
  double NTd = n_log10Tdust_vals;
  double Td0 = 0.0000000; // todo: remove me!
  double dTd = 0.1000000; // todo: remove me!

  // todo: rename gr_Td, dTd since they are related to log10(Tdust) or ln(Tdust)
      my_rates->gr_Td = (double*)malloc(NTd * Nmom * sizeof(double));
      my_rates->gr_Size = NTd * Nmom;
      my_rates->gr_N[0] = Nmom;
      my_rates->gr_N[1] = NTd;
      my_rates->gr_dT   = dTd;
      for(int iTd = 0; iTd < NTd; iTd++) {
        my_rates->gr_Td[iTd] = Td0 + (double)iTd * dTd;
      }

      my_rates->SN0_kpSiM      = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeM      = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMg2SiO4  = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgSiO3   = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFe3O4    = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAC       = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpSiO2D    = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgO      = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeS      = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAl2O3    = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpreforg   = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpvolorg   = (double*)malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpH2Oice   = (double*)malloc(NSN * Nmom * NTd * sizeof(double));


  // zero-out all metal injection yield fractions and dust grain properties
  override_metal_inject_props(inject_pathway_props, 0.0, n_pathways);
  override_dust_inject_props(my_rates, 0.0, n_pathways);

  SetupCallbackCtx ctx = {my_rates, 0};

  int ret = grackle::impl::inj_model_input::input_inject_model_iterate(
      &setup_yield_table_callback, static_cast<void*>(&ctx));
  if (ret != GR_SUCCESS) {
    return GR_FAIL;
  }

  GR_INTERNAL_REQUIRE(ctx.counter == 12, "sanity-checking");

  return GR_SUCCESS;
}

int grackle::impl::free_dust_yields(chemistry_data *my_chemistry,
                                    chemistry_data_storage *my_rates)
{

  if (my_chemistry->metal_chemistry == 0)
    return SUCCESS;

  // reminder: we are holding on to SN0_f<nuclide> as a crude hack (we will
  // delete them later)
  my_rates->SN0_fC = nullptr;
  my_rates->SN0_fO = nullptr;
  my_rates->SN0_fMg = nullptr;
  my_rates->SN0_fAl = nullptr;
  my_rates->SN0_fSi = nullptr;
  my_rates->SN0_fS = nullptr;
  my_rates->SN0_fFe = nullptr;

  my_rates->SN0_r0SiM = nullptr;
  my_rates->SN0_r0FeM = nullptr;
  my_rates->SN0_r0Mg2SiO4 = nullptr;
  my_rates->SN0_r0MgSiO3 = nullptr;
  my_rates->SN0_r0Fe3O4 = nullptr;
  my_rates->SN0_r0AC = nullptr;
  my_rates->SN0_r0SiO2D = nullptr;
  my_rates->SN0_r0MgO = nullptr;
  my_rates->SN0_r0FeS = nullptr;
  my_rates->SN0_r0Al2O3 = nullptr;
  my_rates->SN0_r0reforg = nullptr;
  my_rates->SN0_r0volorg = nullptr;
  my_rates->SN0_r0H2Oice = nullptr;

  GRACKLE_FREE(my_rates->gr_Td);

  GRACKLE_FREE(my_rates->SN0_kpSiM);
  GRACKLE_FREE(my_rates->SN0_kpFeM);
  GRACKLE_FREE(my_rates->SN0_kpMg2SiO4);
  GRACKLE_FREE(my_rates->SN0_kpMgSiO3);
  GRACKLE_FREE(my_rates->SN0_kpFe3O4);
  GRACKLE_FREE(my_rates->SN0_kpAC);
  GRACKLE_FREE(my_rates->SN0_kpSiO2D);
  GRACKLE_FREE(my_rates->SN0_kpMgO);
  GRACKLE_FREE(my_rates->SN0_kpFeS);
  GRACKLE_FREE(my_rates->SN0_kpAl2O3);
  GRACKLE_FREE(my_rates->SN0_kpreforg);
  GRACKLE_FREE(my_rates->SN0_kpvolorg);
  GRACKLE_FREE(my_rates->SN0_kpH2Oice);

  return SUCCESS;
}

