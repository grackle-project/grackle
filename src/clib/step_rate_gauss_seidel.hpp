//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements step_rate_gauss_seidel
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// step_rate_g function from FORTRAN to C++

#ifndef STEP_RATE_GAUSS_SEIDEL_HPP
#define STEP_RATE_GAUSS_SEIDEL_HPP

#include <cmath>
#include <cstdio>

#include "grackle.h"  // gr_float
#include "chemistry_solver_funcs.hpp"
#include "chem_model/infer_species.hpp"
#include "fortran_func_decls.h"  // gr_mask_type
#include "internal_types.hpp"
#include "LUT.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {

namespace step_gs_detail{

/// This function overwrites the species density field values using the values
/// in @p species_tmpdens, which holds the values from the end of the current
/// (sub-)cycle.
inline void update_from_tmpdens(
  const double* dtit, IndexRange idx_range, double* dedot_prev,
  double* HIdot_prev, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  grackle_field_data* my_fields,
  grackle::impl::SpeciesCollection species_tmpdens
) {

  // Get pointers to the various species fields
  // ------------------------------------------
  // -> TODO: stop doing this!
  int ptr_offset_ = layoutleft_3D_index_to_1D_(
      my_fields->grid_dimension, 0, idx_range.j, idx_range.k);

  gr_float* de = my_fields->e_density + ptr_offset_;
  gr_float* HI = my_fields->HI_density + ptr_offset_;
  gr_float* HII = my_fields->HII_density + ptr_offset_;
  gr_float* HeI = my_fields->HeI_density + ptr_offset_;
  gr_float* HeII = my_fields->HeII_density + ptr_offset_;
  gr_float* HeIII = my_fields->HeIII_density + ptr_offset_;
  gr_float* HM = my_fields->HM_density + ptr_offset_;
  gr_float* H2I = my_fields->H2I_density + ptr_offset_;
  gr_float* H2II = my_fields->H2II_density + ptr_offset_;
  gr_float* DI = my_fields->DI_density + ptr_offset_;
  gr_float* DII = my_fields->DII_density + ptr_offset_;
  gr_float* HDI = my_fields->HDI_density + ptr_offset_;
  gr_float* DM = my_fields->DM_density + ptr_offset_;
  gr_float* HDII = my_fields->HDII_density + ptr_offset_;
  gr_float* HeHII = my_fields->HeHII_density + ptr_offset_;
  gr_float* CI = my_fields->CI_density + ptr_offset_;
  gr_float* CII = my_fields->CII_density + ptr_offset_;
  gr_float* COI = my_fields->COI_density + ptr_offset_;
  gr_float* CO2I = my_fields->CO2I_density + ptr_offset_;
  gr_float* OI = my_fields->OI_density + ptr_offset_;
  gr_float* OHI = my_fields->OHI_density + ptr_offset_;
  gr_float* H2OI = my_fields->H2OI_density + ptr_offset_;
  gr_float* O2I = my_fields->O2I_density + ptr_offset_;
  gr_float* SiI = my_fields->SiI_density + ptr_offset_;
  gr_float* SiOI = my_fields->SiOI_density + ptr_offset_;
  gr_float* SiO2I = my_fields->SiO2I_density + ptr_offset_;
  gr_float* CHI = my_fields->CHI_density + ptr_offset_;
  gr_float* CH2I = my_fields->CH2I_density + ptr_offset_;
  gr_float* COII = my_fields->COII_density + ptr_offset_;
  gr_float* OII = my_fields->OII_density + ptr_offset_;
  gr_float* OHII = my_fields->OHII_density + ptr_offset_;
  gr_float* H2OII = my_fields->H2OII_density + ptr_offset_;
  gr_float* H3OII = my_fields->H3OII_density + ptr_offset_;
  gr_float* O2II = my_fields->O2II_density + ptr_offset_;
  gr_float* Mg = my_fields->Mg_density + ptr_offset_;
  gr_float* Al = my_fields->Al_density + ptr_offset_;
  gr_float* S = my_fields->S_density + ptr_offset_;
  gr_float* Fe = my_fields->Fe_density + ptr_offset_;
  gr_float* SiM = my_fields->SiM_dust_density + ptr_offset_;
  gr_float* FeM = my_fields->FeM_dust_density + ptr_offset_;
  gr_float* Mg2SiO4 = my_fields->Mg2SiO4_dust_density + ptr_offset_;
  gr_float* MgSiO3 = my_fields->MgSiO3_dust_density + ptr_offset_;
  gr_float* Fe3O4 = my_fields->Fe3O4_dust_density + ptr_offset_;
  gr_float* AC = my_fields->AC_dust_density + ptr_offset_;
  gr_float* SiO2D = my_fields->SiO2_dust_density + ptr_offset_;
  gr_float* MgO = my_fields->MgO_dust_density + ptr_offset_;
  gr_float* FeS = my_fields->FeS_dust_density + ptr_offset_;
  gr_float* Al2O3 = my_fields->Al2O3_dust_density + ptr_offset_;
  gr_float* reforg = my_fields->ref_org_dust_density + ptr_offset_;
  gr_float* volorg = my_fields->vol_org_dust_density + ptr_offset_;
  gr_float* H2Oice = my_fields->H2O_ice_dust_density + ptr_offset_;

  // --- (E) Set densities from 1D temps to 3D fields ---

  // flag for if Gen Chiaki's dust model is enabled
  // -> to simplify bookkeeping in other parts of the code, we ignore whether
  //    the grain_growth flag was enabled, (the logic for setting the relevant
  //    reaction rates, explicitly sets the rates to 0 if disabled)
  const bool chiaki_model_dust_evolution = my_chemistry->dust_chemistry == 2;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {
      HIdot_prev[i] = std::fabs(HI[i]-species_tmpdens.data[SpLUT::HI][i]) /
              std::fmax((double)(dtit[i] ), tiny8);
      HI[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HI][i] ), tiny_fortran_val);
      HII[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HII][i] ), tiny_fortran_val);
      HeI[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeI][i] ), tiny_fortran_val);
      HeII[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeII][i] ), tiny_fortran_val);
      HeIII[i] = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeIII][i] ), (gr_float)(1e-5)*tiny_fortran_val);

      // temporarily store the electron density from the start of the current
      // subcycle.
      dedot_prev[i] = de[i];

      // Use charge conservation to determine electron fraction
      // -> in other words, we ignore species_tmpdens.data[SpLUT::e] (in
      //    practice, I think that array holds garbage values)
      de[i] = HII[i] + HeII[i]/(gr_float)(4.) +
           HeIII[i]/(gr_float)(2.);
      if (my_chemistry->primordial_chemistry > 1)
           { de[i] = de[i] - HM[i] + H2II[i]/(gr_float)(2.); }

      if (my_chemistry->primordial_chemistry > 2)
           { de[i] = de[i] + DII[i]/(gr_float)(2.); }
      if (my_chemistry->primordial_chemistry > 3)
           { de[i] = de[i] - DM[i]/(gr_float)(2.)
                + HDII[i]/(gr_float)(3.) + HeHII[i]/(gr_float)(5.); }
      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i] != MASK_FALSE) )
           { de[i] = de[i]
                + CII[i]/(gr_float)(12.) + COII[i]/(gr_float)(28.)
                + OII[i]/(gr_float)(16.) + OHII[i]/(gr_float)(17.)
                + H2OII[i]/(gr_float)(18.) + H3OII[i]/(gr_float)(19.)
                + O2II[i]/(gr_float)(32.); }

      // store the time-derivative of the electron-density in dedot
      // (don't forget that we previously stored the value from the start of
      // the current cycle within dedot_prev)
      dedot_prev[i] = std::fabs(de[i]-dedot_prev[i])/
           std::fmax(dtit[i],tiny8);

      if (my_chemistry->primordial_chemistry > 1)  {
        HM[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HM][i] ), tiny_fortran_val);
        H2I[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2I][i]), tiny_fortran_val);
        H2II[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2II][i] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 2)  {
        DI[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DI][i] ), tiny_fortran_val);
        DII[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DII][i] ), tiny_fortran_val);
        HDI[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDI][i] ), tiny_fortran_val);
      }

      if (my_chemistry->primordial_chemistry > 3)  {
        DM[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::DM][i] ), tiny_fortran_val);
        HDII[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HDII][i] ), tiny_fortran_val);
        HeHII[i] = std::fmax((gr_float)(species_tmpdens.data[SpLUT::HeHII][i] ), tiny_fortran_val);
      }

      if ( (my_chemistry->metal_chemistry == 1)  && 
           (itmask_metal[i] != MASK_FALSE) )  {
        CI[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CI][i]      ), tiny_fortran_val);
        CII[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CII][i]     ), tiny_fortran_val);
        COI[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::COI][i]     ), tiny_fortran_val);
        CO2I[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CO2I][i]    ), tiny_fortran_val);
        OI[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OI][i]      ), tiny_fortran_val);
        OHI[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OHI][i]     ), tiny_fortran_val);
        H2OI[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2OI][i]    ), tiny_fortran_val);
        O2I[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2I][i]     ), tiny_fortran_val);
        SiI[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiI][i]     ), tiny_fortran_val);
        SiOI[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiOI][i]    ), tiny_fortran_val);
        SiO2I[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2I][i]   ), tiny_fortran_val);
        CHI[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CHI][i]     ), tiny_fortran_val);
        CH2I[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::CH2I][i]    ), tiny_fortran_val);
        COII[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::COII][i]    ), tiny_fortran_val);
        OII[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OII][i]     ), tiny_fortran_val);
        OHII[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::OHII][i]    ), tiny_fortran_val);
        H2OII[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2OII][i]   ), tiny_fortran_val);
        H3OII[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H3OII][i]   ), tiny_fortran_val);
        O2II[i]    = std::fmax((gr_float)(species_tmpdens.data[SpLUT::O2II][i]    ), tiny_fortran_val);
      }

      if (chiaki_model_dust_evolution && itmask_metal[i] != MASK_FALSE) {

        if (my_chemistry->dust_species > 0)  {
          Mg[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg][i]      ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 1)  {
          Al[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al][i]      ), tiny_fortran_val);
          S[i]       = std::fmax((gr_float)(species_tmpdens.data[SpLUT::S][i]       ), tiny_fortran_val);
          Fe[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe][i]      ), tiny_fortran_val);
        }

        if (my_chemistry->dust_species > 0)  {
          MgSiO3[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgSiO3_dust][i]  ), tiny_fortran_val);
          AC[i]      = std::fmax((gr_float)(species_tmpdens.data[SpLUT::AC_dust][i]      ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 1)  {
          SiM[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiM_dust][i]     ), tiny_fortran_val);
          FeM[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeM_dust][i]     ), tiny_fortran_val);
          Mg2SiO4[i] = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Mg2SiO4_dust][i] ), tiny_fortran_val);
          Fe3O4[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Fe3O4_dust][i]   ), tiny_fortran_val);
          SiO2D[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::SiO2_dust][i]   ), tiny_fortran_val);
          MgO[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::MgO_dust][i]     ), tiny_fortran_val);
          FeS[i]     = std::fmax((gr_float)(species_tmpdens.data[SpLUT::FeS_dust][i]     ), tiny_fortran_val);
          Al2O3[i]   = std::fmax((gr_float)(species_tmpdens.data[SpLUT::Al2O3_dust][i]   ), tiny_fortran_val);
        }
        if (my_chemistry->dust_species > 2)  {
          reforg[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::ref_org_dust][i]   ), tiny_fortran_val);
          volorg[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::vol_org_dust][i]   ), tiny_fortran_val);
          H2Oice[i]  = std::fmax((gr_float)(species_tmpdens.data[SpLUT::H2O_ice_dust][i]   ), tiny_fortran_val);
        }
      }

    }

    if (HI[i] != HI[i])  {
      OMP_PRAGMA_CRITICAL
      {
        std::printf("HUGE HI! ::  %d %d %d %g\n",
                    i, idx_range.j, idx_range.k, HI[i]);
      }
    }
  }

}

}  // namespace step_gs_detail

/// Uses one linearly implicit Gauss-Seidel sweep of a backward-Euler time
/// integrator to advance the rate equations by one (sub-)cycle (dtit).
inline void step_rate_gauss_seidel(
  const double* dtit, IndexRange idx_range, gr_mask_type anydust,
  const double* rhoH, double* dedot_prev,
  double* HIdot_prev, const gr_mask_type* itmask,
  const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
  grackle_field_data* my_fields, SpeciesMultiView<const gr_float> sp_densities,
  SpeciesCollection species_tmpdens, const FullRxnRateBuf rxn_rate_buf,
  const PartMap<SpKind>& sp_kind_map
) {
  View<gr_float**> species_tmpdens_view(species_tmpdens.data[0],
                                        MAX_EVOLVED_SPECIES_FIELDS,
                                        my_fields->grid_dimension[0]);

  // perform the Gauss-Seidel sweep to compute the species densities at the
  // end of the current timestep. The results are saved in species_tmpdens
  chemistry::species_density_updates_gauss_seidel(
      species_tmpdens_view, idx_range, dtit, anydust, rhoH, itmask,
      itmask_metal, my_chemistry, sp_densities, rxn_rate_buf, sp_kind_map);

  // update the entries from my_fields with the values in species_tmpdens
  step_gs_detail::update_from_tmpdens(
      dtit, idx_range, dedot_prev, HIdot_prev, itmask, itmask_metal,
      my_chemistry, my_fields, species_tmpdens);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // STEP_RATE_GAUSS_SEIDEL_HPP
