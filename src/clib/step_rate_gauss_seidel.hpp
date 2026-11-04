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
  const gr_mask_type* itmask_metal, const chemistry_data* my_chemistry,
  SpeciesMultiView<gr_float> spdens,
  grackle::impl::SpeciesCollection species_tmpdens
) {

  // Get pointers to the various species fields
  // ------------------------------------------
  // -> TODO: stop doing this!
  gr_float* de = spdens.contig1d_ptr(SpLUT::e);
  gr_float* HI = spdens.contig1d_ptr(SpLUT::HI);
  gr_float* HII = spdens.contig1d_ptr(SpLUT::HII);
  gr_float* HeI = spdens.contig1d_ptr(SpLUT::HeI);
  gr_float* HeII = spdens.contig1d_ptr(SpLUT::HeII);
  gr_float* HeIII = spdens.contig1d_ptr(SpLUT::HeIII);
  gr_float* HM = spdens.contig1d_ptr(SpLUT::HM);
  gr_float* H2I = spdens.contig1d_ptr(SpLUT::H2I);
  gr_float* H2II = spdens.contig1d_ptr(SpLUT::H2II);
  gr_float* DI = spdens.contig1d_ptr(SpLUT::DI);
  gr_float* DII = spdens.contig1d_ptr(SpLUT::DII);
  gr_float* HDI = spdens.contig1d_ptr(SpLUT::HDI);
  gr_float* DM = spdens.contig1d_ptr(SpLUT::DM);
  gr_float* HDII = spdens.contig1d_ptr(SpLUT::HDII);
  gr_float* HeHII = spdens.contig1d_ptr(SpLUT::HeHII);
  gr_float* CI = spdens.contig1d_ptr(SpLUT::CI);
  gr_float* CII = spdens.contig1d_ptr(SpLUT::CII);
  gr_float* COI = spdens.contig1d_ptr(SpLUT::COI);
  gr_float* CO2I = spdens.contig1d_ptr(SpLUT::CO2I);
  gr_float* OI = spdens.contig1d_ptr(SpLUT::OI);
  gr_float* OHI = spdens.contig1d_ptr(SpLUT::OHI);
  gr_float* H2OI = spdens.contig1d_ptr(SpLUT::H2OI);
  gr_float* O2I = spdens.contig1d_ptr(SpLUT::O2I);
  gr_float* SiI = spdens.contig1d_ptr(SpLUT::SiI);
  gr_float* SiOI = spdens.contig1d_ptr(SpLUT::SiOI);
  gr_float* SiO2I = spdens.contig1d_ptr(SpLUT::SiO2I);
  gr_float* CHI = spdens.contig1d_ptr(SpLUT::CHI);
  gr_float* CH2I = spdens.contig1d_ptr(SpLUT::CH2I);
  gr_float* COII = spdens.contig1d_ptr(SpLUT::COII);
  gr_float* OII = spdens.contig1d_ptr(SpLUT::OII);
  gr_float* OHII = spdens.contig1d_ptr(SpLUT::OHII);
  gr_float* H2OII = spdens.contig1d_ptr(SpLUT::H2OII);
  gr_float* H3OII = spdens.contig1d_ptr(SpLUT::H3OII);
  gr_float* O2II = spdens.contig1d_ptr(SpLUT::O2II);
  gr_float* Mg = spdens.contig1d_ptr(SpLUT::Mg);
  gr_float* Al = spdens.contig1d_ptr(SpLUT::Al);
  gr_float* S = spdens.contig1d_ptr(SpLUT::S);
  gr_float* Fe = spdens.contig1d_ptr(SpLUT::Fe);
  gr_float* SiM = spdens.contig1d_ptr(SpLUT::SiM_dust);
  gr_float* FeM = spdens.contig1d_ptr(SpLUT::FeM_dust);
  gr_float* Mg2SiO4 = spdens.contig1d_ptr(SpLUT::Mg2SiO4_dust);
  gr_float* MgSiO3 = spdens.contig1d_ptr(SpLUT::MgSiO3_dust);
  gr_float* Fe3O4 = spdens.contig1d_ptr(SpLUT::Fe3O4_dust);
  gr_float* AC = spdens.contig1d_ptr(SpLUT::AC_dust);
  gr_float* SiO2D = spdens.contig1d_ptr(SpLUT::SiO2_dust);
  gr_float* MgO = spdens.contig1d_ptr(SpLUT::MgO_dust);
  gr_float* FeS = spdens.contig1d_ptr(SpLUT::FeS_dust);
  gr_float* Al2O3 = spdens.contig1d_ptr(SpLUT::Al2O3_dust);
  gr_float* reforg = spdens.contig1d_ptr(SpLUT::ref_org_dust);
  gr_float* volorg = spdens.contig1d_ptr(SpLUT::vol_org_dust);
  gr_float* H2Oice = spdens.contig1d_ptr(SpLUT::H2O_ice_dust);

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
  grackle_field_data* my_fields, SpeciesMultiView<gr_float> sp_densities,
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
      my_chemistry, sp_densities, species_tmpdens);
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // STEP_RATE_GAUSS_SEIDEL_HPP
