//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the @ref make_consistent function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "field_adaptor.hpp"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "inject_model/misc.hpp"
#include "LUT.hpp"

#include "make_consistent.hpp"

namespace GRIMPL_NAMESPACE_DECL {

void make_consistent(const IndexRange& idx_range, double dom,
                     const chemistry_data* my_chemistry,
                     const GrainMetalInjectPathways* inject_pathway_props,
                     grackle_field_data* my_fields,
                     SpeciesMultiView<gr_float> sp_densities) {
  FortranView<const gr_float***> d(
      const_cast<const gr_float*>(my_fields->density),
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  FortranView<const gr_float***> metal(
      const_cast<const gr_float*>(my_fields->metal_density),
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

  // TODO: when we rewrite this function, we should stop using stop
  // - loading unconditionally accessing pointers (e.g. we should only be
  //   retrieving a pointer if we know we're using it)
  // - replace all occurrences of SpLUT (in practice, we're **probably** going
  //   to convert this whole thing to a functor & prebuild lists of species ids)
  gr_float* de = sp_densities.contig1d_ptr(SpLUT::e);
  gr_float* HI = sp_densities.contig1d_ptr(SpLUT::HI);
  gr_float* HII = sp_densities.contig1d_ptr(SpLUT::HII);
  gr_float* HeI = sp_densities.contig1d_ptr(SpLUT::HeI);
  gr_float* HeII = sp_densities.contig1d_ptr(SpLUT::HeII);
  gr_float* HeIII = sp_densities.contig1d_ptr(SpLUT::HeIII);
  gr_float* HM = sp_densities.contig1d_ptr(SpLUT::HM);
  gr_float* H2I = sp_densities.contig1d_ptr(SpLUT::H2I);
  gr_float* H2II = sp_densities.contig1d_ptr(SpLUT::H2II);
  gr_float* DI = sp_densities.contig1d_ptr(SpLUT::DI);
  gr_float* DII = sp_densities.contig1d_ptr(SpLUT::DII);
  gr_float* HDI = sp_densities.contig1d_ptr(SpLUT::HDI);
  gr_float* DM = sp_densities.contig1d_ptr(SpLUT::DM);
  gr_float* HDII = sp_densities.contig1d_ptr(SpLUT::HDII);
  gr_float* HeHII = sp_densities.contig1d_ptr(SpLUT::HeHII);
  gr_float* CI = sp_densities.contig1d_ptr(SpLUT::CI);
  gr_float* CII = sp_densities.contig1d_ptr(SpLUT::CII);
  gr_float* COI = sp_densities.contig1d_ptr(SpLUT::COI);
  gr_float* CO2I = sp_densities.contig1d_ptr(SpLUT::CO2I);
  gr_float* OI = sp_densities.contig1d_ptr(SpLUT::OI);
  gr_float* OHI = sp_densities.contig1d_ptr(SpLUT::OHI);
  gr_float* H2OI = sp_densities.contig1d_ptr(SpLUT::H2OI);
  gr_float* O2I = sp_densities.contig1d_ptr(SpLUT::O2I);
  gr_float* SiI = sp_densities.contig1d_ptr(SpLUT::SiI);
  gr_float* SiOI = sp_densities.contig1d_ptr(SpLUT::SiOI);
  gr_float* SiO2I = sp_densities.contig1d_ptr(SpLUT::SiO2I);
  gr_float* CHI = sp_densities.contig1d_ptr(SpLUT::CHI);
  gr_float* CH2I = sp_densities.contig1d_ptr(SpLUT::CH2I);
  gr_float* COII = sp_densities.contig1d_ptr(SpLUT::COII);
  gr_float* OII = sp_densities.contig1d_ptr(SpLUT::OII);
  gr_float* OHII = sp_densities.contig1d_ptr(SpLUT::OHII);
  gr_float* H2OII = sp_densities.contig1d_ptr(SpLUT::H2OII);
  gr_float* H3OII = sp_densities.contig1d_ptr(SpLUT::H3OII);
  gr_float* O2II = sp_densities.contig1d_ptr(SpLUT::O2II);
  gr_float* Mg = sp_densities.contig1d_ptr(SpLUT::Mg);
  gr_float* Al = sp_densities.contig1d_ptr(SpLUT::Al);
  gr_float* S = sp_densities.contig1d_ptr(SpLUT::S);
  gr_float* Fe = sp_densities.contig1d_ptr(SpLUT::Fe);
  gr_float* SiM = sp_densities.contig1d_ptr(SpLUT::SiM_dust);
  gr_float* FeM = sp_densities.contig1d_ptr(SpLUT::FeM_dust);
  gr_float* Mg2SiO4 = sp_densities.contig1d_ptr(SpLUT::Mg2SiO4_dust);
  gr_float* MgSiO3 = sp_densities.contig1d_ptr(SpLUT::MgSiO3_dust);
  gr_float* Fe3O4 = sp_densities.contig1d_ptr(SpLUT::Fe3O4_dust);
  gr_float* AC = sp_densities.contig1d_ptr(SpLUT::AC_dust);
  gr_float* SiO2D = sp_densities.contig1d_ptr(SpLUT::SiO2_dust);
  gr_float* MgO = sp_densities.contig1d_ptr(SpLUT::MgO_dust);
  gr_float* FeS = sp_densities.contig1d_ptr(SpLUT::FeS_dust);
  gr_float* Al2O3 = sp_densities.contig1d_ptr(SpLUT::Al2O3_dust);
  gr_float* reforg = sp_densities.contig1d_ptr(SpLUT::ref_org_dust);
  gr_float* volorg = sp_densities.contig1d_ptr(SpLUT::vol_org_dust);
  gr_float* H2Oice = sp_densities.contig1d_ptr(SpLUT::H2O_ice_dust);

  // locals

  // flag for if Gen Chiaki's dust model is enabled with grain growth
  const bool chiaki_model_dust_evolution =
      my_chemistry->dust_chemistry == 2 && my_chemistry->grain_growth == 1;
  double totalD;
  std::vector<double> totalH(my_fields->grid_dimension[0]);
  std::vector<double> totalHe(my_fields->grid_dimension[0]);
  std::vector<double> metalfree(my_fields->grid_dimension[0]);
  gr_float correctH, correctHe, correctD;
  double totalC, totalO, totalMg, totalAl, totalSi, totalS, totalFe;
  double totalCg, totalOg, totalMgg, totalSig, totalFeg;
  double totalCd, totalOd, totalMgd, totalSid, totalFed;
  gr_float correctC, correctO, correctMg, correctAl, correctSi, correctS,
      correctFe;
  gr_float correctCg, correctOg, correctMgg, correctSig, correctFeg;
  gr_float correctCd, correctOd, correctMgd, correctSid, correctFed;

  // when relevant, each view within SN_metal_arr wraps a field specifying the
  // total metal density that corresponds to an injection pathway
  FortranView<const gr_float***>
      SN_metal_arr[inj_model_input::N_Injection_Pathways];

  int n_pathways = 0;
  // construct view of each specified injection pathway metal density field
  if (my_chemistry->metal_chemistry > 0) {
    // note: when (my_chemistry->multi_metals == 0) a view within
    //       SN_metal_arr will wrap my_fields->metal_density. In other words,
    //       that view will be an alias of the `metal` view. This is ok because
    //       my_fields->metal_density is **NOT** mutated by this function.
    n_pathways = inject_pathway_props->n_pathways;
    const gr_float* const* inject_pathway_metal_densities =
        get_inject_pathway_metal_density(my_chemistry, my_fields);

    for (int iSN = 0; iSN < n_pathways; iSN++) {
      SN_metal_arr[iSN] = FortranView<const gr_float***>(
          inject_pathway_metal_densities[iSN], my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    }
  }

  std::vector<double> Ct(my_fields->grid_dimension[0]);
  std::vector<double> Ot(my_fields->grid_dimension[0]);
  std::vector<double> Mgt(my_fields->grid_dimension[0]);
  std::vector<double> Alt(my_fields->grid_dimension[0]);
  std::vector<double> Sit(my_fields->grid_dimension[0]);
  std::vector<double> St(my_fields->grid_dimension[0]);
  std::vector<double> Fet(my_fields->grid_dimension[0]);
  std::vector<double> Cg(my_fields->grid_dimension[0]);
  std::vector<double> Og(my_fields->grid_dimension[0]);
  std::vector<double> Mgg(my_fields->grid_dimension[0]);
  std::vector<double> Alg(my_fields->grid_dimension[0]);
  std::vector<double> Sig(my_fields->grid_dimension[0]);
  std::vector<double> Sg(my_fields->grid_dimension[0]);
  std::vector<double> Feg(my_fields->grid_dimension[0]);
  std::vector<double> Cd(my_fields->grid_dimension[0]);
  std::vector<double> Od(my_fields->grid_dimension[0]);
  std::vector<double> Mgd(my_fields->grid_dimension[0]);
  std::vector<double> Ald(my_fields->grid_dimension[0]);
  std::vector<double> Sid(my_fields->grid_dimension[0]);
  std::vector<double> Sd(my_fields->grid_dimension[0]);
  std::vector<double> Fed(my_fields->grid_dimension[0]);

  // Loop over current index range

  const int k = idx_range.k;
  const int j = idx_range.j;

  // Compute total densities of H and He (ensure non-negativity)

  if ((my_chemistry->metal_cooling) == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      metalfree[i] = d(i, j, k) - metal(i, j, k);
    }
  } else {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      metalfree[i] = d(i, j, k);
    }
  }

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    HI[i] = std::fabs(HI[i]);
    HII[i] = std::fabs(HII[i]);
    HeI[i] = std::fabs(HeI[i]);
    HeII[i] = std::fabs(HeII[i]);
    HeIII[i] = std::fabs(HeIII[i]);
    totalH[i] = HI[i] + HII[i];
    totalHe[i] = HeI[i] + HeII[i] + HeIII[i];
  }

  // include molecular hydrogen

  if (my_chemistry->primordial_chemistry > 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      HM[i] = std::fabs(HM[i]);
      H2II[i] = std::fabs(H2II[i]);
      H2I[i] = std::fabs(H2I[i]);
      totalH[i] = totalH[i] + HM[i] + H2I[i] + H2II[i];
    }
  }

  if (my_chemistry->primordial_chemistry > 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      HDI[i] = std::fabs(HDI[i]);
      totalH[i] = totalH[i] + 1. / 3. * HDI[i];
    }
  }
  // ! GC202005

  if (my_chemistry->primordial_chemistry > 3) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      HDII[i] = std::fabs(HDII[i]);
      HeHII[i] = std::fabs(HeHII[i]);
      totalH[i] = totalH[i] + 1. / 3. * HDII[i] + 1. / 5. * HeHII[i];
      totalHe[i] = totalHe[i] + 4. / 5. * HeHII[i];
    }
  }

  // In the future, we may elect to only consider metal-rich cells by
  // using an ``itmask_metal`` buffer:
  // -> the fortran version of this function actually had some commented out
  //    logic for handling this precise situation. The modern equivalent of
  //    that logic would use the ``mask::fill_itmask_metal`` function to
  //    fill ``itmask_metal``
  // -> however, a better solution may be to make this function only
  //    consider a single index-range of values at a time (rather than
  //    entire fields) and to pass in a prefilled ``itmask_metal`` buffer
  //    as an argument for use within this function
  // -> this may be better if we decide to allow the metal density field to
  //    change in value from dust formation or destruction. If we just
  //    computed ``itmask_metal`` on the fly, then there would be cases
  //    where the computed mask may be different from the mask that was
  //    used to actually decide whether to perform metal/dust chemistry

  if (my_chemistry->metal_chemistry > 0) {
    // compute the expected mass density for each metal nuclide using the
    // yields from every injection pathway
    const grackle::impl::yields::MetalTables& total_metal_yields =
        inject_pathway_props->total_metal_nuclide_yields;
    const grackle::impl::yields::MetalTables& onlygas_metal_yields =
        inject_pathway_props->gas_metal_nuclide_yields;

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      Ct[i] = 0.;
      Cg[i] = 0.;
      Ot[i] = 0.;
      Og[i] = 0.;
      Mgt[i] = 0.;
      Mgg[i] = 0.;
      Alt[i] = 0.;
      Alg[i] = 0.;
      Sit[i] = 0.;
      Sig[i] = 0.;
      St[i] = 0.;
      Sg[i] = 0.;
      Fet[i] = 0.;
      Feg[i] = 0.;
      for (int iSN = 0; iSN < n_pathways; iSN++) {
        gr_float cur_val = SN_metal_arr[iSN](i, j, k);

        Ct[i] = Ct[i] + total_metal_yields.C[iSN] * cur_val;
        Ot[i] = Ot[i] + total_metal_yields.O[iSN] * cur_val;
        Mgt[i] = Mgt[i] + total_metal_yields.Mg[iSN] * cur_val;
        Alt[i] = Alt[i] + total_metal_yields.Al[iSN] * cur_val;
        Sit[i] = Sit[i] + total_metal_yields.Si[iSN] * cur_val;
        St[i] = St[i] + total_metal_yields.S[iSN] * cur_val;
        Fet[i] = Fet[i] + total_metal_yields.Fe[iSN] * cur_val;

        Cg[i] = Cg[i] + onlygas_metal_yields.C[iSN] * cur_val;
        Og[i] = Og[i] + onlygas_metal_yields.O[iSN] * cur_val;
        Mgg[i] = Mgg[i] + onlygas_metal_yields.Mg[iSN] * cur_val;
        Alg[i] = Alg[i] + onlygas_metal_yields.Al[iSN] * cur_val;
        Sig[i] = Sig[i] + onlygas_metal_yields.Si[iSN] * cur_val;
        Sg[i] = Sg[i] + onlygas_metal_yields.S[iSN] * cur_val;
        Feg[i] = Feg[i] + onlygas_metal_yields.Fe[iSN] * cur_val;
      }
    }

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      Cd[i] = Ct[i] - Cg[i];
      Od[i] = Ot[i] - Og[i];
      Mgd[i] = Mgt[i] - Mgg[i];
      Ald[i] = Alt[i] - Alg[i];
      Sid[i] = Sit[i] - Sig[i];
      Sd[i] = St[i] - Sg[i];
      Fed[i] = Fet[i] - Feg[i];
    }

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      // if (itmask_metal(i)) then
      OHI[i] = std::fabs(OHI[i]);
      H2OI[i] = std::fabs(H2OI[i]);
      CHI[i] = std::fabs(CHI[i]);
      CH2I[i] = std::fabs(CH2I[i]);
      OHII[i] = std::fabs(OHII[i]);
      H2OII[i] = std::fabs(H2OII[i]);
      H3OII[i] = std::fabs(H3OII[i]);
      totalH[i] = totalH[i] + OHI[i] / 17. + H2OI[i] / 18. * 2. + CHI[i] / 13. +
                  CH2I[i] / 14. * 2. + OHII[i] / 17. + H2OII[i] / 18. * 2. +
                  H3OII[i] / 19. * 3.;
      // endif
    }
  }

  // Correct densities by keeping fractions the same

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    correctH = (gr_float)(my_chemistry->HydrogenFractionByMass * metalfree[i] /
                          totalH[i]);
    HI[i] = HI[i] * correctH;
    HII[i] = HII[i] * correctH;

    correctHe = (gr_float)((1. - my_chemistry->HydrogenFractionByMass) *
                           metalfree[i] / totalHe[i]);
    HeI[i] = HeI[i] * correctHe;
    HeII[i] = HeII[i] * correctHe;
    HeIII[i] = HeIII[i] * correctHe;

    // Correct molecular hydrogen-related fractions

    if (my_chemistry->primordial_chemistry > 1) {
      HM[i] = HM[i] * correctH;
      H2II[i] = H2II[i] * correctH;
      H2I[i] = H2I[i] * correctH;
    }
    if (my_chemistry->primordial_chemistry > 3) {
      // !          HDII (i,j,k) = HDII (i,j,k)*correctH
      HeHII[i] = HeHII[i] * correctHe;
    }
  }

  // Do the same thing for deuterium (ignore HD) Assumes dtoh is small

  if (my_chemistry->primordial_chemistry > 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      DI[i] = std::fabs(DI[i]);
      DII[i] = std::fabs(DII[i]);
      HDI[i] = std::fabs(HDI[i]);
      totalD = DI[i] + DII[i] + 2. / 3. * HDI[i];
      if (my_chemistry->primordial_chemistry > 3) {
        DM[i] = std::fabs(DM[i]);
        HDII[i] = std::fabs(HDII[i]);
        totalD = totalD + DM[i] + 2. / 3. * HDII[i];
      }
      correctD = (gr_float)(my_chemistry->HydrogenFractionByMass *
                            my_chemistry->DeuteriumToHydrogenRatio *
                            metalfree[i] / totalD);
      DI[i] = DI[i] * correctD;
      DII[i] = DII[i] * correctD;
      HDI[i] = HDI[i] * correctD;
      if (my_chemistry->primordial_chemistry > 3) {
        DM[i] = DM[i] * correctD;
        HDII[i] = HDII[i] * correctD;
      }
    }
  }

  // Do the same thing for metal species

  if (my_chemistry->metal_chemistry == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      CI[i] = std::fabs(CI[i]);
      CII[i] = std::fabs(CII[i]);
      COI[i] = std::fabs(COI[i]);
      CO2I[i] = std::fabs(CO2I[i]);
      OI[i] = std::fabs(OI[i]);
      OHI[i] = std::fabs(OHI[i]);
      H2OI[i] = std::fabs(H2OI[i]);
      O2I[i] = std::fabs(O2I[i]);
      SiI[i] = std::fabs(SiI[i]);
      SiOI[i] = std::fabs(SiOI[i]);
      SiO2I[i] = std::fabs(SiO2I[i]);
      CHI[i] = std::fabs(CHI[i]);
      CH2I[i] = std::fabs(CH2I[i]);
      COII[i] = std::fabs(COII[i]);
      OII[i] = std::fabs(OII[i]);
      OHII[i] = std::fabs(OHII[i]);
      H2OII[i] = std::fabs(H2OII[i]);
      H3OII[i] = std::fabs(H3OII[i]);
      O2II[i] = std::fabs(O2II[i]);
    }
  }

  if (chiaki_model_dust_evolution) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (my_chemistry->dust_species > 0) {
        Mg[i] = std::fabs(Mg[i]);

        MgSiO3[i] = std::fabs(MgSiO3[i]);
        AC[i] = std::fabs(AC[i]);
      }
      if (my_chemistry->dust_species > 1) {
        Al[i] = std::fabs(Al[i]);
        S[i] = std::fabs(S[i]);
        Fe[i] = std::fabs(Fe[i]);

        SiM[i] = std::fabs(SiM[i]);
        FeM[i] = std::fabs(FeM[i]);
        Mg2SiO4[i] = std::fabs(Mg2SiO4[i]);
        Fe3O4[i] = std::fabs(Fe3O4[i]);
        SiO2D[i] = std::fabs(SiO2D[i]);
        MgO[i] = std::fabs(MgO[i]);
        FeS[i] = std::fabs(FeS[i]);
        Al2O3[i] = std::fabs(Al2O3[i]);
      }
      if (my_chemistry->dust_species > 2) {
        reforg[i] = std::fabs(reforg[i]);
        volorg[i] = std::fabs(volorg[i]);
        H2Oice[i] = std::fabs(H2Oice[i]);
      }
    }
  }

  if (my_chemistry->metal_chemistry == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      // if (itmask_metal(i)) then
      //- if (d(i,j,k)*dom .lt. 1.e-2_DKIND) then
      // !       if (d(i,j,k)*dom .lt.
      // !   &    min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
      // !   &       ,1.e6_DKIND)) then
      if (((my_chemistry->metal_cooling == 0) && (d(i, j, k) * dom < 1.e8)) ||
          ((my_chemistry->metal_cooling == 1) &&
           (((metal(i, j, k) <= 1.e-9 * d(i, j, k)) &&
             (d(i, j, k) * dom < 1.e8)) ||
            ((metal(i, j, k) > 1.e-9 * d(i, j, k)) &&
             (d(i, j, k) * dom < 1.e6))))) {
        totalOg = 16. / 28. * COI[i] + 32. / 44. * CO2I[i] + OI[i] +
                  16. / 17. * OHI[i] + 16. / 18. * H2OI[i] + O2I[i] +
                  16. / 44. * SiOI[i] + 32. / 60. * SiO2I[i] +
                  16. / 28. * COII[i] + OII[i] + 16. / 17. * OHII[i] +
                  16. / 18. * H2OII[i] + 16. / 19. * H3OII[i] + O2II[i];
        correctOg = (gr_float)(Og[i] / totalOg);
        COI[i] = COI[i] * correctOg;
        CO2I[i] = CO2I[i] * correctOg;
        OI[i] = OI[i] * correctOg;
        OHI[i] = OHI[i] * correctOg;
        H2OI[i] = H2OI[i] * correctOg;
        O2I[i] = O2I[i] * correctOg;
        SiOI[i] = SiOI[i] * correctOg;
        SiO2I[i] = SiO2I[i] * correctOg;
        COII[i] = COII[i] * correctOg;
        OII[i] = OII[i] * correctOg;
        OHII[i] = OHII[i] * correctOg;
        H2OII[i] = H2OII[i] * correctOg;
        H3OII[i] = H3OII[i] * correctOg;
        O2II[i] = O2II[i] * correctOg;
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalOd = 48. / 100. * MgSiO3[i];
          }
          if (my_chemistry->dust_species > 1) {
            totalOd = totalOd + 64. / 140. * Mg2SiO4[i] +
                      64. / 232. * Fe3O4[i] + 32. / 60. * SiO2D[i] +
                      16. / 40. * MgO[i] + 48. / 102. * Al2O3[i];
          }
          if (my_chemistry->dust_species > 2) {
            totalOd = totalOd + 8. / 22.68 * reforg[i] + 16. / 32. * volorg[i] +
                      16. / 18. * H2Oice[i];
          }
          correctOd = (gr_float)(Od[i] / totalOd);
          if (my_chemistry->dust_species > 0) {
            MgSiO3[i] = MgSiO3[i] * correctOd;
          }
          if (my_chemistry->dust_species > 1) {
            Mg2SiO4[i] = Mg2SiO4[i] * correctOd;
            Fe3O4[i] = Fe3O4[i] * correctOd;
            SiO2D[i] = SiO2D[i] * correctOd;
            MgO[i] = MgO[i] * correctOd;
            Al2O3[i] = Al2O3[i] * correctOd;
          }
          if (my_chemistry->dust_species > 2) {
            reforg[i] = reforg[i] * correctOd;
            volorg[i] = volorg[i] * correctOd;
            H2Oice[i] = H2Oice[i] * correctOd;
          }
        }

        totalCg = CI[i] + CII[i] + 12. / 28. * COI[i] + 12. / 44. * CO2I[i] +
                  12. / 13. * CHI[i] + 12. / 14. * CH2I[i] +
                  12. / 28. * COII[i];
        correctCg = (gr_float)(Cg[i] / totalCg);
        CI[i] = CI[i] * correctCg;
        CII[i] = CII[i] * correctCg;
        COI[i] = COI[i] * correctCg;
        CO2I[i] = CO2I[i] * correctCg;
        CHI[i] = CHI[i] * correctCg;
        CH2I[i] = CH2I[i] * correctCg;
        COII[i] = COII[i] * correctCg;
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalCd = AC[i];
          }
          if (my_chemistry->dust_species > 2) {
            totalCd = totalCd + 12. / 22.68 * reforg[i] + 12. / 32. * volorg[i];
          }
          correctCd = (gr_float)(Cd[i] / totalCd);
          if (my_chemistry->dust_species > 0) {
            AC[i] = AC[i] * correctCd;
          }
          if (my_chemistry->dust_species > 2) {
            reforg[i] = reforg[i] * correctCd;
            volorg[i] = volorg[i] * correctCd;
          }
        }

        totalSig = SiI[i] + 28. / 44. * SiOI[i] + 28. / 60. * SiO2I[i];
        correctSig = (gr_float)(Sig[i] / totalSig);
        SiI[i] = SiI[i] * correctSig;
        SiOI[i] = SiOI[i] * correctSig;
        SiO2I[i] = SiO2I[i] * correctSig;
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalSid = 28. / 100. * MgSiO3[i];
          }
          if (my_chemistry->dust_species > 1) {
            totalSid = totalSid + SiM[i] + 28. / 140. * Mg2SiO4[i] +
                       28. / 60. * SiO2D[i];
          }
          correctSid = (gr_float)(Sid[i] / totalSid);
          if (my_chemistry->dust_species > 0) {
            MgSiO3[i] = MgSiO3[i] * correctSid;
          }
          if (my_chemistry->dust_species > 1) {
            SiM[i] = SiM[i] * correctSid;
            Mg2SiO4[i] = Mg2SiO4[i] * correctSid;
            SiO2D[i] = SiO2D[i] * correctSid;
          }
        }

        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 1) {
            totalFeg = Fe[i];
            correctFeg = (gr_float)(Feg[i] / totalFeg);
            Fe[i] = Fe[i] * correctFeg;

            totalFed = FeM[i] + 168. / 232. * Fe3O4[i] + 56. / 88. * FeS[i];
            correctFed = (gr_float)(Fed[i] / totalFed);
            FeM[i] = FeM[i] * correctFed;
            Fe3O4[i] = Fe3O4[i] * correctFed;
            FeS[i] = FeS[i] * correctFed;
          }

          if (my_chemistry->dust_species > 0) {
            totalMgg = Mg[i];
            correctMgg = (gr_float)(Mgg[i] / totalMgg);
            Mg[i] = Mg[i] * correctMgg;
            totalMgd = 24. / 100. * MgSiO3[i];
            if (my_chemistry->dust_species > 1) {
              totalMgd =
                  totalMgd + 48. / 140. * Mg2SiO4[i] + 24. / 40. * MgO[i];
            }
            correctMgd = (gr_float)(Mgd[i] / totalMgd);
            MgSiO3[i] = MgSiO3[i] * correctMgd;
            if (my_chemistry->dust_species > 1) {
              Mg2SiO4[i] = Mg2SiO4[i] * correctMgd;
              MgO[i] = MgO[i] * correctMgd;
            }
          }

          if (my_chemistry->dust_species > 1) {
            S[i] = Sg[i];
            FeS[i] = 88. / 32. * Sd[i];

            Al[i] = Alg[i];
            Al2O3[i] = 102. / 54. * Ald[i];
          }
        }

      } else {
        totalO = 16. / 28. * COI[i] + 32. / 44. * CO2I[i] + OI[i] +
                 16. / 17. * OHI[i] + 16. / 18. * H2OI[i] + O2I[i] +
                 16. / 44. * SiOI[i] + 32. / 60. * SiO2I[i] +
                 16. / 28. * COII[i] + OII[i] + 16. / 17. * OHII[i] +
                 16. / 18. * H2OII[i] + 16. / 19. * H3OII[i] + O2II[i];
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalO = totalO + 48. / 100. * MgSiO3[i];
          }
          if (my_chemistry->dust_species > 1) {
            totalO = totalO + 64. / 140. * Mg2SiO4[i] + 64. / 232. * Fe3O4[i] +
                     32. / 60. * SiO2D[i] + 16. / 40. * MgO[i] +
                     48. / 102. * Al2O3[i];
          }
          if (my_chemistry->dust_species > 2) {
            totalO = totalO + 8. / 22.68 * reforg[i] + 16. / 32. * volorg[i] +
                     16. / 18. * H2Oice[i];
          }
        }
        if (!chiaki_model_dust_evolution) {
          correctO = (gr_float)(Og[i] / totalO);
          COI[i] = COI[i] * correctO;
          CO2I[i] = CO2I[i] * correctO;
          OI[i] = OI[i] * correctO;
          OHI[i] = OHI[i] * correctO;
          H2OI[i] = H2OI[i] * correctO;
          O2I[i] = O2I[i] * correctO;
          SiOI[i] = SiOI[i] * correctO;
          SiO2I[i] = SiO2I[i] * correctO;
          COII[i] = COII[i] * correctO;
          OII[i] = OII[i] * correctO;
          OHII[i] = OHII[i] * correctO;
          H2OII[i] = H2OII[i] * correctO;
          H3OII[i] = H3OII[i] * correctO;
          O2II[i] = O2II[i] * correctO;
        } else {
          correctO = (gr_float)(Ot[i] / totalO);
          COI[i] = COI[i] * correctO;
          CO2I[i] = CO2I[i] * correctO;
          OI[i] = OI[i] * correctO;
          OHI[i] = OHI[i] * correctO;
          H2OI[i] = H2OI[i] * correctO;
          O2I[i] = O2I[i] * correctO;
          SiOI[i] = SiOI[i] * correctO;
          SiO2I[i] = SiO2I[i] * correctO;
          COII[i] = COII[i] * correctO;
          OII[i] = OII[i] * correctO;
          OHII[i] = OHII[i] * correctO;
          H2OII[i] = H2OII[i] * correctO;
          H3OII[i] = H3OII[i] * correctO;
          O2II[i] = O2II[i] * correctO;
          if (my_chemistry->dust_species > 0) {
            MgSiO3[i] = MgSiO3[i] * correctO;
          }
          if (my_chemistry->dust_species > 1) {
            Mg2SiO4[i] = Mg2SiO4[i] * correctO;
            Fe3O4[i] = Fe3O4[i] * correctO;
            SiO2D[i] = SiO2D[i] * correctO;
            MgO[i] = MgO[i] * correctO;
            Al2O3[i] = Al2O3[i] * correctO;
          }
          if (my_chemistry->dust_species > 2) {
            reforg[i] = reforg[i] * correctO;
            volorg[i] = volorg[i] * correctO;
            H2Oice[i] = H2Oice[i] * correctO;
          }
        }

        totalC = CI[i] + CII[i] + 12. / 28. * COI[i] + 12. / 44. * CO2I[i] +
                 12. / 13. * CHI[i] + 12. / 14. * CH2I[i] + 12. / 28. * COII[i];
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalC = totalC + AC[i];
          }
          if (my_chemistry->dust_species > 2) {
            totalC = totalC + 12. / 22.68 * reforg[i] + 12. / 32. * volorg[i];
          }
        }
        if (!chiaki_model_dust_evolution) {
          correctC = (gr_float)(Cg[i] / totalC);
          CI[i] = CI[i] * correctC;
          CII[i] = CII[i] * correctC;
          COI[i] = COI[i] * correctC;
          CO2I[i] = CO2I[i] * correctC;
          CHI[i] = CHI[i] * correctC;
          CH2I[i] = CH2I[i] * correctC;
          COII[i] = COII[i] * correctC;
        } else {
          correctC = (gr_float)(Ct[i] / totalC);
          CI[i] = CI[i] * correctC;
          CII[i] = CII[i] * correctC;
          COI[i] = COI[i] * correctC;
          CO2I[i] = CO2I[i] * correctC;
          CHI[i] = CHI[i] * correctC;
          CH2I[i] = CH2I[i] * correctC;
          COII[i] = COII[i] * correctC;
          if (my_chemistry->dust_species > 0) {
            AC[i] = AC[i] * correctC;
          }
          if (my_chemistry->dust_species > 2) {
            reforg[i] = reforg[i] * correctC;
            volorg[i] = volorg[i] * correctC;
          }
        }

        totalSi = SiI[i] + 28. / 44. * SiOI[i] + 28. / 60. * SiO2I[i];
        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 0) {
            totalSi = totalSi + 28. / 100. * MgSiO3[i];
          }
          if (my_chemistry->dust_species > 1) {
            totalSi = totalSi + SiM[i] + 28. / 140. * Mg2SiO4[i] +
                      28. / 60. * SiO2D[i];
          }
        }
        if (!chiaki_model_dust_evolution) {
          correctSi = (gr_float)(Sig[i] / totalSi);
          SiI[i] = SiI[i] * correctSi;
          SiOI[i] = SiOI[i] * correctSi;
          SiO2I[i] = SiO2I[i] * correctSi;
        } else {
          correctSi = (gr_float)(Sit[i] / totalSi);
          SiI[i] = SiI[i] * correctSi;
          SiOI[i] = SiOI[i] * correctSi;
          SiO2I[i] = SiO2I[i] * correctSi;
          if (my_chemistry->dust_species > 0) {
            MgSiO3[i] = MgSiO3[i] * correctSi;
          }
          if (my_chemistry->dust_species > 1) {
            SiM[i] = SiM[i] * correctSi;
            Mg2SiO4[i] = Mg2SiO4[i] * correctSi;
            SiO2D[i] = SiO2D[i] * correctSi;
          }
        }

        if (chiaki_model_dust_evolution) {
          if (my_chemistry->dust_species > 1) {
            totalFe =
                Fe[i] + FeM[i] + 168. / 232. * Fe3O4[i] + 56. / 88. * FeS[i];
            correctFe = (gr_float)(Fet[i] / totalFe);
            Fe[i] = Fe[i] * correctFe;
            FeM[i] = FeM[i] * correctFe;
            Fe3O4[i] = Fe3O4[i] * correctFe;
            FeS[i] = FeS[i] * correctFe;
          }

          if (my_chemistry->dust_species > 0) {
            totalMg = Mg[i] + 24. / 100. * MgSiO3[i];
            if (my_chemistry->dust_species > 1) {
              totalMg = totalMg + 48. / 140. * Mg2SiO4[i] + 24. / 40. * MgO[i];
            }
            correctMg = (gr_float)(Mgt[i] / totalMg);
            Mg[i] = Mg[i] * correctMg;
            MgSiO3[i] = MgSiO3[i] * correctMg;
            if (my_chemistry->dust_species > 1) {
              Mg2SiO4[i] = Mg2SiO4[i] * correctMg;
              MgO[i] = MgO[i] * correctMg;
            }
          }

          if (my_chemistry->dust_species > 1) {
            totalS = S[i] + 32. / 88. * FeS[i];
            correctS = (gr_float)(St[i] / totalS);
            S[i] = S[i] * correctS;
            FeS[i] = FeS[i] * correctS;

            totalAl = Al[i] + 54. / 102. * Al2O3[i];
            correctAl = (gr_float)(Alt[i] / totalAl);
            Al[i] = Al[i] * correctAl;
            Al2O3[i] = Al2O3[i] * correctAl;
          }
        }
      }
    }
  }  // my_chemistry->metal_chemistry == 1

  // Set the electron density

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    de[i] = HII[i] + HeII[i] / (gr_float)(4.) + HeIII[i] / (gr_float)(2.);
    if (my_chemistry->primordial_chemistry > 1) {
      de[i] = de[i] - HM[i] + H2II[i] / (gr_float)(2.);
    }
    if (my_chemistry->primordial_chemistry > 3) {
      de[i] = de[i] - DM[i] / (gr_float)(2.) + HDII[i] / (gr_float)(3.) +
              HeHII[i] / (gr_float)(5.);
    }
    if (my_chemistry->metal_chemistry == 1) {
      // if (itmask_metal(i)) then
      de[i] = de[i] + CII[i] / (gr_float)(12.) + COII[i] / (gr_float)(28.) +
              OII[i] / (gr_float)(16.) + OHII[i] / (gr_float)(17.) +
              H2OII[i] / (gr_float)(18.) + H3OII[i] / (gr_float)(19.) +
              O2II[i] / (gr_float)(32.);
      // endif
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL
