//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the `make_consistent_g` function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// make_consistent_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "make_consistent.hpp"

namespace grackle::impl {

void make_consistent(int imetal, double dom, chemistry_data* my_chemistry,
                     chemistry_data_storage* my_rates,
                     grackle_field_data* my_fields) {
  // Arguments

  grackle::impl::View<gr_float***> de(
      my_fields->e_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(
      my_fields->HeI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(
      my_fields->HeII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(
      my_fields->HeIII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(
      my_fields->HM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(
      my_fields->H2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DI(
      my_fields->DI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(
      my_fields->DII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(
      my_fields->HDI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DM(
      my_fields->DM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(
      my_fields->HDII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(
      my_fields->HeHII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(
      my_fields->CI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(
      my_fields->CII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(
      my_fields->CO_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO2(
      my_fields->CO2_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(
      my_fields->OI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(
      my_fields->OH_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(
      my_fields->H2O_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2(
      my_fields->O2_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiI(
      my_fields->SiI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiOI(
      my_fields->SiOI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2I(
      my_fields->SiO2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH(
      my_fields->CH_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH2(
      my_fields->CH2_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> COII(
      my_fields->COII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OII(
      my_fields->OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OHII(
      my_fields->OHII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2OII(
      my_fields->H2OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H3OII(
      my_fields->H3OII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2II(
      my_fields->O2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg(
      my_fields->Mg_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al(
      my_fields->Al_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> S(
      my_fields->S_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe(
      my_fields->Fe_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(
      my_fields->SiM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(
      my_fields->FeM_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(
      my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(
      my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(
      my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(
      my_fields->AC_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(
      my_fields->SiO2_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(
      my_fields->MgO_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(
      my_fields->FeS_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(
      my_fields->Al2O3_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(
      my_fields->ref_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(
      my_fields->vol_org_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(
      my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_loc(
      my_fields->local_ISM_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C13(
      my_fields->ccsn13_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C20(
      my_fields->ccsn20_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C25(
      my_fields->ccsn25_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C30(
      my_fields->ccsn30_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F13(
      my_fields->fsn13_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F15(
      my_fields->fsn15_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F50(
      my_fields->fsn50_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F80(
      my_fields->fsn80_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P170(
      my_fields->pisn170_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P200(
      my_fields->pisn200_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_Y19(
      my_fields->y19_metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // locals

  int i, j, k;
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
  int iSN, nSN, iSN0;
  std::vector<gr_float> SN_metal_data_(my_fields->grid_dimension[0] *
                                       my_rates->SN0_N);
  grackle::impl::View<gr_float**> SN_metal(
      SN_metal_data_.data(), my_fields->grid_dimension[0], my_rates->SN0_N);
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

  // Loop over all zones

  for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
    for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
      // Compute total densities of H and He
      //     (ensure non-negativity)

      if ((imetal) == 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          metalfree[i] = d(i, j, k) - metal(i, j, k);
        }
      } else {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          metalfree[i] = d(i, j, k);
        }
      }

      for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
        HI(i, j, k) = std::fabs(HI(i, j, k));
        HII(i, j, k) = std::fabs(HII(i, j, k));
        HeI(i, j, k) = std::fabs(HeI(i, j, k));
        HeII(i, j, k) = std::fabs(HeII(i, j, k));
        HeIII(i, j, k) = std::fabs(HeIII(i, j, k));
        totalH[i] = HI(i, j, k) + HII(i, j, k);
        totalHe[i] = HeI(i, j, k) + HeII(i, j, k) + HeIII(i, j, k);
      }

      // include molecular hydrogen

      if (my_chemistry->primordial_chemistry > 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          HM(i, j, k) = std::fabs(HM(i, j, k));
          H2II(i, j, k) = std::fabs(H2II(i, j, k));
          H2I(i, j, k) = std::fabs(H2I(i, j, k));
          totalH[i] = totalH[i] + HM(i, j, k) + H2I(i, j, k) + H2II(i, j, k);
        }
      }

      if (my_chemistry->primordial_chemistry > 2) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          HDI(i, j, k) = std::fabs(HDI(i, j, k));
          totalH[i] = totalH[i] + 1. / 3. * HDI(i, j, k);
        }
      }
      // ! GC202005

      if (my_chemistry->primordial_chemistry > 3) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          HDII(i, j, k) = std::fabs(HDII(i, j, k));
          HeHII(i, j, k) = std::fabs(HeHII(i, j, k));
          totalH[i] =
              totalH[i] + 1. / 3. * HDII(i, j, k) + 1. / 5. * HeHII(i, j, k);
          totalHe[i] = totalHe[i] + 4. / 5. * HeHII(i, j, k);
        }
      }

      // Iteration mask for metal-rich cells

      // do i = is+1, ie + 1
      //    itmask_metal(i) = .false.
      // enddo
      // if (imetal .eq. 1) then
      //     do i = is+1, ie + 1
      //        if (metal(i,j,k) .gt. 1.e-9_DKIND * d(i,j,k)) then
      //           itmask_metal(i) = .true.
      //        endif
      //     enddo
      // endif

      if (my_chemistry->metal_chemistry > 0) {
        if (my_chemistry->multi_metals == 0) {
          iSN0 = my_chemistry->metal_abundances;
          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
            Ct[i] = my_rates->SN0_XC[iSN0] * metal(i, j, k);
            Ot[i] = my_rates->SN0_XO[iSN0] * metal(i, j, k);
            Mgt[i] = my_rates->SN0_XMg[iSN0] * metal(i, j, k);
            Alt[i] = my_rates->SN0_XAl[iSN0] * metal(i, j, k);
            Sit[i] = my_rates->SN0_XSi[iSN0] * metal(i, j, k);
            St[i] = my_rates->SN0_XS[iSN0] * metal(i, j, k);
            Fet[i] = my_rates->SN0_XFe[iSN0] * metal(i, j, k);

            Cg[i] = my_rates->SN0_fC[iSN0] * metal(i, j, k);
            Og[i] = my_rates->SN0_fO[iSN0] * metal(i, j, k);
            Mgg[i] = my_rates->SN0_fMg[iSN0] * metal(i, j, k);
            Alg[i] = my_rates->SN0_fAl[iSN0] * metal(i, j, k);
            Sig[i] = my_rates->SN0_fSi[iSN0] * metal(i, j, k);
            Sg[i] = my_rates->SN0_fS[iSN0] * metal(i, j, k);
            Feg[i] = my_rates->SN0_fFe[iSN0] * metal(i, j, k);
          }

        } else {
          //        do i = is+1, ie+1
          //           totalZ = metal_loc(i,j,k)
          // &           + metal_C13(i,j,k) + metal_C20(i,j,k)
          // &           + metal_C25(i,j,k) + metal_C30(i,j,k)
          // &           + metal_F13(i,j,k) + metal_F15(i,j,k)
          // &           + metal_F50(i,j,k) + metal_F80(i,j,k)
          // &           + metal_P170(i,j,k)+ metal_P200(i,j,k)
          // &           + metal_Y19(i,j,k)
          //           correctZ = metal(i,j,k) / totalZ
          //           metal_loc(i,j,k) = metal_loc(i,j,k) * correctZ
          //           metal_C13(i,j,k) = metal_C13(i,j,k) * correctZ
          //           metal_C20(i,j,k) = metal_C20(i,j,k) * correctZ
          //           metal_C25(i,j,k) = metal_C25(i,j,k) * correctZ
          //           metal_C30(i,j,k) = metal_C30(i,j,k) * correctZ
          //           metal_F13(i,j,k) = metal_F13(i,j,k) * correctZ
          //           metal_F15(i,j,k) = metal_F15(i,j,k) * correctZ
          //           metal_F50(i,j,k) = metal_F50(i,j,k) * correctZ
          //           metal_F80(i,j,k) = metal_F80(i,j,k) * correctZ
          //           metal_P170(i,j,k)= metal_P170(i,j,k)* correctZ
          //           metal_P200(i,j,k)= metal_P200(i,j,k)* correctZ
          //           metal_Y19(i,j,k) = metal_Y19(i,j,k) * correctZ
          //        enddo

          nSN = 12;
          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
            SN_metal(i, 0) = metal_loc(i, j, k);
            SN_metal(i, 1) = metal_C13(i, j, k);
            SN_metal(i, 2) = metal_C20(i, j, k);
            SN_metal(i, 3) = metal_C25(i, j, k);
            SN_metal(i, 4) = metal_C30(i, j, k);
            SN_metal(i, 5) = metal_F13(i, j, k);
            SN_metal(i, 6) = metal_F15(i, j, k);
            SN_metal(i, 7) = metal_F50(i, j, k);
            SN_metal(i, 8) = metal_F80(i, j, k);
            SN_metal(i, 9) = metal_P170(i, j, k);
            SN_metal(i, 10) = metal_P200(i, j, k);
            SN_metal(i, 11) = metal_Y19(i, j, k);
          }

          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
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
            for (iSN = 0; iSN < nSN; iSN++) {
              Ct[i] = Ct[i] + my_rates->SN0_XC[iSN] * SN_metal(i, iSN);
              Ot[i] = Ot[i] + my_rates->SN0_XO[iSN] * SN_metal(i, iSN);
              Mgt[i] = Mgt[i] + my_rates->SN0_XMg[iSN] * SN_metal(i, iSN);
              Alt[i] = Alt[i] + my_rates->SN0_XAl[iSN] * SN_metal(i, iSN);
              Sit[i] = Sit[i] + my_rates->SN0_XSi[iSN] * SN_metal(i, iSN);
              St[i] = St[i] + my_rates->SN0_XS[iSN] * SN_metal(i, iSN);
              Fet[i] = Fet[i] + my_rates->SN0_XFe[iSN] * SN_metal(i, iSN);

              Cg[i] = Cg[i] + my_rates->SN0_fC[iSN] * SN_metal(i, iSN);
              Og[i] = Og[i] + my_rates->SN0_fO[iSN] * SN_metal(i, iSN);
              Mgg[i] = Mgg[i] + my_rates->SN0_fMg[iSN] * SN_metal(i, iSN);
              Alg[i] = Alg[i] + my_rates->SN0_fAl[iSN] * SN_metal(i, iSN);
              Sig[i] = Sig[i] + my_rates->SN0_fSi[iSN] * SN_metal(i, iSN);
              Sg[i] = Sg[i] + my_rates->SN0_fS[iSN] * SN_metal(i, iSN);
              Feg[i] = Feg[i] + my_rates->SN0_fFe[iSN] * SN_metal(i, iSN);
            }
          }
        }

        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          Cd[i] = Ct[i] - Cg[i];
          Od[i] = Ot[i] - Og[i];
          Mgd[i] = Mgt[i] - Mgg[i];
          Ald[i] = Alt[i] - Alg[i];
          Sid[i] = Sit[i] - Sig[i];
          Sd[i] = St[i] - Sg[i];
          Fed[i] = Fet[i] - Feg[i];
        }

        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          // if (itmask_metal(i)) then
          OH(i, j, k) = std::fabs(OH(i, j, k));
          H2O(i, j, k) = std::fabs(H2O(i, j, k));
          CH(i, j, k) = std::fabs(CH(i, j, k));
          CH2(i, j, k) = std::fabs(CH2(i, j, k));
          OHII(i, j, k) = std::fabs(OHII(i, j, k));
          H2OII(i, j, k) = std::fabs(H2OII(i, j, k));
          H3OII(i, j, k) = std::fabs(H3OII(i, j, k));
          totalH[i] = totalH[i] + OH(i, j, k) / 17. + H2O(i, j, k) / 18. * 2. +
                      CH(i, j, k) / 13. + CH2(i, j, k) / 14. * 2. +
                      OHII(i, j, k) / 17. + H2OII(i, j, k) / 18. * 2. +
                      H3OII(i, j, k) / 19. * 3.;
          // endif
        }
      }

      // Correct densities by keeping fractions the same

      for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
        correctH = (gr_float)(my_chemistry->HydrogenFractionByMass *
                              metalfree[i] / totalH[i]);
        HI(i, j, k) = HI(i, j, k) * correctH;
        HII(i, j, k) = HII(i, j, k) * correctH;

        correctHe = (gr_float)((1. - my_chemistry->HydrogenFractionByMass) *
                               metalfree[i] / totalHe[i]);
        HeI(i, j, k) = HeI(i, j, k) * correctHe;
        HeII(i, j, k) = HeII(i, j, k) * correctHe;
        HeIII(i, j, k) = HeIII(i, j, k) * correctHe;

        // Correct molecular hydrogen-related fractions

        if (my_chemistry->primordial_chemistry > 1) {
          HM(i, j, k) = HM(i, j, k) * correctH;
          H2II(i, j, k) = H2II(i, j, k) * correctH;
          H2I(i, j, k) = H2I(i, j, k) * correctH;
        }
        if (my_chemistry->primordial_chemistry > 3) {
          // !          HDII (i,j,k) = HDII (i,j,k)*correctH
          HeHII(i, j, k) = HeHII(i, j, k) * correctHe;
        }
      }

      // Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (my_chemistry->primordial_chemistry > 2) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          DI(i, j, k) = std::fabs(DI(i, j, k));
          DII(i, j, k) = std::fabs(DII(i, j, k));
          HDI(i, j, k) = std::fabs(HDI(i, j, k));
          totalD = DI(i, j, k) + DII(i, j, k) + 2. / 3. * HDI(i, j, k);
          if (my_chemistry->primordial_chemistry > 3) {
            DM(i, j, k) = std::fabs(DM(i, j, k));
            HDII(i, j, k) = std::fabs(HDII(i, j, k));
            totalD = totalD + DM(i, j, k) + 2. / 3. * HDII(i, j, k);
          }
          correctD = (gr_float)(my_chemistry->HydrogenFractionByMass *
                                my_chemistry->DeuteriumToHydrogenRatio *
                                metalfree[i] / totalD);
          DI(i, j, k) = DI(i, j, k) * correctD;
          DII(i, j, k) = DII(i, j, k) * correctD;
          HDI(i, j, k) = HDI(i, j, k) * correctD;
          if (my_chemistry->primordial_chemistry > 3) {
            DM(i, j, k) = DM(i, j, k) * correctD;
            HDII(i, j, k) = HDII(i, j, k) * correctD;
          }
        }
      }

      // Do the same thing for metal species

      if (my_chemistry->metal_chemistry == 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          // if (itmask_metal(i)) then
          CI(i, j, k) = std::fabs(CI(i, j, k));
          CII(i, j, k) = std::fabs(CII(i, j, k));
          CO(i, j, k) = std::fabs(CO(i, j, k));
          CO2(i, j, k) = std::fabs(CO2(i, j, k));
          OI(i, j, k) = std::fabs(OI(i, j, k));
          OH(i, j, k) = std::fabs(OH(i, j, k));
          H2O(i, j, k) = std::fabs(H2O(i, j, k));
          O2(i, j, k) = std::fabs(O2(i, j, k));
          SiI(i, j, k) = std::fabs(SiI(i, j, k));
          SiOI(i, j, k) = std::fabs(SiOI(i, j, k));
          SiO2I(i, j, k) = std::fabs(SiO2I(i, j, k));
          CH(i, j, k) = std::fabs(CH(i, j, k));
          CH2(i, j, k) = std::fabs(CH2(i, j, k));
          COII(i, j, k) = std::fabs(COII(i, j, k));
          OII(i, j, k) = std::fabs(OII(i, j, k));
          OHII(i, j, k) = std::fabs(OHII(i, j, k));
          H2OII(i, j, k) = std::fabs(H2OII(i, j, k));
          H3OII(i, j, k) = std::fabs(H3OII(i, j, k));
          O2II(i, j, k) = std::fabs(O2II(i, j, k));
          if ((my_chemistry->grain_growth == 1) ||
              (my_chemistry->dust_sublimation == 1)) {
            if (my_chemistry->dust_species > 0) {
              Mg(i, j, k) = std::fabs(Mg(i, j, k));
            }
            if (my_chemistry->dust_species > 1) {
              Al(i, j, k) = std::fabs(Al(i, j, k));
              S(i, j, k) = std::fabs(S(i, j, k));
              Fe(i, j, k) = std::fabs(Fe(i, j, k));
            }
          }
          // endif
        }
      }

      if ((my_chemistry->grain_growth == 1) ||
          (my_chemistry->dust_sublimation == 1)) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          // if (itmask_metal(i)) then
          if (my_chemistry->dust_species > 0) {
            MgSiO3(i, j, k) = std::fabs(MgSiO3(i, j, k));
            AC(i, j, k) = std::fabs(AC(i, j, k));
          }
          if (my_chemistry->dust_species > 1) {
            SiM(i, j, k) = std::fabs(SiM(i, j, k));
            FeM(i, j, k) = std::fabs(FeM(i, j, k));
            Mg2SiO4(i, j, k) = std::fabs(Mg2SiO4(i, j, k));
            Fe3O4(i, j, k) = std::fabs(Fe3O4(i, j, k));
            SiO2D(i, j, k) = std::fabs(SiO2D(i, j, k));
            MgO(i, j, k) = std::fabs(MgO(i, j, k));
            FeS(i, j, k) = std::fabs(FeS(i, j, k));
            Al2O3(i, j, k) = std::fabs(Al2O3(i, j, k));
          }
          if (my_chemistry->dust_species > 2) {
            reforg(i, j, k) = std::fabs(reforg(i, j, k));
            volorg(i, j, k) = std::fabs(volorg(i, j, k));
            H2Oice(i, j, k) = std::fabs(H2Oice(i, j, k));
          }
          // endif
        }
      }

      if (my_chemistry->metal_chemistry == 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          // if (itmask_metal(i)) then
          //- if (d(i,j,k)*dom .lt. 1.e-2_DKIND) then
          // !       if (d(i,j,k)*dom .lt.
          // !   &    min(1.e6_DKIND/(metal(i,j,k)/d(i,j,k)/0.02d-4)**2
          // !   &       ,1.e6_DKIND)) then
          if (((imetal == 0) && (d(i, j, k) * dom < 1.e8)) ||
              ((imetal == 1) && (((metal(i, j, k) <= 1.e-9 * d(i, j, k)) &&
                                  (d(i, j, k) * dom < 1.e8)) ||
                                 ((metal(i, j, k) > 1.e-9 * d(i, j, k)) &&
                                  (d(i, j, k) * dom < 1.e6))))) {
            totalOg = 16. / 28. * CO(i, j, k) + 32. / 44. * CO2(i, j, k) +
                      OI(i, j, k) + 16. / 17. * OH(i, j, k) +
                      16. / 18. * H2O(i, j, k) + O2(i, j, k) +
                      16. / 44. * SiOI(i, j, k) + 32. / 60. * SiO2I(i, j, k) +
                      16. / 28. * COII(i, j, k) + OII(i, j, k) +
                      16. / 17. * OHII(i, j, k) + 16. / 18. * H2OII(i, j, k) +
                      16. / 19. * H3OII(i, j, k) + O2II(i, j, k);
            correctOg = (gr_float)(Og[i] / totalOg);
            CO(i, j, k) = CO(i, j, k) * correctOg;
            CO2(i, j, k) = CO2(i, j, k) * correctOg;
            OI(i, j, k) = OI(i, j, k) * correctOg;
            OH(i, j, k) = OH(i, j, k) * correctOg;
            H2O(i, j, k) = H2O(i, j, k) * correctOg;
            O2(i, j, k) = O2(i, j, k) * correctOg;
            SiOI(i, j, k) = SiOI(i, j, k) * correctOg;
            SiO2I(i, j, k) = SiO2I(i, j, k) * correctOg;
            COII(i, j, k) = COII(i, j, k) * correctOg;
            OII(i, j, k) = OII(i, j, k) * correctOg;
            OHII(i, j, k) = OHII(i, j, k) * correctOg;
            H2OII(i, j, k) = H2OII(i, j, k) * correctOg;
            H3OII(i, j, k) = H3OII(i, j, k) * correctOg;
            O2II(i, j, k) = O2II(i, j, k) * correctOg;
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalOd = 48. / 100. * MgSiO3(i, j, k);
              }
              if (my_chemistry->dust_species > 1) {
                totalOd =
                    totalOd + 64. / 140. * Mg2SiO4(i, j, k) +
                    64. / 232. * Fe3O4(i, j, k) + 32. / 60. * SiO2D(i, j, k) +
                    16. / 40. * MgO(i, j, k) + 48. / 102. * Al2O3(i, j, k);
              }
              if (my_chemistry->dust_species > 2) {
                totalOd = totalOd + 8. / 22.68 * reforg(i, j, k) +
                          16. / 32. * volorg(i, j, k) +
                          16. / 18. * H2Oice(i, j, k);
              }
              correctOd = (gr_float)(Od[i] / totalOd);
              if (my_chemistry->dust_species > 0) {
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctOd;
              }
              if (my_chemistry->dust_species > 1) {
                Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctOd;
                Fe3O4(i, j, k) = Fe3O4(i, j, k) * correctOd;
                SiO2D(i, j, k) = SiO2D(i, j, k) * correctOd;
                MgO(i, j, k) = MgO(i, j, k) * correctOd;
                Al2O3(i, j, k) = Al2O3(i, j, k) * correctOd;
              }
              if (my_chemistry->dust_species > 2) {
                reforg(i, j, k) = reforg(i, j, k) * correctOd;
                volorg(i, j, k) = volorg(i, j, k) * correctOd;
                H2Oice(i, j, k) = H2Oice(i, j, k) * correctOd;
              }
            }

            totalCg = CI(i, j, k) + CII(i, j, k) + 12. / 28. * CO(i, j, k) +
                      12. / 44. * CO2(i, j, k) + 12. / 13. * CH(i, j, k) +
                      12. / 14. * CH2(i, j, k) + 12. / 28. * COII(i, j, k);
            correctCg = (gr_float)(Cg[i] / totalCg);
            CI(i, j, k) = CI(i, j, k) * correctCg;
            CII(i, j, k) = CII(i, j, k) * correctCg;
            CO(i, j, k) = CO(i, j, k) * correctCg;
            CO2(i, j, k) = CO2(i, j, k) * correctCg;
            CH(i, j, k) = CH(i, j, k) * correctCg;
            CH2(i, j, k) = CH2(i, j, k) * correctCg;
            COII(i, j, k) = COII(i, j, k) * correctCg;
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalCd = AC(i, j, k);
              }
              if (my_chemistry->dust_species > 2) {
                totalCd = totalCd + 12. / 22.68 * reforg(i, j, k) +
                          12. / 32. * volorg(i, j, k);
              }
              correctCd = (gr_float)(Cd[i] / totalCd);
              if (my_chemistry->dust_species > 0) {
                AC(i, j, k) = AC(i, j, k) * correctCd;
              }
              if (my_chemistry->dust_species > 2) {
                reforg(i, j, k) = reforg(i, j, k) * correctCd;
                volorg(i, j, k) = volorg(i, j, k) * correctCd;
              }
            }

            totalSig = SiI(i, j, k) + 28. / 44. * SiOI(i, j, k) +
                       28. / 60. * SiO2I(i, j, k);
            correctSig = (gr_float)(Sig[i] / totalSig);
            SiI(i, j, k) = SiI(i, j, k) * correctSig;
            SiOI(i, j, k) = SiOI(i, j, k) * correctSig;
            SiO2I(i, j, k) = SiO2I(i, j, k) * correctSig;
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalSid = 28. / 100. * MgSiO3(i, j, k);
              }
              if (my_chemistry->dust_species > 1) {
                totalSid = totalSid + SiM(i, j, k) +
                           28. / 140. * Mg2SiO4(i, j, k) +
                           28. / 60. * SiO2D(i, j, k);
              }
              correctSid = (gr_float)(Sid[i] / totalSid);
              if (my_chemistry->dust_species > 0) {
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctSid;
              }
              if (my_chemistry->dust_species > 1) {
                SiM(i, j, k) = SiM(i, j, k) * correctSid;
                Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctSid;
                SiO2D(i, j, k) = SiO2D(i, j, k) * correctSid;
              }
            }

            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 1) {
                totalFeg = Fe(i, j, k);
                correctFeg = (gr_float)(Feg[i] / totalFeg);
                Fe(i, j, k) = Fe(i, j, k) * correctFeg;

                totalFed = FeM(i, j, k) + 168. / 232. * Fe3O4(i, j, k) +
                           56. / 88. * FeS(i, j, k);
                correctFed = (gr_float)(Fed[i] / totalFed);
                FeM(i, j, k) = FeM(i, j, k) * correctFed;
                Fe3O4(i, j, k) = Fe3O4(i, j, k) * correctFed;
                FeS(i, j, k) = FeS(i, j, k) * correctFed;
              }

              if (my_chemistry->dust_species > 0) {
                totalMgg = Mg(i, j, k);
                correctMgg = (gr_float)(Mgg[i] / totalMgg);
                Mg(i, j, k) = Mg(i, j, k) * correctMgg;
                totalMgd = 24. / 100. * MgSiO3(i, j, k);
                if (my_chemistry->dust_species > 1) {
                  totalMgd = totalMgd + 48. / 140. * Mg2SiO4(i, j, k) +
                             24. / 40. * MgO(i, j, k);
                }
                correctMgd = (gr_float)(Mgd[i] / totalMgd);
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctMgd;
                if (my_chemistry->dust_species > 1) {
                  Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctMgd;
                  MgO(i, j, k) = MgO(i, j, k) * correctMgd;
                }
              }

              if (my_chemistry->dust_species > 1) {
                S(i, j, k) = Sg[i];
                FeS(i, j, k) = 88. / 32. * Sd[i];

                Al(i, j, k) = Alg[i];
                Al2O3(i, j, k) = 102. / 54. * Ald[i];
              }
            }

          } else {
            totalO = 16. / 28. * CO(i, j, k) + 32. / 44. * CO2(i, j, k) +
                     OI(i, j, k) + 16. / 17. * OH(i, j, k) +
                     16. / 18. * H2O(i, j, k) + O2(i, j, k) +
                     16. / 44. * SiOI(i, j, k) + 32. / 60. * SiO2I(i, j, k) +
                     16. / 28. * COII(i, j, k) + OII(i, j, k) +
                     16. / 17. * OHII(i, j, k) + 16. / 18. * H2OII(i, j, k) +
                     16. / 19. * H3OII(i, j, k) + O2II(i, j, k);
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalO = totalO + 48. / 100. * MgSiO3(i, j, k);
              }
              if (my_chemistry->dust_species > 1) {
                totalO = totalO + 64. / 140. * Mg2SiO4(i, j, k) +
                         64. / 232. * Fe3O4(i, j, k) +
                         32. / 60. * SiO2D(i, j, k) + 16. / 40. * MgO(i, j, k) +
                         48. / 102. * Al2O3(i, j, k);
              }
              if (my_chemistry->dust_species > 2) {
                totalO = totalO + 8. / 22.68 * reforg(i, j, k) +
                         16. / 32. * volorg(i, j, k) +
                         16. / 18. * H2Oice(i, j, k);
              }
            }
            if ((my_chemistry->grain_growth == 0) &&
                (my_chemistry->dust_sublimation == 0)) {
              correctO = (gr_float)(Og[i] / totalO);
              CO(i, j, k) = CO(i, j, k) * correctO;
              CO2(i, j, k) = CO2(i, j, k) * correctO;
              OI(i, j, k) = OI(i, j, k) * correctO;
              OH(i, j, k) = OH(i, j, k) * correctO;
              H2O(i, j, k) = H2O(i, j, k) * correctO;
              O2(i, j, k) = O2(i, j, k) * correctO;
              SiOI(i, j, k) = SiOI(i, j, k) * correctO;
              SiO2I(i, j, k) = SiO2I(i, j, k) * correctO;
              COII(i, j, k) = COII(i, j, k) * correctO;
              OII(i, j, k) = OII(i, j, k) * correctO;
              OHII(i, j, k) = OHII(i, j, k) * correctO;
              H2OII(i, j, k) = H2OII(i, j, k) * correctO;
              H3OII(i, j, k) = H3OII(i, j, k) * correctO;
              O2II(i, j, k) = O2II(i, j, k) * correctO;
            } else {
              correctO = (gr_float)(Ot[i] / totalO);
              CO(i, j, k) = CO(i, j, k) * correctO;
              CO2(i, j, k) = CO2(i, j, k) * correctO;
              OI(i, j, k) = OI(i, j, k) * correctO;
              OH(i, j, k) = OH(i, j, k) * correctO;
              H2O(i, j, k) = H2O(i, j, k) * correctO;
              O2(i, j, k) = O2(i, j, k) * correctO;
              SiOI(i, j, k) = SiOI(i, j, k) * correctO;
              SiO2I(i, j, k) = SiO2I(i, j, k) * correctO;
              COII(i, j, k) = COII(i, j, k) * correctO;
              OII(i, j, k) = OII(i, j, k) * correctO;
              OHII(i, j, k) = OHII(i, j, k) * correctO;
              H2OII(i, j, k) = H2OII(i, j, k) * correctO;
              H3OII(i, j, k) = H3OII(i, j, k) * correctO;
              O2II(i, j, k) = O2II(i, j, k) * correctO;
              if (my_chemistry->dust_species > 0) {
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctO;
              }
              if (my_chemistry->dust_species > 1) {
                Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctO;
                Fe3O4(i, j, k) = Fe3O4(i, j, k) * correctO;
                SiO2D(i, j, k) = SiO2D(i, j, k) * correctO;
                MgO(i, j, k) = MgO(i, j, k) * correctO;
                Al2O3(i, j, k) = Al2O3(i, j, k) * correctO;
              }
              if (my_chemistry->dust_species > 2) {
                reforg(i, j, k) = reforg(i, j, k) * correctO;
                volorg(i, j, k) = volorg(i, j, k) * correctO;
                H2Oice(i, j, k) = H2Oice(i, j, k) * correctO;
              }
            }

            totalC = CI(i, j, k) + CII(i, j, k) + 12. / 28. * CO(i, j, k) +
                     12. / 44. * CO2(i, j, k) + 12. / 13. * CH(i, j, k) +
                     12. / 14. * CH2(i, j, k) + 12. / 28. * COII(i, j, k);
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalC = totalC + AC(i, j, k);
              }
              if (my_chemistry->dust_species > 2) {
                totalC = totalC + 12. / 22.68 * reforg(i, j, k) +
                         12. / 32. * volorg(i, j, k);
              }
            }
            if ((my_chemistry->grain_growth == 0) &&
                (my_chemistry->dust_sublimation == 0)) {
              correctC = (gr_float)(Cg[i] / totalC);
              CI(i, j, k) = CI(i, j, k) * correctC;
              CII(i, j, k) = CII(i, j, k) * correctC;
              CO(i, j, k) = CO(i, j, k) * correctC;
              CO2(i, j, k) = CO2(i, j, k) * correctC;
              CH(i, j, k) = CH(i, j, k) * correctC;
              CH2(i, j, k) = CH2(i, j, k) * correctC;
              COII(i, j, k) = COII(i, j, k) * correctC;
            } else {
              correctC = (gr_float)(Ct[i] / totalC);
              CI(i, j, k) = CI(i, j, k) * correctC;
              CII(i, j, k) = CII(i, j, k) * correctC;
              CO(i, j, k) = CO(i, j, k) * correctC;
              CO2(i, j, k) = CO2(i, j, k) * correctC;
              CH(i, j, k) = CH(i, j, k) * correctC;
              CH2(i, j, k) = CH2(i, j, k) * correctC;
              COII(i, j, k) = COII(i, j, k) * correctC;
              if (my_chemistry->dust_species > 0) {
                AC(i, j, k) = AC(i, j, k) * correctC;
              }
              if (my_chemistry->dust_species > 2) {
                reforg(i, j, k) = reforg(i, j, k) * correctC;
                volorg(i, j, k) = volorg(i, j, k) * correctC;
              }
            }

            totalSi = SiI(i, j, k) + 28. / 44. * SiOI(i, j, k) +
                      28. / 60. * SiO2I(i, j, k);
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                totalSi = totalSi + 28. / 100. * MgSiO3(i, j, k);
              }
              if (my_chemistry->dust_species > 1) {
                totalSi = totalSi + SiM(i, j, k) +
                          28. / 140. * Mg2SiO4(i, j, k) +
                          28. / 60. * SiO2D(i, j, k);
              }
            }
            if ((my_chemistry->grain_growth == 0) &&
                (my_chemistry->dust_sublimation == 0)) {
              correctSi = (gr_float)(Sig[i] / totalSi);
              SiI(i, j, k) = SiI(i, j, k) * correctSi;
              SiOI(i, j, k) = SiOI(i, j, k) * correctSi;
              SiO2I(i, j, k) = SiO2I(i, j, k) * correctSi;
            } else {
              correctSi = (gr_float)(Sit[i] / totalSi);
              SiI(i, j, k) = SiI(i, j, k) * correctSi;
              SiOI(i, j, k) = SiOI(i, j, k) * correctSi;
              SiO2I(i, j, k) = SiO2I(i, j, k) * correctSi;
              if (my_chemistry->dust_species > 0) {
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctSi;
              }
              if (my_chemistry->dust_species > 1) {
                SiM(i, j, k) = SiM(i, j, k) * correctSi;
                Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctSi;
                SiO2D(i, j, k) = SiO2D(i, j, k) * correctSi;
              }
            }

            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 1) {
                totalFe = Fe(i, j, k) + FeM(i, j, k) +
                          168. / 232. * Fe3O4(i, j, k) +
                          56. / 88. * FeS(i, j, k);
                correctFe = (gr_float)(Fet[i] / totalFe);
                Fe(i, j, k) = Fe(i, j, k) * correctFe;
                FeM(i, j, k) = FeM(i, j, k) * correctFe;
                Fe3O4(i, j, k) = Fe3O4(i, j, k) * correctFe;
                FeS(i, j, k) = FeS(i, j, k) * correctFe;
              }

              if (my_chemistry->dust_species > 0) {
                totalMg = Mg(i, j, k) + 24. / 100. * MgSiO3(i, j, k);
                if (my_chemistry->dust_species > 1) {
                  totalMg = totalMg + 48. / 140. * Mg2SiO4(i, j, k) +
                            24. / 40. * MgO(i, j, k);
                }
                correctMg = (gr_float)(Mgt[i] / totalMg);
                Mg(i, j, k) = Mg(i, j, k) * correctMg;
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * correctMg;
                if (my_chemistry->dust_species > 1) {
                  Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * correctMg;
                  MgO(i, j, k) = MgO(i, j, k) * correctMg;
                }
              }

              if (my_chemistry->dust_species > 1) {
                totalS = S(i, j, k) + 32. / 88. * FeS(i, j, k);
                correctS = (gr_float)(St[i] / totalS);
                S(i, j, k) = S(i, j, k) * correctS;
                FeS(i, j, k) = FeS(i, j, k) * correctS;

                totalAl = Al(i, j, k) + 54. / 102. * Al2O3(i, j, k);
                correctAl = (gr_float)(Alt[i] / totalAl);
                Al(i, j, k) = Al(i, j, k) * correctAl;
                Al2O3(i, j, k) = Al2O3(i, j, k) * correctAl;
              }
            }
          }

          //    CI(i,j,k)      = max(CI(i,j,k), tiny)
          //    CII(i,j,k)     = max(CII(i,j,k), tiny)
          //    CO(i,j,k)      = max(CO(i,j,k), tiny)
          //    CO2(i,j,k)     = max(CO2(i,j,k), tiny)
          //    OI(i,j,k)      = max(OI(i,j,k), tiny)
          //    OH(i,j,k)      = max(OH(i,j,k), tiny)
          //    H2O(i,j,k)     = max(H2O(i,j,k), tiny)
          //    O2(i,j,k)      = max(O2(i,j,k), tiny)
          //    SiI(i,j,k)     = max(SiI(i,j,k), tiny)
          //    SiOI(i,j,k)    = max(SiOI(i,j,k), tiny)
          //    SiO2I(i,j,k)   = max(SiO2I(i,j,k), tiny)
          //    CH(i,j,k)      = max(CH(i,j,k), tiny)
          //    CH2(i,j,k)     = max(CH2(i,j,k), tiny)
          //    COII(i,j,k)    = max(COII(i,j,k), tiny)
          //    OII(i,j,k)     = max(OII(i,j,k), tiny)
          //    OHII(i,j,k)    = max(OHII(i,j,k), tiny)
          //    H2OII(i,j,k)   = max(H2OII(i,j,k), tiny)
          //    H3OII(i,j,k)   = max(H3OII(i,j,k), tiny)
          //    O2II(i,j,k)    = max(O2II(i,j,k), tiny)
          // if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
          // if (idspecies .gt. 0) then
          //    Mg(i,j,k)      = max(Mg(i,j,k), tiny)
          //    MgSiO3(i,j,k)  = max(MgSiO3(i,j,k), tiny)
          //    AC(i,j,k)      = max(AC(i,j,k), tiny)
          // endif
          // if (idspecies .gt. 1) then
          //    Al(i,j,k)      = max(Al(i,j,k), tiny)
          //    S(i,j,k)       = max(S(i,j,k), tiny)
          //    Fe(i,j,k)      = max(Fe(i,j,k), tiny)
          //    SiM(i,j,k)     = max(SiM(i,j,k), tiny)
          //    FeM(i,j,k)     = max(FeM(i,j,k), tiny)
          //    Mg2SiO4(i,j,k) = max(Mg2SiO4(i,j,k), tiny)
          //    Fe3O4(i,j,k)   = max(Fe3O4(i,j,k), tiny)
          //    SiO2D(i,j,k)   = max(SiO2D(i,j,k), tiny)
          //    MgO(i,j,k)     = max(MgO(i,j,k), tiny)
          //    FeS(i,j,k)     = max(FeS(i,j,k), tiny)
          //    Al2O3(i,j,k)   = max(Al2O3(i,j,k), tiny)
          // endif
          // if (idspecies .gt. 2) then
          //    reforg(i,j,k)  = max(reforg(i,j,k), tiny)
          //    volorg(i,j,k)  = max(volorg(i,j,k), tiny)
          //    H2Oice(i,j,k)  = max(H2Oice(i,j,k), tiny)
          // endif
          // endif

          // endif
        }
      }

      //      if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
      //         do i = is+1, ie+1
      // !          if ( itmask_metal(i) ) then
      //            if (idspecies .gt. 0) then
      //               dust(i,j,k) = MgSiO3  (i,j,k)
      //     &                     + AC      (i,j,k)
      //            endif
      //            if (idspecies .gt. 1) then
      //               dust(i,j,k) = dust(i,j,k)
      //     &                     + SiM     (i,j,k)
      //     &                     + FeM     (i,j,k)
      //     &                     + Mg2SiO4 (i,j,k)
      //     &                     + Fe3O4   (i,j,k)
      //     &                     + SiO2D   (i,j,k)
      //     &                     + MgO     (i,j,k)
      //     &                     + FeS     (i,j,k)
      //     &                     + Al2O3   (i,j,k)
      //            endif
      //            if (idspecies .gt. 2) then
      //               dust(i,j,k) = dust(i,j,k)
      //     &                     + reforg  (i,j,k)
      //     &                     + volorg  (i,j,k)
      //     &                     + H2Oice  (i,j,k)
      //            endif
      // !          endif
      //         enddo
      //      endif

      // Set the electron density

      for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
        de(i, j, k) = HII(i, j, k) + HeII(i, j, k) / (gr_float)(4.) +
                      HeIII(i, j, k) / (gr_float)(2.);
        if (my_chemistry->primordial_chemistry > 1) {
          de(i, j, k) =
              de(i, j, k) - HM(i, j, k) + H2II(i, j, k) / (gr_float)(2.);
        }
        if (my_chemistry->primordial_chemistry > 3) {
          de(i, j, k) = de(i, j, k) - DM(i, j, k) / (gr_float)(2.) +
                        HDII(i, j, k) / (gr_float)(3.) +
                        HeHII(i, j, k) / (gr_float)(5.);
        }
        if (my_chemistry->metal_chemistry == 1) {
          // if (itmask_metal(i)) then
          de(i, j, k) = de(i, j, k) + CII(i, j, k) / (gr_float)(12.) +
                        COII(i, j, k) / (gr_float)(28.) +
                        OII(i, j, k) / (gr_float)(16.) +
                        OHII(i, j, k) / (gr_float)(17.) +
                        H2OII(i, j, k) / (gr_float)(18.) +
                        H3OII(i, j, k) / (gr_float)(19.) +
                        O2II(i, j, k) / (gr_float)(32.);
          // endif
        }
      }
    }
  }

  return;
}

}  // namespace grackle::impl
