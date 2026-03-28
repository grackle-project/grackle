//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the scale_fields function
///
//===----------------------------------------------------------------------===//

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "scale_fields.hpp"

namespace grackle::impl {

void scale_fields(int imetal, gr_float factor, chemistry_data* my_chemistry,
                  grackle_field_data* my_fields) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
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
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
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

  int dj, dk;
  dk = my_fields->grid_end[2] - my_fields->grid_start[2] + 1;
  dj = my_fields->grid_end[1] - my_fields->grid_start[1] + 1;

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  OMP_PRAGMA("omp parallel") {
    int i, j, k;

    OMP_PRAGMA("omp for schedule(runtime)")
    for (int t = 0; t <= (dk * dj - 1); t++) {
      k = t / dj + my_fields->grid_start[2];
      j = grackle::impl::mod(t, dj) + my_fields->grid_start[1];

      for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
        d(i, j, k) = d(i, j, k) * factor;
      }

      if (my_chemistry->primordial_chemistry > 0) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          de(i, j, k) = de(i, j, k) * factor;
          HI(i, j, k) = HI(i, j, k) * factor;
          HII(i, j, k) = HII(i, j, k) * factor;
          HeI(i, j, k) = HeI(i, j, k) * factor;
          HeII(i, j, k) = HeII(i, j, k) * factor;
          HeIII(i, j, k) = HeIII(i, j, k) * factor;
        }
      }
      if (my_chemistry->primordial_chemistry > 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          HM(i, j, k) = HM(i, j, k) * factor;
          H2I(i, j, k) = H2I(i, j, k) * factor;
          H2II(i, j, k) = H2II(i, j, k) * factor;
        }
      }
      if (my_chemistry->primordial_chemistry > 2) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          DI(i, j, k) = DI(i, j, k) * factor;
          DII(i, j, k) = DII(i, j, k) * factor;
          HDI(i, j, k) = HDI(i, j, k) * factor;
        }
      }
      if (my_chemistry->primordial_chemistry > 3) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          DM(i, j, k) = DM(i, j, k) * factor;
          HDII(i, j, k) = HDII(i, j, k) * factor;
          HeHII(i, j, k) = HeHII(i, j, k) * factor;
        }
      }

      if (imetal == 1) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          metal(i, j, k) = metal(i, j, k) * factor;
        }

        if (my_chemistry->metal_chemistry == 1) {
          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
            CI(i, j, k) = CI(i, j, k) * factor;
            CII(i, j, k) = CII(i, j, k) * factor;
            CO(i, j, k) = CO(i, j, k) * factor;
            CO2(i, j, k) = CO2(i, j, k) * factor;
            OI(i, j, k) = OI(i, j, k) * factor;
            OH(i, j, k) = OH(i, j, k) * factor;
            H2O(i, j, k) = H2O(i, j, k) * factor;
            O2(i, j, k) = O2(i, j, k) * factor;
            SiI(i, j, k) = SiI(i, j, k) * factor;
            SiOI(i, j, k) = SiOI(i, j, k) * factor;
            SiO2I(i, j, k) = SiO2I(i, j, k) * factor;
            CH(i, j, k) = CH(i, j, k) * factor;
            CH2(i, j, k) = CH2(i, j, k) * factor;
            COII(i, j, k) = COII(i, j, k) * factor;
            OII(i, j, k) = OII(i, j, k) * factor;
            OHII(i, j, k) = OHII(i, j, k) * factor;
            H2OII(i, j, k) = H2OII(i, j, k) * factor;
            H3OII(i, j, k) = H3OII(i, j, k) * factor;
            O2II(i, j, k) = O2II(i, j, k) * factor;
          }
        }

        if (my_chemistry->multi_metals > 0) {
          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
            metal_loc(i, j, k) = metal_loc(i, j, k) * factor;
            metal_C13(i, j, k) = metal_C13(i, j, k) * factor;
            metal_C20(i, j, k) = metal_C20(i, j, k) * factor;
            metal_C25(i, j, k) = metal_C25(i, j, k) * factor;
            metal_C30(i, j, k) = metal_C30(i, j, k) * factor;
            metal_F13(i, j, k) = metal_F13(i, j, k) * factor;
            metal_F15(i, j, k) = metal_F15(i, j, k) * factor;
            metal_F50(i, j, k) = metal_F50(i, j, k) * factor;
            metal_F80(i, j, k) = metal_F80(i, j, k) * factor;
            metal_P170(i, j, k) = metal_P170(i, j, k) * factor;
            metal_P200(i, j, k) = metal_P200(i, j, k) * factor;
            metal_Y19(i, j, k) = metal_Y19(i, j, k) * factor;
          }
        }

        if ((my_chemistry->grain_growth == 1) ||
            (my_chemistry->dust_sublimation == 1)) {
          if (my_chemistry->dust_species > 0) {
            for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0];
                 i++) {
              Mg(i, j, k) = Mg(i, j, k) * factor;
            }
          }
          if (my_chemistry->dust_species > 1) {
            for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0];
                 i++) {
              Al(i, j, k) = Al(i, j, k) * factor;
              S(i, j, k) = S(i, j, k) * factor;
              Fe(i, j, k) = Fe(i, j, k) * factor;
            }
          }
        }

        if (my_chemistry->use_dust_density_field == 1) {
          for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
            dust(i, j, k) = dust(i, j, k) * factor;
          }

          if ((my_chemistry->grain_growth == 1) ||
              (my_chemistry->dust_sublimation == 1)) {
            if (my_chemistry->dust_species > 0) {
              for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0];
                   i++) {
                MgSiO3(i, j, k) = MgSiO3(i, j, k) * factor;
                AC(i, j, k) = AC(i, j, k) * factor;
              }
            }
            if (my_chemistry->dust_species > 1) {
              for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0];
                   i++) {
                SiM(i, j, k) = SiM(i, j, k) * factor;
                FeM(i, j, k) = FeM(i, j, k) * factor;
                Mg2SiO4(i, j, k) = Mg2SiO4(i, j, k) * factor;
                Fe3O4(i, j, k) = Fe3O4(i, j, k) * factor;
                SiO2D(i, j, k) = SiO2D(i, j, k) * factor;
                MgO(i, j, k) = MgO(i, j, k) * factor;
                FeS(i, j, k) = FeS(i, j, k) * factor;
                Al2O3(i, j, k) = Al2O3(i, j, k) * factor;
              }
            }
            if (my_chemistry->dust_species > 2) {
              for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0];
                   i++) {
                reforg(i, j, k) = reforg(i, j, k) * factor;
                volorg(i, j, k) = volorg(i, j, k) * factor;
                H2Oice(i, j, k) = H2Oice(i, j, k) * factor;
              }
            }
          }
        }
      }
    }
  }  // OMP_PRAGMA("omp parallel")

  return;
}

}  // namespace grackle::impl
