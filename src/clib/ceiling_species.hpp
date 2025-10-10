//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file ceiling_species.hpp
/// @brief Implements ceiling_species
/// Applies a floor to all species densities.
///
/// @todo consider renaming this (maybe enforce_species_floor)
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// ceiling_species_g function from FORTRAN to C++

#ifndef CEILING_SPECIES_HPP
#define CEILING_SPECIES_HPP

namespace grackle::impl {

inline void ceiling_species(int imetal, chemistry_data* my_chemistry,
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
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
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

  // locals

  int i, j, k;

  if (my_chemistry->primordial_chemistry > 0) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          de(i, j, k) = std::fmax(de(i, j, k), tiny_fortran_val);
          HI(i, j, k) = std::fmax(HI(i, j, k), tiny_fortran_val);
          HII(i, j, k) = std::fmax(HII(i, j, k), tiny_fortran_val);
          HeI(i, j, k) = std::fmax(HeI(i, j, k), tiny_fortran_val);
          HeII(i, j, k) = std::fmax(HeII(i, j, k), tiny_fortran_val);
          HeIII(i, j, k) =
              std::fmax(HeIII(i, j, k), (gr_float)(1e-5) * tiny_fortran_val);
        }
      }
    }
  }
  if (my_chemistry->primordial_chemistry > 1) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          HM(i, j, k) = std::fmax(HM(i, j, k), tiny_fortran_val);
          H2I(i, j, k) = std::fmax(H2I(i, j, k), tiny_fortran_val);
          H2II(i, j, k) = std::fmax(H2II(i, j, k), tiny_fortran_val);
        }
      }
    }
  }
  if (my_chemistry->primordial_chemistry > 2) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          DI(i, j, k) = std::fmax(DI(i, j, k), tiny_fortran_val);
          DII(i, j, k) = std::fmax(DII(i, j, k), tiny_fortran_val);
          HDI(i, j, k) = std::fmax(HDI(i, j, k), tiny_fortran_val);
        }
      }
    }
  }
  if (my_chemistry->primordial_chemistry > 3) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          DM(i, j, k) = std::fmax(DM(i, j, k), tiny_fortran_val);
          HDII(i, j, k) = std::fmax(HDII(i, j, k), tiny_fortran_val);
          HeHII(i, j, k) = std::fmax(HeHII(i, j, k), tiny_fortran_val);
        }
      }
    }
  }
  if (imetal == 1) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          metal(i, j, k) = std::fmax(metal(i, j, k), tiny_fortran_val);
          if (metal(i, j, k) > d(i, j, k)) {
            eprintf("WARNING: metal density exceeds  total density!\n");
            eprintf("i, j, k, metal, density =  %d %d %d %g %g\n", i, j, k,
                    metal(i, j, k), d(i, j, k));
          }
          // if( immulti .gt. 0 ) then
          //    metal_loc(i,j,k) = max(metal_loc(i,j,k), tiny)
          //    metal_C13(i,j,k) = max(metal_C13(i,j,k), tiny)
          //    metal_C20(i,j,k) = max(metal_C20(i,j,k), tiny)
          //    metal_C25(i,j,k) = max(metal_C25(i,j,k), tiny)
          //    metal_C30(i,j,k) = max(metal_C30(i,j,k), tiny)
          //    metal_F13(i,j,k) = max(metal_F13(i,j,k), tiny)
          //    metal_F15(i,j,k) = max(metal_F15(i,j,k), tiny)
          //    metal_F50(i,j,k) = max(metal_F50(i,j,k), tiny)
          //    metal_F80(i,j,k) = max(metal_F80(i,j,k), tiny)
          //    metal_P170(i,j,k)= max(metal_P170(i,j,k),tiny)
          //    metal_P200(i,j,k)= max(metal_P200(i,j,k),tiny)
          //    metal_Y19(i,j,k) = max(metal_Y19(i,j,k), tiny)
          // endif
          //- !                if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
          if (my_chemistry->metal_chemistry == 1) {
            CI(i, j, k) = std::fmax(CI(i, j, k), tiny_fortran_val);
            CII(i, j, k) = std::fmax(CII(i, j, k), tiny_fortran_val);
            CO(i, j, k) = std::fmax(CO(i, j, k), tiny_fortran_val);
            CO2(i, j, k) = std::fmax(CO2(i, j, k), tiny_fortran_val);
            OI(i, j, k) = std::fmax(OI(i, j, k), tiny_fortran_val);
            OH(i, j, k) = std::fmax(OH(i, j, k), tiny_fortran_val);
            H2O(i, j, k) = std::fmax(H2O(i, j, k), tiny_fortran_val);
            O2(i, j, k) = std::fmax(O2(i, j, k), tiny_fortran_val);
            SiI(i, j, k) = std::fmax(SiI(i, j, k), tiny_fortran_val);
            SiOI(i, j, k) = std::fmax(SiOI(i, j, k), tiny_fortran_val);
            SiO2I(i, j, k) = std::fmax(SiO2I(i, j, k), tiny_fortran_val);
            CH(i, j, k) = std::fmax(CH(i, j, k), tiny_fortran_val);
            CH2(i, j, k) = std::fmax(CH2(i, j, k), tiny_fortran_val);
            COII(i, j, k) = std::fmax(COII(i, j, k), tiny_fortran_val);
            OII(i, j, k) = std::fmax(OII(i, j, k), tiny_fortran_val);
            OHII(i, j, k) = std::fmax(OHII(i, j, k), tiny_fortran_val);
            H2OII(i, j, k) = std::fmax(H2OII(i, j, k), tiny_fortran_val);
            H3OII(i, j, k) = std::fmax(H3OII(i, j, k), tiny_fortran_val);
            O2II(i, j, k) = std::fmax(O2II(i, j, k), tiny_fortran_val);
            if ((my_chemistry->grain_growth == 1) ||
                (my_chemistry->dust_sublimation == 1)) {
              if (my_chemistry->dust_species > 0) {
                Mg(i, j, k) = std::fmax(Mg(i, j, k), tiny_fortran_val);
              }
              if (my_chemistry->dust_species > 1) {
                Al(i, j, k) = std::fmax(Al(i, j, k), tiny_fortran_val);
                S(i, j, k) = std::fmax(S(i, j, k), tiny_fortran_val);
                Fe(i, j, k) = std::fmax(Fe(i, j, k), tiny_fortran_val);
              }
            }
          }
          // !                endif
        }
      }
    }
  }
  if (my_chemistry->use_dust_density_field == 1) {
    for (k = my_fields->grid_start[2]; k <= my_fields->grid_end[2]; k++) {
      for (j = my_fields->grid_start[1]; j <= my_fields->grid_end[1]; j++) {
        for (i = my_fields->grid_start[0]; i <= my_fields->grid_end[0]; i++) {
          dust(i, j, k) = std::fmax(dust(i, j, k), tiny_fortran_val);
          if ((my_chemistry->grain_growth == 1) ||
              (my_chemistry->dust_sublimation == 1)) {
            // !                if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
            if (my_chemistry->dust_species > 0) {
              MgSiO3(i, j, k) = std::fmax(MgSiO3(i, j, k), tiny_fortran_val);
              AC(i, j, k) = std::fmax(AC(i, j, k), tiny_fortran_val);
            }
            if (my_chemistry->dust_species > 1) {
              SiM(i, j, k) = std::fmax(SiM(i, j, k), tiny_fortran_val);
              FeM(i, j, k) = std::fmax(FeM(i, j, k), tiny_fortran_val);
              Mg2SiO4(i, j, k) = std::fmax(Mg2SiO4(i, j, k), tiny_fortran_val);
              Fe3O4(i, j, k) = std::fmax(Fe3O4(i, j, k), tiny_fortran_val);
              SiO2D(i, j, k) = std::fmax(SiO2D(i, j, k), tiny_fortran_val);
              MgO(i, j, k) = std::fmax(MgO(i, j, k), tiny_fortran_val);
              FeS(i, j, k) = std::fmax(FeS(i, j, k), tiny_fortran_val);
              Al2O3(i, j, k) = std::fmax(Al2O3(i, j, k), tiny_fortran_val);
            }
            if (my_chemistry->dust_species > 2) {
              reforg(i, j, k) = std::fmax(reforg(i, j, k), tiny_fortran_val);
              volorg(i, j, k) = std::fmax(volorg(i, j, k), tiny_fortran_val);
              H2Oice(i, j, k) = std::fmax(H2Oice(i, j, k), tiny_fortran_val);
            }
            // !                endif
          }
        }
      }
    }
  }

  return;
}

}  // namespace grackle::impl

#endif /* CEILING_SPECIES_CPP_H */
