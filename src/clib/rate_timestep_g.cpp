//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the rate_timestep_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// rate_timestep_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "LUT.hpp"
#include "utils-cpp.hpp"

#include "rate_timestep_g.hpp"

namespace grackle::impl {

void rate_timestep_g(double* dedot, double* HIdot, gr_mask_type anydust,
                     const double* h2dust, const double* rhoH,
                     const gr_mask_type* itmask, double* edot, double chunit,
                     double dom, chemistry_data* my_chemistry,
                     grackle_field_data* my_fields, IndexRange idx_range,
                     grackle::impl::CollisionalRxnRateCollection kcr_buf,
                     grackle::impl::PhotoRxnRateCollection kshield_buf,
                     grackle::impl::ChemHeatingRates chemheatrates_buf) {
  // Density fields

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
  grackle::impl::View<gr_float***> HM(
      my_fields->HM_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(
      my_fields->H2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Radiative Transfer Fields
  grackle::impl::View<gr_float***> kphHI(
      my_fields->RT_HI_ionization_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeI(
      my_fields->RT_HeI_ionization_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphHeII(
      my_fields->RT_HeII_ionization_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> HDI(
      my_fields->HDI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> CI(
      my_fields->CI_density, my_fields->grid_dimension[0],
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

  grackle::impl::View<gr_float***> kdissHDI(
      my_fields->RT_HDI_dissociation_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphCI(
      my_fields->RT_CI_ionization_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kphOI(
      my_fields->RT_OI_ionization_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissOH(
      my_fields->RT_OH_dissociation_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> kdissH2O(
      my_fields->RT_H2O_dissociation_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // locals

  int i;
  double atten;
  std::vector<double> h2heatfac(my_fields->grid_dimension[0]);
  std::vector<double> H2delta(my_fields->grid_dimension[0]);

  if (my_chemistry->primordial_chemistry == 1) {
    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Compute the electron density rate-of-change

        dedot[i] =
            +kcr_buf.data[CollisionalRxnLUT::k1][i] *
                HI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k3][i] *
                HeI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. +
            kcr_buf.data[CollisionalRxnLUT::k5][i] *
                HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. -
            kcr_buf.data[CollisionalRxnLUT::k2][i] *
                HII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k4][i] *
                HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. -
            kcr_buf.data[CollisionalRxnLUT::k6][i] *
                HeIII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. +
            kcr_buf.data[CollisionalRxnLUT::k57][i] *
                HI(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k58][i] *
                HI(i, idx_range.j, idx_range.k) *
                HeI(i, idx_range.j, idx_range.k) / 4. +
            (kshield_buf.k24[i] * HI(i, idx_range.j, idx_range.k) +
             kshield_buf.k25[i] * HeII(i, idx_range.j, idx_range.k) / 4. +
             kshield_buf.k26[i] * HeI(i, idx_range.j, idx_range.k) / 4.);

        // Compute the HI density rate-of-change

        HIdot[i] = -kcr_buf.data[CollisionalRxnLUT::k1][i] *
                       HI(i, idx_range.j, idx_range.k) *
                       de(i, idx_range.j, idx_range.k) +
                   kcr_buf.data[CollisionalRxnLUT::k2][i] *
                       HII(i, idx_range.j, idx_range.k) *
                       de(i, idx_range.j, idx_range.k) -
                   kcr_buf.data[CollisionalRxnLUT::k57][i] *
                       HI(i, idx_range.j, idx_range.k) *
                       HI(i, idx_range.j, idx_range.k) -
                   kcr_buf.data[CollisionalRxnLUT::k58][i] *
                       HI(i, idx_range.j, idx_range.k) *
                       HeI(i, idx_range.j, idx_range.k) / 4. -
                   kshield_buf.k24[i] * HI(i, idx_range.j, idx_range.k);
      }
    }
  } else {
    // Include molecular hydrogen rates for HIdot

    for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        HIdot[i] =
            -kcr_buf.data[CollisionalRxnLUT::k1][i] *
                de(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k7][i] *
                de(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k8][i] *
                HM(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k9][i] *
                HII(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k10][i] *
                H2II(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) / 2. -
            2. * kcr_buf.data[CollisionalRxnLUT::k22][i] *
                std::pow(HI(i, idx_range.j, idx_range.k), 2) *
                HI(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k2][i] *
                HII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) +
            2. * kcr_buf.data[CollisionalRxnLUT::k13][i] *
                HI(i, idx_range.j, idx_range.k) *
                H2I(i, idx_range.j, idx_range.k) / 2. +
            kcr_buf.data[CollisionalRxnLUT::k11][i] *
                HII(i, idx_range.j, idx_range.k) *
                H2I(i, idx_range.j, idx_range.k) / 2. +
            2. * kcr_buf.data[CollisionalRxnLUT::k12][i] *
                de(i, idx_range.j, idx_range.k) *
                H2I(i, idx_range.j, idx_range.k) / 2. +
            kcr_buf.data[CollisionalRxnLUT::k14][i] *
                HM(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k15][i] *
                HM(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) +
            2. * kcr_buf.data[CollisionalRxnLUT::k16][i] *
                HM(i, idx_range.j, idx_range.k) *
                HII(i, idx_range.j, idx_range.k) +
            2. * kcr_buf.data[CollisionalRxnLUT::k18][i] *
                H2II(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 2. +
            kcr_buf.data[CollisionalRxnLUT::k19][i] *
                H2II(i, idx_range.j, idx_range.k) *
                HM(i, idx_range.j, idx_range.k) / 2. -
            kcr_buf.data[CollisionalRxnLUT::k57][i] *
                HI(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k58][i] *
                HI(i, idx_range.j, idx_range.k) *
                HeI(i, idx_range.j, idx_range.k) / 4. -
            kshield_buf.k24[i] * HI(i, idx_range.j, idx_range.k) +
            2.0 * kshield_buf.k31[i] * H2I(i, idx_range.j, idx_range.k) / 2.0;

        // Add H2 formation on dust grains

        if (anydust != MASK_FALSE) {
          if (metal(i, idx_range.j, idx_range.k) >
              1.e-9 * d(i, idx_range.j, idx_range.k)) {
            HIdot[i] = HIdot[i] - 2. * h2dust[i] * rhoH[i] *
                                      HI(i, idx_range.j, idx_range.k);
          }
        }

        // Compute the electron density rate-of-change

        dedot[i] =
            +kcr_buf.data[CollisionalRxnLUT::k1][i] *
                HI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k3][i] *
                HeI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. +
            kcr_buf.data[CollisionalRxnLUT::k5][i] *
                HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. +
            kcr_buf.data[CollisionalRxnLUT::k8][i] *
                HM(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k15][i] *
                HM(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k17][i] *
                HM(i, idx_range.j, idx_range.k) *
                HII(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k14][i] *
                HM(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k2][i] *
                HII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k4][i] *
                HeII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. -
            kcr_buf.data[CollisionalRxnLUT::k6][i] *
                HeIII(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 4. -
            kcr_buf.data[CollisionalRxnLUT::k7][i] *
                HI(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) -
            kcr_buf.data[CollisionalRxnLUT::k18][i] *
                H2II(i, idx_range.j, idx_range.k) *
                de(i, idx_range.j, idx_range.k) / 2. +
            kcr_buf.data[CollisionalRxnLUT::k57][i] *
                HI(i, idx_range.j, idx_range.k) *
                HI(i, idx_range.j, idx_range.k) +
            kcr_buf.data[CollisionalRxnLUT::k58][i] *
                HI(i, idx_range.j, idx_range.k) *
                HeI(i, idx_range.j, idx_range.k) / 4. +
            (kshield_buf.k24[i] * HI(i, idx_range.j, idx_range.k) +
             kshield_buf.k25[i] * HeII(i, idx_range.j, idx_range.k) / 4. +
             kshield_buf.k26[i] * HeI(i, idx_range.j, idx_range.k) / 4.);

        // HII, HeII, HeIII recombination heating

        edot[i] = edot[i] -
                  chunit * (13.6 * (kcr_buf.data[CollisionalRxnLUT::k1][i] *
                                        HI(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k) -
                                    kcr_buf.data[CollisionalRxnLUT::k2][i] *
                                        HII(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k)) +
                            24.6 * (kcr_buf.data[CollisionalRxnLUT::k3][i] *
                                        HeI(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k) / 4. -
                                    kcr_buf.data[CollisionalRxnLUT::k4][i] *
                                        HeII(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k) / 4.) +
                            79.0 * (kcr_buf.data[CollisionalRxnLUT::k5][i] *
                                        HeII(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k) / 4. -
                                    kcr_buf.data[CollisionalRxnLUT::k6][i] *
                                        HeIII(i, idx_range.j, idx_range.k) *
                                        de(i, idx_range.j, idx_range.k) / 4.));

        // H2 formation heating

        // Equation 23 from Omukai (2000)
        h2heatfac[i] =
            std::pow((1. + (chemheatrates_buf.n_cr_n[i] /
                            (dom * (HI(i, idx_range.j, idx_range.k) *
                                        chemheatrates_buf.n_cr_d1[i] +
                                    H2I(i, idx_range.j, idx_range.k) * 0.5 *
                                        chemheatrates_buf.n_cr_d2[i])))),
                     (-1.));

        // We only want to apply this if the formation dominates, but we
        // need to apply it outside the delta calculation.

        H2delta[i] = HI(i, idx_range.j, idx_range.k) *
                     ((3.53 * kcr_buf.data[CollisionalRxnLUT::k8][i] *
                           HM(i, idx_range.j, idx_range.k) +
                       4.48 * kcr_buf.data[CollisionalRxnLUT::k22][i] *
                           std::pow(HI(i, idx_range.j, idx_range.k), 2.)) *
                          h2heatfac[i] -
                      4.48 * kcr_buf.data[CollisionalRxnLUT::k13][i] *
                          H2I(i, idx_range.j, idx_range.k) / 2.);
        // ! corrected by GC 202002

        // !          if(H2delta(i).gt.0._DKIND) then
        // !            H2delta(i) = H2delta(i) * h2heatfac(i)
        // !          endif

        if (anydust != MASK_FALSE) {
          if (metal(i, idx_range.j, idx_range.k) >
              1.e-9 * d(i, idx_range.j, idx_range.k)) {
            H2delta[i] = H2delta[i] + h2dust[i] *
                                          HI(i, idx_range.j, idx_range.k) *
                                          rhoH[i] * (0.2 + 4.2 * h2heatfac[i]);
          }
        }

        //        H2dmag = abs(H2delta)/(
        // &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2._DKIND
        // &                    + k13(i) * H2I(i,j,k)/2._DKIND))
        //        tau = (H2dmag/1e-5_DKIND)**-1.0_DKIND
        //        tau = max(tau, 1.e-5_DKIND)
        //        atten = min((1.-exp(-tau))/tau,1._DKIND)
        atten = 1.;
        edot[i] = edot[i] + chunit * H2delta[i] * atten;
        //      &       + H2I(i,j,k)*( k21(i) * HI(i,j,k)**2.0_DKIND
        //      &                    - k23(i) * H2I(i,j,k))
        // H * (k22 * H^2 - k13 * H_2) + H_2 * (k21 * H^2 - k23 * H_2) */
      }
    }
  }

  // Add photo-ionization rates if needed

  if (my_chemistry->use_radiative_transfer == 1) {
    if (my_chemistry->radiative_transfer_hydrogen_only == 0) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          HIdot[i] = HIdot[i] - kphHI(i, idx_range.j, idx_range.k) *
                                    HI(i, idx_range.j, idx_range.k);
          dedot[i] = dedot[i] +
                     kphHI(i, idx_range.j, idx_range.k) *
                         HI(i, idx_range.j, idx_range.k) +
                     kphHeI(i, idx_range.j, idx_range.k) *
                         HeI(i, idx_range.j, idx_range.k) / 4. +
                     kphHeII(i, idx_range.j, idx_range.k) *
                         HeII(i, idx_range.j, idx_range.k) / 4.;
        }
      }
    } else {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          HIdot[i] = HIdot[i] - kphHI(i, idx_range.j, idx_range.k) *
                                    HI(i, idx_range.j, idx_range.k);
          dedot[i] = dedot[i] + kphHI(i, idx_range.j, idx_range.k) *
                                    HI(i, idx_range.j, idx_range.k);
        }
      }
    }
    if ((my_chemistry->primordial_chemistry > 2) &&
        (my_chemistry->radiative_transfer_HDI_dissociation > 0)) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          HIdot[i] = HIdot[i] + kdissHDI(i, idx_range.j, idx_range.k) *
                                    HDI(i, idx_range.j, idx_range.k) / 3.0;
        }
      }
    }
    if ((my_chemistry->metal_chemistry > 0) &&
        (my_chemistry->radiative_transfer_metal_ionization > 0)) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          dedot[i] = dedot[i] +
                     kphCI(i, idx_range.j, idx_range.k) *
                         CI(i, idx_range.j, idx_range.k) / 12.0 +
                     kphOI(i, idx_range.j, idx_range.k) *
                         OI(i, idx_range.j, idx_range.k) / 16.0;
        }
      }
    }
    if ((my_chemistry->metal_chemistry > 0) &&
        (my_chemistry->radiative_transfer_metal_dissociation > 0)) {
      for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          HIdot[i] = HIdot[i] +
                     kdissOH(i, idx_range.j, idx_range.k) *
                         OH(i, idx_range.j, idx_range.k) / 17.0 +
                     kdissH2O(i, idx_range.j, idx_range.k) *
                         H2O(i, idx_range.j, idx_range.k) / 18.0;
        }
      }
    }
  }

  return;
}

}  // namespace grackle::impl
