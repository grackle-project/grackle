//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implementing logic for computing gas properties.
///
//===----------------------------------------------------------------------===//
#ifndef GRACKLE_GAS_PROPS_HPP
#define GRACKLE_GAS_PROPS_HPP

#include "calc_temp1d_cloudy_g.hpp"
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.hpp"
#include "support/config.hpp"
#include "utils-cpp.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Approximate mean molecular weight of metals
///
/// @todo
/// Consider moving this to phys_constants.h. It may not be necessary since
/// this may only ultimately get used in 2 places.
inline constexpr double MU_METAL = 16.;

/// Calculate thermal pressure
///
/// @todo find a better home for this function
inline double calc_pressure(double gamma, double density,
                            double specific_eint) {
  return (gamma - 1.) * density * specific_eint;
}

/// calculate basic gas properties for the specified @p idx_range
inline void basic_gas_props(
    int imetal, double* tgas, double* mmw, double* rhoH, gr_mask_type* itmask,
    chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
    grackle_field_data* my_fields, InternalGrUnits internalu,
    IndexRange idx_range, grackle::impl::View<double***> d,
    grackle::impl::View<double***> e, grackle::impl::View<double***> de,
    grackle::impl::View<double***> HI, grackle::impl::View<double***> HII,
    grackle::impl::View<double***> HeI, grackle::impl::View<double***> HeII,
    grackle::impl::View<double***> HeIII, grackle::impl::View<double***> HM,
    grackle::impl::View<double***> H2I, grackle::impl::View<double***> H2II,
    grackle::impl::View<double***> metal, double dom, double zr) {
  // If no chemistry, use a tabulated mean molecular weight
  // and iterate to convergence.

  if (my_chemistry->primordial_chemistry == 0) {
    // fh is H mass fraction in metal-free gas.

    if (imetal == 1) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          rhoH[i] = my_chemistry->HydrogenFractionByMass *
                    (d(i, idx_range.j, idx_range.k) -
                     metal(i, idx_range.j, idx_range.k));
        }
      }
    } else {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          rhoH[i] = my_chemistry->HydrogenFractionByMass *
                    d(i, idx_range.j, idx_range.k);
        }
      }
    }

    grackle::impl::calc_temp1d_cloudy_g(
        rhoH, tgas, mmw, dom, zr, imetal, itmask, my_chemistry,
        my_rates->cloudy_primordial, my_fields, internalu, idx_range);

  } else {
    // Compute mean molecular weight (and temperature) directly

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        mmw[i] = (HeI(i, idx_range.j, idx_range.k) +
                  HeII(i, idx_range.j, idx_range.k) +
                  HeIII(i, idx_range.j, idx_range.k)) /
                     4. +
                 HI(i, idx_range.j, idx_range.k) +
                 HII(i, idx_range.j, idx_range.k) +
                 de(i, idx_range.j, idx_range.k);
        rhoH[i] =
            HI(i, idx_range.j, idx_range.k) + HII(i, idx_range.j, idx_range.k);
      }
    }

    // (include molecular hydrogen, but ignore deuterium)

    if (my_chemistry->primordial_chemistry > 1) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          mmw[i] = mmw[i] + HM(i, idx_range.j, idx_range.k) +
                   (H2I(i, idx_range.j, idx_range.k) +
                    H2II(i, idx_range.j, idx_range.k)) /
                       2.;
          rhoH[i] = rhoH[i] + H2I(i, idx_range.j, idx_range.k) +
                    H2II(i, idx_range.j, idx_range.k);
        }
      }
    }

    // Include metal species

    if (imetal == 1) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          mmw[i] = mmw[i] + metal(i, idx_range.j, idx_range.k) / MU_METAL;
        }
      }
    }

    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        double fixed_adiabat_pressure =
            calc_pressure(my_chemistry->Gamma, d(i, idx_range.j, idx_range.k),
                          e(i, idx_range.j, idx_range.k));
        tgas[i] = std::fmax(fixed_adiabat_pressure * internalu.utem / mmw[i],
                            my_chemistry->TemperatureStart);
        mmw[i] = d(i, idx_range.j, idx_range.k) / mmw[i];
      }
    }

    // Correct temperature for gamma from H2

    if (my_chemistry->primordial_chemistry > 1) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          double nH2 = 0.5 * (H2I(i, idx_range.j, idx_range.k) +
                              H2II(i, idx_range.j, idx_range.k));
          double nother = (HeI(i, idx_range.j, idx_range.k) +
                           HeII(i, idx_range.j, idx_range.k) +
                           HeIII(i, idx_range.j, idx_range.k)) /
                              4. +
                          HI(i, idx_range.j, idx_range.k) +
                          HII(i, idx_range.j, idx_range.k) +
                          de(i, idx_range.j, idx_range.k);

          int iter_tgas = 0;
          double tgas_err = huge8;
          while ((iter_tgas < 100) && (tgas_err > 1.e-3)) {
            double gamma2;
            // tgas0 is used when CALCULATE_TGAS_SELF_CONSISTENTLY is defined
            [[maybe_unused]] double tgas0 = tgas[i];
            if (nH2 / nother > 1.0e-3) {
              double x = 6100. / tgas[i];  // not quite self-consistent
              if (x > 10.) {
                gamma2 = 0.5 * 5.;
              } else {
                gamma2 = 0.5 * (5. + 2. * std::pow(x, 2) * std::exp(x) /
                                         std::pow((std::exp(x) - 1), 2));
              }
            } else {
              gamma2 = 2.5;
            }
            gamma2 =
                1. + (nH2 + nother) /
                         (nH2 * gamma2 + nother / (my_chemistry->Gamma - 1.));
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
            tgas[i] =
                std::fmax((gamma2 - 1.) * mmw[i] *
                              e(i, idx_range.j, idx_range.k) * internalu.utem,
                          my_chemistry->TemperatureStart);
            tgas_err = grackle::impl::dabs(tgas0 - tgas[i]) / tgas0;
            iter_tgas = iter_tgas + 1;
#else
            tgas[i] = tgas[i] * (gamma2 - 1.) / (my_chemistry->Gamma - 1.);
            iter_tgas = 101;
#endif
          }
        }
      }
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // GRACKLE_GAS_PROPS_HPP
