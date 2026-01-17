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

namespace chemistry_T_detail {

/// computes temperature dependent `1 / (γ_H₂ - 1)`
///
/// This is based on equation 6 from Omukai & Nishi (1998)
/// https://ui.adsabs.harvard.edu/abs/1998ApJ...508..141O/abstract
///
/// @param x the value of `(6100 K / Tgas)`
inline double inverse_gm1_for_H2(double x) {
  return 0.5 * (5. + 2. * std::pow(x, 2) * std::exp(x) /
                         std::pow((std::exp(x) - 1), 2));
}

/// calculate the adiabatic index for gas at a given temperature
///
/// This is based on equation 5 from Omukai & Nishi (1998)
/// https://ui.adsabs.harvard.edu/abs/1998ApJ...508..141O/abstract
///
/// @param tgas temperature of the gas
/// @param nH2, n_other number densities of H2 and of all other species. In
///     practice, these can both be multiplied by a constant (typically the
///     Hydrogen mass)
/// @param gamma_other The adiabatic index of all material other than H2
inline double variable_gamma(double tgas, double nH2, double n_other,
                             double gamma_other) {
  double inv_gm1_for_H2;
  if (nH2 / n_other > 1.0e-3) {
    double x = 6100. / tgas;
    if (x > 10.) {
      inv_gm1_for_H2 = 0.5 * 5.;
    } else {
      inv_gm1_for_H2 = inverse_gm1_for_H2(x);
    }
  } else {
    inv_gm1_for_H2 = 2.5;
  }
  double gamma2 = 1. + (nH2 + n_other) / (nH2 * inv_gm1_for_H2 +
                                          n_other / (gamma_other - 1.));
  return gamma2;
}

}  // namespace chemistry_T_detail

/// calculate basic gas properties for the specified @p idx_range
inline void basic_gas_props(int imetal, double* tgas, double* mmw, double* rhoH,
                            gr_mask_type* itmask, chemistry_data* my_chemistry,
                            cloudy_data* primordial_cloudy_data,
                            grackle_field_data* my_fields,
                            InternalGrUnits internalu, IndexRange idx_range,
                            double zr) {
  // construct 3d views
  View<const gr_float***> d(my_fields->density, my_fields->grid_dimension[0],
                            my_fields->grid_dimension[1],
                            my_fields->grid_dimension[2]);
  View<const gr_float***> e(
      my_fields->internal_energy, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  View<const gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // get the appropriate constant
  const double dom = internalu_calc_dom_(internalu);

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
        *primordial_cloudy_data, my_fields, internalu, idx_range);

  } else {
    // get 3D views
    View<const gr_float***> de(
        my_fields->e_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HI(
        my_fields->HI_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HII(
        my_fields->HII_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HeI(
        my_fields->HeI_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HeII(
        my_fields->HeII_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HeIII(
        my_fields->HeIII_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> HM(
        my_fields->HM_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> H2I(
        my_fields->H2I_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    View<const gr_float***> H2II(
        my_fields->H2II_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

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
            // tgas0 is used when CALCULATE_TGAS_SELF_CONSISTENTLY is defined
            [[maybe_unused]] double tgas0 = tgas[i];
            double gamma2 = chemistry_T_detail::variable_gamma(
                tgas[i], nH2, nother, my_chemistry->Gamma);
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
