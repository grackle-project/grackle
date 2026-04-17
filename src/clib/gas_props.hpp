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
  bool exceed_min_T = tgas >= 610.0;  // aka x<=10.0
  bool use_formula = (nH2 > (1.0e-3 * n_other)) && exceed_min_T;
  double inv_gm1_for_H2 = use_formula ? inverse_gm1_for_H2(6100. / tgas) : 2.5;
  return 1. + (nH2 + n_other) /
                  (nH2 * inv_gm1_for_H2 + n_other / (gamma_other - 1.));
}

/// use "fixed point iteration" to solve for a self-consistent gas temperature
///
/// For the uninitiated, "fixed point iteration" is a scheme where we solve for
/// a so-called "fixed-point" of a function `f(x)`. This scheme produces a
/// sequence of values `x₀, x₁, x₂, ...`, where `xₙ₊₁ = f(xₙ)`. Essentially, we
/// "hope" that the sequence converges to a fixed point `xfix` (i.e. where
/// `xfix == f(xfix)`).
///
/// In the context, we try to iterate to a self-consistent gas temperature
/// using the variable adiabatic index implemented in @ref variable_gamma.
///
/// @param tgas0 Initial guess of the gas temperature
/// @param nH2, n_other number densities of H2 and of all other species. In
///     practice, these can both be multiplied by a constant (typically the
///     Hydrogen mass)
/// @param gamma_other The adiabatic index of all material other than H2
/// @param tgas_div_gm1 The gas temperature divided by `(gamma - 1)`. This is
///     always exactly known.
/// @param tgas_floor The floor to apply to the gas temperature
/// @param n_iter Number of iterations
inline double self_consistent_Tgas(double tgas0, double nH2, double n_other,
                                   double gamma_other, double tgas_div_gm1,
                                   double tgas_floor, unsigned int n_iter) {
  double tgas = tgas0;

  for (unsigned int n = 0; n < n_iter; n++) {
    double gamma = variable_gamma(tgas0, nH2, n_other, gamma_other);
    tgas = std::fmax((gamma - 1.) * tgas_div_gm1, tgas_floor);
    if (std::fabs(tgas0 - tgas) <= tgas0 * 1.0e-3) {
      break;
    };
    tgas0 = tgas;
  }
  return tgas;
}

}  // namespace chemistry_T_detail

/// calculate basic gas properties for the specified @p idx_range
///
/// Basic properties include @p tgas, @p mmw, and @p rhoH. For some context,
/// we effectively compute @p tgas and @p mmw at the same time, and @p rhoH is
/// required for the calculation along some execution pathways.
///
/// @param[out] tgas 1D array to hold the computed gas temperatures in the
///     @p idx_range
/// @param[out] mmw 1D array to hold the computed mean molecular weight
///     in the @p idx_range
/// @param[out] rhoH 1D array to hold the computed Hydrogen mass density
///     for the @p idx_range
/// @param[in] imetal Indicates whether metals are evolved
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] primordial_cloudy_data Cloudy cooling table data for primordial
///     species. (Only used when my_chemistry->primordial_chemistry is 0)
/// @param[in] my_fields Specifies the field data.
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[in] idx_range Specifies the current index-range
/// @param[in] zr Current redshift
inline void basic_gas_props(double* tgas, double* mmw, double* rhoH, int imetal,
                            const gr_mask_type* itmask,
                            chemistry_data* my_chemistry,
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

    // If the chemical network includes H2, adjust tgas to account for the fact
    // that adiabatic index actually depends on H2 abundance & on tgas
    //
    // IDEA for the future if we eventually decide that the accuracy of this
    // quantity needs to be improved. We should adjust the initial guess for
    // tgas (in cases where the `Gamma` parameter is 5/3):
    // - the new initial guess should assume that gamma=7/5 for H2 (accounting
    //   for rotational degrees of freedom). Then the adjustment only adjusts
    //   for the relevance of the single vibrational degree of freedom
    // - the current initial guess implicitly assumes H2's gamma is the same
    //   as monatomic species (so it's initially less accurate)
    if (my_chemistry->primordial_chemistry > 1) {
      for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        if (itmask[i] != MASK_FALSE) {
          // compute number densities of H2 and all other primordial species
          // -> technically, we're computing number density time Hydrogen mass,
          //    but that's ok for our purposes
          double nH2 = 0.5 * (H2I(i, idx_range.j, idx_range.k) +
                              H2II(i, idx_range.j, idx_range.k));
          double nother = (HeI(i, idx_range.j, idx_range.k) +
                           HeII(i, idx_range.j, idx_range.k) +
                           HeIII(i, idx_range.j, idx_range.k)) /
                              4. +
                          HI(i, idx_range.j, idx_range.k) +
                          HII(i, idx_range.j, idx_range.k) +
                          de(i, idx_range.j, idx_range.k);

          // the quality of the adjustment varies between branches
#ifdef CALCULATE_TGAS_SELF_CONSISTENTLY
          // iteratively solve for a self-consistent gas temperature
          // -> the current strategy uses fixed point iteration. This isn't
          //    practical (and it's not totally obvious that it converges)
          constexpr unsigned int n_iter = 100u;
          const double tgas_div_gm1 =
              mmw[i] * e(i, idx_range.j, idx_range.k) * internalu.utem;
          const double tgas_floor = my_chemistry->TemperatureStart;
          tgas[i] = chemistry_T_detail::self_consistent_Tgas(
              tgas[i], nH2, nother, my_chemistry->Gamma, tgas_div_gm1,
              tgas_floor, n_iter);
#else
          // in this branch, we roughly approximate the adjustment
          // -> it's equivalent to a single iteration of fixed point iteration
          // -> the resulting tgas is technically not self-consistent (i.e.
          //    applying the adjustment again would yield a different answer)
          double adjusted_gamma = chemistry_T_detail::variable_gamma(
              tgas[i], nH2, nother, my_chemistry->Gamma);
          tgas[i] =
              tgas[i] * (adjusted_gamma - 1.) / (my_chemistry->Gamma - 1.);
#endif
        }
      }
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // GRACKLE_GAS_PROPS_HPP
