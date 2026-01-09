//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// computes contributions to the heating/cooling of gas
///
//===----------------------------------------------------------------------===//
#ifndef DUST_GAS_HEATING_COOLING_HPP
#define DUST_GAS_HEATING_COOLING_HPP

#include "fortran_func_decls.h"  // gr_mask_type
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"

#include <cmath>

namespace grackle::impl::dust {

/// update edot, in place, with contributions from Photo-electric heating by
/// UV-irradiated dust
///
/// Each of the 1D arrays is only valid for the specified @p idx_range
///
/// @param [in,out] edot 1D array being used to accumulate the net rate of
///     change in thermal energy
/// @param[in] tgas 1D array of gas temperatures
/// @param[in] dust2gas 1D array of ratios between dust & gas densities
/// @param[in] rhoH 1D array holding the Hydrogen mass density
/// @param[in] e_density 1D array of electron mass densities (see below)
/// @param[in] isrf 1D array specifying the strength of the interstellar
///     radiation field in Habing units
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] gammah Parameterizes the calculation
/// @param[in] idx_range Specifies the current index-range
/// @param[in] dom_inv
///
/// @todo Perhaps we should extract the relevant parameters from my_chemistry?
///
/// @note
/// For primordial_chemistry > 0, @p e_density is typically just a copy of
/// the electron density field. Otherwise, these values are inferred
inline void update_edot_photoelectric_heat(
    double* edot, const double* tgas, const double* dust2gas,
    const double* rhoH, const double* e_density, const double* isrf,
    const gr_mask_type* itmask, const chemistry_data* my_chemistry,
    double gammah, IndexRange idx_range, double dom_inv) {
  double local_dust_to_gas_ratio = my_chemistry->local_dust_to_gas_ratio;

  // define a lambda function to actually apply the update (once we determine
  // the value of gammaha)
  auto update_edot_fn = [=](int i, double gammaha_eff) {
    edot[i] = edot[i] + gammaha_eff * rhoH[i] * dom_inv * dust2gas[i] /
                            local_dust_to_gas_ratio;
  };

  if (my_chemistry->photoelectric_heating == 1) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        double gammaha_eff = (tgas[i] > 2.e4) ? 0. : gammah;
        update_edot_fn(i, gammaha_eff);
      }
    }

    // Use eqn. 1 of Wolfire et al. (1995)
  } else if (my_chemistry->photoelectric_heating == 2) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Assume constant epsilon = 0.05.
        double gammaha_eff = (tgas[i] > 2.e4) ? 0. : gammah * 0.05 * isrf[i];
        update_edot_fn(i, gammaha_eff);
      }
    }

    // Full calculation of epsilon (eqn. 2 of Wolfire 1995)
  } else if (my_chemistry->photoelectric_heating == 3) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        double pe_X = isrf[i] * dom_inv * std::sqrt(tgas[i]) / e_density[i];
        double pe_eps = (4.9e-2 / (1. + std::pow((pe_X / 1925.), 0.73))) +
                        ((3.7e-2 * std::pow((tgas[i] / 1.e4), 0.7)) /
                         (1. + (pe_X / 5000.)));
        double gammaha_eff = gammah * pe_eps * isrf[i];
        update_edot_fn(i, gammaha_eff);
      }
    }
  }
}

/// update edot, in place, with contributions from electron recombination
/// onto positively charged dust grains.
///
/// Based upon equation 9 of
/// [Wolfire+95](https://ui.adsabs.harvard.edu/abs/1995ApJ...443..152W/abstract).
///
/// Each of the 1D arrays (other than @p regr) is only valid for the
/// specified @p idx_range
///
/// @param [in,out] edot 1D array being used to accumulate the net rate of
///     change in internal energy of the gas
/// @param[in] tgas 1D array of gas temperatures
/// @param[in] dust2gas 1D array of ratios between dust & gas densities
/// @param[in] rhoH 1D array holding the Hydrogen mass density
/// @param[in] e_density 1D array of electron mass densities (see below)
/// @param[in] isrf 1D array specifying the strength of the interstellar
///     radiation field in Habing units
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] local_dust_to_gas_ratio ratio of total dust mass to gas mass
///     in the local Universe
/// @param[in] logTlininterp_buf Hold pre-computed values for each location
///     in @p idx_range with values that are used to linearly interpolate
///     tables of values with respect to the natural log of @p tgas1. Used
///     here to interpolate the @p regr table.
/// @param[in] regr 1D table of rate values (precomputed on the shared gas
///     temperature grid using @ref regr_rate). This is just a cluster of
///     variables taken from the equation.
/// @param[in] idx_range Specifies the current index-range
/// @param[in] dom_inv
///
/// @note
/// For primordial_chemistry > 0, @p e_density is typically just a copy of
/// the electron density field. Otherwise, these values are inferred.
inline void update_edot_dust_recombination(
    double* edot, const double* tgas, const double* dust2gas,
    const double* rhoH, const double* e_density, const double* isrf,
    const gr_mask_type* itmask, double local_dust_to_gas_ratio,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
    const double* regr, IndexRange idx_range, double dom_inv) {
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      double cur_regr =
          regr[logTlininterp_buf.indixe[i] - 1] +
          logTlininterp_buf.tdef[i] * (regr[logTlininterp_buf.indixe[i]] -
                                       regr[logTlininterp_buf.indixe[i] - 1]);

      double grbeta = 0.74 / std::pow(tgas[i], 0.068);
      edot[i] =
          edot[i] -
          cur_regr * std::pow((isrf[i] * dom_inv / e_density[i]), grbeta) *
              e_density[i] * rhoH[i] * dust2gas[i] / local_dust_to_gas_ratio;
    }
  }
}

}  // namespace grackle::impl::dust

#endif  // DUST_GAS_HEATING_COOLING_HPP