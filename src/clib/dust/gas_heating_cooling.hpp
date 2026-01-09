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
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] gammah Parameterizes the calculation
/// @param[in] idx_range Specifies the current index-range
/// @param[in] e_density 1D array of electron mass densities (see below)
/// @param[in] dom_inv
/// @param[in] isrf 1D array specifying the strength of the interstellar
///     radiation field in Habing units
///
/// @todo Perhaps we should extract the relevant parameters from my_chemistry?
///
/// @note
/// For primordial_chemistry > 0, @p e_density is typically just a copy of
/// the electron density field. Otherwise, these values are inferred
inline void update_edot_photoelectric_heat(
    double* edot, const double* tgas, const double* dust2gas,
    const double* rhoH, const gr_mask_type* itmask,
    const chemistry_data* my_chemistry, double gammah, IndexRange idx_range,
    const double* e_density, double dom_inv, const double* isrf) {
  double local_dust_to_gas_ratio = my_chemistry->local_dust_to_gas_ratio;

  // define a lambda function to actually apply the update (once we determine
  // the value of gammaha)
  auto update_edot_fn = [=](int i, double gammaha_eff) {
    edot[i] = edot[i] + gammaha_eff * rhoH[i] * dom_inv * dust2gas[i] /
                            local_dust_to_gas_ratio;
  };

  if (my_chemistry->photoelectric_heating == 1) {
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        double gammaha_eff = (tgas[i] > 2.e4) ? 0. : gammah;
        update_edot_fn(i, gammaha_eff);
      }
    }

    // Use eqn. 1 of Wolfire et al. (1995)
  } else if (my_chemistry->photoelectric_heating == 2) {
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Assume constant epsilon = 0.05.
        double gammaha_eff = (tgas[i] > 2.e4) ? 0. : gammah * 0.05 * isrf[i];
        update_edot_fn(i, gammaha_eff);
      }
    }

    // Full calculation of epsilon (eqn. 2 of Wolfire 1995)
  } else if (my_chemistry->photoelectric_heating == 3) {
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
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

void update_edot_dust_recombination(
    double* edot, const double* tgas, const double* dust2gas,
    const double* rhoH, const gr_mask_type* itmask,
    double local_dust_to_gas_ratio, const chemistry_data_storage* my_rates,
    IndexRange idx_range,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf, double dom_inv,
    const double* isrf) {
  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      cool1dmulti_buf.regr[i] =
          my_rates->regr[logTlininterp_buf.indixe[i] - 1] +
          logTlininterp_buf.tdef[i] *
              (my_rates->regr[logTlininterp_buf.indixe[i] + 1 - 1] -
               my_rates->regr[logTlininterp_buf.indixe[i] - 1]);
    }
  }

  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      double grbeta = 0.74 / std::pow(tgas[i], 0.068);
      edot[i] =
          edot[i] -
          cool1dmulti_buf.regr[i] *
              std::pow((isrf[i] * dom_inv / cool1dmulti_buf.myde[i]), grbeta) *
              cool1dmulti_buf.myde[i] * rhoH[i] * dust2gas[i] /
              local_dust_to_gas_ratio;
    }
  }
}

}  // namespace grackle::impl::dust

#endif  // DUST_GAS_HEATING_COOLING_HPP