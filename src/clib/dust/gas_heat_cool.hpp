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
#ifndef DUST_GAS_HEAT_COOL_HPP
#define DUST_GAS_HEAT_COOL_HPP

#include "fortran_func_decls.h"  // gr_mask_type
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "status_reporting.h"

#include <cmath>

/// this namespace groups logic modeling the various sources of gas
/// heating/cooling from dust-chemistry
///
/// @todo Come up with a better name?
namespace grackle::impl::dust_gas_edot {

/// update edot, in place, with contributions from photo-electric heating by
/// UV-irradiated dust
///
/// Based on eqn 1 of
/// [Wolfire+95](https://ui.adsabs.harvard.edu/abs/1995ApJ...443..152W/abstract):
/// @code{unparsed}
/// Γeff = ε * Γ
/// Γ = (1e-24 * G₀) erg cm⁻³ s⁻¹
/// @endcode
/// where G₀ is the interstellar FUV radiation (in Habing units)
///
/// The value of `photoelectric_heating` config parameter controls the
/// meaning of @p gammah:
/// - when the `photoelectric_heating` is 1, @p gammah directly holds Γ
/// - otherwise, it holds: (Γ/G₀)
///
///  `photoelectric_heating` also controls modeling of ε
///
/// @note
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
/// @param[in] gammah Parameterizes the calculation. Precomputed by
///     @ref gammah_rate.
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

/// update edot, in place, to account for energy transfer between gas and dust
/// grains.
///
/// Each of the 1D arrays is only valid for the specified @p idx_range.
///
/// The rates of energy transfer are computed, from interpolation tables, at
/// the same time that @p tdust (or @p grain_temperatures) is computed. Based
/// on the configuration, the original interpolation table was computed by
/// @ref gasGrain_rate or @ref gasGrain2_rate.
///
/// @todo
/// Consider refactoring so that the updates to @p edot are directly computed
/// at the same time as the dust temperatures. If we opt not to do this, then
/// we should probably refactor to ensure that the definition of the
/// gas-to-grain rate is more consistent between configurations (and so that we
/// don't need to pass d to this function).
///
/// @param [in,out] edot 1D array being used to accumulate the net rate of
///     change in internal energy of the gas
/// @param[in] tgas 1D array of gas temperatures
/// @param[in] tdust, grain_temperatures Specifies dust temperatures. The
///     former is used in some configurations to provide a 1D array of dust
///     temperatures for all dust. The latter is used in other configurations
///     to provide separate 1d dust temperature arrays for each modeled grain
///     species.
/// @param[in] dust2gas 1D array of ratios between dust & gas densities
/// @param[in] rhoH 1D array holding the Hydrogen mass density
/// @param[in] itmask_metal Specifies the iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] my_chemistry Holds various configuration parameters
/// @param[in] idx_range Specifies the current index-range
/// @param[in] d The 3D mass density array
/// @param[in] gasgr, gas_grainsp_heatrate Specifies gas-to-grain heat transfer
///     rates (the precise definition of this "rate" varies a lot!). The
///     former is used in some configurations to provide a 1D array of rates
///     for all dust. The latter is used in other configurations to provide
///     separate 1d dust temperature arrays for each modeled grain species.
void update_edot_dust_cooling_rate(
    double* edot, const double* tgas, const double* tdust,
    const GrainSpeciesCollection& grain_temperatures, const double* dust2gas,
    const double* rhoH, const gr_mask_type* itmask_metal,
    const chemistry_data* my_chemistry, IndexRange idx_range,
    View<gr_float***> d, const double* gasgr,
    const GrainSpeciesCollection& gas_grainsp_heatrate) {
  for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask_metal[i] != MASK_FALSE) {
      double Ldst;
      if (my_chemistry->dust_species == 0) {
        Ldst =
            -gasgr[i] * (tgas[i] - tdust[i]) * dust2gas[i] * rhoH[i] * rhoH[i];

      } else {  // my_chemistry->dust_species > 0
        if (my_chemistry->use_multiple_dust_temperatures == 0) {
          Ldst = -gasgr[i] * (tgas[i] - tdust[i]) *
                 d(i, idx_range.j, idx_range.k) * rhoH[i];
        } else {
          Ldst =
              -(gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] *
                    (tgas[i] -
                     grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i]) +
                gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] *
                    (tgas[i] -
                     grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i])) *
              d(i, idx_range.j, idx_range.k) * rhoH[i];

          if (my_chemistry->dust_species > 1) {
            Ldst =
                Ldst -
                (gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] *
                     (tgas[i] - grain_temperatures
                                    .data[OnlyGrainSpLUT::Mg2SiO4_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] *
                     (tgas[i] -
                      grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i])) *
                    d(i, idx_range.j, idx_range.k) * rhoH[i];
          }

          if (my_chemistry->dust_species > 2) {
            Ldst =
                Ldst -
                (gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] *
                     (tgas[i] - grain_temperatures
                                    .data[OnlyGrainSpLUT::ref_org_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] *
                     (tgas[i] - grain_temperatures
                                    .data[OnlyGrainSpLUT::vol_org_dust][i]) +
                 gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] *
                     (tgas[i] - grain_temperatures
                                    .data[OnlyGrainSpLUT::H2O_ice_dust][i])) *
                    d(i, idx_range.j, idx_range.k) * rhoH[i];
          }
        }
      }

      edot[i] = edot[i] + Ldst;
    }
  }
}

}  // namespace grackle::impl::dust_gas_edot

#endif  // DUST_GAS_HEAT_COOL_HPP