//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the lookup_cool_rates1d_g function.
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// lookup_cool_rates1d_g function from FORTRAN to C++

#ifndef LOOKUP_COOL_RATES1D_HPP
#define LOOKUP_COOL_RATES1D_HPP

#include <vector>

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "internal_types.hpp"
#include "opaque_storage.hpp"
#include "utils-cpp.hpp"

namespace grackle::impl {

/// This function tries to account for secondary electron ionization in
/// high-energy radiation fields (Shull & Steenberg 1985)
///
/// For context, when a photon ionizes a given atom, the resulting photoelectron
/// (that was previously bound) has a kinetic energy given by the difference
/// between the photon's energy and the ionization energy. The photoelectron
/// subsequently undergoes collisions that transfers excess kinetic energy:
/// - collisions with free electrons eventually convert it to heat
/// - collisions can also excite energy states bound electrons and molecular
///   energy states
/// If they have energy, photoelectrons and secondary electrons may produce
/// collisional ionization, which are called secondary ionization.
///
/// Cecondary ionization is most relevant in high-energy radiation fields.
///
/// @note
/// The impetus for factoring out this logic is that it existed behind an ifdef
/// statement that was not being used. By placing the logic in a function, we
/// can confirm that the logic is valid C++ even if the logic isn't used.
void secondary_ionization_adjustments(
    IndexRange idx_range, const gr_mask_type* itmask,
    grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
    InternalGrUnits internalu,
    grackle::impl::PhotoRxnRateCollection kshield_buf) {
  // construct views of HI_density & HII_density fields
  grackle::impl::View<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const double everg = ev2erg_grflt;
  const double e24 = 13.6;
  const double e26 = 24.6;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      double x = std::fmax(
          HII(i, idx_range.j, idx_range.k) / (HI(i, idx_range.j, idx_range.k) +
                                              HII(i, idx_range.j, idx_range.k)),
          1.0e-4);
      double factor = 0.3908 * std::pow((1. - std::pow(x, 0.4092)), 1.7592);
      kshield_buf.k24[i] =
          kshield_buf.k24[i] +
          factor * (my_uvb_rates.piHI + 0.08 * my_uvb_rates.piHeI) /
              (e24 * everg) * internalu.coolunit * internalu.tbase1;
      factor = 0.0554 * std::pow((1. - std::pow(x, 0.4614)), 1.6660);
      kshield_buf.k26[i] =
          kshield_buf.k26[i] +
          factor * (my_uvb_rates.piHI / 0.08 + my_uvb_rates.piHeI) /
              (e26 * everg) * internalu.coolunit * internalu.tbase1;
    }
  }
}

/// interpolate terms used to compute heating due to H2 formation for each each
/// index in the index-range that is also selected by the given itmask
///
/// @param[out] chemheatrates_buf A struct containing the buffers that are
///     filled by this call
/// @param[in]  idx_range Specifies the current index-range
/// @param[in]  my_rates Contains the input interpolation tables
/// @param[in]  itmask Specifies the `idx_range`'s iteration-mask for this
///    calculation
/// @param[in]  logTlininterp_buf Specifies the information related to the
///    position in the logT interpolations (for a number of chemistry zones)
///
/// @todo
/// Ideally, we wouldn't need to pass in the full
/// grackle::impl::ChemHeatingRates and chemistry_data_storage structs, and
/// we could pass in just what we need
inline void interpolate_h2_heating_terms_(
    grackle::impl::ChemHeatingRates chemheatrates_buf, IndexRange idx_range,
    chemistry_data_storage* my_rates, const gr_mask_type* itmask,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf) {
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      chemheatrates_buf.n_cr_n[i] =
          my_rates->n_cr_n[logTlininterp_buf.indixe[i] - 1] +
          (my_rates->n_cr_n[logTlininterp_buf.indixe[i]] -
           my_rates->n_cr_n[logTlininterp_buf.indixe[i] - 1]) *
              logTlininterp_buf.tdef[i];
      chemheatrates_buf.n_cr_d1[i] =
          my_rates->n_cr_d1[logTlininterp_buf.indixe[i] - 1] +
          (my_rates->n_cr_d1[logTlininterp_buf.indixe[i]] -
           my_rates->n_cr_d1[logTlininterp_buf.indixe[i] - 1]) *
              logTlininterp_buf.tdef[i];
      chemheatrates_buf.n_cr_d2[i] =
          my_rates->n_cr_d2[logTlininterp_buf.indixe[i] - 1] +
          (my_rates->n_cr_d2[logTlininterp_buf.indixe[i]] -
           my_rates->n_cr_d2[logTlininterp_buf.indixe[i] - 1]) *
              logTlininterp_buf.tdef[i];
    }
  }
}

/// interpolate the "standard" collisional reaction rates for each each index
/// in the index-range that is also selected by the given itmask
///
/// @param[out] kcol_buf A struct containing the buffers that are filled by
///    this function
/// @param[in] idx_range Specifies the current index-range
/// @param[in] kcol_rate_tables The 1D tables of "standard" collisional
///    reaction rates
/// @param[in] kcol_lut_indices A list of indices corresponding to the reaction
///    rates that need to be interpolated
/// @param[in] n_rates The number of elements in kcol_lut_indices
/// @param[in] itmask Specifies the `idx_range`'s iteration-mask for this
///    calculation
/// @param[in] logTlininterp_buf Specifies the information related to the
///    position in the logT interpolations (for a number of chemistry zones)
inline void interpolate_kcol_rate_tables_(
    grackle::impl::CollisionalRxnRateCollection kcol_buf, IndexRange idx_range,
    grackle::impl::CollisionalRxnRateCollection kcol_rate_tables,
    int* kcol_lut_indices, int n_rates, const gr_mask_type* itmask,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf) {
  // TODO: make this more efficient
  // -> to accomplish this, we probably need to account for the fact that all
  //    of the buffers within grackle::impl::CollisionalRxnRateCollection are
  //    allocated as a single buffer (in other words, they can be treated as a
  //    monolithic 2D array)
  // -> to maximize performance, we may want to change order of values in the
  //    tables and possibly change transpose the data layout. But, we may want
  //    to get GPU support working before we head down this path

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      for (int counter = 0; counter < n_rates; counter++) {
        int idx = kcol_lut_indices[counter];
        const double* table = kcol_rate_tables.data[idx];

        kcol_buf.data[idx][i] = table[logTlininterp_buf.indixe[i] - 1] +
                                (table[logTlininterp_buf.indixe[i]] -
                                 table[logTlininterp_buf.indixe[i] - 1]) *
                                    logTlininterp_buf.tdef[i];
      }
    }
  }
}

/// interpolate collisional reaction rates for each each index in the
/// index-range that is also selected by the given itmask
///
/// @param[out] kcol_buf A struct containing the buffers that are filled by
///    this function
/// @param[in] idx_range Specifies the current index-range
/// @param[in] tgas specifies the gas temperatures for the `idx_range`
/// @param[in] itmask Specifies the `idx_range`'s iteration-mask for this
///    calculation
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] my_rates Contains the input interpolation tables
/// @param[in] my_fields specifies the field data
/// @param[in] logTlininterp_buf Specifies the information related to the
///    position in the logT interpolations (for a number of chemistry zones)
inline void interpolate_collisional_rxn_rates_(
    grackle::impl::CollisionalRxnRateCollection kcol_buf, IndexRange idx_range,
    const double* tgas1d, const gr_mask_type* itmask, double dom,
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    chemistry_data_storage* my_rates,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf) {
  // There are 2 parts to this function
  // ----------------------------------

  // Part 1: Handle interpolations for most collisional rxn rates
  interpolate_kcol_rate_tables_(
      kcol_buf, idx_range, *(my_rates->opaque_storage->kcol_rate_tables),
      my_rates->opaque_storage->used_kcol_rate_indices,
      my_rates->opaque_storage->n_kcol_rate_indices, itmask, logTlininterp_buf);

  // Part 2: possibly override k13 using density dependent values
  if (my_chemistry->primordial_chemistry > 1) {
    grackle::impl::View<gr_float***> HI(
        my_fields->HI_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    // construct the view of the k13 table
    grackle::impl::View<double**> k13dda(
        my_rates->k13dd, my_chemistry->NumberOfTemperatureBins, 14);

    // define an inline function to fill an array with the 14 interpolated
    // k13dda values for an arbitrary `i`
    // -> important: `i` must satisfy: `itmask[i] != MASK_FALSE` (so that
    //    logTlininterp_buf.indixe[i] corresponds to a valid index)
    auto fill_k13dd_buf_fn = [k13dda, logTlininterp_buf](int i, double* out) {
      for (int n1 = 0; n1 < 14; n1++) {
        out[n1] = k13dda(logTlininterp_buf.indixe[i] - 1, n1) +
                  (k13dda(logTlininterp_buf.indixe[i], n1) -
                   k13dda(logTlininterp_buf.indixe[i] - 1, n1)) *
                      logTlininterp_buf.tdef[i];
      }
    };

    //   If using H2, and using the density-dependent collisional
    //     H2 dissociation rate, then replace the the density-independant
    //        k13 rate with the new one.
    // May/00: there appears to be a problem with the density-dependent
    //     collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
    if (my_chemistry->three_body_rate == 0) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          double nh = std::fmin(HI(i, idx_range.j, idx_range.k) * dom, 1.0e9);
          kcol_buf.data[CollisionalRxnLUT::k13][i] = tiny8;
          if (tgas1d[i] >= 500. && tgas1d[i] < 1.0e6) {
            // define a local buffer and fill it with values
            double k13dd[14];
            fill_k13dd_buf_fn(i, k13dd);

            // Direct collisional dissociation
            double k13_CID =
                k13dd[0] -
                k13dd[1] / (1. + std::pow((nh / k13dd[4]), k13dd[6])) +
                k13dd[2] -
                k13dd[3] / (1. + std::pow((nh / k13dd[5]), k13dd[6]));
            k13_CID = std::fmax(std::pow(10., k13_CID), tiny8);
            // Dissociative tunnelling
            double k13_DT =
                k13dd[7] -
                k13dd[8] / (1. + std::pow((nh / k13dd[11]), k13dd[13])) +
                k13dd[9] -
                k13dd[10] / (1. + std::pow((nh / k13dd[12]), k13dd[13]));
            k13_DT = std::fmax(std::pow(10., k13_DT), tiny8);
            //
            kcol_buf.data[CollisionalRxnLUT::k13][i] = k13_DT + k13_CID;
          }
        }
      }
    }
    // #define USE_PALLA_SALPETER_STAHLER1983
    // #ifdef USE_PALLA_SALPETER_STAHLER1983
    //            if (ispecies .gt. 1 .and. ithreebody .eq. 1) then
    //               do i = is+1, ie+1
    //                  if (itmask(i)) then
    //                  nh = (HI(i,j,k) + H2I(i,j,k)/2._DKIND)*dom
    //                  k13ind = 1._DKIND / (1._DKIND + nh / k13dd(i,3))
    //                  k13(i) = 10._DKIND**(
    //     &                     (1._DKIND-k13ind) * k13dd(i,2)
    //     &                             + k13ind  * k13dd(i,1) )
    //                  endif
    //               enddo
    //            endif
    // #endif
#endif  //  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
  }
}

struct ShieldFactorCalculator {
  const double* tgas1d;
  double k24_div_tbase1;
  double k26_div_tbase1;
  double crsHI;
  double crsHeI;
  double dom;
  int primordial_chemistry;
  IndexRange idx_range;

  grackle::impl::View<const gr_float***> HI;
  grackle::impl::View<const gr_float***> HII;
  grackle::impl::View<const gr_float***> HeI;
  grackle::impl::View<const gr_float***> HeII;
  grackle::impl::View<const gr_float***> HeIII;
  grackle::impl::View<const gr_float***> HM;
  grackle::impl::View<const gr_float***> H2I;
  grackle::impl::View<const gr_float***> H2II;
  grackle::impl::View<const gr_float***> DI;
  grackle::impl::View<const gr_float***> DII;
  grackle::impl::View<const gr_float***> HDI;
};

/// construct a ShieldFactorCalculator instance
inline ShieldFactorCalculator setup_shield_factor_calculator(
    const double* tgas1d, IndexRange idx_range, double dom,
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    photo_rate_storage my_uvb_rates, InternalGrUnits internalu) {
  ShieldFactorCalculator calc;
  calc.tgas1d = tgas1d;
  calc.k24_div_tbase1 = my_uvb_rates.k24 / internalu.tbase1;
  calc.k26_div_tbase1 = my_uvb_rates.k26 / internalu.tbase1;
  calc.crsHI = my_uvb_rates.crsHI;
  calc.crsHeI = my_uvb_rates.crsHeI;
  calc.dom = dom;
  calc.primordial_chemistry = my_chemistry->primordial_chemistry;
  calc.idx_range = idx_range;

  calc.HI = grackle::impl::View<const gr_float***>(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  calc.HII = grackle::impl::View<const gr_float***>(
      my_fields->HII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  calc.HeI = grackle::impl::View<const gr_float***>(
      my_fields->HeI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  calc.HeII = grackle::impl::View<const gr_float***>(
      my_fields->HeII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  calc.HeIII = grackle::impl::View<const gr_float***>(
      my_fields->HeIII_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  if (my_chemistry->primordial_chemistry > 1) {
    calc.HM = grackle::impl::View<const gr_float***>(
        my_fields->HM_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    calc.H2I = grackle::impl::View<const gr_float***>(
        my_fields->H2I_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    calc.H2II = grackle::impl::View<const gr_float***>(
        my_fields->H2II_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    if (my_chemistry->primordial_chemistry > 2) {
      calc.DI = grackle::impl::View<const gr_float***>(
          my_fields->DI_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      calc.DII = grackle::impl::View<const gr_float***>(
          my_fields->DII_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      calc.HDI = grackle::impl::View<const gr_float***>(
          my_fields->HDI_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    }
  }
  return calc;
}

struct ShieldFactor {
  double f_shield_H;
  double f_shield_He;
};

inline ShieldFactor calc_shield_factor(const ShieldFactorCalculator* calc,
                                       int i) {
  // Compute shielding factor for H
  double nSSh = 6.73e-3 * std::pow((calc->crsHI / 2.49e-18), (-2. / 3.)) *
                std::pow((calc->tgas1d[i] / 1.0e4), (0.17)) *
                std::pow((calc->k24_div_tbase1 / 1.0e-12), (2.0 / 3.0));

  // Compute the total Hydrogen number density
  double nratio = (calc->HI(i, calc->idx_range.j, calc->idx_range.k) +
                   calc->HII(i, calc->idx_range.j, calc->idx_range.k));
  if (calc->primordial_chemistry > 1) {
    nratio = nratio + calc->HM(i, calc->idx_range.j, calc->idx_range.k) +
             calc->H2I(i, calc->idx_range.j, calc->idx_range.k) +
             calc->H2II(i, calc->idx_range.j, calc->idx_range.k);

    if (calc->primordial_chemistry > 2) {
      nratio = nratio +
               0.5 * (calc->DI(i, calc->idx_range.j, calc->idx_range.k) +
                      calc->DII(i, calc->idx_range.j, calc->idx_range.k)) +
               2.0 * calc->HDI(i, calc->idx_range.j, calc->idx_range.k) / 3.0;
    }
  }

  nratio = nratio * calc->dom / nSSh;

  double f_shield_H = (0.98 * std::pow((1.0 + std::pow(nratio, 1.64)), -2.28) +
                       0.02 * std::pow(1.0 + nratio, -0.84));

  // Compute shielding factor for He

  nSSh = 6.73e-3 * std::pow((calc->crsHeI / 2.49e-18), (-2. / 3.)) *
         std::pow((calc->tgas1d[i] / 1.0e4), (0.17)) *
         std::pow((calc->k26_div_tbase1 / 1.0e-12), (2.0 / 3.0));

  nratio = 0.25 *
           (calc->HeI(i, calc->idx_range.j, calc->idx_range.k) +
            calc->HeII(i, calc->idx_range.j, calc->idx_range.k) +
            calc->HeIII(i, calc->idx_range.j, calc->idx_range.k)) *
           calc->dom / nSSh;

  double f_shield_He = (0.98 * std::pow(1.0 + std::pow(nratio, 1.64), -2.28) +
                        0.02 * std::pow(1.0 + nratio, -0.84));
  return ShieldFactor{f_shield_H, f_shield_He};
}

/// Use self-shielding factors to adjust miscellaneous photorates
///
/// @note
/// The caller should ensure that the value of `idx_range` passed to this
/// function matches the value that was used to construct `calculator`
void apply_misc_shield_factors(
    grackle::impl::PhotoRxnRateCollection kshield_buf, IndexRange idx_range,
    const gr_mask_type* itmask, int self_shielding_method,
    photo_rate_storage my_uvb_rates, const ShieldFactorCalculator* calculator) {
  if (self_shielding_method == 1) {
    // approximate self shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // to shield HI, while leaving HeI and HeII optically thin

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling. (rate k29)
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        ShieldFactor tmp = calc_shield_factor(calculator, i);
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * tmp.f_shield_H;
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * tmp.f_shield_H;
        }

        kshield_buf.k25[i] = my_uvb_rates.k25;
        kshield_buf.k26[i] = my_uvb_rates.k26;
      }
    }

  } else if (self_shielding_method == 2) {
    // Better self-shielding in HI using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // approximate self shielding in HeI and HeII

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling as HI. (rate k29)

    //   Attenuate radiation rates for H2+ dissociation (30 eV)
    //   using same scaling as HeII. (rate k28 and k30)

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        ShieldFactor tmp = calc_shield_factor(calculator, i);
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * tmp.f_shield_H;
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * tmp.f_shield_H;
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i] = 0.;
        } else {
          kshield_buf.k26[i] = kshield_buf.k26[i] * tmp.f_shield_He;
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i] = 0.0;
        } else {
          kshield_buf.k28[i] = kshield_buf.k28[i] * tmp.f_shield_He;
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i] = 0.0;
        } else {
          kshield_buf.k30[i] = kshield_buf.k30[i] * tmp.f_shield_He;
        }

        kshield_buf.k25[i] = my_uvb_rates.k25;
      }
    }

  } else if (self_shielding_method == 3) {
    // shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // in HI and HeI, but ignoring HeII heating entirely
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        ShieldFactor tmp = calc_shield_factor(calculator, i);
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * tmp.f_shield_H;
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * tmp.f_shield_H;
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i] = 0.;
        } else {
          kshield_buf.k26[i] = kshield_buf.k26[i] * tmp.f_shield_He;
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i] = 0.0;
        } else {
          kshield_buf.k28[i] = kshield_buf.k28[i] * tmp.f_shield_He;
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i] = 0.0;
        } else {
          kshield_buf.k30[i] = kshield_buf.k30[i] * tmp.f_shield_He;
        }

        kshield_buf.k25[i] = 0.0;
      }
    }
  }
}

/// This routine uses the temperature to look up the chemical rates which are
/// tabulated in a log table as a function of temperature.
///
/// > [!important]
/// > TODO: The role of the `dt` argument **MUST** be clarified! It is passed
/// > different values in different areas of the codebase!!!!
/// > - `solve_rate_cool_g` passes in the value of the total timestep that the
/// >   chemistry is evolved. This is the traditional meaning of `dt`
/// > - the time derivative calculation within `step_rate_newton_raphson`
/// >   passes the timestep of the current subcycle (effectively the whole
/// >   function is only being called for a single element idx_range)
/// >
/// > Internally, this arg only appears to be used to determine dust grain
/// > destruction rate.
/// > - the dust destruction rate is 0 for all temperatures below some
/// >   threshold (the threshold depends on the grain species)
/// > - above the threshold, the destruction rate is essentially the current
/// >   grain density divided by the value of the `dt` argument
/// >
/// > If you think about it:
/// > - I'd argue that setting `dt` to the whole timestep that we are evolving
/// >   the zone over is blatantly wrong. It violates the principle that you
/// >   should get consistent results whether you invoke grackle 100 separate
/// >   times or just 1 time. (The amount of dust heating would change)
/// > - setting `dt` to the current subcycle timestep makes a lot more sense
/// >   (and is the only logical choice)
/// >   - It is roughly equivalent to saying that dust is immediately destroyed
/// >     once the gas reaches a threshold temperature.
/// >   - the model is overly simplistic since dust grains can survive for
/// >     quite in ionized gas (see for example
/// >     https://ui.adsabs.harvard.edu/abs/2024ApJ...974...81R/abstract)
/// >
/// > If we stick with this instantaneous destruction model, then all
/// > dust-grain related heating and cooling should probably assume that the
/// > dust-grain density is already 0.
inline void lookup_cool_rates1d(
    IndexRange idx_range, gr_mask_type anydust, const double* tgas1d,
    const double* mmw, const double* tdust, const double* dust2gas,
    double* h2dust, double dom, double dx_cgs, double c_ljeans,
    const gr_mask_type* itmask, const gr_mask_type* itmask_metal, double dt,
    chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
    grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
    InternalGrUnits internalu,
    grackle::impl::GrainSpeciesCollection grain_growth_rates,
    grackle::impl::GrainSpeciesCollection grain_temperatures,
    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
    grackle::impl::CollisionalRxnRateCollection kcol_buf,
    grackle::impl::PhotoRxnRateCollection kshield_buf,
    grackle::impl::ChemHeatingRates chemheatrates_buf) {
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  // Allocate temporary buffers
  // --------------------------
  // TODO: we should create a struct that holds pre-allocated buffers that we
  //       reference here (to avoid heap allocations in this function)

  // collection of buffers for intermediate quantities used in dust-routines
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf =
      grackle::impl::new_InternalDustPropBuf(my_fields->grid_dimension[0],
                                             my_rates->gr_N[1]);

  // Construct views of fields referenced in several parts of this function.
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> H2I, H2II, kdissH2I;
  if (my_chemistry->primordial_chemistry > 1) {
    H2I = grackle::impl::View<gr_float***>(
        my_fields->H2I_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    H2II = grackle::impl::View<gr_float***>(
        my_fields->H2II_density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    if (my_chemistry->use_radiative_transfer == 1) {
      kdissH2I = grackle::impl::View<gr_float***>(
          my_fields->RT_H2_dissociation_rate, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    }
  }

  // Linearly Interpolate the Collisional Rxn Rates
  // ----------------------------------------------

  // Set log values of start and end of lookup tables
  const double logtem_start = std::log(my_chemistry->TemperatureStart);
  const double logtem_end = std::log(my_chemistry->TemperatureEnd);
  const double dlogtem = (std::log(my_chemistry->TemperatureEnd) -
                          std::log(my_chemistry->TemperatureStart)) /
                         (double)(my_chemistry->NumberOfTemperatureBins - 1);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Compute temp-centered temperature (and log)

      // logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
      logTlininterp_buf.logtem[i] = std::log(tgas1d[i]);
      logTlininterp_buf.logtem[i] = grackle::impl::clamp(
          logTlininterp_buf.logtem[i], logtem_start, logtem_end);

      // Find index into tble and precompute interpolation values

      logTlininterp_buf.indixe[i] = grackle::impl::clamp(
          (long long)((logTlininterp_buf.logtem[i] - logtem_start) / dlogtem) +
              1LL,
          1LL, (long long)my_chemistry->NumberOfTemperatureBins - 1LL);
      logTlininterp_buf.t1[i] =
          (logtem_start + (logTlininterp_buf.indixe[i] - 1) * dlogtem);
      logTlininterp_buf.t2[i] =
          (logtem_start + (logTlininterp_buf.indixe[i]) * dlogtem);
      logTlininterp_buf.tdef[i] =
          (logTlininterp_buf.logtem[i] - logTlininterp_buf.t1[i]) /
          (logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
    }
  }

  // interpolate all collisional reaction rates
  interpolate_collisional_rxn_rates_(kcol_buf, idx_range, tgas1d, itmask, dom,
                                     my_chemistry, my_fields, my_rates,
                                     logTlininterp_buf);

  // interpolate terms used to compute H2 formation heating terms.
  // (this is honestly a little out of place in this function)
  if (my_chemistry->primordial_chemistry > 1) {
    interpolate_h2_heating_terms_(chemheatrates_buf, idx_range, my_rates,
                                  itmask, logTlininterp_buf);
  }

  // Compute grain size increment

  if ((anydust != MASK_FALSE) && (my_chemistry->dust_species > 0)) {
    f_wrap::calc_grain_size_increment_1d(dom, idx_range, itmask_metal,
                                         my_chemistry, my_rates, my_fields,
                                         internal_dust_prop_buf);
  }

  // Look-up rate for H2 formation on dust & (when relevant) grain growth rates

  if (anydust != MASK_FALSE) {
    // these are some legacy variables that referene allocations now tracked
    // within internal_dust_prop_buf.
    //
    // TODO: directly access members of internal_dust_prop_buf (and get rid of
    // these demporary variables)
    double* sgSiM = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::SiM_dust];
    double* sgFeM = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::FeM_dust];
    double* sgMg2SiO4 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                            .data[OnlyGrainSpLUT::Mg2SiO4_dust];
    double* sgMgSiO3 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                           .data[OnlyGrainSpLUT::MgSiO3_dust];
    double* sgFe3O4 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                          .data[OnlyGrainSpLUT::Fe3O4_dust];
    double* sgAC = internal_dust_prop_buf.grain_sigma_per_gas_mass
                       .data[OnlyGrainSpLUT::AC_dust];
    double* sgSiO2D = internal_dust_prop_buf.grain_sigma_per_gas_mass
                          .data[OnlyGrainSpLUT::SiO2_dust];
    double* sgMgO = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::MgO_dust];
    double* sgFeS = internal_dust_prop_buf.grain_sigma_per_gas_mass
                        .data[OnlyGrainSpLUT::FeS_dust];
    double* sgAl2O3 = internal_dust_prop_buf.grain_sigma_per_gas_mass
                          .data[OnlyGrainSpLUT::Al2O3_dust];
    double* sgreforg = internal_dust_prop_buf.grain_sigma_per_gas_mass
                           .data[OnlyGrainSpLUT::ref_org_dust];
    double* sgvolorg = internal_dust_prop_buf.grain_sigma_per_gas_mass
                           .data[OnlyGrainSpLUT::vol_org_dust];
    double* sgH2Oice = internal_dust_prop_buf.grain_sigma_per_gas_mass
                           .data[OnlyGrainSpLUT::H2O_ice_dust];

    const double logTdust_start = std::log(my_chemistry->DustTemperatureStart);
    const double logTdust_end = std::log(my_chemistry->DustTemperatureEnd);
    const double dlogTdust =
        (std::log(my_chemistry->DustTemperatureEnd) -
         std::log(my_chemistry->DustTemperatureStart)) /
        (double)(my_chemistry->NumberOfDustTemperatureBins - 1);

    // we should probably enforce the following at initialization!
    GRIMPL_REQUIRE(my_chemistry->dust_species >= 0, "sanity-check!");

    if (my_chemistry->dust_species == 0) {
      // in this branch, we are just tracking a single generic dust field

      // construct a view the h2dust interpolation table
      grackle::impl::View<double**> h2dusta(
          my_rates->h2dust, my_chemistry->NumberOfTemperatureBins,
          my_chemistry->NumberOfDustTemperatureBins);

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          // Assume dust melts at Tdust > DustTemperatureEnd, in the context of
          // computing the H2 formation rate
          //
          // important: at the time of writing, whem using a generic dust
          // density field (my_chemistry->use_dust_density_field > 0), I'm
          // 99% sure that we don't mutate that density field. This contrasts
          // with Grackle's behavior when tracking dust species fields.

          if (tdust[i] > my_chemistry->DustTemperatureEnd) {
            h2dust[i] = tiny8;
          } else {
            // Get log dust temperature

            double logTdust = grackle::impl::clamp(
                std::log(tdust[i]), logTdust_start, logTdust_end);

            // Find index into table and precompute interpolation values

            long long d_indixe = grackle::impl::clamp(
                (long long)((logTdust - logTdust_start) / dlogTdust) + 1LL, 1LL,
                (long long)(my_chemistry->NumberOfDustTemperatureBins) - 1LL);
            double d_t1 = (logTdust_start + (d_indixe - 1) * dlogTdust);
            double d_t2 = (logTdust_start + d_indixe * dlogTdust);
            double d_tdef = (logTdust - d_t1) / (d_t2 - d_t1);

            // Get rate from 2D interpolation

            double dusti1 =
                h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe - 1) +
                (h2dusta(logTlininterp_buf.indixe[i], d_indixe - 1) -
                 h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe - 1)) *
                    logTlininterp_buf.tdef[i];
            double dusti2 =
                h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe) +
                (h2dusta(logTlininterp_buf.indixe[i], d_indixe) -
                 h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe)) *
                    logTlininterp_buf.tdef[i];
            h2dust[i] = dusti1 + (dusti2 - dusti1) * d_tdef;

            // Multiply by dust to gas ratio

            h2dust[i] = h2dust[i] * dust2gas[i];
          }
        }
      }

    } else {  // my_chemistry->dust_species > 0

      // before we do anything else, let's construct views of some species
      // density fields (we only need it for part of the calculation, but there
      // is enough boilerplate here that breaks everything up, a lot!)

      // The Fortran version of this function implicitly assumed that the
      // following condition was satisfied when my_chemistry->dust_species > 0
      // (we are just making it more explicit)
      GRIMPL_REQUIRE(my_chemistry->metal_chemistry == 1, "sanity-check!");

      // load views of some metal species and molecular species
      grackle::impl::View<gr_float***> CI(
          my_fields->CI_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      grackle::impl::View<gr_float***> H2O(
          my_fields->H2O_density, my_fields->grid_dimension[0],
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

      // construct some views of dust grain densities (only load in species
      // that are explicitly enabled by my_chemistry->dust_species)
      grackle::impl::View<gr_float***> MgSiO3, AC, SiM, FeM, Mg2SiO4, Fe3O4,
          SiO2D, MgO, FeS, Al2O3, reforg, volorg, H2Oice;

      if (my_chemistry->dust_species > 0) {
        MgSiO3 = grackle::impl::View<gr_float***>(
            my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        AC = grackle::impl::View<gr_float***>(
            my_fields->AC_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      }

      if (my_chemistry->dust_species > 1) {
        SiM = grackle::impl::View<gr_float***>(
            my_fields->SiM_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        FeM = grackle::impl::View<gr_float***>(
            my_fields->FeM_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        Mg2SiO4 = grackle::impl::View<gr_float***>(
            my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        Fe3O4 = grackle::impl::View<gr_float***>(
            my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        SiO2D = grackle::impl::View<gr_float***>(
            my_fields->SiO2_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        MgO = grackle::impl::View<gr_float***>(
            my_fields->MgO_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        FeS = grackle::impl::View<gr_float***>(
            my_fields->FeS_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        Al2O3 = grackle::impl::View<gr_float***>(
            my_fields->Al2O3_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      }
      if (my_chemistry->dust_species > 2) {
        reforg = grackle::impl::View<gr_float***>(
            my_fields->ref_org_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        volorg = grackle::impl::View<gr_float***>(
            my_fields->vol_org_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        H2Oice = grackle::impl::View<gr_float***>(
            my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      }

      // Look-up rate for H2 formation on dust
      // -------------------------------------
      // in this branch, we are effectively computing:
      //   h2dust[i]
      //     = ∑ₛ coef_fnₛ(Tgas[i], Tdustₛ[i]) * grain_sigma_per_gas_massₛ[i]
      // where, the "s" subscript corresponds to a dust species.
      //
      // In practice, `coef_fnₛ` is implemented through 2D interpolation
      // of tabulated values

      // load properties of the interpolation table
      gr_interp_grid_props& interp_props =
          my_rates->opaque_storage->h2dust_grain_interp_props;

      // load the tables that we are interpolating over
      //
      // TODO: remove these local variables
      // - strictly speaking, these local variables are **NOT** necessary even
      //   now. They just exist for expositionary purposes
      // - we should only delete these local variables **AFTER** the struct
      //   members holding the values have better names
      const double* h2rate_silicate_coef_table = my_rates->h2dustS;
      const double* h2rate_carbonaceous_coef_table = my_rates->h2dustC;

      // At the time of writing:
      // - we **only** use h2rate_carbonaceous_coef_table for the AC_dust
      //   species (amorphous carbon grains)
      // - we use h2rate_silicate_coef_table for **ALL** other grain species
      //   (including species that are obviously NOT silicates)
      //
      // It is not obvious at all why this is the case...
      //
      // TODO: we should add documentation clearly addressing each of the
      // following questions and replace this todo-item with a reference to the
      // part of the documentation that discusses this! The documentation
      // should make sure to address:
      // - if this choice physically motivated or a common convention? (or if
      //   it was a quick & dirty choice because Gen didn't know what to do)
      // - do we have a sense for how "wrong" the choice is? (e.g. will it
      //   generally overpredict/underpredict?)
      // - why don't we use the generic my_rates->h2dusta rate as a generic
      //   fallback for non-silicate grains instead of using h2dustS? I realize
      //   that h2dusta has very different units...

      // now its time to actually perform the calculation
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          // don't forget, we are only executing this logic if we already know
          // that my_chemistry->dust_species > 0
          if (my_chemistry->use_multiple_dust_temperatures == 0) {
            // in this branch, all grains share a single dust temperature

            double logTdust = std::log(tdust[i]);
            double h2dust_silicate_coef = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                interp_props.parameters[0], dlogTdust,
                interp_props.parameters[1], dlogtem, interp_props.data_size,
                h2rate_silicate_coef_table);

            double h2AC = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                interp_props.parameters[0], dlogTdust,
                interp_props.parameters[1], dlogtem, interp_props.data_size,
                h2rate_carbonaceous_coef_table);

            h2dust[i] = h2dust_silicate_coef * sgMgSiO3[i] + h2AC * sgAC[i];
            if (my_chemistry->dust_species > 1) {
              h2dust[i] = h2dust[i] + h2dust_silicate_coef * sgSiM[i] +
                          h2dust_silicate_coef * sgFeM[i] +
                          h2dust_silicate_coef * sgMg2SiO4[i] +
                          h2dust_silicate_coef * sgFe3O4[i] +
                          h2dust_silicate_coef * sgSiO2D[i] +
                          h2dust_silicate_coef * sgMgO[i] +
                          h2dust_silicate_coef * sgFeS[i] +
                          h2dust_silicate_coef * sgAl2O3[i];
            }
            if (my_chemistry->dust_species > 2) {
              h2dust[i] = h2dust[i] + h2dust_silicate_coef * sgreforg[i] +
                          h2dust_silicate_coef * sgvolorg[i] +
                          h2dust_silicate_coef * sgH2Oice[i];
            }

          } else {
            double logTdust = std::log(
                grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i]);
            double h2MgSiO3 = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                interp_props.parameters[0], dlogTdust,
                interp_props.parameters[1], dlogtem, interp_props.data_size,
                h2rate_silicate_coef_table);

            logTdust =
                std::log(grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i]);
            double h2AC = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                interp_props.parameters[0], dlogTdust,
                interp_props.parameters[1], dlogtem, interp_props.data_size,
                h2rate_carbonaceous_coef_table);

            h2dust[i] = h2MgSiO3 * sgMgSiO3[i] + h2AC * sgAC[i];

            if (my_chemistry->dust_species > 1) {
              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i]);
              double h2SiM = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i]);
              double h2FeM = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i]);
              double h2Mg2SiO4 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i]);
              double h2Fe3O4 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i]);
              double h2SiO2D = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i]);
              double h2MgO = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i]);
              double h2FeS = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i]);
              double h2Al2O3 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              h2dust[i] = h2dust[i] + h2SiM * sgSiM[i] + h2FeM * sgFeM[i] +
                          h2Mg2SiO4 * sgMg2SiO4[i] + h2Fe3O4 * sgFe3O4[i] +
                          h2SiO2D * sgSiO2D[i] + h2MgO * sgMgO[i] +
                          h2FeS * sgFeS[i] + h2Al2O3 * sgAl2O3[i];
            }

            if (my_chemistry->dust_species > 2) {
              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][i]);
              double h2reforg = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i]);
              double h2volorg = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i]);
              double h2H2Oice = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
                  interp_props.parameters[0], dlogTdust,
                  interp_props.parameters[1], dlogtem, interp_props.data_size,
                  h2rate_silicate_coef_table);

              h2dust[i] = h2dust[i] + h2reforg * sgreforg[i] +
                          h2volorg * sgvolorg[i] + h2H2Oice * sgH2Oice[i];
            }
          }
        }
      }

      // Compute net grain growth rates
      // ------------------------------

      long long nratec_single_elem_arr[1] = {
          (long long)(my_chemistry->NumberOfTemperatureBins)};

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          if (my_chemistry->grain_growth == 1) {
            double kd;
            if (my_chemistry->dust_species > 0) {
              kd = f_wrap::interpolate_1d_g(
                  logTlininterp_buf.logtem[i], nratec_single_elem_arr,
                  interp_props.parameters[1], dlogtem,
                  nratec_single_elem_arr[0], my_rates->grain_growth_rate);

              grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i] =
                  kd * sgMgSiO3[i] * d(i, idx_range.j, idx_range.k) *
                  grackle::impl::fmin(
                      Mg(i, idx_range.j, idx_range.k) / std::pow(24., 1.5),
                      SiOI(i, idx_range.j, idx_range.k) / std::pow(44., 1.5),
                      H2O(i, idx_range.j, idx_range.k) / std::pow(18., 1.5) /
                          2.);
              // !             if ( idsub .eq. 1 )
              // !   &            kdMgSiO3  (i) = kdMgSiO3  (i) * ( 1.d0 -
              // !   &           sqrt(tMgSiO3 (i) / tgas1d(i))
              // !   &         * exp(-min(37.2400d4/tgas1d(i) -
              // 104.872d0, 5.d1)) / ( !   &           (     Mg
              // (i,j,k)*dom/24._DKIND * kboltz*tgas1d(i)) !   &         * (
              // SiOI(i,j,k)*dom/44._DKIND * kboltz*tgas1d(i)) !   &         *
              // (2.d0*H2O (i,j,k)*dom/18._DKIND * kboltz*tgas1d(i)) !   & /
              // (1.d0*H2I (i,j,k)*dom/ 2._DKIND * kboltz*tgas1d(i)) !   & ) )
              // !             Formulation from Nozawa et al. (2003, 2012)

              grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i] =
                  kd * sgAC[i] * d(i, idx_range.j, idx_range.k) *
                  CI(i, idx_range.j, idx_range.k) / std::pow(12., 1.5);
            }

            if (my_chemistry->dust_species > 1) {
              grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i] =
                  kd * sgSiM[i] * d(i, idx_range.j, idx_range.k) *
                  SiI(i, idx_range.j, idx_range.k) / std::pow(28., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i] =
                  kd * sgFeM[i] * d(i, idx_range.j, idx_range.k) *
                  Fe(i, idx_range.j, idx_range.k) / std::pow(56., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] =
                  kd * sgMg2SiO4[i] * d(i, idx_range.j, idx_range.k) *
                  grackle::impl::fmin(
                      Mg(i, idx_range.j, idx_range.k) / std::pow(24., 1.5) / 2.,
                      SiOI(i, idx_range.j, idx_range.k) / std::pow(44., 1.5),
                      H2O(i, idx_range.j, idx_range.k) / std::pow(18., 1.5) /
                          3.);

              grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i] =
                  kd * sgFe3O4[i] * d(i, idx_range.j, idx_range.k) *
                  std::fmin(
                      Fe(i, idx_range.j, idx_range.k) / std::pow(56., 1.5) / 3.,
                      H2O(i, idx_range.j, idx_range.k) / std::pow(18., 1.5) /
                          4.);

              grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i] =
                  kd * sgSiO2D[i] * d(i, idx_range.j, idx_range.k) *
                  SiO2I(i, idx_range.j, idx_range.k) / std::pow(60., 1.5);

              grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i] =
                  kd * sgMgO[i] * d(i, idx_range.j, idx_range.k) *
                  std::fmin(
                      Mg(i, idx_range.j, idx_range.k) / std::pow(24., 1.5),
                      H2O(i, idx_range.j, idx_range.k) / std::pow(18., 1.5));

              grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i] =
                  kd * sgFeS[i] * d(i, idx_range.j, idx_range.k) *
                  std::fmin(
                      S(i, idx_range.j, idx_range.k) / std::pow(32., 1.5),
                      Fe(i, idx_range.j, idx_range.k) / std::pow(56., 1.5));

              grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i] =
                  kd * sgAl2O3[i] * d(i, idx_range.j, idx_range.k) *
                  std::fmin(
                      Al(i, idx_range.j, idx_range.k) / std::pow(27., 1.5) / 2.,
                      H2O(i, idx_range.j, idx_range.k) / std::pow(18., 1.5) /
                          3.);
            }

            // We do not consider the growth of refractory organics, volatile
            // organics, and water ice because their sublimation temperatures
            // are low (100-600 K). They sublimate before the growth occurs.
            if (my_chemistry->dust_species > 2) {
              grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i] = 0.e0;
            }
          }

          if (my_chemistry->dust_sublimation == 1) {
            if (my_chemistry->dust_species > 0) {
              grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i] = 0.e0;
            }
            if (my_chemistry->dust_species > 1) {
              grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i] = 0.e0;
            }
            if (my_chemistry->dust_species > 2) {
              grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i] = 0.e0;
              grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i] = 0.e0;
            }

            if (my_chemistry->use_multiple_dust_temperatures == 0) {
              if (my_chemistry->dust_species > 0) {
                if (tdust[i] > 1222.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i] =
                      (tiny8 - MgSiO3(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1800.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i] =
                      (tiny8 - AC(i, idx_range.j, idx_range.k)) / dt;
                }
              }

              if (my_chemistry->dust_species > 1) {
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i] =
                      (tiny8 - SiM(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i] =
                      (tiny8 - FeM(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1277.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] =
                      (tiny8 - Mg2SiO4(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i] =
                      (tiny8 - Fe3O4(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i] =
                      (tiny8 - SiO2D(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i] =
                      (tiny8 - MgO(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 680.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i] =
                      (tiny8 - FeS(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i] =
                      (tiny8 - Al2O3(i, idx_range.j, idx_range.k)) / dt;
                }
              }

              if (my_chemistry->dust_species > 2) {
                if (tdust[i] > 575.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i] =
                      (tiny8 - reforg(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 375.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i] =
                      (tiny8 - volorg(i, idx_range.j, idx_range.k)) / dt;
                }
                if (tdust[i] > 153.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i] =
                      (tiny8 - H2Oice(i, idx_range.j, idx_range.k)) / dt;
                }
              }

            } else {
              if (my_chemistry->dust_species > 0) {
                if (grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i] >
                    1222.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][i] =
                      (tiny8 - MgSiO3(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i] >
                    1800.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][i] =
                      (tiny8 - AC(i, idx_range.j, idx_range.k)) / dt;
                }
              }

              if (my_chemistry->dust_species > 1) {
                if (grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][i] =
                      (tiny8 - SiM(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][i] =
                      (tiny8 - FeM(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] >
                    1277.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] =
                      (tiny8 - Mg2SiO4(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][i] =
                      (tiny8 - Fe3O4(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][i] =
                      (tiny8 - SiO2D(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][i] =
                      (tiny8 - MgO(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i] >
                    680.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][i] =
                      (tiny8 - FeS(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i] >
                    1500.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][i] =
                      (tiny8 - Al2O3(i, idx_range.j, idx_range.k)) / dt;
                }
              }

              if (my_chemistry->dust_species > 2) {
                if (grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][i] >
                    575.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][i] =
                      (tiny8 - reforg(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i] >
                    375.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][i] =
                      (tiny8 - volorg(i, idx_range.j, idx_range.k)) / dt;
                }
                if (grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i] >
                    153.e0) {
                  grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][i] =
                      (tiny8 - H2Oice(i, idx_range.j, idx_range.k)) / dt;
                }
              }
            }
          }
        }
      }
    }
  }

  // Deal with the photo reaction rates
  // ----------------------------------
  // In the simplest configuration, these rates are same everywhere. Things get
  // messier if Grackle is configured to approximate self-shielding or use
  // custom external radiation fields.

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      kshield_buf.k24[i] = my_uvb_rates.k24;
      kshield_buf.k25[i] = my_uvb_rates.k25;
      kshield_buf.k26[i] = my_uvb_rates.k26;
      kshield_buf.k28[i] = my_uvb_rates.k28;
      kshield_buf.k29[i] = my_uvb_rates.k29;
      kshield_buf.k30[i] = my_uvb_rates.k30;
    }
  }

  // H2 self-shielding (Sobolev-like, spherically averaged, Wolcott-Green+ 2011)

  if (my_chemistry->primordial_chemistry > 1) {
    // If radiative transfer for LW photons have been already solved in
    // your hydro code, add kdissH2I later
    if (my_chemistry->use_radiative_transfer == 0 ||
        my_chemistry->radiative_transfer_use_H2_shielding == 1) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          kshield_buf.k31[i] = my_uvb_rates.k31;
        }
      }
    } else {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          kshield_buf.k31[i] =
              my_uvb_rates.k31 + kdissH2I(i, idx_range.j, idx_range.k);
        }
      }
    }

    if (my_chemistry->H2_self_shielding > 0) {
      // conditionally construct a view
      grackle::impl::View<gr_float***> xH2shield;
      if (my_chemistry->H2_self_shielding == 2) {
        xH2shield = grackle::impl::View<gr_float***>(
            my_fields->H2_self_shielding_length, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      }

      // now compute the self-shielding
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          double l_H2shield;
          // Calculate a Sobolev-like length assuming a 3D grid.
          if (my_chemistry->H2_self_shielding == 1) {
            int j = idx_range.j;
            int k = idx_range.k;
            double divrhoa[6] = {
                d(i + 1, j, k) - d(i, j, k), d(i - 1, j, k) - d(i, j, k),
                d(i, j + 1, k) - d(i, j, k), d(i, j - 1, k) - d(i, j, k),
                d(i, j, k + 1) - d(i, j, k), d(i, j, k - 1) - d(i, j, k)};
            double divrho = tiny_fortran_val;
            // Exclude directions with (drho/ds > 0)
            for (int n1 = 0; n1 < 6; n1++) {
              if (divrhoa[n1] < 0.) {
                divrho = divrho + divrhoa[n1];
              }
            }
            // (rho / divrho) is the Sobolev-like length in cell widths
            l_H2shield = std::fmin(dx_cgs * d(i, idx_range.j, idx_range.k) /
                                       (2. * std::fabs(divrho)),
                                   internalu.xbase1);

            // User-supplied length-scale field.
          } else if (my_chemistry->H2_self_shielding == 2) {
            l_H2shield =
                xH2shield(i, idx_range.j, idx_range.k) * internalu.uxyz;

            // Jeans Length
          } else if (my_chemistry->H2_self_shielding == 3) {
            l_H2shield =
                c_ljeans * std::sqrt(tgas1d[i] /
                                     (d(i, idx_range.j, idx_range.k) * mmw[i]));

          } else {
            l_H2shield = (gr_float)(0.);
          }

          double N_H2 = dom * H2I(i, idx_range.j, idx_range.k) * l_H2shield;

          // update: self-shielding following Wolcott-Green & Haiman (2019)
          // range of validity: T=100-8000 K, n<=1e7 cm^-3

          double tgas_touse = grackle::impl::clamp(tgas1d[i], 1e2, 8e3);
          double ngas_touse =
              std::fmin(d(i, idx_range.j, idx_range.k) * dom / mmw[i], 1e7);

          double aWG2019 = (0.8711 * std::log10(tgas_touse) - 1.928) *
                               std::exp(-0.2856 * std::log10(ngas_touse)) +
                           (-0.9639 * std::log10(tgas_touse) + 3.892);

          double x = 2.0e-15 * N_H2;
          double b_doppler =
              1e-5 * std::sqrt(2. * kboltz_grflt * tgas1d[i] / (2. * mh_grflt));
          double f_shield =
              0.965 / std::pow((1. + x / b_doppler), aWG2019) +
              0.035 * std::exp(-8.5e-4 * std::sqrt(1. + x)) / std::sqrt(1. + x);

          // avoid f>1
          f_shield = std::fmin(f_shield, 1.);

          kshield_buf.k31[i] = f_shield * kshield_buf.k31[i];
        }
      }
    }
    // If radiative transfer for LW photons have been already solved in
    // your hydro code, add kdissH2I here
    if (my_chemistry->use_radiative_transfer == 1 &&
        my_chemistry->radiative_transfer_use_H2_shielding == 1) {
      // write(*,*) 'kdissH2I included'
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          kshield_buf.k31[i] =
              kshield_buf.k31[i] + kdissH2I(i, idx_range.j, idx_range.k);
        }
      }
    }

    // Custom H2 shielding
    if (my_chemistry->H2_custom_shielding > 0) {
      // create a view of the field of custom shielding values
      grackle::impl::View<gr_float***> f_shield_custom(
          my_fields->H2_custom_shielding_factor, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          kshield_buf.k31[i] =
              f_shield_custom(i, idx_range.j, idx_range.k) * kshield_buf.k31[i];
        }
      }
    }
  }

  // apply some miscellaneous self-shielding adjustments
  if (my_chemistry->self_shielding_method > 0) {
    ShieldFactorCalculator calculator =
        setup_shield_factor_calculator(tgas1d, idx_range, dom, my_chemistry,
                                       my_fields, my_uvb_rates, internalu);
    apply_misc_shield_factors(
        kshield_buf, idx_range, itmask, my_chemistry->self_shielding_method,
        my_uvb_rates, &calculator);
  }

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
  // If using a high-energy radiation field, then account for
  //   effects of secondary electrons (Shull * Steenberg 1985)
  //   (see calc_rate.src)
  secondary_ionization_adjustments(idx_range, itmask, my_fields, my_uvb_rates,
                                   internalu, kshield_buf);
#endif

  drop_InternalDustPropBuf(&internal_dust_prop_buf);

  return;
}

}  // namespace grackle::impl

#endif  // LOOKUP_COOL_RATES1D_HPP
