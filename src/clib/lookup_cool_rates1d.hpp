//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the lookup_cool_rates1d function.
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// lookup_cool_rates1d_g function from FORTRAN to C++

#ifndef LOOKUP_COOL_RATES1D_HPP
#define LOOKUP_COOL_RATES1D_HPP

#include "grackle.h"
#include "dust_props.hpp"
#include "dust/lookup_dust_rates1d.hpp"
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
/// @param[in] tgas1d specifies the gas temperatures for the `idx_range`
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

/// adjust the rate of neutral H2 photodissoication by modelling self-shielding
inline void model_H2I_dissociation_shielding(
    grackle::impl::PhotoRxnRateCollection kshield_buf, IndexRange idx_range,
    const double* tgas1d, const double* mmw, double dom, double dx_cgs,
    double c_ljeans, const gr_mask_type* itmask, chemistry_data* my_chemistry,
    grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
    InternalGrUnits internalu) {
  if (my_chemistry->primordial_chemistry <= 1) {
    return;
  }

  // Construct views of fields referenced in several parts of this function.
  grackle::impl::View<const gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<const gr_float***> HI(
      my_fields->HI_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<const gr_float***> H2I(
      my_fields->H2I_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<const gr_float***> H2II(
      my_fields->H2II_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<const gr_float***> kdissH2I;
  if (my_chemistry->use_radiative_transfer == 1) {
    kdissH2I = grackle::impl::View<const gr_float***>(
        my_fields->RT_H2_dissociation_rate, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  }

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
    grackle::impl::View<const gr_float***> xH2shield;
    if (my_chemistry->H2_self_shielding == 2) {
      xH2shield = grackle::impl::View<const gr_float***>(
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
          l_H2shield = xH2shield(i, idx_range.j, idx_range.k) * internalu.uxyz;

          // Jeans Length
        } else if (my_chemistry->H2_self_shielding == 3) {
          l_H2shield =
              c_ljeans *
              std::sqrt(tgas1d[i] / (d(i, idx_range.j, idx_range.k) * mmw[i]));

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
    grackle::impl::View<const gr_float***> f_shield_custom(
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

// ShieldFactorCalculator should be treated like a "callable class." In
// idiomatic C++:
// - setup_shield_factor_calculator would be the constructor
// - calc_shield_factor would be `ShieldFactorCalculator::operator()` method
//   (which is like __call__ in python)

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
inline void apply_misc_shield_factors(
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

/// This routine uses the gas temperature to calculate rate at each location
/// in the specified index range.
///
/// In more detail, this function does a lot (probably too much):
/// - it computes collisional reaction rates
/// - shielding-adjusted photo-rates (related to the UV background)
/// - a few heating/cooling rates
/// - dust-related rates (details depend on the dust model)
/// - logTlininterp_buf is considered an output too (at the time of writing, it
///   is used for some subsequent calculations)
///
/// > [!note]
/// > A case could be made to handle the dust-related rates in a separate
/// > function
///
/// @param[in] idx_range Specifies the current index-range
/// @param[in] anydust Whether to model any dust
/// @param[in] tgas1d specifies the gas temperatures for the @p idx_range
/// @param[in] mmw specifies the mean molecular weight for the @p idx_range
/// @param[in] tdust Precomputed dust temperatures at each location in the
///     index range. This **ONLY** holds meaningful values when using variants
///     of the classic 1-field dust-model or using variant of the
///     multi-grain-species model where all grains are configured to share a
///     single temperature.
/// @param[in] dust2gas Holds the dust-to-gas ratio at each location in the
///     index range. In other words, this holds the dust mass per unit gas mass
///     (only used in certain configuration)
/// @param[out] h2dust Buffer that gets filled with the rate for forming
///     molecular hydrogen on dust grains. (This is filled whenever @p anydust
///     holds `MASK_TRUE`.
/// @param[in] dom a standard quantity used throughout the codebase
/// @param[in] dx_cgs The width of a cell in comoving cm (I think). Used in
///     certain self-shielding calculations.
/// @param[in] c_ljeans Coefficient used for computing the Jeans length
/// @param[in] itmask Specifies the general iteration-mask of the @p idx_range
///     for this calculation.
/// @param[in] itmask_metal Specifies the iteration-mask of the @p idx_range for
///     performing metal and dust calculations.
/// @param[in] dt See the warning at the end of the docstring
/// @param[in] my_chemistry holds a number of configuration parameters.
/// @param[in] my_rates Holds assorted rate data and other internal
///     configuration info.
/// @param[in] my_fields Specifies the field data.
/// @param[in] my_uvb_rates Holds precomputed photorates that depend on the UV
///     background. These rates do not include the effects of self-shielding.
/// @param[in] internalu Specifies Grackle's internal unit-system
/// @param[out] grain_growth_rates output buffers that are used to hold the
///     net grain growth rates at each @p idx_range (only used in certain
///     configurations)
/// @param[in] grain_temperatures individual grain species temperatures. This
///     is only used in certain configurations (i.e. when we aren't using the
///     tdust argument)
/// @param[out] logTlininterp_buf Buffers that are filled with values for each
///     location in @p idx_range with valuea that are used to linearly
///     interpolate tables with respect to the natural log of @p tgas1d
/// @param[out] kcol_buf Buffers filled with the collisional reaction rates for
///     each location in @p idx_range
/// @param[out] kshield_buf Buffers filled with shielding-adjusted photo
///     reaction rates for @p idx_range
/// @param[out] chemheatrates_buf Buffers that are filled with interpolated
///     values that are used to compute heating from certain chemical reactions.
/// @param[inout] internal_dust_prop_scratch_buf Scratch space used to hold
///     temporary grain species properties (only used in certain configurations)
///
/// > [!important]
/// > TODO: The role of the `dt` argument **MUST** be clarified! See the
/// > docstring of @ref grackle::impl::lookup_dust_rates1d for more details
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
    grackle::impl::ChemHeatingRates chemheatrates_buf,
    grackle::impl::InternalDustPropBuf internal_dust_prop_scratch_buf) {
  // Construct views of fields referenced in several parts of this function.

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

  // Look-up rate for H2 formation on dust & (when relevant) grain growth rates

  if (anydust != MASK_FALSE) {
    lookup_dust_rates1d(idx_range, dlogtem, tdust, dust2gas, h2dust, dom,
                        itmask_metal, dt, my_chemistry, my_rates, my_fields,
                        grain_growth_rates, grain_temperatures,
                        logTlininterp_buf, internal_dust_prop_scratch_buf);
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

  // model the effects of H2 self-shielding
  if (my_chemistry->primordial_chemistry > 1) {
    model_H2I_dissociation_shielding(kshield_buf, idx_range, tgas1d, mmw, dom,
                                     dx_cgs, c_ljeans, itmask, my_chemistry,
                                     my_fields, my_uvb_rates, internalu);
  }

  // apply some miscellaneous self-shielding adjustments
  if (my_chemistry->self_shielding_method > 0) {
    ShieldFactorCalculator calculator =
        setup_shield_factor_calculator(tgas1d, idx_range, dom, my_chemistry,
                                       my_fields, my_uvb_rates, internalu);
    apply_misc_shield_factors(kshield_buf, idx_range, itmask,
                              my_chemistry->self_shielding_method, my_uvb_rates,
                              &calculator);
  }

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
  // If using a high-energy radiation field, then account for
  //   effects of secondary electrons (Shull * Steenberg 1985)
  //   (see calc_rate.src)
  secondary_ionization_adjustments(idx_range, itmask, my_fields, my_uvb_rates,
                                   internalu, kshield_buf);
#endif
}

}  // namespace grackle::impl

#endif  // LOOKUP_COOL_RATES1D_HPP
