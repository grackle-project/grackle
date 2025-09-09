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

inline void simple_interp_lnT_rate(
  double* out, const double* table, const gr_mask_type* itmask, int i_start,
  int i_stop, grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf
) {
  // todo: we REALLY want to get rid of references to itmask in this loop,
  //       (for better performance). To do that, we need to ensure that
  //       logTlininterp_buf has reasonable values at every index!
  for (int i = i_start; i < i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      out[i] = table[logTlininterp_buf.indixe[i] - 1] +
          (table[logTlininterp_buf.indixe[i]] -
           table[logTlininterp_buf.indixe[i] - 1]) * logTlininterp_buf.tdef[i];
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
    double* k13dd_data_, double* h2dust, double dom, double dx_cgs,
    double c_ljeans, const gr_mask_type* itmask,
    const gr_mask_type* itmask_metal, double dt, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
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
  // TODO: we should create a struct that holds the k13dd_data_ scratch
  //       buffer and all other pre-allocated buffers that we reference
  //       here (to avoid heap allocations in this function)

  // collection of buffers for intermediate quantities used in dust-routines
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf =
      grackle::impl::new_InternalDustPropBuf(my_fields->grid_dimension[0],
                                             my_rates->gr_N[1]);
  // with some refactoring, we may be able to avoid allocating these buffers
  std::vector<double> f_shield_H(my_fields->grid_dimension[0]);
  std::vector<double> f_shield_He(my_fields->grid_dimension[0]);
  // TODO: these buffers should really be pre-initialized and tracked somewhere
  //       within my_rates (maybe in a gr_interp_grid_props instance). The fact
  //       that we reinitialize these every single time is EXTREMELY inefficient
  std::vector<double> d_Td(my_chemistry->NumberOfDustTemperatureBins);
  std::vector<double> d_Tg(my_chemistry->NumberOfTemperatureBins);

  // Construct some Views
  // --------------------
  // with some refactoring, we could entirely eliminate the k13dd_data_ buffer
  grackle::impl::View<double**> k13dd(k13dd_data_, my_fields->grid_dimension[0],
                                      14);

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

  // access the kcol_rate_tables from my_rates
  grackle::impl::CollisionalRxnRateCollection kcol_rate_tables =
      *(my_rates->opaque_storage->kcol_rate_tables);

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

  // Do linear table lookup (in log temperature)
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k1], kcol_rate_tables.data[CollisionalRxnLUT::k1], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k2], kcol_rate_tables.data[CollisionalRxnLUT::k2], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k3], kcol_rate_tables.data[CollisionalRxnLUT::k3], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k4], kcol_rate_tables.data[CollisionalRxnLUT::k4], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k5], kcol_rate_tables.data[CollisionalRxnLUT::k5], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k6], kcol_rate_tables.data[CollisionalRxnLUT::k6], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k57], kcol_rate_tables.data[CollisionalRxnLUT::k57], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k58], kcol_rate_tables.data[CollisionalRxnLUT::k58], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);

  // Look-up for 9-species model
  if (my_chemistry->primordial_chemistry > 1) {
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k7], kcol_rate_tables.data[CollisionalRxnLUT::k7], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k8], kcol_rate_tables.data[CollisionalRxnLUT::k8], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k9], kcol_rate_tables.data[CollisionalRxnLUT::k9], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k10], kcol_rate_tables.data[CollisionalRxnLUT::k10], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k11], kcol_rate_tables.data[CollisionalRxnLUT::k11], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k12], kcol_rate_tables.data[CollisionalRxnLUT::k12], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k13], kcol_rate_tables.data[CollisionalRxnLUT::k13], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k14], kcol_rate_tables.data[CollisionalRxnLUT::k14], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k15], kcol_rate_tables.data[CollisionalRxnLUT::k15], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k16], kcol_rate_tables.data[CollisionalRxnLUT::k16], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k17], kcol_rate_tables.data[CollisionalRxnLUT::k17], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k18], kcol_rate_tables.data[CollisionalRxnLUT::k18], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k19], kcol_rate_tables.data[CollisionalRxnLUT::k19], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k22], kcol_rate_tables.data[CollisionalRxnLUT::k22], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  }


  // Look-up for 12-species model
  if (my_chemistry->primordial_chemistry > 2) {
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k50], kcol_rate_tables.data[CollisionalRxnLUT::k50], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k51], kcol_rate_tables.data[CollisionalRxnLUT::k51], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k52], kcol_rate_tables.data[CollisionalRxnLUT::k52], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k53], kcol_rate_tables.data[CollisionalRxnLUT::k53], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k54], kcol_rate_tables.data[CollisionalRxnLUT::k54], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k55], kcol_rate_tables.data[CollisionalRxnLUT::k55], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k56], kcol_rate_tables.data[CollisionalRxnLUT::k56], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  }

  // Look-up for 15-species model

  if (my_chemistry->primordial_chemistry > 3) {
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k125], kcol_rate_tables.data[CollisionalRxnLUT::k125], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k129], kcol_rate_tables.data[CollisionalRxnLUT::k129], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k130], kcol_rate_tables.data[CollisionalRxnLUT::k130], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k131], kcol_rate_tables.data[CollisionalRxnLUT::k131], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k132], kcol_rate_tables.data[CollisionalRxnLUT::k132], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k133], kcol_rate_tables.data[CollisionalRxnLUT::k133], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k134], kcol_rate_tables.data[CollisionalRxnLUT::k134], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k135], kcol_rate_tables.data[CollisionalRxnLUT::k135], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k136], kcol_rate_tables.data[CollisionalRxnLUT::k136], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k137], kcol_rate_tables.data[CollisionalRxnLUT::k137], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k148], kcol_rate_tables.data[CollisionalRxnLUT::k148], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k149], kcol_rate_tables.data[CollisionalRxnLUT::k149], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k150], kcol_rate_tables.data[CollisionalRxnLUT::k150], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k151], kcol_rate_tables.data[CollisionalRxnLUT::k151], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k152], kcol_rate_tables.data[CollisionalRxnLUT::k152], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::k153], kcol_rate_tables.data[CollisionalRxnLUT::k153], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  }

  // Look-up for metal species model

  if (my_chemistry->metal_chemistry == 1) {
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz15], kcol_rate_tables.data[CollisionalRxnLUT::kz15], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz16], kcol_rate_tables.data[CollisionalRxnLUT::kz16], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz17], kcol_rate_tables.data[CollisionalRxnLUT::kz17], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz18], kcol_rate_tables.data[CollisionalRxnLUT::kz18], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz19], kcol_rate_tables.data[CollisionalRxnLUT::kz19], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz20], kcol_rate_tables.data[CollisionalRxnLUT::kz20], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz21], kcol_rate_tables.data[CollisionalRxnLUT::kz21], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz22], kcol_rate_tables.data[CollisionalRxnLUT::kz22], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz23], kcol_rate_tables.data[CollisionalRxnLUT::kz23], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz24], kcol_rate_tables.data[CollisionalRxnLUT::kz24], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz25], kcol_rate_tables.data[CollisionalRxnLUT::kz25], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz26], kcol_rate_tables.data[CollisionalRxnLUT::kz26], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz27], kcol_rate_tables.data[CollisionalRxnLUT::kz27], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz28], kcol_rate_tables.data[CollisionalRxnLUT::kz28], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz29], kcol_rate_tables.data[CollisionalRxnLUT::kz29], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz30], kcol_rate_tables.data[CollisionalRxnLUT::kz30], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz31], kcol_rate_tables.data[CollisionalRxnLUT::kz31], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz32], kcol_rate_tables.data[CollisionalRxnLUT::kz32], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz33], kcol_rate_tables.data[CollisionalRxnLUT::kz33], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz34], kcol_rate_tables.data[CollisionalRxnLUT::kz34], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz35], kcol_rate_tables.data[CollisionalRxnLUT::kz35], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz36], kcol_rate_tables.data[CollisionalRxnLUT::kz36], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz37], kcol_rate_tables.data[CollisionalRxnLUT::kz37], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz38], kcol_rate_tables.data[CollisionalRxnLUT::kz38], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz39], kcol_rate_tables.data[CollisionalRxnLUT::kz39], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz40], kcol_rate_tables.data[CollisionalRxnLUT::kz40], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz41], kcol_rate_tables.data[CollisionalRxnLUT::kz41], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz42], kcol_rate_tables.data[CollisionalRxnLUT::kz42], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz43], kcol_rate_tables.data[CollisionalRxnLUT::kz43], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz44], kcol_rate_tables.data[CollisionalRxnLUT::kz44], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz45], kcol_rate_tables.data[CollisionalRxnLUT::kz45], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz46], kcol_rate_tables.data[CollisionalRxnLUT::kz46], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz47], kcol_rate_tables.data[CollisionalRxnLUT::kz47], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz48], kcol_rate_tables.data[CollisionalRxnLUT::kz48], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz49], kcol_rate_tables.data[CollisionalRxnLUT::kz49], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz50], kcol_rate_tables.data[CollisionalRxnLUT::kz50], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz51], kcol_rate_tables.data[CollisionalRxnLUT::kz51], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz52], kcol_rate_tables.data[CollisionalRxnLUT::kz52], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz53], kcol_rate_tables.data[CollisionalRxnLUT::kz53], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
    simple_interp_lnT_rate(kcol_buf.data[CollisionalRxnLUT::kz54], kcol_rate_tables.data[CollisionalRxnLUT::kz54], itmask, idx_range.i_start, idx_range.i_stop, logTlininterp_buf);
  }

  // interpolate a few more rate tables
  if (my_chemistry->primordial_chemistry > 1) {
        // H2 formation heating terms.

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

    // construct the view of the k13 table
    grackle::impl::View<double**> k13dda(
        my_rates->k13dd, my_chemistry->NumberOfTemperatureBins, 14);

    for (int n1 = 0; n1 < 14; n1++) {
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask[i] != MASK_FALSE) {
          k13dd(i, n1) = k13dda(logTlininterp_buf.indixe[i] - 1, n1) +
                         (k13dda(logTlininterp_buf.indixe[i], n1) -
                          k13dda(logTlininterp_buf.indixe[i] - 1, n1)) *
                             logTlininterp_buf.tdef[i];
        }
      }
    }
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
      long long d_N[2] = {
          (long long)(my_chemistry->NumberOfDustTemperatureBins),
          (long long)(my_chemistry->NumberOfTemperatureBins)};
      long long d_Size = d_N[0] * d_N[1];
      // note: it is inefficient to repeatedly reinitialize d_Td & d_Tg
      for (int idx = 0; idx < my_chemistry->NumberOfDustTemperatureBins;
           idx++) {
        d_Td[idx] = logTdust_start + (double)idx * dlogTdust;
      }
      for (int idx = 0; idx < my_chemistry->NumberOfTemperatureBins; idx++) {
        d_Tg[idx] = logtem_start + (double)idx * dlogtem;
      }

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
                logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                dlogTdust, d_Tg.data(), dlogtem, d_Size,
                h2rate_silicate_coef_table);

            double h2AC = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                dlogTdust, d_Tg.data(), dlogtem, d_Size,
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
                logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                dlogTdust, d_Tg.data(), dlogtem, d_Size,
                h2rate_silicate_coef_table);

            logTdust =
                std::log(grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i]);
            double h2AC = f_wrap::interpolate_2d_g(
                logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                dlogTdust, d_Tg.data(), dlogtem, d_Size,
                h2rate_carbonaceous_coef_table);

            h2dust[i] = h2MgSiO3 * sgMgSiO3[i] + h2AC * sgAC[i];

            if (my_chemistry->dust_species > 1) {
              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i]);
              double h2SiM = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i]);
              double h2FeM = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i]);
              double h2Mg2SiO4 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i]);
              double h2Fe3O4 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i]);
              double h2SiO2D = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i]);
              double h2MgO = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i]);
              double h2FeS = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i]);
              double h2Al2O3 = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
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
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i]);
              double h2volorg = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
                  h2rate_silicate_coef_table);

              logTdust = std::log(
                  grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i]);
              double h2H2Oice = f_wrap::interpolate_2d_g(
                  logTdust, logTlininterp_buf.logtem[i], d_N, d_Td.data(),
                  dlogTdust, d_Tg.data(), dlogtem, d_Size,
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
              kd = f_wrap::interpolate_1d_g(logTlininterp_buf.logtem[i],
                                            nratec_single_elem_arr, d_Tg.data(),
                                            dlogtem, nratec_single_elem_arr[0],
                                            my_rates->grain_growth_rate);

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

  if (my_chemistry->self_shielding_method > 0) {
    // Compute shielding factors

    // conditionally construct views of some primordial species density fields
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

    grackle::impl::View<gr_float***> HM, DI, DII, HDI;
    if (my_chemistry->primordial_chemistry > 1) {
      HM = grackle::impl::View<gr_float***>(
          my_fields->HM_density, my_fields->grid_dimension[0],
          my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      if (my_chemistry->primordial_chemistry > 2) {
        DI = grackle::impl::View<gr_float***>(
            my_fields->DI_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        DII = grackle::impl::View<gr_float***>(
            my_fields->DII_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
        HDI = grackle::impl::View<gr_float***>(
            my_fields->HDI_density, my_fields->grid_dimension[0],
            my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
      }
    }

    // now do the actual calculation
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        // Compute shielding factor for H
        double nSSh = 6.73e-3 *
                      std::pow((my_uvb_rates.crsHI / 2.49e-18), (-2. / 3.)) *
                      std::pow((tgas1d[i] / 1.0e4), (0.17)) *
                      std::pow((my_uvb_rates.k24 / internalu.tbase1 / 1.0e-12),
                               (2.0 / 3.0));

        // Compute the total Hydrogen number density
        double nratio = (HI(i, idx_range.j, idx_range.k) +
                         HII(i, idx_range.j, idx_range.k));
        if (my_chemistry->primordial_chemistry > 1) {
          nratio = nratio + HM(i, idx_range.j, idx_range.k) +
                   H2I(i, idx_range.j, idx_range.k) +
                   H2II(i, idx_range.j, idx_range.k);

          if (my_chemistry->primordial_chemistry > 2) {
            nratio = nratio +
                     0.5 * (DI(i, idx_range.j, idx_range.k) +
                            DII(i, idx_range.j, idx_range.k)) +
                     2.0 * HDI(i, idx_range.j, idx_range.k) / 3.0;
          }
        }

        nratio = nratio * dom / nSSh;

        f_shield_H[i] =
            (0.98 * std::pow((1.0 + std::pow(nratio, (1.64))), (-2.28)) +
             0.02 * std::pow((1.0 + nratio), (-0.84)));

        // Compute shielding factor for He

        nSSh = 6.73e-3 *
               std::pow((my_uvb_rates.crsHeI / 2.49e-18), (-2. / 3.)) *
               std::pow((tgas1d[i] / 1.0e4), (0.17)) *
               std::pow((my_uvb_rates.k26 / internalu.tbase1 / 1.0e-12),
                        (2.0 / 3.0));

        nratio = 0.25 *
                 (HeI(i, idx_range.j, idx_range.k) +
                  HeII(i, idx_range.j, idx_range.k) +
                  HeIII(i, idx_range.j, idx_range.k)) *
                 dom / nSSh;

        f_shield_He[i] =
            (0.98 * std::pow((1.0 + std::pow(nratio, (1.64))), (-2.28)) +
             0.02 * std::pow((1.0 + nratio), (-0.84)));
      }
    }
  }

  if (my_chemistry->self_shielding_method == 1) {
    // approximate self shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // to shield HI, while leaving HeI and HeII optically thin

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling. (rate k29)
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * f_shield_H[i];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * f_shield_H[i];
        }

        kshield_buf.k25[i] = my_uvb_rates.k25;
        kshield_buf.k26[i] = my_uvb_rates.k26;
      }
    }

  } else if (my_chemistry->self_shielding_method == 2) {
    // Better self-shielding in HI using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // approximate self shielding in HeI and HeII

    //   Attenuate radiation rates for direct H2 ionization (15.4 eV)
    //   using same scaling as HI. (rate k29)

    //   Attenuate radiation rates for H2+ dissociation (30 eV)
    //   using same scaling as HeII. (rate k28 and k30)

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * f_shield_H[i];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * f_shield_H[i];
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i] = 0.;
        } else {
          kshield_buf.k26[i] = kshield_buf.k26[i] * f_shield_He[i];
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i] = 0.0;
        } else {
          kshield_buf.k28[i] = kshield_buf.k28[i] * f_shield_He[i];
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i] = 0.0;
        } else {
          kshield_buf.k30[i] = kshield_buf.k30[i] * f_shield_He[i];
        }

        kshield_buf.k25[i] = my_uvb_rates.k25;
      }
    }

  } else if (my_chemistry->self_shielding_method == 3) {
    // shielding using Eq. 13 and 14 from
    // Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
    // in HI and HeI, but ignoring HeII heating entirely
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        if (my_uvb_rates.k24 < tiny8) {
          kshield_buf.k24[i] = 0.;
        } else {
          kshield_buf.k24[i] = kshield_buf.k24[i] * f_shield_H[i];
        }

        // Scale H2 direct ionization radiation
        if (my_uvb_rates.k29 < tiny8) {
          kshield_buf.k29[i] = 0.;
        } else {
          kshield_buf.k29[i] = kshield_buf.k29[i] * f_shield_H[i];
        }

        // Apply same equations to HeI (assumes HeI closely follows HI)

        if (my_uvb_rates.k26 < tiny8) {
          kshield_buf.k26[i] = 0.;
        } else {
          kshield_buf.k26[i] = kshield_buf.k26[i] * f_shield_He[i];
        }

        // Scale H2+ dissociation radiation
        if (my_uvb_rates.k28 < tiny8) {
          kshield_buf.k28[i] = 0.0;
        } else {
          kshield_buf.k28[i] = kshield_buf.k28[i] * f_shield_He[i];
        }

        if (my_uvb_rates.k30 < tiny8) {
          kshield_buf.k30[i] = 0.0;
        } else {
          kshield_buf.k30[i] = kshield_buf.k30[i] * f_shield_He[i];
        }

        kshield_buf.k25[i] = 0.0;
      }
    }
  }

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
  // If using a high-energy radiation field, then account for
  //   effects of secondary electrons (Shull * Steenberg 1985)
  //   (see calc_rate.src)
  secondary_ionization_adjustments(idx_range, itmask, my_fields, my_uvb_rates,
                                   internalu, kshield_buf);
#endif

  //   If using H2, and using the density-dependent collisional
  //     H2 dissociation rate, then replace the the density-independant
  //        k13 rate with the new one.
  // May/00: there appears to be a problem with the density-dependent
  //     collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
  if (my_chemistry->primordial_chemistry > 1 &&
      my_chemistry->three_body_rate == 0) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        double nh = std::fmin(HI(i, idx_range.j, idx_range.k) * dom, 1.0e9);
        kcol_buf.data[CollisionalRxnLUT::k13][i] = tiny8;
        if (tgas1d[i] >= 500. && tgas1d[i] < 1.0e6) {
          // Direct collisional dissociation
          double k13_CID =
              k13dd(i, 0) -
              k13dd(i, 1) / (1. + std::pow((nh / k13dd(i, 4)), k13dd(i, 6))) +
              k13dd(i, 2) -
              k13dd(i, 3) / (1. + std::pow((nh / k13dd(i, 5)), k13dd(i, 6)));
          k13_CID = std::fmax(std::pow(10., k13_CID), tiny8);
          // Dissociative tunnelling
          double k13_DT =
              k13dd(i, 7) -
              k13dd(i, 8) / (1. + std::pow((nh / k13dd(i, 11)), k13dd(i, 13))) +
              k13dd(i, 9) -
              k13dd(i, 10) / (1. + std::pow((nh / k13dd(i, 12)), k13dd(i, 13)));
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

  drop_InternalDustPropBuf(&internal_dust_prop_buf);

  return;
}

}  // namespace grackle::impl

#endif  // LOOKUP_COOL_RATES1D_HPP
