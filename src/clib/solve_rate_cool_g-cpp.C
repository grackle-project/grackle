// See LICENSE file for license and copyright information

/// @file solve_rate_cool_g-cpp.C
/// @brief Declares signature of solve_rate_cool_g

// This file was initially generated automatically during conversion of the
// solve_rate_cool_g function from FORTRAN to C++

#include <cstdio>
#include <cstdlib> // std::malloc, std::free
#include <cstring> // std::memcpy
#include <vector>

#include "grackle.h"
#include "fortran_func_wrappers.hpp"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "step_rate_newton_raphson.hpp"
#include "utils-cpp.hpp"
#include "visitor/common.hpp"
#include "visitor/memory.hpp"

#include "ceiling_species.hpp"
#include "solve_rate_cool_g-cpp.h"

/// overrides the subcycle timestep (for each index in the index-range that is
/// selected by the given itmask) with the maximum allowed heating/cooling
/// timestep when the current value is larger.
///
/// @param[out]    dtit buffer tracking the current subcycle timestep for each
///    index in the index-range. If the current value exceeds the maximum
///    allowed heating/cooling timestep, the values will overwritten
/// @param[in]     idx_range Specifies the current index-range
/// @param[in]     dt tracks the full timestep that all the subcycles will
///    eventually add up to
/// @param[in]     ttot tracks the total time that has already elapsed from
///    previous subcycles for each location in `idx_range`
/// @param[in]     itmask Specifies the `idx_range`'s iteration-mask for this
///    calculation
/// @param[in]     tgas specifies the gas temperatures for the `idx_range`
/// @param[in]     p2d specifies the pressures for the `idx_range`. This is
///    computed user-specified nominal adiabatic index value (i.e. no attempts
///    are made to correct for presence of H2)
/// @param[in,out] edot specifies the time derivative of internal energy
///    density for each location in `idx_range`. This may be overwritten to
///    enforce the floor.
static void enforce_max_heatcool_subcycle_dt_(
  double* dtit, IndexRange idx_range, double dt, const double* ttot,
  const gr_mask_type* itmask, const double* tgas, const double* p2d,
  double* edot, const chemistry_data* my_chemistry
) {

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {

    if (itmask[i] != MASK_FALSE) {
      // Set energy per unit volume of this cell based in the pressure
      // (the gamma used here is the right one even for H2 since p2d
      //  is calculated with this gamma).
      double energy = std::fmax(p2d[i]/(my_chemistry->Gamma-1.), tiny8);

      // If the temperature is at the bottom of the temperature look-up
      // table and edot < 0, then shut off the cooling.
      if (tgas[i] <= 1.01*my_chemistry->TemperatureStart && edot[i] < 0.) {
        edot[i] = tiny8;
      }

      // enforce the floor
      if (std::fabs(edot[i]) < tiny8) { edot[i] = tiny8; }

      // Compute timestep for 10% change
      dtit[i] = grackle::impl::fmin(
        (double)(std::fabs(0.1 * energy / edot[i])), dt - ttot[i], dtit[i]
      );

      if (dtit[i] != dtit[i]) {
        OMP_PRAGMA_CRITICAL
        {
          eprintf("HUGE dtit ::  %g %g %g %g %g %g %g\n",
                  energy, edot[i], dtit[i], dt, ttot[i],
                  std::fabs(0.1 * energy / edot[i]),
                  (double)(std::fabs(0.1 * energy / edot [i])) );
        }
      }

    }
  }

}

// -------------------------------------------------------------

/// Set up the masks (and `imp_eng`) that identify the schemes that will be
/// used to evolve the chemistry network
///
/// There are 2 schemes:
///   1. Gauss-Seidel (for low-density zones)
///   2. Newton-Raphson (for high-density zones) in 2 modes:
///      - internal energy evolution is operator-split (handled separately from
///        chemistry network)
///      - internal energy is coupled with the rest of the chemistry network
///
/// @param[in]  idx_range Specifies the current index-range
/// @param[in]  itmask Specifies all locations to be evolved
///     during the current subcycle (in `idx_range`).
/// @param[out] itmask_gs Buffer for `idx_range` that is used to specify
///     locations where we will apply Gauss-Seidel scheme
/// @param[out] itmask_nr Buffer for `idx_range` that is used to specify
///     locations where we will apply Newton-Raphson scheme
/// @param[out] imp_eng Buffer for `idx_range` where the choice of
///     energy-evolution handling is recorded for the Newton-Raphson scheme
/// @param[in]  mask_len the length of the iteration masks
/// @param[in]  imetal specifies whether or not the caller provided a metal
///     density field
/// @param[in]  min_metallicity specifies the minimum metallicity where we
///     consider metal chemistry/cooling
/// @param[in]  ddom specifies precomputed product of mass density and the
///    `dom` quantity for each location in `idx_range`
/// @param[in]  tgas specifies the gas temperatures for the `idx_range`
/// @param[in]  metallicity specifies the metallicity for the `idx_range`
static void setup_chem_scheme_masks_(
  IndexRange idx_range, const gr_mask_type* itmask, gr_mask_type* itmask_gs,
  gr_mask_type* itmask_nr, int* imp_eng, int mask_len, int imetal,
  double min_metallicity, const double* ddom, const double* tgas,
  const double* metallicity, const chemistry_data* my_chemistry
) {

  std::memcpy(itmask_gs, itmask, sizeof(gr_mask_type)*mask_len);
  std::memcpy(itmask_nr, itmask, sizeof(gr_mask_type)*mask_len);

  // would it be more robust to use my_chemistry->metal_cooling than imetal?

  // netwon-raphson solver can only coevolves energy when this is true
  const bool has_nr_coevolve_eint_prereqs = (
    (my_chemistry->with_radiative_cooling == 1) &&
    (my_chemistry->primordial_chemistry > 1)
  );

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if ( itmask[i] != MASK_FALSE )  {
      bool usemetal = (imetal == 1) && (metallicity[i] > min_metallicity);
      bool is_hi_dens = (ddom[i] >= 1.e8) || (usemetal && (ddom[i] >= 1.0e6));

      itmask_gs[i] = (!is_hi_dens) ? MASK_TRUE : MASK_FALSE;
      itmask_nr[i] = (is_hi_dens) ? MASK_TRUE : MASK_FALSE;

      // when true, the newton-raphson scheme coevolves internal-energy
      // alongside chemistry
      bool nr_coevolve_energy = (
        has_nr_coevolve_eint_prereqs && (ddom[i] > 1.e7) && (tgas[i] > 1650.0)
      );
      // we don't care about imp_eng[i] where itmask_nr[i] == MASK_FALSE
      imp_eng[i] = (nr_coevolve_energy) ? 1 : 0;
    }
  }

}

// -------------------------------------------------------------

/// Computes the timescale given by `ndens_Heq / (d ndens_Heq / d t)`
///
/// This is used to help compute the subcycle timestep when using a primordial
/// chemistry solver
///
/// @param my_chemistry holds a number of configuration parameters
/// @param my_rates holds assorted rate data. In this function, this is being
///    used to specify some the interpolation tables of some relevant reaction
///    rates (they are tabulated with respect to logT)
/// @param dlogtem Specifies the constant spacing shared by the relevant rate
///    interpolation tables
/// @param logTlininterp_buf Specifies the information related to the position
///    in the logT interpolations (for a number of chemistry zones)
/// @param k13, k22 1D arrays specifying the previously looked up, local values
///    of the k13 and k22 rates.
/// @param local_rho specifies the local (total) mass density
/// @param tgas 1D array specifying the temperature
/// @param p2d 1D array specifying the pressures values. This is computed from
///    the user-specified nominal adiabatic index value (i.e. no attempts
///    are made to correct for presence of H2)
/// @param edot 1D array specifying the time derivative of the internal energy
///    density
/// @param i Specifies the index of the relevant zone in the 1D array. (**BE
///    AWARE:** this is a 0-based index)
///
/// @note
/// The `static` annotation indicates that this function is only visible to the
/// current translation unit
static double calc_Heq_div_dHeqdt_(
  const chemistry_data* my_chemistry,
  const chemistry_data_storage* my_rates,
  double dlogtem,
  const grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  const double* k13,
  const double* k22,
  double local_rho,
  const double* tgas,
  const double* p2d,
  const double* edot,
  int i
) {

  // Equilibrium value for H is:
  // Heq = (-1._DKIND / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
  // We want to know dH_eq/dt.
  // - We can trivially get dH_eq/dT.
  // - We have de/dt.
  // - We need dT/de.
  //
  // T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
  // dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
  // dH_eq / dT (see above; we can calculate the derivative here)
  // dT / de = utem * (gamma - 1._DKIND) / N == tgas / p2d
  // de / dt = edot
  // Now we use our estimate of dT/de to get the estimated
  // difference in the equilibrium
  double eqt2 = std::fmin(std::log(tgas[i]) + 0.1*dlogtem, logTlininterp_buf.t2[i]);
  double eqtdef = (eqt2 - logTlininterp_buf.t1[i])/(logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
  double eqk222 = my_rates->k22[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k22[logTlininterp_buf.indixe[i]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double eqk132 = my_rates->k13[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k13[logTlininterp_buf.indixe[i]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double heq2 = (-1. / (4.*eqk222)) * (eqk132-
    std::sqrt(8.*eqk132*eqk222*
              my_chemistry->HydrogenFractionByMass*local_rho+
              std::pow(eqk132,2.)));

  double eqt1 = std::fmax(std::log(tgas[i]) - 0.1*dlogtem, logTlininterp_buf.t1[i]);
  eqtdef = (eqt1 - logTlininterp_buf.t1[i])/(logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
  double eqk221 = my_rates->k22[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k22[logTlininterp_buf.indixe[i]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double eqk131 = my_rates->k13[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k13[logTlininterp_buf.indixe[i]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double heq1 = (-1. / (4.*eqk221)) * (eqk131-
    std::sqrt(8.*eqk131*eqk221*
              my_chemistry->HydrogenFractionByMass*local_rho+std::pow(eqk131,2.)));

  double dheq = (std::fabs(heq2-heq1)/(std::exp(eqt2) - std::exp(eqt1)))
    * (tgas[i]/p2d[i]) * edot[i];
  double heq = (-1. / (4.*k22[i])) * (k13[i]-
    std::sqrt(8.*k13[i]*k22[i]*
              my_chemistry->HydrogenFractionByMass*local_rho+std::pow(k13[i],2.)));

  return heq / dheq;
}

// -------------------------------------------------------------

/// Sets the current subcycle timestep for each index in the index-range
/// if it exceeds maximum the allowed chemistry-rate timestep.
///
/// @param[out] dtit buffer tracking the current subcycle timestep for each
///    index in the index-range. Values will be modified in place.
/// @param[in] idx_range Specifies the current index-range
/// @param[in] iter current subcycle iteration
/// @param[in] dt tracks the full timestep that all the subcycles will
///    eventually add up to
/// @param[in] ttot tracks the total time that has already elapsed from
///    previous subcycles for each location in `idx_range`
/// @param[in] itmask_gs Specifies the `idx_range`'s iteration-mask for the
///    Gauss-Seidel scheme
/// @param[in] itmask_nr Specifies the `idx_range`'s iteration-mask for the
///    Newton-Raphson scheme
/// @param[in] imp_eng Specifies how Newton-Raphson scheme handles energy
///    evolution at each `idx_range` location
/// @param[in] dedot, HIdot respectively specify the time derivative of the
///    free electrons and HI for the `idx_range`
/// @param[in] dedot_prev, HIdot_prev respectively specify the time derivative
///    of the free electron density and HI density for the `idx_range` from the
///    previous subcycle (they're allowed to hold garbage data in 1st subcycle)
/// @param[in] ddom specifies precomputed product of mass density and the
///    `dom` quantity for each location in `idx_range`
/// @param[in] tgas specifies the gas temperatures for the `idx_range`
/// @param[in] p2d specifies the pressures for the `idx_range`. This is
///    computed user-specified nominal adiabatic index value (i.e. no attempts
///    are made to correct for presence of H2)
/// @param[in] edot specifies the time derivative of the internal energy
///    density for the `idx_range`.
/// @param[in] my_chemistry holds a number of configuration parameters
/// @param[in] my_rates holds assorted rate data. In this function, this is
///    being used to specify the interpolation tables of some relevant reaction
///    rates (they are tabulated with respect to logT)
/// @param[in] dlogtem Specifies the constant spacing shared by the relevant
///    rate interpolation tables
/// @param[in] logTlininterp_buf Specifies the information related to the
///    position in the logT interpolations (for a number of chemistry zones)
/// @param[in] my_fields specifies the field data
/// @param[in] kcr_buf holds various pre-computed chemical reaction rates for
///    each location in `idx_range`.
///
/// @todo
/// Consider breaking this into 2 functions that separately determines dtit for
/// Gauss-Seidel and Newton-Raphson. (At the time of writing, the included
/// logic for Newton-Raphson doesn't care about the chemistry-rates, instead
/// it sets the timestep based on the energy evolution)
static void set_subcycle_dt_from_chemistry_scheme_(
  double* dtit, IndexRange idx_range, int iter, double dt, const double* ttot,
  const gr_mask_type* itmask_gs, const gr_mask_type* itmask_nr,
  const int* imp_eng, double* dedot, double* HIdot,
  const double* dedot_prev, const double* HIdot_prev,
  const double* ddom, const double* tgas, const double* p2d, const double* edot,
  const chemistry_data* my_chemistry, const chemistry_data_storage* my_rates,
  double dlogtem,
  const grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle_field_data* my_fields,
  grackle::impl::CollisionalRxnRateCollection kcr_buf
) {
  const int j = idx_range.j;
  const int k = idx_range.k;

  grackle::impl::View<gr_float***> de(my_fields->e_density,
                                      my_fields->grid_dimension[0],
                                      my_fields->grid_dimension[1],
                                      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density,
                                      my_fields->grid_dimension[0],
                                      my_fields->grid_dimension[1],
                                      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density,
                                       my_fields->grid_dimension[0],
                                       my_fields->grid_dimension[1],
                                       my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> d(my_fields->density,
                                     my_fields->grid_dimension[0],
                                     my_fields->grid_dimension[1],
                                     my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy,
                                     my_fields->grid_dimension[0],
                                     my_fields->grid_dimension[1],
                                     my_fields->grid_dimension[2]);


  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask_gs[i] != MASK_FALSE) {
      // in this case, the chemical network will be evolved with Gauss-Seidel

      // Part 1 of 2: adjust values of dedot and HIdot
      // ---------------------------------------------

      // Bound from below to prevent numerical errors
      if (std::fabs(dedot[i]) < tiny8) {
        dedot[i] = std::fmin(tiny_fortran_val, de(i,j,k));
      }
      if (std::fabs(HIdot[i]) < tiny8){
        HIdot[i] = std::fmin(tiny_fortran_val, HI(i,j,k));
      }

      // If the net rate is almost perfectly balanced then set
      //     it to zero (since it is zero to available precision)
      {
        double ion_rate = std::fabs(kcr_buf.data[CollisionalRxnLUT::k1][i] *
                                    de(i,j,k) * HI(i,j,k));
        double recomb_rate = std::fabs(kcr_buf.data[CollisionalRxnLUT::k2][i] *
                                       HII(i,j,k) * de(i,j,k));
        double ratio = (std::fmin(ion_rate, recomb_rate) /
                        std::fmax(std::fabs(dedot[i]), std::fabs(HIdot[i])));
        if (ratio > 1.0e6) {
          dedot[i] = tiny8;
          HIdot[i] = tiny8;
        }
      }

      // If the iteration count is high then take the smaller of
      //   the calculated dedot and last time step's actual dedot.
      //   This is intended to get around the problem of a low
      //   electron or HI fraction which is in equilibrium with high
      //   individual terms (which all nearly cancel).
      if (iter > 50)  {
        dedot[i] = std::fmin(std::fabs(dedot[i]), std::fabs(dedot_prev[i]));
        HIdot[i] = std::fmin(std::fabs(HIdot[i]), std::fabs(HIdot_prev[i]));
      }

      // Part 2 of 2: compute minimum rate timestep
      // ------------------------------------------

      double olddtit = dtit[i];
      dtit[i] = grackle::impl::fmin(std::fabs(0.1*de(i,j,k)/dedot[i]),
                                    std::fabs(0.1*HI(i,j,k)/HIdot[i]),
                                    dt-ttot[i],
                                    0.5*dt);

      if (ddom[i] > 1.e8  && edot[i] > 0. &&
          my_chemistry->primordial_chemistry > 1)  {
        // here, we ensure that that the equilibrium mass density of
        // Hydrogen changes by 10% or less
        double Heq_div_dHeqdt = calc_Heq_div_dHeqdt_(
          my_chemistry, my_rates, dlogtem, logTlininterp_buf,
          kcr_buf.data[CollisionalRxnLUT::k13],
          kcr_buf.data[CollisionalRxnLUT::k22],
          d(i,j,k), tgas, p2d, edot, i
        );

        dtit[i] = std::fmin(dtit[i], 0.1*Heq_div_dHeqdt);
      }

      if (iter > 10) {
        dtit[i] = std::fmin(olddtit*1.5, dtit[i]);
      }

    } else if ((itmask_nr[i]!=MASK_FALSE) && (imp_eng[i]==0))  {
      // we may want to handle this case and the next case in a separate
      // function (they determine the timestep using very different logic than
      // in the above case)
      dtit[i] = grackle::impl::fmin(std::fabs(0.1*e(i,j,k)/edot[i]*d(i,j,k)),
                                    dt-ttot[i],
                                    0.5*dt);

    } else if ((itmask_nr[i]!=MASK_FALSE) && (imp_eng[i]==1))  {
      dtit[i] = dt - ttot[i];

    } else {
      dtit[i] = dt;
    }
  }
}

// -------------------------------------------------------------

/// Updates the iteration mask in the case where the user has specified that we
/// are using grackle as a part of a coupled radiative transfer calculation
///
/// @param[out] itmask the mask that will be overriden
/// @param[in] idx_range specifies the index-range
/// @param[in] my_chemistry specifies grackle settings (we probably don't need
///     to pass in everything)
/// @param[in] my_fields used to access HI photo-ionization rate field

static inline void coupled_rt_modify_itmask_(
  gr_mask_type* itmask,
  IndexRange idx_range,
  const chemistry_data* my_chemistry,
  grackle_field_data* my_fields
)
{
  grackle::impl::View<const gr_float***> kphHI(my_fields->RT_HI_ionization_rate,
                                               my_fields->grid_dimension[0],
                                               my_fields->grid_dimension[1],
                                               my_fields->grid_dimension[2]);

  // adjust iteration mask if the caller indicates that they're using
  // Grackle in a coupled radiative-transfer/chemistry-energy calculation
  // (that has intermediate steps)
  if (my_chemistry->use_radiative_transfer == 1 &&
      my_chemistry->radiative_transfer_coupled_rate_solver == 1)  {
    // we only define behavior for radiative_transfer_intermediate_step
    // values of 0 or 1

    const int j = idx_range.j;
    const int k = idx_range.k;

    if (my_chemistry->radiative_transfer_intermediate_step == 1) {
      // the caller has invoked this chemistry-energy solver as an
      // intermediate step of a coupled radiative-transfer/chemistry-energy
      // calculation and they only want the solver consider cells where
      // the radiation is non-zero
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = (kphHI(i,j,k) > 0) ? MASK_TRUE : MASK_FALSE;
      }
    } else if (my_chemistry->radiative_transfer_intermediate_step == 0) {
      // the caller has invoked this chemistry-energy solver outside
      // of their coupled radiative-transfer/chemistry-energy calculation.
      // They want to apply the solver to cells where radiation is 0 (i.e.
      // locations where skipped by the coupled calculation)
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = (kphHI(i,j,k) > 0) ? MASK_FALSE : MASK_TRUE;
      }
    }
  }
}

// -------------------------------------------------------------

namespace grackle::impl {

/// Aggregates buffers used as scratch space in rate-related calculations
///
/// This exists to encapsulate the logic for all of the local buffers used in
/// rate calculations within solve_rate_cool_g (this is useful given the size
/// of the function). It is not expected to be used outside of the function.
///
/// @note
/// The purpose of this type is similar in spirit to the purpose of the
/// ScratchBuf data structrues described in the inernal_types.hpp header and
/// we consequently observe those conventions relating to (con|de)structors.
///
/// @note
/// The majority of the time, functions only need a subset of data stored by
/// this struct. In these cases, you should strongly prefer to only pass the
/// members that are required as function arguments (if we unnecessarily pass
/// this whole struct, that can make it difficult to visualize the data flow)
struct SpeciesRateSolverScratchBuf {

  /// specifies precomputed product of mass density and the `dom` quantity for
  /// each location in the index-range. Used to pick the scheme for solving the
  /// rate equations and in setting the max allowed chemistry-timstep
  double* ddom;

  /// buffers to hold time derivatives of free electron and HI mass densities
  /// for index_range
  double *dedot, *HIdot;
  /// buffers to hold time derivatives of free electron and HI mass densities
  /// for index_range computed during the previous cycle
  double *dedot_prev, *HIdot_prev;

  /// buffer used to track the rate of H2 formation on dust grains
  double* h2dust;

  /// scratch space used only within lookup_cool_rates1d_g. This is 14 times
  /// larger than most of the other buffers.
  ///
  /// (with minimal refactoring, this buffer could probably be removed)
  double *k13dd;

  /// iteration mask denoting where the Gauss-Seidel scheme will be used
  gr_mask_type* itmask_gs;

  /// iteration mask denoting where the Newton-Raphson scheme will be used
  gr_mask_type* itmask_nr;

  /// buffer specifying how the Newton-Raphson scheme handles energy evolution
  int* imp_eng;

  // buffers in the following data structure are used to temporarily hold
  // the evolved density of various species as we evolve over a subcycle
  // (currently only used by step_rate_g)
  grackle::impl::SpeciesCollection species_tmpdens;

  // buffers in the following data structure are used to temporarily hold
  // the interpolated Collisional Rxn Rates that have been
  // interpolated using the standard 1D log temperature table.
  grackle::impl::CollisionalRxnRateCollection kcr_buf;

  // buffers in the following data structure are used to temporarily hold
  // the computed radiative reaction rates
  grackle::impl::PhotoRxnRateCollection kshield_buf;

  // buffers in the following data structure are used to temporarily hold
  // the interpolated chemistry-heating rates at each index-range location
  grackle::impl::ChemHeatingRates chemheatrates_buf;

  // holds computed grain growth/destruction rates:
  grackle::impl::GrainSpeciesCollection grain_growth_rates;

};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template <class BinaryVisitor>
void visit_member_pair(SpeciesRateSolverScratchBuf& obj0,
                       SpeciesRateSolverScratchBuf& obj1,
                       BinaryVisitor f) {
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("SpeciesRateSolverScratchBuf", f);
  f(VIS_MEMBER_NAME("ddom"), obj0.ddom, obj1.ddom, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("dedot"), obj0.dedot, obj1.dedot, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("HIdot"), obj0.HIdot, obj1.HIdot, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("dedot_prev"), obj0.dedot_prev, obj1.dedot_prev, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("HIdot_prev"), obj0.HIdot_prev, obj1.HIdot_prev, vis::idx_range_len_multiple(1));
  // the next line is NOT a typo
  f(VIS_MEMBER_NAME("k13dd"), obj0.k13dd, obj1.k13dd, vis::idx_range_len_multiple(14));
  f(VIS_MEMBER_NAME("h2dust"), obj0.h2dust, obj1.h2dust, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("itmask_gs"), obj0.itmask_gs, obj1.itmask_gs, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("itmask_nr"), obj0.itmask_nr, obj1.itmask_nr, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("imp_eng"), obj0.imp_eng, obj1.imp_eng, vis::idx_range_len_multiple(1));

  vis::previsit_struct_member(VIS_MEMBER_NAME("species_tmpdens"), f);
  visit_member_pair(obj0.species_tmpdens, obj1.species_tmpdens, f);

  vis::previsit_struct_member(VIS_MEMBER_NAME("kcr_buf"), f);
  visit_member_pair(obj0.kcr_buf, obj1.kcr_buf, f);

  vis::previsit_struct_member(VIS_MEMBER_NAME("kshield_buf"), f);
  visit_member_pair(obj0.kshield_buf, obj1.kshield_buf, f);

  vis::previsit_struct_member(VIS_MEMBER_NAME("chemheatrates_buf"), f);
  visit_member_pair(obj0.chemheatrates_buf, obj1.chemheatrates_buf, f);

  vis::previsit_struct_member(VIS_MEMBER_NAME("grain_growth_rates"), f);
  visit_member_pair(obj0.grain_growth_rates, obj1.grain_growth_rates, f);

  vis::end_visit(f);
}

template <class UnaryFn>
void visit_member(SpeciesRateSolverScratchBuf* obj, UnaryFn fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, SpeciesRateSolverScratchBuf, obj, fn)
}

/// allocates the contents of a new SpeciesRateSolverScratchBuf
///
/// @param nelem The number of elements a buffer is expected to have in order
///    to store values for the standard sized index-range
SpeciesRateSolverScratchBuf new_SpeciesRateSolverScratchBuf(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  SpeciesRateSolverScratchBuf out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out, grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

/// performs cleanup of the contents of SpeciesRateSolverScratchBuf
///
/// This effectively invokes a destructor
void drop_SpeciesRateSolverScratchBuf(SpeciesRateSolverScratchBuf* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}


} // namespace grackle::impl

// -------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int solve_rate_cool_g(
  int imetal, double dt, InternalGrUnits internalu,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage* my_uvb_rates
)
{
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

#ifdef GRACKLE_FLOAT_4
  const gr_float tolerance = (gr_float)(1.0e-05);
#else
  const gr_float tolerance = (gr_float)(1.0e-10);
#endif

  // Set error indicator (we will return this value)
  int ierr = GR_SUCCESS;

  // Set flag for dust-related options
  const gr_mask_type anydust =
    ((my_chemistry->h2_on_dust > 0)  ||  (my_chemistry->dust_chemistry > 0))
    ? MASK_TRUE
    : MASK_FALSE;

  // ignore metal chemistry/cooling below this metallicity
  const double min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;
      
  // Set units
  const double dom      = internalu_calc_dom_(internalu);
  const double chunit   = internalu_get_chunit_(internalu);

  const double dx_cgs = my_fields->grid_dx * internalu.xbase1;
  const double c_ljeans = internalu_calc_coef_ljeans_(internalu,
                                                      my_chemistry->Gamma);

  const double dlogtem = (
    (std::log(my_chemistry->TemperatureEnd) -
     std::log(my_chemistry->TemperatureStart)) /
    (double)(my_chemistry->NumberOfTemperatureBins-1 )
  );

  // Convert densities from comoving to proper

  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,(-3)) );
    f_wrap::scale_fields_g(imetal, factor, my_chemistry, my_fields);
  }

  grackle::impl::ceiling_species_g(imetal, my_chemistry, my_fields);

  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  OMP_PRAGMA("omp parallel")
  {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

    // holds computed grain temperatures:
    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf =
      grackle::impl::new_Cool1DMultiScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::CoolHeatScratchBuf coolingheating_buf =
      grackle::impl::new_CoolHeatScratchBuf(my_fields->grid_dimension[0]);

    // holds buffers exclusively used for solving species rate equations
    // (i.e. in the future, we could have the constructor skip allocations of
    // all contained data structures when using primordial_chemistry == 0)
    grackle::impl::SpeciesRateSolverScratchBuf spsolvbuf =
      grackle::impl::new_SpeciesRateSolverScratchBuf(
        my_fields->grid_dimension[0]
      );

    // the following variables aren't embedded in structs because they are used
    // in a number of different internal routines. Sorting these into
    // additional structs (or leaving them free-standing) will become more
    // obvious as we transcribe more routines.
    std::vector<double> dtit(my_fields->grid_dimension[0]);
    std::vector<double> ttot(my_fields->grid_dimension[0]);
    std::vector<double> p2d(my_fields->grid_dimension[0]);
    std::vector<double> tgas(my_fields->grid_dimension[0]);
    std::vector<double> tdust(my_fields->grid_dimension[0]);
    std::vector<double> metallicity(my_fields->grid_dimension[0]);
    std::vector<double> dust2gas(my_fields->grid_dimension[0]);
    std::vector<double> rhoH(my_fields->grid_dimension[0]);
    std::vector<double> mmw(my_fields->grid_dimension[0]);
    std::vector<double> edot(my_fields->grid_dimension[0]);

    // iteration masks
    std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
    std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);

    // create views of density and internal energy fields to support 3D access
    grackle::impl::View<gr_float***> d(my_fields->density,
                                       my_fields->grid_dimension[0],
                                       my_fields->grid_dimension[1],
                                       my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> e(my_fields->internal_energy,
                                       my_fields->grid_dimension[0],
                                       my_fields->grid_dimension[1],
                                       my_fields->grid_dimension[2]);

    // The following for-loop is a flattened loop over every k,j combination.
    // OpenMP divides this loop between all threads. Within the loop, we
    // complete calculations for the constructed index-range construct
    // (an index range corresponds to an "i-slice")
    OMP_PRAGMA("omp for schedule(runtime)")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);
      const int k = idx_range.k; // use 0-based index
      const int j = idx_range.j; // use 0-based index

      // `tolerance = 1.0e-06_DKIND * dt` was some commented logic in the
      // original fortran subroutine in this location

      // Initialize iteration mask to true for all cells.
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;
      }

      // adjust iteration mask (but only if using Grackle in a coupled
      // radiative-transfer calculation)
      coupled_rt_modify_itmask_(itmask.data(), idx_range, my_chemistry,
                                my_fields);

      // Set time elapsed to zero for each cell in 1D section

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        ttot[i] = 0.;
      }

      // A useful slice variable since we do this a lot
      // -> we don't need it for primordial_chemistry==0
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        spsolvbuf.ddom[i] = d(i,j,k) * dom;
      }

      // declare 2 variables (primarily used for subcycling, but also used in
      // error reporting)
      int iter;
      double ttmin;

      // ------------------ Loop over subcycles ----------------

      for (iter = 1; iter<=(my_chemistry->max_iterations); iter++) {

        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (itmask[i] != MASK_FALSE)  {
            dtit[i] = huge8;
          }
        }

        // Compute the cooling rate, tgas, tdust, and metallicity for this row
        f_wrap::cool1d_multi_g(
          imetal, idx_range, iter, edot.data(), tgas.data(),
          mmw.data(), p2d.data(), tdust.data(), metallicity.data(),
          dust2gas.data(), rhoH.data(), itmask.data(), itmask_metal.data(),
          my_chemistry, my_rates, my_fields, *my_uvb_rates, internalu,
          grain_temperatures, logTlininterp_buf, cool1dmulti_buf,
          coolingheating_buf
        );

        if (my_chemistry->primordial_chemistry > 0)  {

          // Look-up rates as a function of temperature for 1D set of zones
          //  (maybe should add itmask to this call)
          //
          // -> TODO: passing dt to this function is probably incorrect. See
          //    the C++ docstring for a longer discussion
          f_wrap::lookup_cool_rates1d_g(
            idx_range, anydust, tgas.data(), mmw.data(), tdust.data(),
            dust2gas.data(), spsolvbuf.k13dd, spsolvbuf.h2dust,
            dom, dx_cgs, c_ljeans, itmask.data(), itmask_metal.data(),
            imetal, rhoH.data(), dt, my_chemistry, my_rates, my_fields,
            *my_uvb_rates, internalu, spsolvbuf.grain_growth_rates,
            grain_temperatures, logTlininterp_buf, spsolvbuf.kcr_buf,
            spsolvbuf.kshield_buf, spsolvbuf.chemheatrates_buf
          );

          // Compute dedot and HIdot, the rates of change of de and HI
          //   (should add itmask to this call)

          f_wrap::rate_timestep_g(
            spsolvbuf.dedot, spsolvbuf.HIdot, anydust, idx_range,
            spsolvbuf.h2dust, rhoH.data(), itmask.data(), edot.data(),
            chunit, dom, my_chemistry, my_fields, *my_uvb_rates,
            spsolvbuf.kcr_buf, spsolvbuf.kshield_buf,
            spsolvbuf.chemheatrates_buf
          );

          // Setup masks to identify which chemistry schemes to use. We split
          // cells by density:
          //    => low-density: Gauss-Seidel scheme, tracked by itmask_gs
          //    => high-density: Newton-Raphson scheme, tracked by itmask_nr
          setup_chem_scheme_masks_(
            idx_range, itmask.data(), spsolvbuf.itmask_gs, spsolvbuf.itmask_nr,
            spsolvbuf.imp_eng, my_fields->grid_dimension[0], imetal,
            min_metallicity, spsolvbuf.ddom, tgas.data(), metallicity.data(),
            my_chemistry
          );

          // Set the max timestep for the current subcycle based on our scheme
          // for updating the chemical network:
          // - for Gauss-Seidel, pick a timestep that keeps relative chemical
          //   changes below 10%
          // - do something else for Newton-Raphson

          set_subcycle_dt_from_chemistry_scheme_(
            dtit.data(), idx_range, iter, dt, ttot.data(), spsolvbuf.itmask_gs,
            spsolvbuf.itmask_nr, spsolvbuf.imp_eng,
            spsolvbuf.dedot, spsolvbuf.HIdot,
            spsolvbuf.dedot_prev, spsolvbuf.HIdot_prev,
            spsolvbuf.ddom, tgas.data(), p2d.data(), edot.data(),
            my_chemistry, my_rates, dlogtem, logTlininterp_buf, my_fields,
            spsolvbuf.kcr_buf
          );
        }

        const gr_mask_type* energy_itmask =
          (my_chemistry->primordial_chemistry == 0)
          ? itmask.data() : spsolvbuf.itmask_gs;

        // Update dtit (the current subcycle timestep) to ensure it doesn't
        // exceed the max timestep for cooling/heating
        // -> zones that will use Newton-Raphson scheme are ignored
        enforce_max_heatcool_subcycle_dt_(
          dtit.data(), idx_range, dt, ttot.data(), energy_itmask, tgas.data(),
          p2d.data(), edot.data(), my_chemistry
        );

        // Update total and gas energy
        // -> zones that will use Newton-Raphson scheme are ignored
        if (my_chemistry->with_radiative_cooling == 1)  {
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            if (energy_itmask[i] != MASK_FALSE) {
              e(i,j,k) = e(i,j,k) + (gr_float)(edot[i]/d(i,j,k)*dtit[i]);
            }
          }
        }

        if (my_chemistry->primordial_chemistry > 0)  {

          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask_gs)
          f_wrap::step_rate_g(
            dtit.data(), idx_range, anydust, spsolvbuf.h2dust, rhoH.data(),
            spsolvbuf.dedot_prev, spsolvbuf.HIdot_prev, spsolvbuf.itmask_gs,
            itmask_metal.data(), imetal, my_chemistry, my_fields,
            *my_uvb_rates, spsolvbuf.grain_growth_rates,
            spsolvbuf.species_tmpdens, spsolvbuf.kcr_buf, spsolvbuf.kshield_buf
          );

          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask_nr)
          grackle::impl::step_rate_newton_raphson(
            imetal, idx_range, iter, dom, chunit, dx_cgs, c_ljeans,
            dtit.data(), p2d.data(), tgas.data(), tdust.data(),
            metallicity.data(), dust2gas.data(), rhoH.data(), mmw.data(),
            spsolvbuf.h2dust, edot.data(), anydust, spsolvbuf.itmask_nr,
            itmask_metal.data(), spsolvbuf.imp_eng, my_chemistry, my_rates,
            my_fields, *my_uvb_rates, internalu, grain_temperatures,
            logTlininterp_buf, cool1dmulti_buf, coolingheating_buf,
            spsolvbuf.chemheatrates_buf
          );

        }

        // Add the timestep to the elapsed time for each cell and find
        //  minimum elapsed time step in this row
        ttmin = huge8;
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          ttot[i] = std::fmin(ttot[i] + dtit[i], dt);

          if (std::fabs(dt-ttot[i]) < tolerance*dt) { itmask[i] = MASK_FALSE; }

          if (ttot[i]<ttmin) { ttmin = ttot[i]; }
        }

        // If all cells are done (in idx_range), break out of subcycle loop
        if (std::fabs(dt-ttmin) < tolerance*dt) { break; }

      }  // subcycle iteration loop (for current idx_range)

      // review number of iterations that were spent in the subcycle loop

      if (iter > my_chemistry->max_iterations)  {
        OMP_PRAGMA_CRITICAL
        {
          printf("inside if statement solve rate cool: %d %d\n",
                 my_fields->grid_start[0],
                 my_fields->grid_end[0]);
          eprintf("MULTI_COOL iter >  %d  at j_0based,k_0based = %d %d\n",
                  my_chemistry->max_iterations, idx_range.j, idx_range.k);
          printf("FATAL error (2) in MULTI_COOL\n");
          printf("( dt = %.17e ttmin = %.17e )", dt, ttmin);
          grackle::impl::print_contiguous_row_(
            dtit.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            ttot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            edot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            itmask.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );

          if (my_chemistry->exit_after_iterations_exceeded == 1)  {
            ierr = GR_FAIL;
          }
        }  // OMP_PRAGMA_CRITICAL
      }

      if (iter > my_chemistry->max_iterations/2) { // WARNING_MESSAGE
        OMP_PRAGMA_CRITICAL
        {
          eprintf("MULTI_COOL iter,j_0based,k_0based = %d %d %d\n",
                  iter, idx_range.j, idx_range.k);
        }
      }

    }  // outer-loop (index t) - each of these correspond to j,k pairs

    // cleanup manually allocated temporaries
    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_Cool1DMultiScratchBuf(&cool1dmulti_buf);
    grackle::impl::drop_CoolHeatScratchBuf(&coolingheating_buf);

    grackle::impl::drop_SpeciesRateSolverScratchBuf(&spsolvbuf);

  }  // OMP_PRAGMA("omp parallel")

  // If an error has been produced, return now.

  if (ierr != GR_SUCCESS)  {
    return ierr;
  }

  // Convert densities back to comoving from proper

  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,3) );
    f_wrap::scale_fields_g(imetal, factor, my_chemistry, my_fields);
  }

  if (my_chemistry->primordial_chemistry > 0)  {

    // Correct the species to ensure consistency (i.e. type conservation)

    f_wrap::make_consistent_g(imetal, dom, my_chemistry, my_rates, my_fields);

  }

  return ierr;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
