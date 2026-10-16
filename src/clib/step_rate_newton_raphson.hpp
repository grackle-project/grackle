//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines/declares the step_rate_newton_raphson function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// step_rate_newton_raphson function from FORTRAN to C++

#ifndef STEP_RATE_NEWTON_RAPHSON_HPP
#define STEP_RATE_NEWTON_RAPHSON_HPP

#include <vector>

#include "grackle.h"             // gr_float
#include "field_adaptor.hpp"
#include "fortran_func_decls.h"  // gr_mask_int
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "internal_types.hpp"
#include "internal_units.hpp"
#include "math/integrate.hpp"
#include "opaque_storage.hpp"
#include "support/View.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "time_deriv_0d.hpp"
#include "utils-cpp.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// acts as a temporary adaptor layer for time_deriv_0d::derivatives:
///
/// @note
/// Some additional context
/// - during transcription, time_deriv_0d::derivatives (formerly called
///   lookup_cool_rates0d) needed to be adapted to make use of data structures
///   already adopted in other parts of the transcribed code
///   - (recall that these data structures were introduced to remove HUNDREDS
///     of function arguments and HUNDREDS of local variables)
/// - since the newton-raphson solver already does a bunch of "translation"
///   between data representations, it should directly map into these preferred
///   data structures
/// - as a short-term (and extremely temporary) solution, we this function
///   to copy between formats & wrap the time_deriv_0d::derivatives function.
///   This allows to finish polishing time_deriv_0d::derivatives without
///   refactoring the newton-raphson scheme)
inline void wrapped_calc_derivatives(
  const double* dsp,
  double* dspdot,
  grackle::impl::time_deriv_0d::ContextPack& pack,
  std::vector<gr_float>& rhosp_grflt,         // <- preallocated temp buffer
  grackle::impl::SpeciesCollection& rhosp_dot // <- preallocated temp buffer
) {

  gr_float local_eint[1];
  double eint_dot_specific[1];

  // 1. copy values out of dsp into rhosp_grflt and local_eint
  //    - at some point, we need to add some logic to deal with edge cases
  //      that can arise when `sizeof(gr_float) < sizeof(double)
  for (int sp_idx = 0; sp_idx < MAX_EVOLVED_SPECIES_FIELDS; sp_idx++) {
    rhosp_grflt[sp_idx] = static_cast<gr_float>(dsp[sp_idx]);
  }
  local_eint[0] = static_cast<gr_float>(dsp[MAX_EVOLVED_SPECIES_FIELDS]);

  // 2. call the wrapped function
  grackle::impl::time_deriv_0d::derivatives(
    rhosp_grflt.data(), rhosp_dot, local_eint, eint_dot_specific, pack
  );

  // 3. copy the computed derivatives into dspdot
  for (int sp_idx = 0; sp_idx < MAX_EVOLVED_SPECIES_FIELDS; sp_idx++) {
    dspdot[sp_idx] = rhosp_dot.data[sp_idx][0];
  }
  dspdot[MAX_EVOLVED_SPECIES_FIELDS] = eint_dot_specific[0];
}


/// An alternative to step_rate_g for evolving the species rate equations that
/// employs the Newton-Raphson scheme rather than Gauss-Seidel
///
/// Whereas the Gauss-Seidel scheme always operator-splits evolution of the
/// energy equation, this function internally evolves the energy equation
///
/// @note
/// The values in the various buffers (tgas, tdust, metallicity, dust2gas,
/// rhoH, mmw, edot, grain_temperatures, logTlininterp_buf,
/// coolingheating_buf, chemheatrates_buf), hold undefined
/// values after this function call wherever `itmask_nr` indicates that this
/// function runs.
/// - this has **ALWAYS** been the case. Historically these buffers have been
///   reused when computing finite differences
inline void step_rate_newton_raphson(
    IndexRange idx_range, double dom, double chunit, double dx_cgs,
    double c_ljeans, double* dtit, double* tgas, double* tdust,
    double* metallicity, double* dust2gas, double* rhoH, double* mmw,
    double* nelec_times_mH, double* edot, gr_mask_type anydust,
    const gr_mask_type* itmask_nr, const gr_mask_type* itmask_metal,
    const int* imp_eng, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    SpeciesMultiView<gr_float> sp_densities,
    photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
    grackle::impl::GrainSpeciesCollection grain_temperatures,
    grackle::impl::LnTLinInterpBuf logTlininterp_buf,
    grackle::impl::CoolHeatScratchBuf coolingheating_buf,
    grackle::impl::ChemHeatingRates chemheatrates_buf) {
  // shorten `grackle::impl::time_deriv_0d` to `t_deriv` within this function
  namespace t_deriv = ::grackle::impl::time_deriv_0d;

  // Density, energy and velocity fields fields

  FortranView<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  FortranView<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Local variable
  int isp;

  // the following specifies the historical 1-based index that we would use to
  // hold energy
  const int i_eng = MAX_EVOLVED_SPECIES_FIELDS + 1;
  // There may be an argument for allocating the following at a higher
  // level function, but we will leave that for after transcription
  std::vector<double> dsp(i_eng);
  std::vector<double> reduced_dsp(i_eng);
  std::vector<double> full_dsp_buf(i_eng);
  std::vector<double> full_dspdot_buf(i_eng);
  std::vector<double> full_dspdot_buf1(i_eng);

  // initialize the ode_solver
  // -> a solution is considered converged when the max relative difference
  //    magnitude of any vector component does not exceed this value
  double max_rtol = 1.e-8;
  // -> maximum number of iterations before giving up on a provided timestep
  int maxiter_per_dt = 20;
  integrate::StiffNewtonRaphson ode_solver(i_eng, max_rtol, maxiter_per_dt);

  std::vector<double> ode_scratch_buf(ode_solver.num_scratch_buf_elements());

  // (In the future, we may want to reconsider when/how we allocate
  // the following 3 variables)
  std::vector<int> idsp;
  // Another parameter
  const double eps = 1.e-4;

  // allocate some scratch buffers. In the future, when we refactor to support
  // calculations of the time derivatives of multiple values at once, we should
  // forward the buffers passed into this routine as arguments rather than
  // allocating separate buffers
  t_deriv::MainScratchBuf main_scratch_buf =
    t_deriv::new_MainScratchBuf(
        grackle::impl::GrainMetalInjectPathways_get_n_log10Tdust_vals(
            my_rates->opaque_storage->inject_pathway_props
        )
    );

  // collect args that are forwarded to the time-derivative calculation and are
  // effectively frozen between various calls
  t_deriv::FrozenSimpleArgs frozen_tderiv_args = {
    dom, chunit, dx_cgs, c_ljeans, anydust, my_chemistry,
    my_rates, my_uvb_rates, internalu
  };

  // partially initialize the struct we will use for the time derivative calc
  t_deriv::ContextPack pack = t_deriv::new_ContextPack(
    frozen_tderiv_args, main_scratch_buf
  );

  std::vector<gr_float> rhosp_grflt(MAX_EVOLVED_SPECIES_FIELDS);
  grackle::impl::SpeciesCollection rhosp_dot =
    grackle::impl::new_SpeciesCollection(1);

  const PartMap<SpKind>& sp_kind_map = my_rates->opaque_storage->sp_kind_map;
  // if SpKind::METAL or SpKind::DUST don't correspond to any species in the
  // current configuration, then the corresponding IndexInterval1D object's
  // start & stop data members will be equal to each other
  const IndexInterval1D species_idx_bounds_arr[3] = {
    sp_kind_map.part_bounds(SpKind::PRIMORDIAL),
    sp_kind_map.part_bounds(SpKind::METAL),
    sp_kind_map.part_bounds(SpKind::DUST)
  };

  int n_species_idx_bounds_nometal = 1;
  int n_species_idx_bounds_all = 3;

  // the following check was inspired by a compiler warning indicating that
  // nsp won't be initialized if this condition isn't met
  GRIMPL_REQUIRE((my_chemistry->primordial_chemistry > 0),
    "this function can't support primordial_chemistry == %d",
    my_chemistry->primordial_chemistry
  );

  // The following was extracted from another subroutine
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    const int j = idx_range.j;
    const int k = idx_range.k;

    // remap the 3D index to a 1D index
    const int field_idx1d = layoutleft_3D_index_to_1D_(
        my_fields->grid_dimension, i, j, k
    );

    if (itmask_nr[i] != MASK_FALSE)  {

      // If density and temperature are low, update gas energy explicitly

      if (my_chemistry->with_radiative_cooling == 1)  {
        if (imp_eng[i] == 0)  {
          e(i,j,k)  = e(i,j,k) +
               (gr_float)(edot[i]/d(i,j,k)*dtit[i] );
        }
      }

      const int n_species_idx_bounds = (itmask_metal[i] != MASK_FALSE)
          ? n_species_idx_bounds_all : n_species_idx_bounds_nometal;

      // initialize arrays
      int nsp = imp_eng[i];
      for (int bound_i = 0; bound_i < n_species_idx_bounds; bound_i++) {
        const IndexInterval1D& bounds = species_idx_bounds_arr[bound_i];
        nsp += (bounds.stop - bounds.start);
      }
      idsp.resize(nsp);

      // copy values into dsp from sp_densities
      for (int bound_i = 0; bound_i < n_species_idx_bounds; bound_i++) {
        const IndexInterval1D& bounds = species_idx_bounds_arr[bound_i];
        for (int sp_idx = bounds.start; sp_idx < bounds.stop; sp_idx++) {
          dsp[sp_idx] = sp_densities(sp_idx, i);
        }
      }
      dsp[i_eng-1] = e(i,j,k);

      // here, we fill in the idsp array
      // -> the idsp array is used for mapping between the compressed vector
      //    form (that only has nsp elements) and dsp (that always has
      //    `MAX_EVOLVED_SPECIES_FIELDS + 1` entries)
      // -> in the future, we should construct this array first (before doing
      //    anything else):
      //    -> it gives us the value of nsp for free!
      //    -> using this array with SpeciesMultiView lets us cut out
      //       (duplicated) logic when we copy values from my_fields to dsp and
      //       the reverse (we can cut ~70 lines in each place)
      {
        int offset = 0;
        for (int bound_i = 0; bound_i < n_species_idx_bounds; bound_i++) {
          const IndexInterval1D& bounds = species_idx_bounds_arr[bound_i];
          int n_to_add = bounds.stop - bounds.start;
          for (int index = 0; index < n_to_add; index++) {
            idsp[index + offset] = index + bounds.start;
          }
          offset += n_to_add;
        }
        if ( imp_eng[i] ==1 )  {
          idsp[offset++] = i_eng-1;
        }
        GRIMPL_REQUIRE((offset == nsp), "SANITY CHECK FAILED");
      }

      // configure pack for use during the current index
      t_deriv::configure_ContextPack(
        &pack, my_fields, field_idx1d, itmask_metal[i], imp_eng[i]
      );

      // copy values at the current index from the external scratch buffers
      // passed into the current function as args into the corresponding
      // buffers managed by pack
      // -> see the called function's docstring for an explanation as to why
      //    this is not currently combined with the prior function call
      // -> in a lot of cases, the need to copy values is dictated by the value
      //    of imp_eng
      // -> in the future, we should really refactor so that this isn't
      //    necessary
      t_deriv::scratchbufs_copy_into_pack(
        i, &pack, tgas, tdust, metallicity, dust2gas, rhoH, mmw, nelec_times_mH,
        edot, grain_temperatures, logTlininterp_buf,
        coolingheating_buf, chemheatrates_buf
      );

      // we want to integrate a system of coupled stiff ODEs
      //
      // From a purely mathematical perspective:
      // - let `y(t)` specifies the (vector-valued) state of this system at a
      //   time `t`
      // - let `dy/dt` be equivalent to a (vector-valued) function `f(y)`.
      // - If the current time is t₀, then we want to compute y(t₀ + Δt) from
      //   y(t₀) for a timestep Δt.
      // - To do this, we must define a callback function that computes
      //   f(y) and the Jacobian matrix for f(y).
      //
      // From a more practical perspective:
      // - the `dsp` buffer holds all variables actively evolved by the current
      //   grackle configuration.
      // - for the current timestep, only a subset of the variables may be
      //    evolved. The selection of evolved variables is tracked by the
      //   `idsp` buffer.
      // - this reduced set of variables corresponds to `y(t)`

      // initialize full_dsp_buf with copies of all entries in dsp
      // -> this buffer is used internally by the callback function. It has
      //    space for every variable that the current grackle configuration
      //    evolves.
      // -> the basic idea, is that the callback will overwrite the elements
      //    corresponding to the subset of variables being evolved in the
      //    current timestep, won't mutate other elements (if any), and will
      //    pass this to the buffer that does the heavy lifting of computing
      //    derivatives.
      std::memcpy(full_dsp_buf.data(), dsp.data(), sizeof(double)*i_eng);

      // actually define the callback function:
      // - reduced_dsp corresponds to the current value of the variable y
      // - reduced_dspdot is an output buffer for holding f(y)
      // - jacobian is an output buffer for holding the Jacobian of f
      //
      // Important: index `i` in these buffers correspond to the `i`th variable
      //            being evolved in the current timestep. The corresponding
      //            entry in `dsp` is at the index `idsp[isp]`.
      const auto calc_f_and_jacobian = [&](
          const double* reduced_dsp, double* reduced_dspdot,
          View<double**>& jacobian) -> void {

        // copy values out of reduced_dsp into full_dsp_buf
        for (int isp = 0; isp < nsp; isp++) {
          full_dsp_buf[idsp[isp]] = reduced_dsp[isp];
        }

        wrapped_calc_derivatives(full_dsp_buf.data(),
                                 full_dspdot_buf.data(), pack,
                                 rhosp_grflt, rhosp_dot);

        // copy the reduced set of values that we care about from
        // full_dspdot_buf into dspdot
        for (int isp = 0; isp < nsp; isp++) {
          reduced_dspdot[isp] = full_dspdot_buf[idsp[isp]];
        }

        // fill in the jacobian matrix for the time derivative
        // -> to accomplish this, we use finite differences to estimate
        //    partial derivative for each evolved variable (i.e. the species
        //    densities and possibly the total energy)
        for (int jsp = 0; jsp < nsp; jsp++) {
          double dspj = eps * reduced_dsp[jsp];
          for (int isp = 0; isp < nsp; isp++) {
            if (isp == jsp) {
              full_dsp_buf[idsp[isp]] = reduced_dsp[isp] + dspj;
            } else {
              full_dsp_buf[idsp[isp]] = reduced_dsp[isp];
            }
          }

          wrapped_calc_derivatives(full_dsp_buf.data(),
                                   full_dspdot_buf1.data(),
                                   pack, rhosp_grflt, rhosp_dot);

          for (int isp = 0; isp < nsp; isp++) {
            double tderiv = full_dspdot_buf[idsp[isp]];
            double offset_tderiv = full_dspdot_buf1[idsp[isp]];
            if ((reduced_dsp[isp] == 0.0) && (offset_tderiv == tderiv)) {
              jacobian(isp, jsp) = 0.0;
            } else {
              jacobian(isp, jsp) = (offset_tderiv - tderiv) / dspj;
            }
          }
        }
      };

      // iteratively try to evolve chemistry equations until answer converges
      // -> First, we try to use the current timestep value tracked by dtit[i].
      //    If the answer doesn't converge, we will try again with shorter and
      //    shorter timesteps until chemistry converges
      // -> when we exit this loop, reduced_dsp will hold the values of all of
      //    the evolved variables
      bool is_converged = false;
      while (!is_converged) {

        // copy values into reduced_dsp for variables that are being evolved
        // (we need to do this each time we enter the loop since reduced_dsp
        // is mutated in place as we try to evolve things)
        for (int isp = 0; isp < nsp; isp++) {
          reduced_dsp[isp] = dsp[idsp[isp]];
        }

        // Iteration to solve ODEs



        const bool enforce_positive_non_NaN =
      (imp_eng[i] == 1) && (my_chemistry->primordial_chemistry > 0) &&
      (my_chemistry->with_radiative_cooling == 1);

        integrate::IntegrateResult rslt = ode_solver.step(
            reduced_dsp.data(), dtit[i], calc_f_and_jacobian, nsp,
            ode_scratch_buf.data(), enforce_positive_non_NaN, d(i, j, k));
        is_converged = rslt.converged();

        // Check if the fractions are valid after an iteration
        if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
          for (isp = 0; isp<nsp; isp++) {
            if ( std::isnan(reduced_dsp[isp]) || (reduced_dsp[isp] <= 0.) )  {
              is_converged = false;
            }
          }
        }

        // cut the timestep in half for the next loop
        // todo: handle worst case scenario (i.e. it never converging or no
        //       converging unless timestep is 0)
        if(!is_converged)  {
          dtit[i] = 0.5 * dtit[i];
        }
      }

      // copy values from reduced_dsp back to dsp
      for (int isp = 0; isp < nsp; isp++) {
        dsp[idsp[isp]] = reduced_dsp[isp];
      } 


      // overwrite the scratch-buffer values
      // -> It makes no logical sense for us to do this, but we do it for the
      //    sake of backwards compatability
      // -> we should totally delete this function (we may already be able to
      //    do so)
      t_deriv::scratchbufs_copy_from_pack(
        i, &pack, tgas, tdust, metallicity, dust2gas, rhoH, mmw, nelec_times_mH,
        edot, grain_temperatures, logTlininterp_buf,
        coolingheating_buf, chemheatrates_buf
      );


      // overwrite the fields tracked by my_fields with the evolved values
      for (int bound_i = 0; bound_i < n_species_idx_bounds; bound_i++) {
        const IndexInterval1D& bounds = species_idx_bounds_arr[bound_i];
        for (int sp_idx = bounds.start; sp_idx < bounds.stop; sp_idx++) {
          sp_densities(sp_idx, i) = dsp[sp_idx];
        }
      }
      e(i,j,k)     = dsp[i_eng-1];

      idsp.clear();
    }
  }

  t_deriv::drop_MainScratchBuf(&main_scratch_buf);
  grackle::impl::drop_SpeciesCollection(&rhosp_dot);
}

} // namespace GRIMPL_NAMESPACE_DECL


#endif /* STEP_RATE_NEWTON_RAPHSON_HPP */
