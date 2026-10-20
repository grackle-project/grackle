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

#include <optional>
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
#include "support/status_reporting.hpp"
#include "time_deriv_0d.hpp"
#include "utils-cpp.hpp"

namespace GRIMPL_NAMESPACE_DECL {

namespace step_nr_detail {

/// Describes the layout of the state vector evolved by newton-raphson solver
/// and if energy-evolution is decoupled from chemistry, specifies the fixed
/// energy
///
/// The first @ref n_contained_species of a state vector hold densities of
/// dynamically evolved species. If internal energy is in the state-vector
/// (it is stored as the last entry of the vector). If energy is not tracked
/// as part of the state-vector, an object of this type should then track
/// the fixed_specific_energy for the current timestep
///
/// Implementation Note
/// -------------------
/// @ref StateLayoutAndE::species_ids is an array-like object that maps the
/// `i`th index of the state-vector (for `i < n_contained_species`) to a
/// species_id.
///
/// This needs to be an array-like object because grackle currently associates
/// species ids with species at compile-time even if the species is not used.
/// Thus, it's possible for grackle to be configured in a way where the
/// species_ids of actively evolved species are not contiguous. For example:
/// - if `primordial_chemistry = 3` & `metal_chemistry = 1`, the species_id
///   associated with certain primordial species (e.g. `DM`) aren't used
/// - if `metal_chemistry = 1`, `dust_chemistry = 2` & `dust_species = 1`:
///   the species_ids of certain metal species (e.g. `Fe`) aren't used
///
/// After we eliminate @ref SpLUT from the codebase, some light refactoring
/// will be possible to configure the precise mapping between species ids and
/// species at configuration time (this while preserving relative orderings of
/// species that are encoded in objects like @ref PrimoridalSpLUT,
/// @ref MetalSpLUT, and @ref DustSpLUT). This will make it possible to ensure
/// that all active species ids have contiguous values. At that point, we'll be
/// able to convert the @ref StateLayoutAndE::species_ids will no longer need
/// to be an array-like object
struct StateLayoutAndE {

  /// the number of indices of the state-vector that encode the density of a
  /// dynamically evolved species
  int n_contained_species;
  /// `species_ids[idx]` gives the `species_id` stored at index `idx` of the
  /// state-vector for the first @ref n_contained_species indices of the state
  /// vector
  const int* species_ids;
  /// If energy evolution is decoupled from energy evolution, this specifies
  /// the fixed internal energy
  std::optional<double> fixed_specific_energy;

  /// returns whether the specific internal energy is tracked as part of the
  /// dynamically evolved state vector (i.e. whether energy evolution is
  /// coupled to chemistry evolution)
  bool coupled_e_evolution() const {return !fixed_specific_energy.has_value();}

  /// return state-vector index holding the specific internal energy
  int energy_index() const {
    if (!coupled_e_evolution()) {
      GR_INTERNAL_ERROR("SANITY CHECK FAILED!");
    } else {
      return n_contained_species;
    }
  }

  int state_vec_len() const {
    return n_contained_species + coupled_e_evolution();
  }
};

/// @brief max length of a state vector evolved with the newton raphson solver
///
/// For context:
/// - the state vector is typically referred to by the variable name `y`
/// - the state vector has distinct components for each evolved species and may
///   have a component to track the specific internal energy
inline constexpr int MAX_Y_LEN = MAX_EVOLVED_SPECIES_FIELDS + 1;

/// @brief Creates a @ref StateLayoutAndE object describing the state vector for
///    the current calculation.
///
/// At the time of writing, a given grackle configuration generally has 4
/// possible layouts for the state vector:
/// 1. exclude metal-related chemistry with operator-split energy evolution
/// 2. exclude metal-related chemistry with coupled energy evolution
/// 3. include metal-related chemistry with operator-split energy evolution
/// 4. include metal-related chemistry with coupled energy evolution
///
/// Obviously the fourth case is most general-purpose, but the state vector
/// is definitely the longest. In the current configuration, coupling energy
/// evolution also significantly increases the cost of the calculation (this
/// is definitely fixable -- especially if we compute analytic derivatives)
///
/// Implementation Note
/// -------------------
/// For the initial implementation this is intended to be an ephemeral object
/// that is temporarily constructed during each invocation of
/// @ref step_rate_newton_raphson. If the cost of rebuilding this object is a
/// concern, we could prebuild this object and store it inside of opaque
/// storage.
///
/// It's worth mentioning that once we finish eliminating all references of
/// @ref SpLUT from the codebase, the optimizations described in the docstring
/// for @ref StateLayoutAndE will become possible. If we implement those
/// optimizations, constructing a @ref StateLayoutFactory object will become
/// EXTREMELY cheap
class StateLayoutFactory {
  std::vector<int> idsp_;
  int n_sp_prim_;
  int n_sp_all_;

  public:

  StateLayoutFactory(const StateLayoutFactory&)=delete;
  StateLayoutFactory& operator=(const StateLayoutFactory&)=delete;

  explicit StateLayoutFactory(const PartMap<SpKind>& sp_kind_map) {
    idsp_.resize(MAX_Y_LEN);
    // if SpKind::METAL or SpKind::DUST don't correspond to any species in the
    // current configuration, then the corresponding IndexInterval1D object's
    // start & stop data members will be equal to each other
    const IndexInterval1D species_idx_bounds_arr[3] = {
      sp_kind_map.part_bounds(SpKind::PRIMORDIAL),
      sp_kind_map.part_bounds(SpKind::METAL),
      sp_kind_map.part_bounds(SpKind::DUST)
    };
    n_sp_prim_ = species_idx_bounds_arr[0].stop - species_idx_bounds_arr[0].start;

    int offset = 0;
    for (const IndexInterval1D& bounds: species_idx_bounds_arr) {
      int n_to_add = bounds.stop - bounds.start;
      for (int index = 0; index < n_to_add; index++) {
        idsp_[index + offset] = index + bounds.start;
      }
      offset += n_to_add;
    }
    n_sp_all_ = offset;
  }

  // it's the caller's responsibility to ensure that the StateLayoutAndE objects
  // constructed by the following methods don't outlive the factory (otherwise,
  // the species_ids pointer will be a dangling)

  StateLayoutAndE create_coupled_e_evo(bool include_metal) const {
    return StateLayoutAndE{
        .n_contained_species = include_metal ? n_sp_all_ : n_sp_prim_,
        .species_ids = idsp_.data(),
        .fixed_specific_energy = std::nullopt
      };
  }

  StateLayoutAndE create_fixed_e(bool include_metal,
                                 double fixed_specific_energy) const {
    return StateLayoutAndE{
        .n_contained_species = include_metal ? n_sp_all_ : n_sp_prim_,
        .species_ids = idsp_.data(),
        .fixed_specific_energy = std::optional<double>(fixed_specific_energy)
      };
  }
};

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
  const double* y,
  double* ydot,
  const StateLayoutAndE& state_layout,
  grackle::impl::time_deriv_0d::ContextPack& pack,
  std::vector<gr_float>& rhosp_grflt,         // <- preallocated temp buffer
  grackle::impl::SpeciesCollection& rhosp_dot // <- preallocated temp buffer
) {

  gr_float local_eint[1];
  double eint_dot_specific[1];

  for (int sp_idx = 0; sp_idx < MAX_EVOLVED_SPECIES_FIELDS; sp_idx++) {
    rhosp_grflt[sp_idx] = static_cast<gr_float>(0.0);
  }

  // 1. copy values out of y into rhosp_grflt & (if applicable) local_eint
  //    - at some point, we need to add some logic to deal with edge cases
  //      that can arise when `sizeof(gr_float) < sizeof(double)
  for (int vec_i = 0; vec_i < state_layout.n_contained_species; vec_i++) {
    int species_id = state_layout.species_ids[vec_i];
    rhosp_grflt[species_id] = static_cast<gr_float>(y[vec_i]);
  }
  if (state_layout.fixed_specific_energy.has_value()) {
    local_eint[0] =
        static_cast<gr_float>(state_layout.fixed_specific_energy.value());
  } else {
    local_eint[0] = static_cast<gr_float>(y[state_layout.energy_index()]);
  }

  // 2. call the wrapped function
  grackle::impl::time_deriv_0d::derivatives(
    rhosp_grflt.data(), rhosp_dot, local_eint, eint_dot_specific, pack
  );

  // 3. copy the computed derivatives into ydot
  for (int vec_i = 0; vec_i < state_layout.n_contained_species; vec_i++) {
    int species_id = state_layout.species_ids[vec_i];
    ydot[vec_i] = rhosp_dot.data[species_id][0];
  }
  if (!state_layout.fixed_specific_energy.has_value()) {
    ydot[state_layout.energy_index()] = eint_dot_specific[0];
  }
}

}  // namespace step_nr_detail

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

  // There may be an argument for allocating the following at a higher
  // level function, but we will leave that for after transcription
  std::vector<double> y_initial(step_nr_detail::MAX_Y_LEN);
  std::vector<double> y(step_nr_detail::MAX_Y_LEN);
  std::vector<double> y_buf1(step_nr_detail::MAX_Y_LEN);
  std::vector<double> ydot_buf1(step_nr_detail::MAX_Y_LEN);

  // initialize the ode_solver
  // -> a solution is considered converged when the max relative difference
  //    magnitude of any vector component does not exceed this value
  double max_rtol = 1.e-8;
  // -> maximum number of iterations before giving up on a provided timestep
  int maxiter_per_dt = 20;
  integrate::StiffNewtonRaphson ode_solver(step_nr_detail::MAX_Y_LEN, max_rtol,
                                           maxiter_per_dt);

  std::vector<double> ode_scratch_buf(ode_solver.num_scratch_buf_elements());

  const step_nr_detail::StateLayoutFactory ylayout_factory(
    my_rates->opaque_storage->sp_kind_map);

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

  // the following check was inspired by a compiler warning indicating that
  // y_len won't be initialized if this condition isn't met
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
      bool coupled_e_evolution = (imp_eng[i] == 1);

      // If density and temperature are low, update gas energy explicitly
      if (my_chemistry->with_radiative_cooling == 1 && !coupled_e_evolution) {
        e(i,j,k)  = e(i,j,k) + (gr_float)(edot[i]/d(i,j,k)*dtit[i] );
      }

      bool include_metal = (itmask_metal[i] != MASK_FALSE);
      step_nr_detail::StateLayoutAndE state_layout = (coupled_e_evolution)
        ? ylayout_factory.create_coupled_e_evo(include_metal)
        : ylayout_factory.create_fixed_e(include_metal, e(i,j,k));
      int y_len = state_layout.state_vec_len();

      // copy values into y_initial
      for (int isp = 0; isp < state_layout.n_contained_species; isp++) {
        int species_id = state_layout.species_ids[isp];
        y_initial[isp] = sp_densities(species_id, i);
      }
      if (coupled_e_evolution) {
        y_initial[state_layout.energy_index()] = e(i,j,k);
      }

      // configure pack for use during the current index
      t_deriv::configure_ContextPack(
        &pack, my_fields, field_idx1d, itmask_metal[i], coupled_e_evolution
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
      // - the `y` buffer holds all variables actively evolved by the current
      //   grackle configuration.
      // - for the current timestep, only a subset of the variables may be
      //    evolved. The selection of evolved variables is tracked by the
      //   `idsp` buffer.
      // - this reduced set of variables corresponds to `y(t)`

      // actually define the callback function:
      // - `y` corresponds to the current value of the variable y
      // - `ydot` is an output buffer for holding f(y)
      // - `jacobian` is an output buffer for holding the Jacobian of f
      //
      // Important: index `i` in these buffers correspond to the `i`th variable
      //            being evolved in the current timestep. The corresponding
      //            entry in `dsp` is at the index `idsp[isp]`.
      const auto calc_f_and_jacobian = [&](
          const double* y, double* ydot,
          View<double**>& jacobian) -> void {

        wrapped_calc_derivatives(y, ydot, state_layout,
                                 pack, rhosp_grflt, rhosp_dot);

        // fill in the jacobian matrix for the time derivative
        // -> to accomplish this, we use finite differences to estimate
        //    partial derivative for each evolved variable (i.e. the species
        //    densities and possibly the total energy)
        for (int jsp = 0; jsp < y_len; jsp++) {
          double dspj = eps * y[jsp];
          for (int isp = 0; isp < y_len; isp++) {
            if (isp == jsp) {
              y_buf1[isp] = y[isp] + dspj;
            } else {
              y_buf1[isp] = y[isp];
            }
          }

          wrapped_calc_derivatives(y_buf1.data(),
                                   ydot_buf1.data(), state_layout,
                                   pack, rhosp_grflt, rhosp_dot);

          for (int isp = 0; isp < y_len; isp++) {
            double tderiv = ydot[isp];
            double offset_tderiv = ydot_buf1[isp];
            if ((y[isp] == 0.0) && (offset_tderiv == tderiv)) {
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
      // -> when we exit this loop, y will hold the values of all of
      //    the evolved variables
      bool is_converged = false;
      while (!is_converged) {

        // copy values into y for variables that are being evolved
        // (we need to do this each time we enter the loop since y
        // is mutated in place as we try to evolve things)
        for (int isp = 0; isp < state_layout.state_vec_len(); isp++) {
          y[isp] = y_initial[isp];
        }

        // Iteration to solve ODEs



        const bool enforce_positive_non_NaN =
      (imp_eng[i] == 1) && (my_chemistry->primordial_chemistry > 0) &&
      (my_chemistry->with_radiative_cooling == 1);

        integrate::IntegrateResult rslt = ode_solver.step(
            y.data(), dtit[i], calc_f_and_jacobian, y_len,
            ode_scratch_buf.data(), enforce_positive_non_NaN, d(i, j, k));
        is_converged = rslt.converged();

        // Check if the fractions are valid after an iteration
        if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
          for (isp = 0; isp<y_len; isp++) {
            if ( std::isnan(y[isp]) || (y[isp] <= 0.) )  {
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

      // overwrite the dynamically evolved species fields (tracked by my_fields
      // and accessed thru sp_densities) with the evolved values
      for (int isp = 0; isp < state_layout.n_contained_species; isp++) {
        sp_densities(state_layout.species_ids[isp], i) = y[isp];
      }

      if (coupled_e_evolution) {
        e(i,j,k) = y[state_layout.energy_index()];
      }
    }
  }

  t_deriv::drop_MainScratchBuf(&main_scratch_buf);
  grackle::impl::drop_SpeciesCollection(&rhosp_dot);
}

} // namespace GRIMPL_NAMESPACE_DECL


#endif /* STEP_RATE_NEWTON_RAPHSON_HPP */
