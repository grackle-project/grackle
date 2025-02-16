// See LICENSE file for license and copyright information

/// @file time_deriv_0d.hpp
/// @brief Defines machinery to calculate the time derivative for a single zone

#ifndef TIME_DERIV_0D_HPP
#define TIME_DERIV_0D_HPP

#include "chemistry_solver_funcs.hpp"
#include "fortran_func_wrappers.hpp"
#include "grackle.h"
#include "grackle_macros.h" // GRACKLE_FREE
#include "index_helper.h"
#include "internal_types.hpp"
#include "utils-field.hpp"

// we choose to adopt a longer, more descriptive namespace here so that the
// handful of functions defined in this file can have shorter names (in the
// future, if we are willing to define methods on a struct, we can definitely
// shorten the namespace name)
namespace grackle::impl::time_deriv_0d {

/// this is a collection of the arguments that won't change between successive
/// time derivative calculations
struct FrozenSimpleArgs {
  // the following batch of args are all forwarded
  int imetal;
  int iter;
  double dom;
  double chunit;
  double dx_cgs;
  double c_ljeans;
  gr_mask_type anydust;
  chemistry_data* my_chemistry;
  chemistry_data_storage* my_rates;
  photo_rate_storage my_uvb_rates;
  InternalGrUnits internalu;
};

/// these are collections of scratch buffers for solving the 0d time derivative
///
/// @note
/// Since the time derivative calculation only operates on a single zone, we
/// choose to freshly allocate data for each of these buffers (this is a small
/// fraction of the required memory). In the future, if we refactor the
/// calculation to operate on multiple elements at once, then this should reuse
/// preallocated buffers.
struct MainScratchBuf {
  GrainSpeciesCollection grain_temperatures;
  LogTLinInterpScratchBuf logTlininterp_buf;
  Cool1DMultiScratchBuf cool1dmulti_buf;
  CoolHeatScratchBuf coolingheating_buf;
  ChemHeatingRates chemheatrates_buf;

  // the remaining buffers were originally reallocated (mostly on the stack)
  // every time calculated the time derivatives were computed
  ColRecRxnRateCollection kcr_buf;
  PhotoRxnRateCollection kshield_buf;
  GrainSpeciesCollection grain_growth_rates;
  double* k13dd; // <- only used within lookup_cool_rates1d_g
};

MainScratchBuf new_MainScratchBuf(void) {
  int nelem = 1;
  MainScratchBuf out;
  out.grain_temperatures = new_GrainSpeciesCollection(nelem);
  out.logTlininterp_buf = new_LogTLinInterpScratchBuf(nelem);
  out.cool1dmulti_buf = new_Cool1DMultiScratchBuf(nelem);
  out.coolingheating_buf = new_CoolHeatScratchBuf(nelem);
  out.chemheatrates_buf = new_ChemHeatingRates(nelem);

  out.kcr_buf = new_ColRecRxnRateCollection(nelem);
  out.kshield_buf = new_PhotoRxnRateCollection(nelem);
  out.grain_growth_rates = new_GrainSpeciesCollection(nelem);
  out.k13dd = (double*)malloc(sizeof(double)*14*nelem);
  return out;
}

void drop_MainScratchBuf(MainScratchBuf* ptr) {
  drop_GrainSpeciesCollection(&ptr->grain_temperatures);
  drop_LogTLinInterpScratchBuf(&ptr->logTlininterp_buf);
  drop_Cool1DMultiScratchBuf(&ptr->cool1dmulti_buf);
  drop_CoolHeatScratchBuf(&ptr->coolingheating_buf);
  drop_ChemHeatingRates(&ptr->chemheatrates_buf);

  drop_ColRecRxnRateCollection(&ptr->kcr_buf);
  drop_PhotoRxnRateCollection(&ptr->kshield_buf);
  drop_GrainSpeciesCollection(&ptr->grain_growth_rates);
  GRACKLE_FREE(ptr->k13dd);
}

/// this is a collections of values intended to act as 1-element arrays and
/// that don't need to be dynamically allocated
///
/// @note
/// When we ultimately refactor the time derivative calculation to operate on
/// multiple elements at once, we will need to merge this object with
/// MainScratchBuf and track pointers to previously allocated memory buffer
/// for all cases
struct Assorted1ElemBuf {
  double p2d[1];
  double tgas[1];
  double tdust[1];
  double metallicity[1];
  double dust2gas[1];
  double rhoH[1];
  double mmw[1];
  double h2dust[1];
  double edot[1];

  // the remaining buffers were originally reallocated (on the stack)
  // every time calculated the time derivatives were computed.
  gr_mask_type itmask[1];
  // These are used to compute values that we totally ignore in this context
  double dedot[1];
  double HIdot[1];
};

/// this struct is used to organize some temporary data that is used for
/// computing time derivatives of species data (and possibly, internal energy)
/// in a single zone
///
/// most of the data here acts a little like an adaptor layer
/// - we effectively adapt a representation of all species (and possibly internal
///   energy) from a vector form to the standard data structures to do typical
///   calculations and then we adapt back to the vector format
/// - to facillitate this, we effectively create an instance of
///   grackle_field_data that acts like a 1-elemt slice of the grackle_field_data
///   instance that the user passed in.
/// - this is highly inefficient, but it is logically consistent with the
///   original fortran code from before transcription. (We should refactor this
///   in the future after we finish transcription)
struct ContextPack {
  /// holds the local value of itmask_metal
  gr_mask_type local_itmask_metal;
  /// specifies how we handling edot in the calculation. This comes from
  /// an array historically known as imp_eng
  int local_edot_handling;

  /// this represents the precomputed index-range
  /// - if we transition to an implementation where we compute the time
  ///   derivatives for multiple sets of physical conditions at a time, this
  ///   will need to be computed on the fly
  IndexRange idx_range_1_element;

  /// the idea is that this will hold data for a single zone
  grackle_field_data fields;

  /// @defgroup field_data_metadata
  /// The arrays in this group are used to hold the various 3-element arrays
  /// that are used as members of the `fields` member
  /** @{ */
  int grid_dimension[3];
  int grid_start[3];
  int grid_end[3];
  /** @} */

  /// @defgroup general_time_deriv_packs
  /// The members in this group store packs of data that are used within the
  /// time derivative calculation. There are usually analogues to the data used
  /// in othere parts of Grackle
  /** @{ */

  /// the following batch of arguments essentially forward several arguments
  /// onto the calculation. These are "frozen"; they don't change across calls
  /// to the calculation for different spatial locations.
  FrozenSimpleArgs fwd_args;

  /// the struct containing the primary scratch buffers (the lifetimes of the
  /// buffers in this member are intended to all outlive the lifetime of this
  /// struct)
  MainScratchBuf main_scratch_buf;

  /// collection of other assorted scratch buffers.
  ///
  /// @note
  /// While main_scratch_buf holds collections of buffers, this holds loose
  /// buffers. Because this is entirely composed of loose buffers and the
  /// calculation currently only operates on a single zone, these buffers are
  /// all specified as statically sized arrays (i.e. there's no need to
  /// manually allocate and deallocate the buffers)
  Assorted1ElemBuf other_scratch_buf;

  /** @} */

};

typedef struct ContextPack ContextPack;

/// partially set up a new time_deriv_pack
///
/// Ideally, we would initialize everything all at once in full generality, but
/// we reuse this object multiple times (maybe revisit this in the future?)
inline ContextPack new_ContextPack(
  FrozenSimpleArgs fwd_args,
  MainScratchBuf scratch_buf
) {
  ContextPack pack;
  for (int i = 0; i < 3; i++) {
    pack.grid_dimension[i] = 1;
    pack.grid_start[i] = 0;
    pack.grid_end[i] = 0;
  }
  gr_initialize_field_data(&pack.fields);
  pack.fields.grid_rank=3;
  pack.fields.grid_dimension = pack.grid_dimension;
  pack.fields.grid_start = pack.grid_start;
  pack.fields.grid_end = pack.grid_end;

  // precompute the 1-element index-range
  // - we explicitly follow the standard idiom for constructing an IndexRange
  //   and avoid directly constructing it (if the internals change we don't
  //   want to fix it here).
  const grackle_index_helper idx_helper = build_index_helper_(&pack.fields);
  pack.idx_range_1_element = make_idx_range_(0, &idx_helper);

  // initialize other members
  pack.fwd_args = fwd_args;
  pack.main_scratch_buf = scratch_buf;

  // nothing needs to be done for pack.other_scratch_buf
  return pack;
}

/// configure a ContextPack to be used at a particular index
///
/// @param[in,out] ptr pointer to an already initialized time_deriv_0d_pack
/// @param[in] my_fields the multidimensional field instance
/// @param[in] field_idx1d the index where we will be performing the calculation
///    (it has already been remapped from multiple dimensions to 1D).
/// @param[in] local_itmask_metal specifies the local values of itmask_metal
/// @param[in] local_edot_handling Specifies how we handle energy evolution.
///    This is read from an array historically known as `imp_eng`
inline void configure_ContextPack(
  ContextPack* pack, const grackle_field_data* my_fields,
  int field_idx1d, gr_mask_type local_itmask_metal, int local_edot_handling
) {
  pack->local_itmask_metal = local_itmask_metal;
  pack->local_edot_handling = local_edot_handling;
  pack->fields.grid_dx = my_fields->grid_dx;

  // here, we overwrite each field in pack.fields_1zone with pointers from
  // each field in my_fields corresponding to the current location (i,j,k)
  copy_offset_fieldmember_ptrs_(&pack->fields, my_fields, field_idx1d);

  // in the near future, we will overwrite the pointers of each
  // species-field (& the internal_energy) in pack.fields_1zone with a
  // pointer corresponding to the appropriate entry in `dsp`
}

/// here we copy the values from the scratch buffers (used by grackle's main
/// loop) at a single index into the corresponding buffers tracked by the
/// ContextPack.
///
/// @todo
/// We are copying more values than are strictly necessary (we are doing this
/// for historical consistency). We should be copying as few values as
/// possible. Honestly, creating this routine has made it clear that some of
/// the logic of the 0d time-derivative time calculation should probably be
/// reconsidered.
///
/// @todo
/// Once we remove scratchbuf_copy_from_pack (which does the inverse of this
/// function), and we have removed the unnecessary logic from this function,
/// this should be combined with configure_ContextPack
inline void scratchbufs_copy_into_pack(
  int index, ContextPack* pack, const double* p2d, const double* tgas,
  const double* tdust, const double* metallicity, const double* dust2gas,
  const double* rhoH, const double* mmw, const double* h2dust,
  const double* edot, grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
) {

  // it's a little unclear how many of the following copy-operations are
  // required by the current implementation
  // -> the need for copying depending may depend on whether internal_energy
  //    evolved with the species densities
  // -> we comment in special cases where there is an obvious choice

  // first, we copy the values into the buffers of pack->main_scratch_buf
  {
    // to help out, we define a lambda function (it captures index by value)
    auto copy_fn = [index](
      MemberInfo member_info, auto*& pack_buf, auto*& external_buf
    ) { pack_buf[0] = external_buf[index]; };

    // unclear if the current implementation depends on the following copy, but
    // it SHOULD never be necessary
    visit_member_pair(
      pack->main_scratch_buf.logTlininterp_buf, logTlininterp_buf, copy_fn
    );

    // I think the following case is necessary (if we aren't co-evolving
    // internal energy). But, it seems like we should be recomputing this in
    // all cases (these depend on the number density of dust, right?)
    visit_member_pair(
      pack->main_scratch_buf.grain_temperatures, grain_temperatures, copy_fn
    );

    // the tgasold buffer is definitely needed when the current implementation
    // co-evolves the internal energy
    visit_member_pair(
      pack->main_scratch_buf.cool1dmulti_buf, cool1dmulti_buf, copy_fn
    );

    // unclear whether the following cases need to be copied
    visit_member_pair(
      pack->main_scratch_buf.coolingheating_buf, coolingheating_buf, copy_fn
    );
    visit_member_pair(
      pack->main_scratch_buf.chemheatrates_buf, chemheatrates_buf, copy_fn
    );
  }

  // second, we copy the remaining values
  pack->other_scratch_buf.p2d[0] = p2d[index];
  // we may want to recalculate this regardless of whether we are co-evolving
  // internal-energy (since temperature is dependent on the species number
  // densities
  pack->other_scratch_buf.tgas[0] = tgas[index];
  pack->other_scratch_buf.tdust[0] = tdust[index];
  pack->other_scratch_buf.metallicity[0] = metallicity[index];
  pack->other_scratch_buf.dust2gas[0] = dust2gas[index];
  pack->other_scratch_buf.rhoH[0] = rhoH[index];
  pack->other_scratch_buf.mmw[0] = mmw[index];
  pack->other_scratch_buf.h2dust[0] = h2dust[index];
  pack->other_scratch_buf.edot[0] = edot[index];
}

/// here we copy the values into the scratch buffers (used by grackle's main
/// loop) at a single index from the corresponding buffers tracked by the
/// ContextPack.
///
/// @note
/// This function should **NOT** exist and should be removed (it only exists
/// right now as we pursue transcription). In particular, it makes no logical
/// sense to overwrite the value of cool1dmulti_buf.tgasold
inline void scratchbufs_copy_from_pack(
  int index, ContextPack* pack, double* p2d, double* tgas,
  double* tdust, double* metallicity, double* dust2gas, double* rhoH,
  double* mmw, double* h2dust, double* edot,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
) {

  // first, we copy the values from the buffers of pack->main_scratch_buf
  {
    // to help out, we define a lambda function (it captures index by value)
    auto copy_fn = [index](
      MemberInfo member_info, auto*& pack_buf, auto*& external_buf
    ) { external_buf[index] = pack_buf[0]; };

    visit_member_pair(
      pack->main_scratch_buf.logTlininterp_buf, logTlininterp_buf, copy_fn
    );
    visit_member_pair(
      pack->main_scratch_buf.grain_temperatures, grain_temperatures, copy_fn
    );
    visit_member_pair(
      pack->main_scratch_buf.cool1dmulti_buf, cool1dmulti_buf, copy_fn
    );
    visit_member_pair(
      pack->main_scratch_buf.coolingheating_buf, coolingheating_buf, copy_fn
    );
    visit_member_pair(
      pack->main_scratch_buf.chemheatrates_buf, chemheatrates_buf, copy_fn
    );
  }

  // second, we copy the remaining values
  p2d[index] = pack->other_scratch_buf.p2d[0];
  tgas[index] = pack->other_scratch_buf.tgas[0];
  tdust[index] = pack->other_scratch_buf.tdust[0];
  metallicity[index] = pack->other_scratch_buf.metallicity[0];
  dust2gas[index] = pack->other_scratch_buf.dust2gas[0];
  rhoH[index] = pack->other_scratch_buf.rhoH[0];
  mmw[index] = pack->other_scratch_buf.mmw[0];
  h2dust[index] = pack->other_scratch_buf.h2dust[0];
  edot[index] = pack->other_scratch_buf.edot[0];

}


/// calculate the time derivatives
///
/// @param[in] dt_FIXME Specifies the timestep passed to the
///   lookup_cool_rates1d_g function. See the C++ docstring for that function
///   for more details (this needs to be revisited)
/// @param[in] ycur A buffer representing a mathematical vector that holds the
///   initial values. This always has enough space for each species and the
///   internal energy (even if some entries are not evolved)
/// @param[out] ydot A buffer representing a mathematical vector that holds the
///   computed derivatives. This always has enough space for each species and
///   the internal energy (even if some entries are not evolved)
/// @param[in] pack Specifies extra buffers and information to be used in the
///   calculation
///
/// @note
/// In the future, we probably want to refactor in order to:
///   1. directly compute the jacobian while computing the derivative
///   2. accept an array of inputs and compute the derivatives at all values at
///      once (one might pick multiple simultaneous inputs for use in the
///      finite derivatives that are used to estimate the jacobian matrix)
///
/// @todo
/// we should replace `ycur` and `ydot` with instances of `SpeciesCollection`
/// and then pass `internal_energy` and `edot` as separate arguments.
///   - This will allow make the code easier to understand at a glance (and
///     will be important for reducing code duplication between this function
///     and `step_rate_g`)
///   - Furthermore, `step_rate_newton_raphson` already is responsible for
///     remapping for remapping to the `ycur` format. Effectively, we're only
///     tweaking the format and not producing any extra work.
void derivatives(
  double dt_FIXME, double* ycur, double* ydot, ContextPack& pack
) {

  chemistry_data* my_chemistry = pack.fwd_args.my_chemistry;
  chemistry_data_storage* my_rates = pack.fwd_args.my_rates;
  photo_rate_storage my_uvb_rates = pack.fwd_args.my_uvb_rates;
  InternalGrUnits internalu = pack.fwd_args.internalu;

  pack.other_scratch_buf.itmask[0] = MASK_TRUE;

  // todo: remove this variable
  int nsp = -1; // <- dummy value! (it isn't actually used)

  // some aliases to temporarily use
  double* dtit = &dt_FIXME;
  double* dsp = ycur;
  double* dspdot = ydot;

  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  // shorten `grackle::impl::chemistry` to `gr_chem` within this function
  namespace gr_chem = ::grackle::impl::chemistry;

  const int i_eng = 52;

  double comp1, comp2;  // in the future, these won't need to be passed to
                        // cool1d_multi_g

  // locals

  double atten, H2delta, h2heatfac, min_metallicity;

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // here, we are injecting some logic to help with refactoring

  // first, we declare some storage for our local copies of the species density
  // and internal energy AND we overwrite the pointers tracked by pack.fields
  // corresponding to these quantities
  // -> it would probably be better to do all of this ahead of time
  // -> it would also probably be better to not stack-allocate the storage,
  //    but it isn't much worse than what was here before
  gr_float rho_species[SpLUT::NUM_ENTRIES];
  gr_float e_internal[1];

  copy_contigSpTable_fieldmember_ptrs_(&pack.fields, rho_species, 1);
  pack.fields.internal_energy = &e_internal[0];

  // Now we store a local copy of all of the species (after casting).
  // -> we are remaining consistent with historical behavior in terms of
  //    copying. When gr_float == double, maybe we could reuse storage?
  // -> we aren't particular consistent with the historical behavior when
  //    gr_float == float. That's because the historical behavior was blatently
  //    wrong in that scenario (I doubt it was ever tested)
  SpeciesLUTFieldAdaptor field_lut_adaptor{pack.fields};
  for (int sp_idx = 0; sp_idx < SpLUT::NUM_ENTRIES; sp_idx++) {
    // TODO: somewhere, we need to have some logic to deal with edge cases when
    //       sizeof(gr_float) != sizeof(double)
    field_lut_adaptor.get_ptr_dynamic(sp_idx)[0] = static_cast<gr_float>(
      dsp[sp_idx]
    );
  }
  pack.fields.internal_energy[0] = static_cast<gr_float>(
    dsp[SpLUT::NUM_ENTRIES]
  );


  grackle_field_data* my_fields = &pack.fields;
  GRIMPL_REQUIRE(
    (
      (my_fields->grid_dimension[0] == 1) &&
      (my_fields->grid_dimension[1] == 1) &&
      (my_fields->grid_dimension[2] == 1)
    ),
    "sanity check!"
  );

  // Compute the cooling rate, tgas, tdust, and metallicity for this row

  if (pack.local_edot_handling == 1) {
    f_wrap::cool1d_multi_g(
      pack.fwd_args.imetal, pack.idx_range_1_element, pack.fwd_args.iter,
      pack.other_scratch_buf.edot, pack.other_scratch_buf.tgas,
      pack.other_scratch_buf.mmw, pack.other_scratch_buf.p2d,
      pack.other_scratch_buf.tdust, pack.other_scratch_buf.metallicity,
      pack.other_scratch_buf.dust2gas, pack.other_scratch_buf.rhoH,
      pack.other_scratch_buf.itmask, &pack.local_itmask_metal, my_chemistry,
      my_rates, &pack.fields,
      my_uvb_rates, internalu, pack.main_scratch_buf.grain_temperatures,
      pack.main_scratch_buf.logTlininterp_buf,
      pack.main_scratch_buf.cool1dmulti_buf,
      pack.main_scratch_buf.coolingheating_buf
    );
  }

  // uses the temperature to look up the chemical rates (they are interpolated
  // with respect to log temperature from input tables)
  f_wrap::lookup_cool_rates1d_g(
    pack.idx_range_1_element, pack.fwd_args.anydust,
    pack.other_scratch_buf.tgas, pack.other_scratch_buf.mmw,
    pack.other_scratch_buf.tdust, pack.other_scratch_buf.dust2gas,
    pack.main_scratch_buf.k13dd, pack.other_scratch_buf.h2dust,
    pack.fwd_args.dom, pack.fwd_args.dx_cgs,
    pack.fwd_args.c_ljeans, pack.other_scratch_buf.itmask,
    &pack.local_itmask_metal, pack.fwd_args.imetal,
    pack.other_scratch_buf.rhoH, dtit[0],
    my_chemistry, my_rates, &pack.fields, my_uvb_rates, internalu,
    pack.main_scratch_buf.grain_growth_rates,
    pack.main_scratch_buf.grain_temperatures,
    pack.main_scratch_buf.logTlininterp_buf,
    pack.main_scratch_buf.kcr_buf, pack.main_scratch_buf.kshield_buf,
    pack.main_scratch_buf.chemheatrates_buf
  );


  // The following function nominally computes dedot and HIdot (the time
  // derivatives in the mass densities of electrons and HI)
  // -> we don't care about these quantities (we recompute them later)
  // -> I'm pretty sure we care about how the following function also modifies
  //    edot
  if (pack.local_edot_handling == 1)  {

    f_wrap::rate_timestep_g(
      pack.other_scratch_buf.dedot, pack.other_scratch_buf.HIdot,
      pack.fwd_args.anydust, pack.idx_range_1_element,
      pack.other_scratch_buf.h2dust, pack.other_scratch_buf.rhoH,
      pack.other_scratch_buf.itmask, pack.other_scratch_buf.edot,
      pack.fwd_args.chunit, pack.fwd_args.dom, my_chemistry, &pack.fields,
      my_uvb_rates, pack.main_scratch_buf.kcr_buf,
      pack.main_scratch_buf.kshield_buf,
      pack.main_scratch_buf.chemheatrates_buf
    );
  }

  // Initialize

  std::memset(dspdot, 0, sizeof(double)*i_eng);


  // Heating/cooling rate (per unit volume -> gas mass)

  dspdot[i_eng-1] = *(pack.other_scratch_buf.edot) / *(pack.fields.density);

  // we are missing some logic!
  
  gr_chem::species_density_derivatives_0d(
    dspdot, pack.fwd_args.anydust, pack.other_scratch_buf.h2dust,
    pack.other_scratch_buf.rhoH, &pack.local_itmask_metal, my_chemistry,
    &pack.fields, my_uvb_rates, pack.main_scratch_buf.grain_growth_rates,
    pack.main_scratch_buf.kcr_buf, pack.main_scratch_buf.kshield_buf
  );
}


} // namespace grackle::impl::time_deriv_0d

#endif /* TIME_DERIV_0D_HPP */
