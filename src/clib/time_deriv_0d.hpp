// See LICENSE file for license and copyright information

/// @file time_deriv_0d.hpp
/// @brief Defines machinery to calculate the time derivative for a single zone
///
/// When we transcribe `lookup_cool_rates0d`, the plan is directly embed its
/// logic into the `grackle::impl::time_deriv_0d::derivatives` function.

#ifndef TIME_DERIV_0D_HPP
#define TIME_DERIV_0D_HPP

#include "fortran_func_wrappers.hpp"
#include "grackle.h"
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
  // todo: we can delete `iter`. This is only here for historical reasons
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
};

MainScratchBuf new_MainScratchBuf(void) {
  MainScratchBuf out;
  out.grain_temperatures = new_GrainSpeciesCollection(1);
  out.logTlininterp_buf = new_LogTLinInterpScratchBuf(1);
  out.cool1dmulti_buf = new_Cool1DMultiScratchBuf(1);
  out.coolingheating_buf = new_CoolHeatScratchBuf(1);
  out.chemheatrates_buf = new_ChemHeatingRates(1);
  return out;
}

void drop_MainScratchBuf(MainScratchBuf* ptr) {
  drop_GrainSpeciesCollection(&ptr->grain_temperatures);
  drop_LogTLinInterpScratchBuf(&ptr->logTlininterp_buf);
  drop_Cool1DMultiScratchBuf(&ptr->cool1dmulti_buf);
  drop_CoolHeatScratchBuf(&ptr->coolingheating_buf);
  drop_ChemHeatingRates(&ptr->chemheatrates_buf);
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
};

/// this struct is used to organize some temporary data that is used for
/// computing time derivatives of species data (and possibly, internal energy)
/// in a single zone
///
/// most of the data here acts a little like an adaptor layer
/// - we effectively adapt a representation of all species (and possibly
///   internal energy) from a vector form to the standard data structures to do
///   typical calculations and then we adapt back to the vector format
/// - to facillitate this, we effectively create an instance of
///   grackle_field_data that acts like a 1-element slice of the
///   grackle_field_data instance that the user passed in.
/// - this is highly inefficient, but it is logically consistent with the
///   original fortran code from before transcription. (We should refactor this
///   in the future after we finish transcription)
struct ContextPack {
  /// holds the local value of itmask_metal
  gr_mask_type local_itmask_metal;
  /// specifies how we handling edot in the calculation. This comes from
  /// an array historically known as imp_eng
  int local_edot_handling;

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
      const char* name, auto*& pack_buf, auto*& external_buf,
      const grackle::impl::visitor::BufLenSpec& spec
    ) {
      pack_buf[0] = external_buf[index];
    };

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
      const char* name, auto*& pack_buf, auto*& external_buf,
      const grackle::impl::visitor::BufLenSpec& spec
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

// we temporarily forward declare the following signature
// (in the future, the following function will be consolidated with the
// derivatives function)
void lookup_cool_rates0d(
  double* dtit, int* nsp, double* dsp, double* dspdot, gr_mask_type* itmask,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
  grackle::impl::time_deriv_0d::ContextPack pack
);

/// calculate the time derivatives
///
/// @param[in] dt_FIXME Specifies the timestep passed to the
///   lookup_cool_rates0d_g function. See the C++ docstring for that function
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
/// when we transcribe `lookup_cool_rates0d`, we should replace `ycur` and
/// `ydot` with instances of `SpeciesCollection` and then pass
/// `internal_energy` and `edot` as separate arguments.
///   - This will allow make the code easier to understand at a glance (and
///     will be important for reducing code duplication between this function
///     and `step_rate_g`)
///   - Furthermore, `step_rate_newton_raphson` already is responsible for
///     remapping for remapping to the `ycur` format. Effectively, we're only
///     tweaking the format and not producing any extra work.
void derivatives(
  double dt_FIXME, double* ycur, double* ydot, ContextPack& pack
) {

  // once we transcribe lookup_cool_rates0d, we will need to update the
  // appropriate members of pack.fields to point to entries of ycur

  chemistry_data* my_chemistry = pack.fwd_args.my_chemistry;
  chemistry_data_storage* my_rates = pack.fwd_args.my_rates;
  photo_rate_storage my_uvb_rates = pack.fwd_args.my_uvb_rates;
  InternalGrUnits internalu = pack.fwd_args.internalu;

  // todo: remove this variable
  gr_mask_type local_itmask_nr = MASK_TRUE;

  // todo: remove this variable
  int nsp = -1; // <- dummy value! (it isn't actually used)

  lookup_cool_rates0d(
    &dt_FIXME, &nsp, ycur, ydot, &local_itmask_nr, my_chemistry, my_rates,
    my_uvb_rates, internalu, pack
  );
}

void lookup_cool_rates0d(
  double* dtit, int* nsp, double* dsp, double* dspdot, gr_mask_type* itmask,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
  grackle::impl::time_deriv_0d::ContextPack pack
)
{

  // -------------------------------------------------------------------

  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  const int i_eng = 52;

  // these should not be re-allocated every time we enter this function...
  // (they should all be preallocated ahead of time!)
  double dedot[1];
  double HIdot[1];
  double k13dd[14];
  grackle::impl::CollisionalRxnRateCollection kcr_buf =
    grackle::impl::new_CollisionalRxnRateCollection(1);
  grackle::impl::PhotoRxnRateCollection kshield_buf =
    grackle::impl::new_PhotoRxnRateCollection(1);
  grackle::impl::GrainSpeciesCollection grain_growth_rates =
    grackle::impl::new_GrainSpeciesCollection(1);

  // locals

  double scoef, acoef;

  // we want to avoid directly constructing an IndexRange (if the internals
  // change we don't want to fix it here). But we should probably preconstruct
  // it ahead of time
  const grackle_index_helper idx_helper = build_index_helper_(&pack.fields);
  IndexRange idx_range = make_idx_range_(0, &idx_helper);
  
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

  // define some local variables carried over from the fortran version:
  // - the goal is to remove all of these by time we are done with cleanup
  // - originally we only conditionally defined some of these variables, but
  //   there honestly isn't any benefit to doing that (the memory is allocated
  //   already)
  // - originally, each of these variables was a stack allocated variable that
  //   held a copy of the corresponding entry in dsp. Now, these variables are
  //   references to casted copies taken from dsp
  // - we aren't being that consistent with the historic implementation when
  //   gr_float isn't the same as double (but I don't thi

  gr_float& de      = pack.fields.e_density[0];
  gr_float& HI      = pack.fields.HI_density[0];
  gr_float& HII     = pack.fields.HII_density[0];
  gr_float& HeI     = pack.fields.HeI_density[0];
  gr_float& HeII    = pack.fields.HeII_density[0];
  gr_float& HeIII   = pack.fields.HeIII_density[0];
  gr_float& HM      = pack.fields.HM_density[0];
  gr_float& H2I     = pack.fields.H2I_density[0];
  gr_float& H2II    = pack.fields.H2II_density[0];
  gr_float& DI      = pack.fields.DI_density[0];
  gr_float& DII     = pack.fields.DII_density[0];
  gr_float& HDI     = pack.fields.HDI_density[0];
  gr_float& DM      = pack.fields.DM_density[0];
  gr_float& HDII    = pack.fields.HDII_density[0];
  gr_float& HeHII   = pack.fields.HeHII_density[0];
  gr_float& CI      = pack.fields.CI_density[0];
  gr_float& CII     = pack.fields.CII_density[0];
  gr_float& CO      = pack.fields.CO_density[0];
  gr_float& CO2     = pack.fields.CO2_density[0];
  gr_float& OI      = pack.fields.OI_density[0];
  gr_float& OH      = pack.fields.OH_density[0];
  gr_float& H2O     = pack.fields.H2O_density[0];
  gr_float& O2      = pack.fields.O2_density[0];
  gr_float& SiI     = pack.fields.SiI_density[0];
  gr_float& SiOI    = pack.fields.SiOI_density[0];
  gr_float& SiO2I   = pack.fields.SiO2I_density[0];
  gr_float& CH      = pack.fields.CH_density[0];
  gr_float& CH2     = pack.fields.CH2_density[0];
  gr_float& COII    = pack.fields.COII_density[0];
  gr_float& OII     = pack.fields.OII_density[0];
  gr_float& OHII    = pack.fields.OHII_density[0];
  gr_float& H2OII   = pack.fields.H2OII_density[0];
  gr_float& H3OII   = pack.fields.H3OII_density[0];
  gr_float& O2II    = pack.fields.O2II_density[0];
  gr_float& Mg      = pack.fields.Mg_density[0];
  gr_float& Al      = pack.fields.Al_density[0];
  gr_float& S       = pack.fields.S_density[0];
  gr_float& Fe      = pack.fields.Fe_density[0];
  gr_float& MgSiO3  = pack.fields.MgSiO3_dust_density[0];
  gr_float& AC      = pack.fields.AC_dust_density[0];
  gr_float& SiM     = pack.fields.SiM_dust_density[0];
  gr_float& FeM     = pack.fields.FeM_dust_density[0];
  gr_float& Mg2SiO4 = pack.fields.Mg2SiO4_dust_density[0];
  gr_float& Fe3O4   = pack.fields.Fe3O4_dust_density[0];
  gr_float& SiO2D   = pack.fields.SiO2_dust_density[0];
  gr_float& MgO     = pack.fields.MgO_dust_density[0];
  gr_float& FeS     = pack.fields.FeS_dust_density[0];
  gr_float& Al2O3   = pack.fields.Al2O3_dust_density[0];
  gr_float& reforg  = pack.fields.ref_org_dust_density[0];
  gr_float& volorg  = pack.fields.vol_org_dust_density[0];
  gr_float& H2Oice  = pack.fields.H2O_ice_dust_density[0];
  // the next line isn't used (we can delete it)
  //gr_float& e       = pack.fields.internal_energy[0];

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
    // this is a hacky bugfix
    // -> we need my_local_iter to be 1. If it has any other value,
    //    `cool1d_multi_g` will assume that `tgasold` was previously
    //    initialized (it's not!) and try to use its contained value
    int my_local_iter = 1;
    f_wrap::cool1d_multi_g(
      pack.fwd_args.imetal, idx_range, my_local_iter,
      pack.other_scratch_buf.edot, pack.other_scratch_buf.tgas,
      pack.other_scratch_buf.mmw, pack.other_scratch_buf.p2d,
      pack.other_scratch_buf.tdust, pack.other_scratch_buf.metallicity,
      pack.other_scratch_buf.dust2gas, pack.other_scratch_buf.rhoH,
      itmask, &pack.local_itmask_metal, my_chemistry, my_rates, &pack.fields,
      my_uvb_rates, internalu, pack.main_scratch_buf.grain_temperatures,
      pack.main_scratch_buf.logTlininterp_buf,
      pack.main_scratch_buf.cool1dmulti_buf,
      pack.main_scratch_buf.coolingheating_buf
    );
  }

  // uses the temperature to look up the chemical rates (they are interpolated
  // with respect to log temperature from input tables)
  f_wrap::lookup_cool_rates1d_g(
    idx_range, pack.fwd_args.anydust, pack.other_scratch_buf.tgas,
    pack.other_scratch_buf.mmw, pack.other_scratch_buf.tdust,
    pack.other_scratch_buf.dust2gas, k13dd, pack.other_scratch_buf.h2dust,
    pack.fwd_args.dom, pack.fwd_args.dx_cgs, pack.fwd_args.c_ljeans,
    itmask, &pack.local_itmask_metal,
    pack.fwd_args.imetal, pack.other_scratch_buf.rhoH, dtit[0],
    my_chemistry, my_rates, &pack.fields, my_uvb_rates, internalu,
    grain_growth_rates, pack.main_scratch_buf.grain_temperatures,
    pack.main_scratch_buf.logTlininterp_buf, kcr_buf, kshield_buf,
    pack.main_scratch_buf.chemheatrates_buf
  );


  // The following function nominally computes dedot and HIdot (the time
  // derivatives in the mass densities of electrons and HI)
  // -> we don't care about these quantities (we recompute them later)
  // -> I'm pretty sure we care about how the following function also modifies
  //    edot
  if (pack.local_edot_handling == 1)  {

    f_wrap::rate_timestep_g(
      dedot, HIdot, pack.fwd_args.anydust, idx_range,
      pack.other_scratch_buf.h2dust, pack.other_scratch_buf.rhoH,
      itmask, pack.other_scratch_buf.edot, pack.fwd_args.chunit,
      pack.fwd_args.dom, my_chemistry, &pack.fields, my_uvb_rates,
      kcr_buf, kshield_buf, pack.main_scratch_buf.chemheatrates_buf
    );
  }

  // Initialize

  std::memset(dspdot, 0, sizeof(double)*i_eng);


  // Heating/cooling rate (per unit volume -> gas mass)

  dspdot[i_eng-1] = *(pack.other_scratch_buf.edot) / *(pack.fields.density);


  // TODO: Deal with the remaining logic in this function (~1200 lines of code)
  //   1. First, it should be moved into a different function!
  //   2. Second, we should replace all usage of dspdot (and hardcoded
  //      integers) to use grackle::impl::SpeciesCollection and SpLUT
  //      -> this will finish the migration of solve_rate_newton_raphson
  //         away from hardcoded integers to SpLUT
  //      -> this will aide with the next step
  //   3. Third, we should consider how to deduplicate this logic and most of
  //      the logic in step_rate_g
  //      -> I am pretty sure Gen performed a copy-and-paste when he wrote this
  //      -> it is totally unsustainable for us to maintain both implementations
  //      -> A common obstacle for unifying this logic is:
  //         - step_rate_g uses these calculations with a backward difference
  //           formula (see 3.2 of the grackle paper)
  //         - Here, we're computing the net time derivative.
  //         We can DEFINITELY overcome this obstacle!
  //      -> there's another obstacle that applies for unifying another part of
  //         the logic:
  //         - Parts (A), (B), and (C) in the original `step_rate_g` all apply
  //           partial updates (again, see 3.2 of the grackle paper)
  //         - notably, part (D) does not do this (but I'm not sure that there
  //           is a reasonable motivation for that choice)

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {




    // 1) HI

    scoef  = kcr_buf.data[CollisionalRxnLUT::k2][0]   *HII       *de;
    acoef  = kcr_buf.data[CollisionalRxnLUT::k1][0]   *de
           + kcr_buf.data[CollisionalRxnLUT::k57][0]   *HI
           + kcr_buf.data[CollisionalRxnLUT::k58][0]   *HeI       /4.
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(pack.fields.RT_HI_ionization_rate); }
    dspdot[2-1] = dspdot[2-1] + (scoef - acoef * HI);














    // 2) HII
    scoef  = kcr_buf.data[CollisionalRxnLUT::k1][0]   *HI    *de
           + kcr_buf.data[CollisionalRxnLUT::k57][0]   *HI    *HI
           + kcr_buf.data[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
           + kshield_buf.k24[0]   *HI;
    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)       *HI; }
    acoef  = kcr_buf.data[CollisionalRxnLUT::k2][0]   *de;
    dspdot[3-1] = dspdot[3-1] + (scoef - acoef * HII);
















    // 3) Electron density

    scoef = 0.
               + kcr_buf.data[CollisionalRxnLUT::k57][0]   *HI    *HI
               + kcr_buf.data[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
               + kshield_buf.k24[0]   *HI
               + kshield_buf.k25[0]   *HeII       /4.
               + kshield_buf.k26[0]   *HeI       /4.;

    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 0) )
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)        * HI
              + *(pack.fields.RT_HeI_ionization_rate)         * HeI         / 4.
              + *(pack.fields.RT_HeII_ionization_rate)        * HeII        / 4.; }
    if ( (my_chemistry->use_radiative_transfer == 1)  &&  ( my_chemistry->radiative_transfer_hydrogen_only == 1) )
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)        * HI; }



    acoef = -(kcr_buf.data[CollisionalRxnLUT::k1][0]   *HI             - kcr_buf.data[CollisionalRxnLUT::k2][0]   *HII
            + kcr_buf.data[CollisionalRxnLUT::k3][0]   *HeI       /4. -
         kcr_buf.data[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            + kcr_buf.data[CollisionalRxnLUT::k5][0]   *HeII       /4. -
         kcr_buf.data[CollisionalRxnLUT::k4][0]   *HeII       /4.);
    dspdot[1-1] = dspdot[1-1] + (scoef - acoef * de);





  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---




  // 4) HeI

  scoef  = kcr_buf.data[CollisionalRxnLUT::k4][0]   *HeII       *de;
  acoef  = kcr_buf.data[CollisionalRxnLUT::k3][0]   *de
               + kshield_buf.k26[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(pack.fields.RT_HeI_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    scoef = scoef +  4. * ( 0.
        + kcr_buf.data[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
        + kcr_buf.data[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.
       );
    acoef = acoef
        + kcr_buf.data[CollisionalRxnLUT::k148][0]    *   HII
        + kcr_buf.data[CollisionalRxnLUT::k149][0]    *   HII
        + kcr_buf.data[CollisionalRxnLUT::k150][0]    *  H2II        /  2.;
  }
  dspdot[4-1] = dspdot[4-1] + (scoef - acoef * HeI);


  // 5) HeII

  scoef  = kcr_buf.data[CollisionalRxnLUT::k3][0]   *HeI    *de
         + kcr_buf.data[CollisionalRxnLUT::k6][0]   *HeIII       *de
         + kshield_buf.k26[0]   *HeI;

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(pack.fields.RT_HeI_ionization_rate)       *HeI; }

  acoef  = kcr_buf.data[CollisionalRxnLUT::k4][0]   *de        + kcr_buf.data[CollisionalRxnLUT::k5][0]   *de
         + kshield_buf.k25[0];

  if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { acoef = acoef + *(pack.fields.RT_HeII_ionization_rate); }
  if (my_chemistry->primordial_chemistry > 3)  {
    acoef = acoef
        + kcr_buf.data[CollisionalRxnLUT::k151][0]    *    HI;
  }
  dspdot[5-1] = dspdot[5-1] + (scoef - acoef * HeII);


  // 6) HeIII

  scoef   = kcr_buf.data[CollisionalRxnLUT::k5][0]   *HeII    *de
          + kshield_buf.k25[0]   *HeII;
  if ((my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0))
      { scoef = scoef + *(pack.fields.RT_HeII_ionization_rate)        * HeII; }
  acoef   = kcr_buf.data[CollisionalRxnLUT::k6][0]   *de;
  dspdot[6-1] = dspdot[6-1] + (scoef - acoef * HeIII);





  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms




    // 1) HI
    scoef  =      kcr_buf.data[CollisionalRxnLUT::k2][0]    * HII        * de
           + 2.*kcr_buf.data[CollisionalRxnLUT::k13][0]   * HI         * H2I       /2.
           +      kcr_buf.data[CollisionalRxnLUT::k11][0]   * HII        * H2I       /2.
           + 2.*kcr_buf.data[CollisionalRxnLUT::k12][0]   * de         * H2I       /2.
           +      kcr_buf.data[CollisionalRxnLUT::k14][0]   * HM         * de
           +      kcr_buf.data[CollisionalRxnLUT::k15][0]   * HM         * HI
           + 2.*kcr_buf.data[CollisionalRxnLUT::k16][0]   * HM         * HII
           + 2.*kcr_buf.data[CollisionalRxnLUT::k18][0]   * H2II       * de       /2.
           +      kcr_buf.data[CollisionalRxnLUT::k19][0]   * H2II       * HM       /2.
           + 2.*kshield_buf.k31[0]      * H2I       /2.;

    acoef  =      kcr_buf.data[CollisionalRxnLUT::k1][0]    * de
           +      kcr_buf.data[CollisionalRxnLUT::k7][0]    * de
           +      kcr_buf.data[CollisionalRxnLUT::k8][0]    * HM
           +      kcr_buf.data[CollisionalRxnLUT::k9][0]    * HII
           +      kcr_buf.data[CollisionalRxnLUT::k10][0]   * H2II       /2.
           + 2.*kcr_buf.data[CollisionalRxnLUT::k22][0]   * std::pow(HI       ,2)
           +      kcr_buf.data[CollisionalRxnLUT::k57][0]   * HI
           +      kcr_buf.data[CollisionalRxnLUT::k58][0]   * HeI       /4.
           + kshield_buf.k24[0];

    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(pack.fields.RT_HI_ionization_rate); }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if ((my_chemistry->primordial_chemistry > 2) && (my_chemistry->radiative_transfer_HDI_dissociation > 0))  {
        scoef = scoef
          + *(pack.fields.RT_HDI_dissociation_rate)        * HDI       /3.0;
      }
      if ((my_chemistry->metal_chemistry == 1)  && 
          (pack.local_itmask_metal != MASK_FALSE))  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef
            + *(pack.fields.RT_OH_dissociation_rate)         * OH        /17.0
            + *(pack.fields.RT_H2O_dissociation_rate)        * H2O       /18.0;
        }
      }
    }

    if (pack.fwd_args.anydust != MASK_FALSE)  {
      if(pack.local_itmask_metal != MASK_FALSE   )  {
        acoef = acoef + 2. * *(pack.other_scratch_buf.h2dust)    * *(pack.other_scratch_buf.rhoH);
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcr_buf.data[CollisionalRxnLUT::k50][0]    * HII        * DI         / 2.
            + kcr_buf.data[CollisionalRxnLUT::k54][0]    * H2I        * DI         / 4.;
      acoef = acoef
            + kcr_buf.data[CollisionalRxnLUT::k51][0]    * DII        / 2.
            + kcr_buf.data[CollisionalRxnLUT::k55][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcr_buf.data[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k135][0]    *    HM        *    DI        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
          + kcr_buf.data[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k125][0]    *  HDII        /  3.
          + kcr_buf.data[CollisionalRxnLUT::k130][0]    *   DII        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k136][0]    *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k137][0]    *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k151][0]    *  HeII        /  4.
          + kcr_buf.data[CollisionalRxnLUT::k152][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (pack.local_itmask_metal != MASK_FALSE))  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcr_buf.data[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcr_buf.data[CollisionalRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcr_buf.data[CollisionalRxnLUT::kz26][0]    *    OH        *    CO        / 476.
          + kcr_buf.data[CollisionalRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcr_buf.data[CollisionalRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcr_buf.data[CollisionalRxnLUT::kz33][0]    *    OI        *    OH        / 272.
          + kcr_buf.data[CollisionalRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz36][0]    *   HII        *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz37][0]    *   CII        *    OH        / 204.
          + kcr_buf.data[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
          + kcr_buf.data[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
          + kcr_buf.data[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
          + kcr_buf.data[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 9.5
          + kcr_buf.data[CollisionalRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcr_buf.data[CollisionalRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::kz15][0]    *    CH        / 13.
          + kcr_buf.data[CollisionalRxnLUT::kz16][0]    *   CH2        / 14.
          + kcr_buf.data[CollisionalRxnLUT::kz17][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz18][0]    *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz19][0]    *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz27][0]    *    CI        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz30][0]    *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz39][0]    *   OII        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz43][0]    *  COII        / 28.;
    }
    dspdot[2-1] = dspdot[2-1] + (scoef - acoef * HI);













    // 2) HII

    scoef  =    kcr_buf.data[CollisionalRxnLUT::k1][0]     * HI        * de
           +    kcr_buf.data[CollisionalRxnLUT::k10][0]    * H2II       *HI       /2.
           +    kcr_buf.data[CollisionalRxnLUT::k57][0]    * HI        * HI
           +    kcr_buf.data[CollisionalRxnLUT::k58][0]    * HI        * HeI       /4.
           + kshield_buf.k24[0]   *HI;

    if (my_chemistry->use_radiative_transfer == 1)
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)        * HI; }

    acoef  =    kcr_buf.data[CollisionalRxnLUT::k2][0]     * de
           +    kcr_buf.data[CollisionalRxnLUT::k9][0]     * HI
           +    kcr_buf.data[CollisionalRxnLUT::k11][0]    * H2I       /2.
           +    kcr_buf.data[CollisionalRxnLUT::k16][0]    * HM
           +    kcr_buf.data[CollisionalRxnLUT::k17][0]    * HM;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcr_buf.data[CollisionalRxnLUT::k51][0]    * HI         * DII        / 2.
            + kcr_buf.data[CollisionalRxnLUT::k52][0]    * H2I        * DII        / 4.;
      acoef = acoef
            + kcr_buf.data[CollisionalRxnLUT::k50][0]    * DI         / 2.
            + kcr_buf.data[CollisionalRxnLUT::k53][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k129][0]    *    DI        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k134][0]    *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k148][0]    *   HeI        /  4.
          + kcr_buf.data[CollisionalRxnLUT::k149][0]    *   HeI        /  4.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (pack.local_itmask_metal != MASK_FALSE))  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::kz22][0]    *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz34][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz35][0]    *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz36][0]    *    O2        / 32.;
    }
    dspdot[3-1] = dspdot[3-1] + (scoef - acoef * HII);

    
    // 3) electrons:

    scoef =   kcr_buf.data[CollisionalRxnLUT::k8][0]    * HM        * HI
           +  kcr_buf.data[CollisionalRxnLUT::k15][0]   * HM        * HI
           +  kcr_buf.data[CollisionalRxnLUT::k17][0]   * HM        * HII
           +  kcr_buf.data[CollisionalRxnLUT::k57][0]   * HI        * HI
           +  kcr_buf.data[CollisionalRxnLUT::k58][0]   * HI        * HeI       /4.
    // 
           + kshield_buf.k24[0]   *HI
           + kshield_buf.k25[0]   *HeII    /4.
           + kshield_buf.k26[0]   *HeI    /4.;

    if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 0) )
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)        * HI
              + *(pack.fields.RT_HeI_ionization_rate)         * HeI      / 4.
              + *(pack.fields.RT_HeII_ionization_rate)        * HeII     / 4.; }
    if ( (my_chemistry->use_radiative_transfer == 1)  &&  (my_chemistry->radiative_transfer_hydrogen_only == 1) )
        { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)        * HI; }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if ((my_chemistry->metal_chemistry == 1)  && 
          (pack.local_itmask_metal != MASK_FALSE))  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(pack.fields.RT_CI_ionization_rate)        * CI       /12.0
            + *(pack.fields.RT_OI_ionization_rate)        * OI       /16.0;
        }
      }
    }

    acoef = - (kcr_buf.data[CollisionalRxnLUT::k1][0]    *HI           - kcr_buf.data[CollisionalRxnLUT::k2][0]   *HII
            +  kcr_buf.data[CollisionalRxnLUT::k3][0]    *HeI       /4. -
         kcr_buf.data[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            +  kcr_buf.data[CollisionalRxnLUT::k5][0]    *HeII       /4. -
         kcr_buf.data[CollisionalRxnLUT::k4][0]   *HeII       /4.
            +  kcr_buf.data[CollisionalRxnLUT::k14][0]   *HM
            -  kcr_buf.data[CollisionalRxnLUT::k7][0]    *HI
            -  kcr_buf.data[CollisionalRxnLUT::k18][0]   *H2II       /2.);
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcr_buf.data[CollisionalRxnLUT::k56][0]    * DI         * HM        / 2.;
      acoef = acoef
            - kcr_buf.data[CollisionalRxnLUT::k1][0]     * DI         / 2.
            + kcr_buf.data[CollisionalRxnLUT::k2][0]     * DII        / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k131][0]    *  HDII        /  3.
          + kcr_buf.data[CollisionalRxnLUT::k132][0]    *    DI        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k153][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (pack.local_itmask_metal != MASK_FALSE))  {
      scoef = scoef;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::kz44][0]    *   CII        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz45][0]    *   OII        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz46][0]    * H2OII        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz47][0]    * H2OII        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz48][0]    * H3OII        / 19.
          + kcr_buf.data[CollisionalRxnLUT::kz49][0]    * H3OII        / 19.
          + kcr_buf.data[CollisionalRxnLUT::kz50][0]    *  O2II        / 32.;
    }
    dspdot[1-1] = dspdot[1-1] + (scoef - acoef * de);


    // 7) H2

    scoef = 2.*(kcr_buf.data[CollisionalRxnLUT::k8][0]     * HM          * HI
          +       kcr_buf.data[CollisionalRxnLUT::k10][0]    * H2II        * HI       /2.
          +       kcr_buf.data[CollisionalRxnLUT::k19][0]    * H2II        * HM       /2.
          +       kcr_buf.data[CollisionalRxnLUT::k22][0]    * HI        * std::pow((HI       ),2.));
    acoef = ( kcr_buf.data[CollisionalRxnLUT::k13][0]   *HI        + kcr_buf.data[CollisionalRxnLUT::k11][0]   *HII
            + kcr_buf.data[CollisionalRxnLUT::k12][0]   *de        )
            + kshield_buf.k29[0]    + kshield_buf.k31[0];

    if (pack.fwd_args.anydust != MASK_FALSE)  {
      if(pack.local_itmask_metal != MASK_FALSE   )  {
        scoef = scoef + 2. * *(pack.other_scratch_buf.h2dust)    *
             HI        * *(pack.other_scratch_buf.rhoH);
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef + 2. * (
              kcr_buf.data[CollisionalRxnLUT::k53][0]    * HDI        * HII        / 3.
            + kcr_buf.data[CollisionalRxnLUT::k55][0]    * HDI        * HI         / 3.
               );
      acoef = acoef
            + kcr_buf.data[CollisionalRxnLUT::k52][0]    * DII        / 2.
            + kcr_buf.data[CollisionalRxnLUT::k54][0]    * DI         / 2.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (pack.local_itmask_metal != MASK_FALSE))  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcr_buf.data[CollisionalRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcr_buf.data[CollisionalRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
         );
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::kz20][0]    *    CI        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz21][0]    *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz23][0]    *    CH        / 13.
          + kcr_buf.data[CollisionalRxnLUT::kz24][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz40][0]    *   OII        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz41][0]    *  OHII        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz42][0]    * H2OII        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz51][0]    *    CI        / 12.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          scoef = scoef + 2. *
                grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      * 2.;

        }
        if (my_chemistry->dust_species > 1)  {
          scoef = scoef + 2. * (
                grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 3.
              + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       * 4.
              + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]
              + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       * 3.
            );
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      / H2I        * 2. * 2.;
        }
      }
    }
    dspdot[8-1] = dspdot[8-1] + (scoef - acoef * H2I);


    // 8) H-

    scoef = kcr_buf.data[CollisionalRxnLUT::k7][0]    * HI        * de;
    acoef = (kcr_buf.data[CollisionalRxnLUT::k8][0]     + kcr_buf.data[CollisionalRxnLUT::k15][0]   )  * HI        +
            (kcr_buf.data[CollisionalRxnLUT::k16][0]    + kcr_buf.data[CollisionalRxnLUT::k17][0]   )  * HII        +
            kcr_buf.data[CollisionalRxnLUT::k14][0]    * de        + kcr_buf.data[CollisionalRxnLUT::k19][0]    * H2II       /2.0f +
            my_uvb_rates.k27;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      acoef = acoef
            + kcr_buf.data[CollisionalRxnLUT::k56][0]    * DI         / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcr_buf.data[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k135][0]    *    DI        /  2.;
    }
    dspdot[7-1] = dspdot[7-1] + (scoef - acoef * HM);



    // 9) H2+

    scoef =    2.*( kcr_buf.data[CollisionalRxnLUT::k9][0]    *HI    *HII
                  +   kcr_buf.data[CollisionalRxnLUT::k11][0]   *H2I    /2.*HII
                  +   kcr_buf.data[CollisionalRxnLUT::k17][0]   *HM    *HII
                  + kshield_buf.k29[0]   *H2I    /2.
                  );
    acoef =         kcr_buf.data[CollisionalRxnLUT::k10][0]   *HI     + kcr_buf.data[CollisionalRxnLUT::k18][0]   *de
                  + kcr_buf.data[CollisionalRxnLUT::k19][0]   *HM
                  + (kshield_buf.k28[0]   +kshield_buf.k30[0]   );
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
         );
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k150][0]    *   HeI        /  4.;
    }
    dspdot[9-1] = dspdot[9-1] + (scoef - acoef * H2II);









  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {


    
    // 1) DI
    scoef =   (       kcr_buf.data[CollisionalRxnLUT::k2][0]    * DII        * de
               +      kcr_buf.data[CollisionalRxnLUT::k51][0]   * DII        * HI
               + 2.*kcr_buf.data[CollisionalRxnLUT::k55][0]   * HDI        *
            HI       /3.
               );
    acoef  =    kcr_buf.data[CollisionalRxnLUT::k1][0]    * de
           +    kcr_buf.data[CollisionalRxnLUT::k50][0]    * HII
           +    kcr_buf.data[CollisionalRxnLUT::k54][0]    * H2I       /2.
           +    kcr_buf.data[CollisionalRxnLUT::k56][0]    * HM
           + kshield_buf.k24[0];
    if (my_chemistry->use_radiative_transfer == 1) { acoef = acoef + *(pack.fields.RT_HI_ionization_rate); }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcr_buf.data[CollisionalRxnLUT::k133][0]    *   DII        *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcr_buf.data[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.
          );
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k129][0]    *   HII
          + kcr_buf.data[CollisionalRxnLUT::k132][0]    *    de
          + kcr_buf.data[CollisionalRxnLUT::k135][0]    *    HM;
    }
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        scoef = scoef
          + 2. * *(pack.fields.RT_HDI_dissociation_rate)        * HDI       /3.0;
      }
    }
    dspdot[10-1] = dspdot[10-1] + (scoef - acoef * DI);
                                                    

    // 2) DII
    scoef =   (   kcr_buf.data[CollisionalRxnLUT::k1][0]     * DI        * de
          +       kcr_buf.data[CollisionalRxnLUT::k50][0]    * HII       * DI
          +  2.*kcr_buf.data[CollisionalRxnLUT::k53][0]    * HII       * HDI       /3.
          )
          + kshield_buf.k24[0]   *DI;
    acoef = 0.;
    // ! initialize GC202002
    if (my_chemistry->use_radiative_transfer == 1) { scoef = scoef + *(pack.fields.RT_HI_ionization_rate)       *DI; }
    acoef =    kcr_buf.data[CollisionalRxnLUT::k2][0]     * de
          +    kcr_buf.data[CollisionalRxnLUT::k51][0]    * HI
          +    kcr_buf.data[CollisionalRxnLUT::k52][0]    * H2I       /2.;
    if (my_chemistry->primordial_chemistry > 3)  {
      acoef = acoef
          + kcr_buf.data[CollisionalRxnLUT::k130][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::k133][0]    *    DM        /  2.;
    }
    dspdot[11-1] = dspdot[11-1] + (scoef - acoef * DII);


    // 3) HDI
    scoef = 3.*(kcr_buf.data[CollisionalRxnLUT::k52][0]    * DII       *
         H2I       /2./2.
         + kcr_buf.data[CollisionalRxnLUT::k54][0]    * DI        * H2I       /2./2.
    // !   &           + 2._DKIND*kcr_buf.data[CollisionalRxnLUT::k56][0]    * DI        * HM       /2._DKIND
    //- ! corrected by GC202005
         +          kcr_buf.data[CollisionalRxnLUT::k56][0]    * DI        * HM       /2.
               );
    acoef  =    kcr_buf.data[CollisionalRxnLUT::k53][0]    * HII
           +    kcr_buf.data[CollisionalRxnLUT::k55][0]    * HI;
    if (my_chemistry->use_radiative_transfer == 1)  {
      if (my_chemistry->radiative_transfer_HDI_dissociation > 0)  {
        acoef = acoef
          + *(pack.fields.RT_HDI_dissociation_rate);
      }
    }
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  3. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.
          + kcr_buf.data[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.
          );
    }
    dspdot[12-1] = dspdot[12-1] + (scoef - acoef * HDI);




  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {



    // 1) DM

    scoef =
          kcr_buf.data[CollisionalRxnLUT::k132][0]    *    DI        *    de
        + kcr_buf.data[CollisionalRxnLUT::k135][0]    *    HM        *    DI;
    acoef =
          kcr_buf.data[CollisionalRxnLUT::k133][0]    *   DII        /  2.
        + kcr_buf.data[CollisionalRxnLUT::k134][0]    *   HII
        + kcr_buf.data[CollisionalRxnLUT::k136][0]    *    HI
        + kcr_buf.data[CollisionalRxnLUT::k137][0]    *    HI;

    dspdot[13-1] = dspdot[13-1] + (scoef - acoef * DM);


    // 2) HDII

    scoef = 3. * (
          kcr_buf.data[CollisionalRxnLUT::k129][0]    *    DI        *   HII        /  2.
        + kcr_buf.data[CollisionalRxnLUT::k130][0]    *   DII        *    HI        /  2.
       );
    acoef =
          kcr_buf.data[CollisionalRxnLUT::k125][0]    *    HI
        + kcr_buf.data[CollisionalRxnLUT::k131][0]    *    de;

    dspdot[14-1] = dspdot[14-1] + (scoef - acoef * HDII);


    // 3) HeHII

    scoef = 5. * (
          kcr_buf.data[CollisionalRxnLUT::k148][0]    *   HeI        *   HII        /  4.
        + kcr_buf.data[CollisionalRxnLUT::k149][0]    *   HeI        *   HII        /  4.
        + kcr_buf.data[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
        + kcr_buf.data[CollisionalRxnLUT::k151][0]    *  HeII        *    HI        /  4.
       );
    acoef =
          kcr_buf.data[CollisionalRxnLUT::k152][0]    *    HI
        + kcr_buf.data[CollisionalRxnLUT::k153][0]    *    de;

    dspdot[15-1] = dspdot[15-1] + (scoef - acoef * HeHII);




  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {

    if (pack.local_itmask_metal != MASK_FALSE   )  {

      // ***** CI **********
      scoef = 0. + 12. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz15][0]    *    HI        *    CH        / 13.
          + kcr_buf.data[CollisionalRxnLUT::kz44][0]    *   CII        *    de        / 12.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz20][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz27][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz28][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz29][0]    *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz51][0]    *   H2I        /  2.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][0]          / CI        * 12.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          acoef = acoef
            + *(pack.fields.RT_CI_ionization_rate);
        }
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef + 12. *
              *(pack.fields.RT_CO_dissociation_rate)         * CO        /28.0;
        }
      }

      dspdot[16-1] = dspdot[16-1] + (scoef - acoef * CI);



      // ***** CII **********
      scoef = 0. + 12. * ( 0.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz37][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz38][0]    *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz44][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(pack.fields.RT_CI_ionization_rate)        * CI;
        }
      }

      dspdot[17-1] = dspdot[17-1] + (scoef - acoef * CII);



      // ***** CO **********
      scoef = 0. + 28. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz28][0]    *    CI        *    OH        / 204.
          + kcr_buf.data[CollisionalRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcr_buf.data[CollisionalRxnLUT::kz32][0]    *    OI        *    CH        / 208.
          + kcr_buf.data[CollisionalRxnLUT::kz38][0]    *   CII        *    O2        / 384.
          + kcr_buf.data[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz26][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      / CO        * 17. * 0.5
          + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      / CO        * 17.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(pack.fields.RT_CO_dissociation_rate);
        }
      }

      dspdot[18-1] = dspdot[18-1] + (scoef - acoef * CO);



      // ***** CO2 **********
      scoef = 0. + 44. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz26][0]    *    OH        *    CO        / 476.
         );
      acoef = 0.;

      dspdot[19-1] = dspdot[19-1] + (scoef - acoef * CO2);



      // ***** OI **********
      scoef = 0. + 16. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz17][0]    *    HI        *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcr_buf.data[CollisionalRxnLUT::kz29][0]    *    CI        *    O2        / 384.
          + kcr_buf.data[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz45][0]    *   OII        *    de        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz50][0]    *  O2II        *    de        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz21][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz22][0]    *   HII
          + kcr_buf.data[CollisionalRxnLUT::kz30][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz31][0]    *    OI        / 8.
          + kcr_buf.data[CollisionalRxnLUT::kz32][0]    *    CH        / 13.
          + kcr_buf.data[CollisionalRxnLUT::kz33][0]    *    OH        / 17.;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          acoef = acoef
            + *(pack.fields.RT_OI_ionization_rate);
        }
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          scoef = scoef + 16. *
            ( *(pack.fields.RT_OH_dissociation_rate)         * OH        /17.0
            + *(pack.fields.RT_CO_dissociation_rate)         * CO        /28.0);
        }
      }

      dspdot[20-1] = dspdot[20-1] + (scoef - acoef * OI);



      // ***** OH **********
      scoef = 0. + 17. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz18][0]    *    HI        *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz19][0]    *    HI        *    O2        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcr_buf.data[CollisionalRxnLUT::kz30][0]    *    OI        *    HI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz17][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz24][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz25][0]    *    OH        / 8.5
          + kcr_buf.data[CollisionalRxnLUT::kz26][0]    *    CO        / 28.
          + kcr_buf.data[CollisionalRxnLUT::kz28][0]    *    CI        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz33][0]    *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz34][0]    *   HII
          + kcr_buf.data[CollisionalRxnLUT::kz37][0]    *   CII        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz52][0]    *   SiI        / 28.
          + kcr_buf.data[CollisionalRxnLUT::kz54][0]    *  SiOI        / 44.;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(pack.fields.RT_OH_dissociation_rate);
          scoef = scoef + 17. *
              *(pack.fields.RT_H2O_dissociation_rate)        * H2O       /18.0;
        }
      }

      dspdot[21-1] = dspdot[21-1] + (scoef - acoef * OH);



      // ***** H2O **********
      scoef = 0. + 18. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz24][0]    *   H2I        *    OH        / 34.
          + kcr_buf.data[CollisionalRxnLUT::kz25][0]    *    OH        *    OH        / 289.
          + kcr_buf.data[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz18][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz35][0]    *   HII;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / H2O        * 18. * 2.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / H2O        * 18. * 3.
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       / H2O        * 18. * 4.
          + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         / H2O        * 18.
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       / H2O        * 18. * 3.;
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][0]      / H2O        * 18.;
        }
      }
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_dissociation > 0)  {
          acoef = acoef
            + *(pack.fields.RT_H2O_dissociation_rate);
        }
      }

      dspdot[22-1] = dspdot[22-1] + (scoef - acoef * H2O);



      // ***** O2 **********
      scoef = 0. + 32. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz31][0]    *    OI        *    OI        / 256.
          + kcr_buf.data[CollisionalRxnLUT::kz33][0]    *    OI        *    OH        / 272.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz19][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz29][0]    *    CI        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz36][0]    *   HII
          + kcr_buf.data[CollisionalRxnLUT::kz38][0]    *   CII        / 12.
          + kcr_buf.data[CollisionalRxnLUT::kz53][0]    *   SiI        / 28.;

      dspdot[23-1] = dspdot[23-1] + (scoef - acoef * O2);



      // ***** SiI **********
      scoef = 0. + 28. * ( 0.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz52][0]    *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz53][0]    *    O2        / 32.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][0]         / SiI        * 28.;
        }
      }

      dspdot[24-1] = dspdot[24-1] + (scoef - acoef * SiI);



      // ***** SiOI **********
      scoef = 0. + 44. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz52][0]    *   SiI        *    OH        / 476.
          + kcr_buf.data[CollisionalRxnLUT::kz53][0]    *   SiI        *    O2        / 896.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz54][0]    *    OH        / 17.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / SiOI        * 44.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / SiOI        * 44.;
        }
      }

      dspdot[25-1] = dspdot[25-1] + (scoef - acoef * SiOI);



      // ***** SiO2I **********
      scoef = 0. + 60. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz54][0]    *  SiOI        *    OH        / 748.
         );
      acoef = 0.;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][0]       / SiO2I        * 60.;
        }
      }

      dspdot[26-1] = dspdot[26-1] + (scoef - acoef * SiO2I);



      // ***** CH **********
      scoef = 0. + 13. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz16][0]    *    HI        *   CH2        / 14.
          + kcr_buf.data[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcr_buf.data[CollisionalRxnLUT::kz27][0]    *    CI        *    HI        / 12.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz15][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz23][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz32][0]    *    OI        / 16.;

      dspdot[27-1] = dspdot[27-1] + (scoef - acoef * CH);



      // ***** CH2 **********
      scoef = 0. + 14. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz23][0]    *   H2I        *    CH        / 26.
          + kcr_buf.data[CollisionalRxnLUT::kz51][0]    *   H2I        *    CI        / 24.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz16][0]    *    HI;
      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      / CH2        * 14. * 0.5;
        }
      }

      dspdot[28-1] = dspdot[28-1] + (scoef - acoef * CH2);



      // ***** COII **********
      scoef = 0. + 28. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz37][0]    *   CII        *    OH        / 204.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz43][0]    *    HI;

      dspdot[29-1] = dspdot[29-1] + (scoef - acoef * COII);



      // ***** OII **********
      scoef = 0. + 16. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcr_buf.data[CollisionalRxnLUT::kz38][0]    *   CII        *    O2        / 384.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz39][0]    *    HI
          + kcr_buf.data[CollisionalRxnLUT::kz40][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz45][0]    *    de;
      if (my_chemistry->use_radiative_transfer == 1)  {
        if (my_chemistry->radiative_transfer_metal_ionization > 0)  {
          scoef = scoef
            + *(pack.fields.RT_OI_ionization_rate)        * OI;
        }
      }

      dspdot[30-1] = dspdot[30-1] + (scoef - acoef * OII);



      // ***** OHII **********
      scoef = 0. + 17. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz34][0]    *   HII        *    OH        / 17.
          + kcr_buf.data[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz41][0]    *   H2I        /  2.;

      dspdot[31-1] = dspdot[31-1] + (scoef - acoef * OHII);



      // ***** H2OII **********
      scoef = 0. + 18. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz35][0]    *   HII        *   H2O        / 18.
          + kcr_buf.data[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz42][0]    *   H2I        /  2.
          + kcr_buf.data[CollisionalRxnLUT::kz46][0]    *    de
          + kcr_buf.data[CollisionalRxnLUT::kz47][0]    *    de;

      dspdot[32-1] = dspdot[32-1] + (scoef - acoef * H2OII);



      // ***** H3OII **********
      scoef = 0. + 19. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz48][0]    *    de
          + kcr_buf.data[CollisionalRxnLUT::kz49][0]    *    de;

      dspdot[33-1] = dspdot[33-1] + (scoef - acoef * H3OII);



      // ***** O2II **********
      scoef = 0. + 32. * ( 0.
          + kcr_buf.data[CollisionalRxnLUT::kz36][0]    *   HII        *    O2        / 32.
         );
      acoef = 0.
          + kcr_buf.data[CollisionalRxnLUT::kz50][0]    *    de;

      dspdot[34-1] = dspdot[34-1] + (scoef - acoef * O2II);



      if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
        if (my_chemistry->dust_species > 0)  {
          // ***** Mg **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      / Mg        * 24.;
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / Mg        * 24. * 2.
            + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         / Mg        * 24.;
          }

          dspdot[35-1] = dspdot[35-1] + (scoef - acoef * Mg);


        }

        if (my_chemistry->dust_species > 1)  {
          // ***** Al **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       / Al        * 27. * 2.;

          dspdot[36-1] = dspdot[36-1] + (scoef - acoef * Al);



          // ***** S  **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         / S        * 32.;

          dspdot[37-1] = dspdot[37-1] + (scoef - acoef * S);



          // ***** Fe **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][0]         / Fe        * 56.
          + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       / Fe        * 56. * 3.
          + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         / Fe        * 56.;

          dspdot[38-1] = dspdot[38-1] + (scoef - acoef * Fe);


        }
      }

    }

  }

  // --- (D4) Now do dust species ---
  if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {

    if (pack.local_itmask_metal != MASK_FALSE   )  {

      if (my_chemistry->dust_species > 0)  {
        // ***** MgSiO3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::MgSiO3_dust][0]      * 100.;
        acoef = 0.;

        dspdot[39-1] = dspdot[39-1] + (scoef - acoef * MgSiO3);



        // ***** AC **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::AC_dust][0]          * 12.;
        acoef = 0.;

        dspdot[40-1] = dspdot[40-1] + (scoef - acoef * AC);


      }

      if (my_chemistry->dust_species > 1)  {
        // ***** SiM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::SiM_dust][0]         * 28.;
        acoef = 0.;

        dspdot[41-1] = dspdot[41-1] + (scoef - acoef * SiM);



        // ***** FeM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::FeM_dust][0]         * 56.;
        acoef = 0.;

        dspdot[42-1] = dspdot[42-1] + (scoef - acoef * FeM);



        // ***** Mg2SiO4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 140.;
        acoef = 0.;

        dspdot[43-1] = dspdot[43-1] + (scoef - acoef * Mg2SiO4);



        // ***** Fe3O4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Fe3O4_dust][0]       * 232.;
        acoef = 0.;

        dspdot[44-1] = dspdot[44-1] + (scoef - acoef * Fe3O4);



        // ***** SiO2D **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::SiO2_dust][0]       * 60.;
        acoef = 0.;

        dspdot[45-1] = dspdot[45-1] + (scoef - acoef * SiO2D);



        // ***** MgO **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::MgO_dust][0]         * 40.;
        acoef = 0.;

        dspdot[46-1] = dspdot[46-1] + (scoef - acoef * MgO);



        // ***** FeS **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::FeS_dust][0]         * 88.;
        acoef = 0.;

        dspdot[47-1] = dspdot[47-1] + (scoef - acoef * FeS);



        // ***** Al2O3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::Al2O3_dust][0]       * 102.;
        acoef = 0.;

        dspdot[48-1] = dspdot[48-1] + (scoef - acoef * Al2O3);


      }

      if (my_chemistry->dust_species > 2)  {
        // ***** reforg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::ref_org_dust][0]      * 22.68;
        acoef = 0.;

        dspdot[49-1] = dspdot[49-1] + (scoef - acoef * reforg);



        // ***** volorg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::vol_org_dust][0]      * 32.;
        acoef = 0.;

        dspdot[50-1] = dspdot[50-1] + (scoef - acoef * volorg);



        // ***** H2Oice **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates.data[OnlyGrainSpLUT::H2O_ice_dust][0]      * 18.;
        acoef = 0.;

        dspdot[51-1] = dspdot[51-1] + (scoef - acoef * H2Oice);


      }

    }

  }

  grackle::impl::drop_CollisionalRxnRateCollection(&kcr_buf);
  grackle::impl::drop_PhotoRxnRateCollection(&kshield_buf);
  grackle::impl::drop_GrainSpeciesCollection(&grain_growth_rates);

  return;
}


} // namespace grackle::impl::time_deriv_0d

#endif /* TIME_DERIV_0D_HPP */
