// See LICENSE file for license and copyright information

/// @file time_deriv_0d.hpp
/// @brief Defines machinery to calculate the time derivative for a single zone
///
/// When we transcribe `lookup_cool_rates0d`, the plan is directly embed its
/// logic into the `grackle::impl::time_deriv_0d::derivatives` function.

#ifndef TIME_DERIV_0D_HPP
#define TIME_DERIV_0D_HPP

#include "grackle.h"
#include "utils-field.hpp"
#include "internal_types.hpp"

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
    pack.grid_end[i] = 1;
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

  FORTRAN_NAME(lookup_cool_rates0d)(&dt_FIXME,
    pack.fields.density, pack.fields.x_velocity, pack.fields.y_velocity, pack.fields.z_velocity,
    &nsp, ycur, ydot, &my_chemistry->NumberOfTemperatureBins,
    &internalu.extfields_in_comoving, &my_chemistry->primordial_chemistry, &pack.fwd_args.imetal, &my_chemistry->metal_cooling,
    &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
    &my_chemistry->ih2co, &my_chemistry->ipiht, &pack.fwd_args.iter, &my_chemistry->photoelectric_heating,
    &internalu.a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
    &internalu.utem, &internalu.uxyz, &internalu.a_units, &internalu.urho, &internalu.tbase1,
    &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
    my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
    my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
    my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
    &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates.comp_xray, &my_uvb_rates.temp_xray,
    &my_uvb_rates.piHI, &my_uvb_rates.piHeI, &my_uvb_rates.piHeII,
    pack.fields.metal_density, pack.fields.dust_density,
    my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
    pack.main_scratch_buf.coolingheating_buf.hyd01k, pack.main_scratch_buf.coolingheating_buf.h2k01, pack.main_scratch_buf.coolingheating_buf.vibh, pack.main_scratch_buf.coolingheating_buf.roth, pack.main_scratch_buf.coolingheating_buf.rotl,
    my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, pack.main_scratch_buf.coolingheating_buf.gpldl, pack.main_scratch_buf.coolingheating_buf.gphdl,
    my_rates->HDlte, my_rates->HDlow, pack.main_scratch_buf.coolingheating_buf.hdlte, pack.main_scratch_buf.coolingheating_buf.hdlow,
    my_rates->GAHI, my_rates->GAH2, my_rates->GAHe,
    my_rates->GAHp, my_rates->GAel,
    my_rates->H2LTE, my_rates->gas_grain,
    pack.main_scratch_buf.coolingheating_buf.ceHI, pack.main_scratch_buf.coolingheating_buf.ceHeI, pack.main_scratch_buf.coolingheating_buf.ceHeII, pack.main_scratch_buf.coolingheating_buf.ciHI, pack.main_scratch_buf.coolingheating_buf.ciHeI,
    pack.main_scratch_buf.coolingheating_buf.ciHeIS, pack.main_scratch_buf.coolingheating_buf.ciHeII,
    pack.main_scratch_buf.coolingheating_buf.reHII, pack.main_scratch_buf.coolingheating_buf.reHeII1, pack.main_scratch_buf.coolingheating_buf.reHeII2, pack.main_scratch_buf.coolingheating_buf.reHeIII, pack.main_scratch_buf.coolingheating_buf.brem,
    pack.main_scratch_buf.logTlininterp_buf.indixe, pack.main_scratch_buf.logTlininterp_buf.t1, pack.main_scratch_buf.logTlininterp_buf.t2, pack.main_scratch_buf.logTlininterp_buf.logtem, pack.main_scratch_buf.logTlininterp_buf.tdef, pack.other_scratch_buf.edot,
    pack.other_scratch_buf.tgas, pack.main_scratch_buf.cool1dmulti_buf.tgasold, pack.other_scratch_buf.mmw, pack.other_scratch_buf.p2d, pack.other_scratch_buf.tdust, pack.other_scratch_buf.metallicity,
    pack.other_scratch_buf.dust2gas, pack.other_scratch_buf.rhoH, pack.main_scratch_buf.cool1dmulti_buf.mynh, pack.main_scratch_buf.cool1dmulti_buf.myde,
    pack.main_scratch_buf.cool1dmulti_buf.gammaha_eff, pack.main_scratch_buf.cool1dmulti_buf.gasgr_tdust, pack.main_scratch_buf.cool1dmulti_buf.regr,
    &my_chemistry->self_shielding_method, &my_uvb_rates.crsHI, &my_uvb_rates.crsHeI, &my_uvb_rates.crsHeII,
    &my_chemistry->use_radiative_transfer, &my_chemistry->radiative_transfer_hydrogen_only,
    &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, pack.main_scratch_buf.coolingheating_buf.cieco,
    &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
    &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
    my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
    &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
    &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
    my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
    &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
    &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, pack.fields.volumetric_heating_rate, pack.fields.specific_heating_rate,
    &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, pack.fields.temperature_floor,
    &my_chemistry->use_isrf_field, pack.fields.isrf_habing,
    &my_chemistry->three_body_rate, &pack.fwd_args.anydust, &my_chemistry->H2_self_shielding,
    my_rates->k1, my_rates->k2, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k7, my_rates->k8, my_rates->k9, my_rates->k10,
    my_rates->k11, my_rates->k12, my_rates->k13, my_rates->k13dd, my_rates->k14, my_rates->k15, my_rates->k16,
    my_rates->k17, my_rates->k18, my_rates->k19, my_rates->k22,
    &my_uvb_rates.k24, &my_uvb_rates.k25, &my_uvb_rates.k26, &my_uvb_rates.k27, &my_uvb_rates.k28, &my_uvb_rates.k29, &my_uvb_rates.k30, &my_uvb_rates.k31,
    my_rates->k50, my_rates->k51, my_rates->k52, my_rates->k53, my_rates->k54, my_rates->k55, my_rates->k56,
    my_rates->k57, my_rates->k58, &my_chemistry->NumberOfDustTemperatureBins, &my_chemistry->DustTemperatureStart, &my_chemistry->DustTemperatureEnd, my_rates->h2dust,
    my_rates->n_cr_n, my_rates->n_cr_d1, my_rates->n_cr_d2,
    pack.other_scratch_buf.h2dust, pack.main_scratch_buf.chemheatrates_buf.n_cr_n, pack.main_scratch_buf.chemheatrates_buf.n_cr_d1, pack.main_scratch_buf.chemheatrates_buf.n_cr_d2,
    &pack.fwd_args.dom, &internalu.coolunit, &internalu.tbase1, &internalu.xbase1, &pack.fwd_args.dx_cgs, &pack.fwd_args.c_ljeans,
    pack.fields.RT_HI_ionization_rate, pack.fields.RT_HeI_ionization_rate, pack.fields.RT_HeII_ionization_rate, pack.fields.RT_H2_dissociation_rate,
    pack.fields.RT_heating_rate, pack.fields.H2_self_shielding_length, &pack.fwd_args.chunit, &local_itmask_nr,
    &pack.local_itmask_metal,
     &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature,
     my_rates->k125, my_rates->k129, my_rates->k130, my_rates->k131, my_rates->k132,
     my_rates->k133, my_rates->k134, my_rates->k135, my_rates->k136, my_rates->k137,
     my_rates->k148, my_rates->k149, my_rates->k150, my_rates->k151, my_rates->k152,
     my_rates->k153,
     my_rates->kz15, my_rates->kz16, my_rates->kz17, my_rates->kz18, my_rates->kz19,
     my_rates->kz20, my_rates->kz21, my_rates->kz22, my_rates->kz23, my_rates->kz24,
     my_rates->kz25, my_rates->kz26, my_rates->kz27, my_rates->kz28, my_rates->kz29,
     my_rates->kz30, my_rates->kz31, my_rates->kz32, my_rates->kz33, my_rates->kz34,
     my_rates->kz35, my_rates->kz36, my_rates->kz37, my_rates->kz38, my_rates->kz39,
     my_rates->kz40, my_rates->kz41, my_rates->kz42, my_rates->kz43, my_rates->kz44,
     my_rates->kz45, my_rates->kz46, my_rates->kz47, my_rates->kz48, my_rates->kz49,
     my_rates->kz50, my_rates->kz51, my_rates->kz52, my_rates->kz53, my_rates->kz54,
     my_rates->cieY06,
     my_rates->LH2.props.dimension, &my_rates->LH2.props.data_size,
     my_rates->LH2.props.parameters[0], my_rates->LH2.props.parameters[1], my_rates->LH2.props.parameters[2],
     &my_rates->LH2.props.parameter_spacing[0], &my_rates->LH2.props.parameter_spacing[1], &my_rates->LH2.props.parameter_spacing[2], my_rates->LH2.data,
     my_rates->LHD.props.dimension, &my_rates->LHD.props.data_size,
     my_rates->LHD.props.parameters[0], my_rates->LHD.props.parameters[1], my_rates->LHD.props.parameters[2],
     &my_rates->LHD.props.parameter_spacing[0], &my_rates->LHD.props.parameter_spacing[1], &my_rates->LHD.props.parameter_spacing[2], my_rates->LHD.data,
     my_rates->LCI.props.dimension, &my_rates->LCI.props.data_size,
     my_rates->LCI.props.parameters[0], my_rates->LCI.props.parameters[1], my_rates->LCI.props.parameters[2],
     &my_rates->LCI.props.parameter_spacing[0], &my_rates->LCI.props.parameter_spacing[1], &my_rates->LCI.props.parameter_spacing[2], my_rates->LCI.data,
     my_rates->LCII.props.dimension, &my_rates->LCII.props.data_size,
     my_rates->LCII.props.parameters[0], my_rates->LCII.props.parameters[1], my_rates->LCII.props.parameters[2],
     &my_rates->LCII.props.parameter_spacing[0], &my_rates->LCII.props.parameter_spacing[1], &my_rates->LCII.props.parameter_spacing[2], my_rates->LCII.data,
     my_rates->LOI.props.dimension, &my_rates->LOI.props.data_size,
     my_rates->LOI.props.parameters[0], my_rates->LOI.props.parameters[1], my_rates->LOI.props.parameters[2],
     &my_rates->LOI.props.parameter_spacing[0], &my_rates->LOI.props.parameter_spacing[1], &my_rates->LOI.props.parameter_spacing[2], my_rates->LOI.data,
     my_rates->LCO.props.dimension, &my_rates->LCO.props.data_size,
     my_rates->LCO.props.parameters[0], my_rates->LCO.props.parameters[1], my_rates->LCO.props.parameters[2],
     &my_rates->LCO.props.parameter_spacing[0], &my_rates->LCO.props.parameter_spacing[1], &my_rates->LCO.props.parameter_spacing[2], my_rates->LCO.data,
     my_rates->LOH.props.dimension, &my_rates->LOH.props.data_size,
     my_rates->LOH.props.parameters[0], my_rates->LOH.props.parameters[1], my_rates->LOH.props.parameters[2],
     &my_rates->LOH.props.parameter_spacing[0], &my_rates->LOH.props.parameter_spacing[1], &my_rates->LOH.props.parameter_spacing[2], my_rates->LOH.data,
     my_rates->LH2O.props.dimension, &my_rates->LH2O.props.data_size,
     my_rates->LH2O.props.parameters[0], my_rates->LH2O.props.parameters[1], my_rates->LH2O.props.parameters[2],
     &my_rates->LH2O.props.parameter_spacing[0], &my_rates->LH2O.props.parameter_spacing[1], &my_rates->LH2O.props.parameter_spacing[2], my_rates->LH2O.data,
     my_rates->alphap.props.dimension, &my_rates->alphap.props.data_size,
     my_rates->alphap.props.parameters[0], my_rates->alphap.props.parameters[1], &my_rates->alphap.props.parameter_spacing[0], &my_rates->alphap.props.parameter_spacing[1],
     my_rates->alphap.data,
     &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
     pack.fields.local_ISM_metal_density, pack.fields.ccsn13_metal_density, pack.fields.ccsn20_metal_density,
     pack.fields.ccsn25_metal_density, pack.fields.ccsn30_metal_density, pack.fields.fsn13_metal_density,
     pack.fields.fsn15_metal_density, pack.fields.fsn50_metal_density, pack.fields.fsn80_metal_density,
     pack.fields.pisn170_metal_density, pack.fields.pisn200_metal_density, pack.fields.y19_metal_density,
     &my_rates->SN0_N,
     my_rates->SN0_fSiM, my_rates->SN0_fFeM, my_rates->SN0_fMg2SiO4, my_rates->SN0_fMgSiO3,
     my_rates->SN0_fFe3O4, my_rates->SN0_fAC, my_rates->SN0_fSiO2D, my_rates->SN0_fMgO,
     my_rates->SN0_fFeS, my_rates->SN0_fAl2O3,
     my_rates->SN0_freforg, my_rates->SN0_fvolorg, my_rates->SN0_fH2Oice,
     my_rates->SN0_r0SiM, my_rates->SN0_r0FeM, my_rates->SN0_r0Mg2SiO4, my_rates->SN0_r0MgSiO3,
     my_rates->SN0_r0Fe3O4, my_rates->SN0_r0AC, my_rates->SN0_r0SiO2D, my_rates->SN0_r0MgO,
     my_rates->SN0_r0FeS, my_rates->SN0_r0Al2O3,
     my_rates->SN0_r0reforg, my_rates->SN0_r0volorg, my_rates->SN0_r0H2Oice,
     my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td,
     my_rates->SN0_kpSiM, my_rates->SN0_kpFeM, my_rates->SN0_kpMg2SiO4, my_rates->SN0_kpMgSiO3,
     my_rates->SN0_kpFe3O4, my_rates->SN0_kpAC, my_rates->SN0_kpSiO2D, my_rates->SN0_kpMgO,
     my_rates->SN0_kpFeS, my_rates->SN0_kpAl2O3,
     my_rates->SN0_kpreforg, my_rates->SN0_kpvolorg, my_rates->SN0_kpH2Oice,
     my_rates->h2dustS, my_rates->h2dustC, my_rates->grain_growth_rate,
     pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
     pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::AC_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust],
     pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], pack.main_scratch_buf.grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust],
     my_rates->gas_grain2, &my_rates->gamma_isrf2,
     &pack.local_edot_handling,
     &my_chemistry->radiative_transfer_HDI_dissociation, pack.fields.RT_HDI_dissociation_rate, &my_chemistry->radiative_transfer_metal_ionization, pack.fields.RT_CI_ionization_rate, pack.fields.RT_OI_ionization_rate,
     &my_chemistry->radiative_transfer_metal_dissociation, pack.fields.RT_CO_dissociation_rate, pack.fields.RT_OH_dissociation_rate, pack.fields.RT_H2O_dissociation_rate,
     &my_chemistry->radiative_transfer_use_H2_shielding, &my_chemistry->H2_custom_shielding, pack.fields.H2_custom_shielding_factor
    );

}

} // namespace grackle::impl::time_deriv_0d

#endif /* TIME_DERIV_0D_HPP */
