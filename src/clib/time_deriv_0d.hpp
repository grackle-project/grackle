// See LICENSE file for license and copyright information

/// @file time_deriv_0d.hpp
/// @brief Defines machinery to calculate the time derivative for a single zone
///
/// The plan is to put lookup_cool_rates0d in this file (and to eventually
/// give it a name more consistent with the other functionality)

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
  return pack;
}

/// configure a ContextPack to be used at a particular index
///
/// @param[in,out] ptr pointer to an already initialized time_deriv_0d_pack
/// @param[in] my_fields the multidimensional field instance
/// @param[in] field_idx1d the index where we will be performing the calculation
///    (it has already been remapped from multiple dimensions to 1D).
inline void configure_ContextPack(ContextPack* pack,
                                  const grackle_field_data* my_fields,
                                  int field_idx1d) {
  pack->fields.grid_dx = my_fields->grid_dx;

  // here, we overwrite each field in pack.fields_1zone with pointers from
  // each field in my_fields corresponding to the current location (i,j,k)
  copy_offset_fieldmember_ptrs_(&pack->fields, my_fields, field_idx1d);

  // in the near future, we will overwrite the pointers of each
  // species-field (& the internal_energy) in pack.fields_1zone with a
  // pointer corresponding to the appropriate entry in `dsp`
}

} // namespace grackle::impl::time_deriv_0d

#endif /* TIME_DERIV_0D_HPP */
