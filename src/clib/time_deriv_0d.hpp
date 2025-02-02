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

// we choose to adopt a longer, more descriptive namespace here so that the
// handful of functions defined in this file can have shorter names (in the
// future, if we are willing to define methods on a struct, we can definitely
// shorten the namespace name)
namespace grackle::impl::time_deriv_0d {

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

  /// the following 3 arrays are used within fields
  int grid_dimension[3];
  int grid_start[3];
  int grid_end[3];
};

typedef struct ContextPack ContextPack;

/// partially set up a new time_deriv_pack
///
/// Ideally, we would initialize everything all at once in full generality, but
/// we reuse this object multiple times (maybe revisit this in the future?)
inline ContextPack new_ContextPack(void) {
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
