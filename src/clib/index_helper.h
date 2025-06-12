/***********************************************************************
/
/ Declare Grackle's index_helper type and associated functions. This is only
/ intended to be used internally.
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_PRIVATE_H_
#define __GRACKLE_PRIVATE_H_

#include "grackle_types.h" // grackle_field_data

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

/// internal type to assist with iterating over 3D index-space.
///
/// An instance, `obj`, is commonly used in a 2-level nested for-loop:
/// - the outer loop iterates over index `t` for `0 <= t < obj.outer_ind_size`.
///   In 3D, this loop corresponds to `j`,`k` index pairs.
/// - the inner loop iterate over the "index-range" (constructed from `obj`
///   and `t`). In 3D, this loop corresponds to the `i` axis.
typedef struct
{
  int i_start;
  int i_end;
  int i_dim;

  int j_start;
  int j_dim;
  /* j_end isn't needed */

  int k_start;
  /* k_end & k_dim aren't needed */

  int num_j_inds;
  int outer_ind_size;

} grackle_index_helper;

/// Specifies a range of indices for grackle's 3D fields, for use when you 
/// treat the fields as flattened 1d arrays
typedef struct
{
  int start;
  int end;
} field_flat_index_range;

/***********************************************************************
/  
/ FUNCTION DECLARATIONS
/
************************************************************************/

// to help the compiler optimize the associated for-loops, this function:
//   - is implemented inline
//   - returns results as a struct rather than by modifying pointer arguments
static inline field_flat_index_range inner_flat_range_(
  int outer_index, const grackle_index_helper* ind_helper
)
{
  int k = (outer_index / ind_helper->num_j_inds) + ind_helper->k_start;
  int j = (outer_index % ind_helper->num_j_inds) + ind_helper->j_start;
  int outer_offset = ind_helper->i_dim * (j + ind_helper->j_dim * k);
  field_flat_index_range out = {ind_helper->i_start + outer_offset,
                                ind_helper->i_end + outer_offset};
  return out;
}

grackle_index_helper build_index_helper_(const grackle_field_data *my_fields);

#ifdef __cplusplus
} // extern "C"
#endif /* __cplusplus */

#endif
