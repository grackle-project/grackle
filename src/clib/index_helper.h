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

/// Specifies the range of indices for grackle's 3D fields.
///
/// Essentially, you use this struct to easily get the 3D indices for elements
/// in an "islice".
/// - Grackle commonly uses buffers that are used to hold temporary variable
///   that you just use the i-indices to access
/// - If you want to use a 3D index to get an element from one of grackle's 3D
///   field, special logic is required to remap the 3 indices to a single index
///   in the underlying flat 1d buffer that holds a field's data. In some cases
///   (e.g. Fortran Arrays) the remapping is done behind the scenes and in
///   cases, it is more explicit (but the logic is always somewhere)
/// - For context, the `field_flat_index_range` can be used to specify the same
///   range for 3D fields, but the remapping logic is pre-applied (this makes
///   it much less useful when you have these custom buffers)
///
/// @par Role in Transcription
/// This is intended to aide the transcription process of fortran subroutines
/// - prior to transcription, the 1-based `j` and `k` index values are passed
///   throughout the call-stack. Furthermore, the 0-based `i_start` & `i_end`
///   values are passed around.
/// - the automatic transcription process ends up making the indexing look a
///   little messy when transcribing from Fortran to C++ (in order to account
///   for the difference between 1-based & 0-based indexing). Thus gradual,
///   manual intervention is required to make indexing appear idiomatic.
/// - the manual process introduces a lot of opportunities for error if we
///   aren't extremely careful with keeping track of whether a local variable
///   still refers to a 1-based index or been converted to a 0-based index.
///   By using this struct, we can be very explicit and certain about whether
///   an expressions refers to a 1-based index or a 0-based index'
///
/// @note
/// Naively, this all might sound like it involves a bunch of extra work, but
/// the following it will work really well with the transcription machinery:
/// - Consider a fortran subroutine, `my_subroutine`, where one of the
///   arguments is the `j` index (which is expected to use 1-based indexing)
///   and the signature looks like `void my_subroutine(..., int* j, ...)`
/// - Now suppose we have a C/C++ function that calls `my_subroutine` and that
///   there is an `IndexRange` in the function called `idx_range`.
/// - If the call to `my_subroutine` in the C/C++ function looks like
///   `my_subroutine(..., &idx_range.jp1, ...)`, then when applying the
///   transcription tools to `my_subroutine` we can have the transcription
///   machinery replace every reference to `j` in the transcribed function
///   with `idx_range.jp1` (it will also modify the argument list properly).
///
/// @note
/// ultimately, it will be more idiomatic to make for-loops use `i_stop` rather
/// than `i_end`. Obviously, this is a **very low** priority (but this is why
/// we put both variables in this struct)
///
/// @par Future Usage
/// If we continue using this type after we complete transcription, we can have
/// it take on a prominent role in adding GPU-support.
typedef struct IndexRange
{
  // specifies the fixed (0-based) j and k indices to be used with the range
  int j;
  int k;

  // j & k indices for 1-based indexing (remove when transcription is done)
  int jp1;
  int kp1;

  // specifies the starting and stopping (0-based) indices specifying the range
  // along the i-axis (the stopping index isn't included in the range)
  int i_start;
  int i_stop;

  // holds the value of `i_stop-1` (remove when transcription is done)
  int i_end;
} IndexRange;

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

/// constructs an IndexRange, which holds the 3D index information for an
/// "islice."
static inline IndexRange make_idx_range_(
  int outer_index, const grackle_index_helper* idx_helper
)
{
  IndexRange out;
  out.k = (outer_index / idx_helper->num_j_inds) + idx_helper->k_start;
  out.j = (outer_index % idx_helper->num_j_inds) + idx_helper->j_start;
  out.jp1 = out.j+1;
  out.kp1 = out.k+1;
  out.i_start = idx_helper->i_start;
  out.i_end = idx_helper->i_end;
  out.i_stop = idx_helper->i_end+1;
  return out;
}

grackle_index_helper build_index_helper_(const grackle_field_data *my_fields);

#ifdef __cplusplus
} // extern "C"
#endif /* __cplusplus */

#endif
