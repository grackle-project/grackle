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

/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

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

typedef struct
{
  int start;
  int end;
} grackle_index_range;

/***********************************************************************
/  
/ FUNCTION DECLARATIONS
/
************************************************************************/

// to help the compiler optimize the associated for-loops, this function:
//   - is implemented inline
//   - returns results as a struct rather than by modifying pointer arguments
static inline grackle_index_range _inner_range(int outer_index, 
                                               const grackle_index_helper* ind_helper)
{
  int k = (outer_index / ind_helper->num_j_inds) + ind_helper->k_start;
  int j = (outer_index % ind_helper->num_j_inds) + ind_helper->j_start;
  int outer_offset = ind_helper->i_dim * (j + ind_helper->j_dim * k);
  grackle_index_range out = {ind_helper->i_start + outer_offset,
                             ind_helper->i_end + outer_offset};
  return out;
}

grackle_index_helper _build_index_helper(const grackle_field_data *my_fields);

#endif
