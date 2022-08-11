/***********************************************************************
/
/ Define routines useful for indexing
/
/
/ Copyright (c) Enzo/Grackle Development Team. All rights reserved.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include "grackle_types.h"
#include "index_helper.h"

// to help the compiler optimize the associated for-loops, this function:
//   - is implemented inline
//   - returns results as a struct rather than by modifying pointer arguments
grackle_index_range _inner_range(int outer_index,
                                        const grackle_index_helper* ind_helper)
{
  int k = (outer_index / ind_helper->num_j_inds) + ind_helper->k_start;
  int j = (outer_index % ind_helper->num_j_inds) + ind_helper->j_start;
  int outer_offset = ind_helper->i_dim * (j + ind_helper->j_dim * k);
  grackle_index_range out = {ind_helper->i_start + outer_offset,
                             ind_helper->i_end + outer_offset};
  return out;
}

grackle_index_helper _build_index_helper(const grackle_field_data *my_fields)
{
  grackle_index_helper out;
  const int rank = my_fields->grid_rank;

  /* handle i indices */
  out.i_dim   = my_fields->grid_dimension[0];
  out.i_start = my_fields->grid_start[0];
  out.i_end   = my_fields->grid_end[0];

  /* handle j indices (j_end isn't tracked by grackle_index_helper) */
  out.j_dim   = (rank >= 2) ? my_fields->grid_dimension[1] : 1;
  out.j_start = (rank >= 2) ? my_fields->grid_start[1]     : 0;
  int j_end   = (rank >= 2) ? my_fields->grid_end[1]       : 0;
  out.num_j_inds = (j_end - out.j_start) + 1;

  /* handle k indices (k_end & k_dim aren't tracked by grackle_index_helper) */
  out.k_start = (rank >= 3) ? my_fields->grid_start[2]     : 0;
  int k_end   = (rank >= 3) ? my_fields->grid_end[2]       : 0;
  int num_k_inds = (k_end - out.k_start) + 1;

  out.outer_ind_size = num_k_inds * out.num_j_inds;
  return out;
}
