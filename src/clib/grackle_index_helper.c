/***********************************************************************
/
/ Initialize chemistry and cooling rate data
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

grackle_index_helper _build_index_helper(const grackle_field_data *my_fields)
{
  grackle_index_helper out;
  const int rank = my_fields->grid_rank;

  /* handle i indices */
  out.i_dim = my_fields->grid_dimension[0];
  out.i_start = my_fields->grid_start[0];
  out.i_end = my_fields->grid_end[0];

  /* handle j indices */
  out.j_dim   = (rank >= 2) ? my_fields->grid_dimension[1] : 1;
  out.j_start = (rank >= 2) ? my_fields->grid_start[1]     : 0;
  out.j_end   = (rank >= 2) ? my_fields->grid_end[1]       : 0;

  /* handle k indices */
  out.k_dim   = (rank >= 3) ? my_fields->grid_dimension[2] : 1;
  out.k_start = (rank >= 3) ? my_fields->grid_start[2]     : 0;
  out.k_end   = (rank >= 3) ? my_fields->grid_end[2]       : 0;

  out.num_j_inds = (out.j_end - out.j_start) + 1;
  int num_k_inds = (out.k_end - out.k_start) + 1;

  out.outer_ind_size = num_k_inds * out.num_j_inds;
  return out;
}
