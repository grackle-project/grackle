/***********************************************************************
/
/ Declare utility functions related to the interpolation tables
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef INTERP_TABLE_UTILS_H
#define INTERP_TABLE_UTILS_H

#include "grackle.h" // gr_interp_grid, GRACKLE_CLOUDY_TABLE_MAX_DIMENSION
#include "grackle_macros.h" // GRACKLE_FREE

#ifdef __cplusplus
extern "C" {
#endif

/// Free memory associated with a #gr_interp_grid_props instance
static inline void free_interp_grid_props_(gr_interp_grid_props* props)
{
  for (int i = 0; i < GRACKLE_CLOUDY_TABLE_MAX_DIMENSION; i++) {
    GRACKLE_FREE(props->parameters[i]);
  }
}

/// Free memory associated with a #gr_interp_grid
static inline void free_interp_grid_(gr_interp_grid* grid)
{
  free_interp_grid_props_(&(grid->props));
  GRACKLE_FREE(grid->data);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* INTERP_TABLE_UTILS_H */
