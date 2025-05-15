// See LICENSE file for license and copyright information

/// @file solve_rate_cool_g-cpp.h
/// @brief Declares signature of ceiling_species_g

// This file was initially generated automatically during conversion of the
// ceiling_species_g function from FORTRAN to C++

#ifndef CEILING_SPECIES_CPP_H
#define CEILING_SPECIES_CPP_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

void ceiling_species_g(
  int* imetal, chemistry_data* my_chemistry, grackle_field_data* my_fields
);

#endif /* CEILING_SPECIES_CPP_H */
