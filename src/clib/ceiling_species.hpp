//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file ceiling_species.hpp
/// @brief Declares signature of ceiling_species
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// ceiling_species_g function from FORTRAN to C++

#ifndef CEILING_SPECIES_CPP_H
#define CEILING_SPECIES_CPP_H

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int

namespace grackle::impl {

void ceiling_species(int imetal, chemistry_data* my_chemistry,
                     grackle_field_data* my_fields);

} // namespace grackle::impl

#endif /* CEILING_SPECIES_CPP_H */
