// See LICENSE file for license and copyright information

/// @file scale_fields_g-cpp.h
/// @brief Declares signature of scale_fields_g

// This file was initially generated automatically during conversion of the
// scale_fields_g function from FORTRAN to C++

#ifndef SCALE_FIELDS_G_CPP_H
#define SCALE_FIELDS_G_CPP_H

#include "grackle.h"  // gr_float

namespace grackle::impl {

void scale_fields_g(int imetal, gr_float* factor, chemistry_data* my_chemistry,
                    grackle_field_data* my_fields);

} // namespace grackle::impl

#endif /* SCALE_FIELDS_G_CPP_H */
