//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the signature of the scale_fields_g function
///
//===----------------------------------------------------------------------===//

#ifndef SCALE_FIELDS_G_CPP_H
#define SCALE_FIELDS_G_CPP_H

#include "grackle.h"  // gr_float

namespace grackle::impl {

void scale_fields_g(int imetal, gr_float factor, chemistry_data* my_chemistry,
                    grackle_field_data* my_fields);

} // namespace grackle::impl

#endif /* SCALE_FIELDS_G_CPP_H */
