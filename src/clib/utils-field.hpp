// See LICENSE file for license and copyright information

/// @file utils-field.hpp
/// @brief Declares utility routines related to grackle_field_data

#ifndef UTILS_FIELD_HPP
#define UTILS_FIELD_HPP

#include "grackle.h"
#include "utils-cpp.hpp"

namespace grackle::impl {

/// maps a 3D index to the underlying 1D index of a contiguous memory buffer
/// "using a layout_left mapping policy." We assume 0-based indexing
///
/// "layout_left" is terminology adapted from the C++ standard:
/// - it refers to the idea that the leftmost extent/index has a stride of 1
///   and the sizes of strides increase as you go left to right
/// - this is the mapping used by default in fortran (note - we sill use
///   0-based indexing). It is sometimes called "column-major order"
///
/// @note
/// This may not be the best spot for this definition
inline int layoutleft_3D_index_to_1D_(const int* extent, int i, int j, int k) {
  return i + extent[0] * (j + extent[1] * k);
}

/// This helper function essentially is used to make a kind of slice of an
/// instance of grackle_field_data
///
/// In more detail, each data-member of `dest` that hold spatial 3D data (i.e.
/// everything other than the "grid_*" members) is overwritten with:
/// - NULL if the corresponding pointer from src is NULL
/// - the pointer-addresses from src plus a constant offset.
inline void copy_offset_fieldmember_ptrs_(grackle_field_data* dest,
                                          const grackle_field_data* src,
                                          long long offset) {
// we use X-Macros to do this in 2 parts
// - both of them will use the following logic helper macro
#define GRIMPL_OFFSET_PTR_CPY(MEMBER_NAME)                                     \
  dest->MEMBER_NAME =                                                          \
      (src->MEMBER_NAME == NULL) ? NULL : src->MEMBER_NAME + offset;

// part 1: handle species-field members that grackle can evolve
#define ENTRY(SPECIES_NAME) GRIMPL_OFFSET_PTR_CPY(SPECIES_NAME##_density)
#include "field_data_evolved_species.def"
#undef ENTRY

// part 2: modify all other field members
#define ENTRY(MEMBER_NAME) GRIMPL_OFFSET_PTR_CPY(MEMBER_NAME)
#include "field_data_misc_fdatamembers.def"
#undef ENTRY

#undef GRIMPL_OFFSET_PTR_CPY
}

/// this is an adaptor to support using SpLUT with grackle_field_data
struct SpeciesLUTFieldAdaptor {
  grackle_field_data data;

  /// lookup the pointer corresponding to the field-index
  gr_float* get_ptr_dynamic(int lut_idx) {
    switch (lut_idx) {
#define ENTRY(SPECIES_NAME)                                                    \
  case SpLUT::SPECIES_NAME:                                                    \
    return data.SPECIES_NAME##_density;
#include "field_data_evolved_species.def"
#undef ENTRY

      default:
        GRIMPL_ERROR("%d is not a valid SpLUT index", lut_idx);
    }
  }

  /// an alternative version of get_ptr_dynamic that has no runtime branching
  ///
  /// @note
  /// we are able to perform all branching at compile-time thanks to use of a
  /// template and if-constexpr
  template <int lut_idx>
  gr_float* get_ptr_static() {
    if constexpr (lut_idx < 0) {
      GRIMPL_ERROR("lut_idx can't be negative");
#define ENTRY(SPECIES_NAME)                                                    \
  }                                                                            \
  else if constexpr (lut_idx == SpLUT::SPECIES_NAME) {                         \
    return data.SPECIES_NAME##_density;
#include "field_data_evolved_species.def"
#undef ENTRY
    } else {
      GRIMPL_ERROR("lut_idx must be smaller than %s", SpLUT::NUM_ENTRIES);
    }
  }
};

}  // namespace grackle::impl

#endif /* UTILS_FIELD_HPP */
