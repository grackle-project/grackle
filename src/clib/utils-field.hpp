// See LICENSE file for license and copyright information

/// @file utils-field.hpp
/// @brief Declares utility routines related to grackle_field_data


#ifndef UTILS_FIELD_HPP
#define UTILS_FIELD_HPP

#include "grackle.h"

namespace grackle::impl {

/// This helper function essentially is used to make a kind of slice of an
/// instance of grackle_field_data
///
/// In more detail, each data-member of `dest` that hold spatial 3D data (i.e.
/// everything other than the "grid_*" members) is overwritten with:
/// - NULL if the corresponding pointer from src is NULL
/// - the pointer-addresses from src plus a constant offset.
inline void copy_offset_fieldmember_ptrs_(grackle_field_data *dest,
                                          const grackle_field_data *src,
                                          long long offset)
{
  // we use X-Macros to do this in 2 parts
  // - both of them will use the following logic helper macro
  #define GRIMPL_OFFSET_PTR_CPY(MEMBER_NAME)                          \
      dest->MEMBER_NAME = (src->MEMBER_NAME == NULL)                  \
      ? NULL : src->MEMBER_NAME + offset;

  // part 1: handle species-field members that grackle can evolve
  #define ENTRY(SPECIES_NAME) GRIMPL_OFFSET_PTR_CPY(SPECIES_NAME ## _density)
  #include "field_data_evolved_species.def"
  #undef ENTRY

  // part 2: modify all other field members
  #define ENTRY(MEMBER_NAME) GRIMPL_OFFSET_PTR_CPY(MEMBER_NAME)
  #include "field_data_misc_fdatamembers.def"
  #undef ENTRY

  #undef GRIMPL_OFFSET_PTR_CPY
}

} // namespace grackle::impl

#endif /* UTILS_FIELD_HPP */
