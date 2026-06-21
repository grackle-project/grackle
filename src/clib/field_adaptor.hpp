//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares/Defines logic pertaining to adapting grackle_field_data to a more
/// convenient internal format
///
//===----------------------------------------------------------------------===//

#ifndef FIELD_ADAPTOR_HPP
#define FIELD_ADAPTOR_HPP

#include <cstddef>
#include <vector>
#include "grackle.h"
#include "LUT.hpp"
#include "support/index_helper.hpp"
#include "support/config.hpp"
#include "support/status_reporting.hpp"
#include "support/View.hpp"

namespace GRIMPL_NAMESPACE_DECL {

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

  // Part 3: modify inject pathway density field slots
  for (int i = 0; i < GRIMPL_MAX_INJ_PATHWAYS; i++) {
    gr_float* p = src->inject_pathway_metal_density[i];
    dest->inject_pathway_metal_density[i] =
        (p == nullptr) ? nullptr : p + offset;
  }

#undef GRIMPL_OFFSET_PTR_CPY
}

/// This helper function is used to store store pointers corresponding to
/// contiguous chunks of a species density table within pointers of a
/// grackle_field_data struct.
///
/// The underlying assumption is that the species table is organized as a
/// (nelem_elem_per_species, SpLUT::NUM_ENTRIES) array, where the values of
/// a given species are all contiguous
///
/// @note
/// This should only be used as a temporary measure during initial
/// transcription. In the long-term, the full grackle_field_data struct is not
/// well suited to be the common container that all chemistry functions operate
/// on.
inline void copy_contigSpTable_fieldmember_ptrs_(grackle_field_data* my_fields,
                                                 gr_float* species_table,
                                                 long long nelem_per_species) {
  GRIMPL_REQUIRE(nelem_per_species > 0,
                 "The number of elements per species must exceed 0");

#define ENTRY(SPECIES_NAME)                                                    \
  my_fields->SPECIES_NAME##_density =                                          \
      (&species_table[nelem_per_species * SpLUT::SPECIES_NAME]);
#include "field_data_evolved_species.def"
#undef ENTRY
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

/// @brief Acts as 2D view of all dynamically evolved species data.
///
/// In more detail:
/// - The first axis (i.e. the contiguous, fast axis), corresponds to the
///   spatial i index from an associated @ref IndexRange object. In grackle
///   configurations that handling multidimensional grids, the index offset for
///   the other axes (i.e. the j-axis and k-axis) should already be baked into
///   an object.
/// - The second axis (i.e. the slow axis) has an index for each dynamically
///   evolved Species.
///
/// Motivation
/// ----------
/// The plan is to transition to accessing dynamically evolved species data in
/// Grackle's core calculation routines through instances of this type (rather
/// than directly using @ref grackle_field_data).
///
/// We have explicitly chosen to create a type alias, rather than use
/// @brief Multi1DView to make it easier to update the underlying type in the
/// future. This might happen if we decide to copy the species data into a
/// contiguous 2D array (there are tradeoffs to such a choice, but it seeks
/// like a net positive).
///
/// Other Thoughts
/// --------------
/// In the future we may want to consider renaming this type-alias and consider
/// tracking other dynamically evolved quantities within instances of this
/// type. At this time, specific internal energy is the only other quantity
/// (in the near future, it may contain density as we add better dust networks)
template <typename T>
using SpeciesMultiView = Multi1DView<T>;

/// @brief This manages logic (and data) associated with adapting the
///        information within @ref grackle_field_data to a more convenient
///        format
///
/// Right now, a @ref FieldAdaptorManager instance is constructed so that it
/// tracks a pointer to a @ref grackle_field_data pointer. For the duration of
/// the instance's lifetime, that pointer must remain valid and none of its
/// data should be directly mutated (other than values within the pointers of
/// dynamically evolved fields). For that reason, instances of this type can't
/// persist between Grackle API calls (i.e. if we create an instance during an
/// API call, we destroy the instance after completing the associated work).
///
/// The impetus for creating this type was to have a single place where we
/// could stick the memory management associated with and the logic for
/// creating a @ref SpeciesMultiView.
///
/// The role of this type will probably evolve with time. For example, one
/// could imagine creating a @ref Multi1DView dedicated to tracking RT fields
/// or perhaps a @ref Multi1DView dedicated to any/all non-dynamically evolved
/// fields. In the future, (e.g. in a 4.0 release) we might restructure the
/// internals of @ref grackle_field_data and this type may need to do less work
class FieldAdaptorManager {
  // the field_data that is wrapped
  const grackle_field_data* field_data_;

  // we'll probably need to tailor the precise memory management strategies to
  // the backend, so we'll start out simple
  std::vector<gr_float*> sp_arr_of_ptrs_;

public:
  /// Override the wrapped @ref grackle_field_data pointer
  ///
  /// This is pretty dangerous to use:
  /// - challenging, difficult-to-understand bugs will arise if you try to
  ///   access a @ref SpeciesMultiView that was created before you called this
  ///   method.
  /// - this **ONLY** exists as an optimization for the case where you would
  ///   overwrite one instance with a new instance (i.e. it lets avoid a heap
  ///   allocation). We should remove this (or make it private) when we stop
  ///   tracking an instance within @ref time_deriv_0d::ContextPack.
  void unsafe_reset_wrapped_ptr(const grackle_field_data* field_data) {
    // we are going to need to use a chemistry_data* pointer to carry out this
    // work in the near future
    GRIMPL_REQUIRE(field_data != nullptr, "field_data is a nullptr");
    field_data_ = field_data;

    // the following call does nothing if the size doesn't change
    sp_arr_of_ptrs_.resize(MAX_EVOLVED_SPECIES_FIELDS);

#define ENTRY(SPECIES_NAME)                                                    \
  sp_arr_of_ptrs_[SpLUT::SPECIES_NAME] = field_data->SPECIES_NAME##_density;
#include "field_data_evolved_species.def"
#undef ENTRY
  }

  /// @brief Construct a new instance
  explicit FieldAdaptorManager(const grackle_field_data* field_data)
      : field_data_(nullptr) {
    unsafe_reset_wrapped_ptr(field_data);
  }

  FieldAdaptorManager() = default;  // this exists to support move-operations
  FieldAdaptorManager(FieldAdaptorManager&&) = default;
  FieldAdaptorManager& operator=(FieldAdaptorManager&&) = default;

  // use of copy-constructor or copy-assignment is probably out of error
  // (we can always define them later)
  FieldAdaptorManager(const FieldAdaptorManager&) = delete;
  FieldAdaptorManager& operator=(const FieldAdaptorManager&) = delete;

  /// @brief retrieve views of the dynamically evolved species data
  SpeciesMultiView<gr_float> get_species_data(IndexRange idx_range) const {
    int i_extent = field_data_->grid_dimension[0];
    SpeciesMultiView<gr_float> out(sp_arr_of_ptrs_.data(), i_extent,
                                   sp_arr_of_ptrs_.size());

    // we intentionally use i = 0 so that we properly support loops from
    // idx_range.i_start to idx_range.i_stop
    int ptr_offset = layoutleft_3D_index_to_1D_(field_data_->grid_dimension, 0,
                                                idx_range.j, idx_range.k);
    out.override_ptr_offset_and_ilen(ptr_offset, i_extent);
    return out;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // FIELD_ADAPTOR_HPP
