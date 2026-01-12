//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares logic for constructing a runtime species lookup table
///
//===----------------------------------------------------------------------===//
#ifndef RUNTIME_SPLUT_HPP
#define RUNTIME_SPLUT_HPP

#include "chemistry/species.hpp"
#include "dust/grain_species_info.hpp"
#include "support/config.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/PartMap.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// Describes the kinds of dynamically evolved species
struct SpKind {
  enum {
    CHEMICAL = 0,  ///< primordial or metal species (includes electrons)
    DUST = 1,  ///< an individual kind of grain (or possibly a grain ensemble)
  };
};

/// the result of @ref make_runtime_splut
struct PartitionedKeyIdxBiMapRslt {
  bool success;
  FrozenKeyIdxBiMap bimap;
  PartMap part_map;
};

/// construct a runtime species lookup table
///
/// The resulting bimap encodes a runtime counterpart to the compile-time
/// SpLUT. It maps the name of each dynamically evolved chemical species or
/// grain species to an index and vice versa. The keys of the table are
/// partitioned based on the species kind (see @ref SpKind).
///
/// @par Future Thoughts
/// It's a little weird to pass GrainSpeciesInfo directly into this function.
/// We may want to revisit that choice in the future, (especially if other
/// dust models are implemented that affect the dynamical species fields)
PartitionedKeyIdxBiMapRslt make_runtime_splut(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info);

namespace runtime_splut_detail {

/// carries callback context
struct ChemicalSpCtx {
  int count;
  const char** names;
};

inline void chemsp_counter(const char* name, ChemicalSpKind kind, void* ctx) {
  static_cast<ChemicalSpCtx*>(ctx)->count++;
}

inline void chemsp_name_cpy(const char* name, ChemicalSpKind kind, void* ctx) {
  int& total_count = static_cast<ChemicalSpCtx*>(ctx)->count;
  static_cast<ChemicalSpCtx*>(ctx)->names[total_count] = name;
  total_count++;
}

}  // namespace runtime_splut_detail

inline PartitionedKeyIdxBiMapRslt make_runtime_splut(
    int primordial_chemistry, int metal_chemistry,
    const GrainSpeciesInfo* grain_info) {
  MetalNuclideMassSpChoice choice =
      (grain_info == nullptr)
          ? MetalNuclideMassSpChoice::NONE
          : from_dust_species_param(grain_info->dust_species_parameter);
  int n_grain_sp = (grain_info == nullptr) ? 0 : grain_info->n_species;

  // count up the names of the chemical species
  runtime_splut_detail::ChemicalSpCtx counter_ctx{0, nullptr};
  visit_canonical_chem_species(primordial_chemistry, metal_chemistry, choice,
                               runtime_splut_detail::chemsp_counter,
                               static_cast<void*>(&counter_ctx));

  int n_chemical_sp = counter_ctx.count;
  int total_size = n_chemical_sp + n_grain_sp;
  const char** names = nullptr;

  if (total_size > 0) {
    names = new const char*[total_size];

    // copy every chemical species name into names
    runtime_splut_detail::ChemicalSpCtx name_cpy_ctx{0, names};
    visit_canonical_chem_species(primordial_chemistry, metal_chemistry, choice,
                                 runtime_splut_detail::chemsp_name_cpy,
                                 static_cast<void*>(&name_cpy_ctx));

    // copy every dust species name into names
    if (grain_info != nullptr) {
      for (int i = 0; i < grain_info->n_species; i++) {
        names[n_chemical_sp + i] =
            FrozenKeyIdxBiMap_inverse_find(&grain_info->name_map, i);
      }
    }
  }

  // now it's time to construct the objects
  PartitionedKeyIdxBiMapRslt out;
  int species_kinds[2] = {SpKind::CHEMICAL, SpKind::DUST};
  int sizes[2] = {n_chemical_sp, n_grain_sp};
  out.part_map = new_PartMap(species_kinds, sizes, 2);
  if (!PartMap_is_ok(&out.part_map)) {
    out.success = false;
    out.bimap = mk_invalid_FrozenKeyIdxBiMap();
  } else {
    out.bimap =
        new_FrozenKeyIdxBiMap(names, total_size, BiMapMode::COPIES_KEYDATA);
    out.success = FrozenKeyIdxBiMap_is_ok(&out.bimap);
  }

  if (names != nullptr) {
    delete[] names;
  }
  return out;
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // RUNTIME_SPLUT_HPP