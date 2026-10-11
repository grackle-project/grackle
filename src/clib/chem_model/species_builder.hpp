//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare the @ref SpeciesBuilder type
///
//===----------------------------------------------------------------------===//
#ifndef CHEM_MODEL_SPECIES_BUILDER_HPP
#define CHEM_MODEL_SPECIES_BUILDER_HPP

#include <string>
#include <utility>  // std::move
#include <variant>
#include <vector>

#include "LUT.hpp"
#include "chem_model/infer_species.hpp"
#include "support/config.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief denotes chemical species that are special cases
///
/// Each enumerator are explicitly named after their canonical species names.
/// It is used to configure a @ref SpeciesBuilder object about the
/// existence about the existence of the represented species.
///
/// Future Thoughts
/// ---------------
/// The existence of this type only makes sense while the number of special
/// cases is small. This is a bit of a kludge that we should reconsider if
/// more special cases are added.
///
/// Entries belong to one of 3 categories:
/// 1. free electrons (i.e. the "e" species)
/// 2. negatively charged species composed of nuclides ("HM" and "DM")
/// 3. metal species that aren't actually associated with a charge at all.
///    Currently, they're just used for accounting purposes for dust growth.
///    In practice, they must have a charge of 0 in order for Grackle to
///    conserve charge (because they participate in dust reactions that
///    consume/create ordinary chemical species that explicitly have charge)
///
/// The first category remains special. But the other 2 categories could be
/// viewed as chemical species with special names.
///
/// Thoughts about Metal Accounting Species
/// ---------------------------------------
/// As noted above, we currently treat metal accounting species as a regular
/// metal species with a special name and a charge of 0.
///
/// Going forward we may want to either:
/// - treat them as ordinary metal species (and rename them accordingly)
/// - treat them as a distinct category of species (like a dust species).
///   This might mean
///   - they aren't considered a @ref SpKind::METAL kind of species (e.g.
///     we might define a new species kind)
///   - they don't have an associated @ref ChemSpeciesProps object (or perhaps
///     the associated objects should reflect the distinct category).
enum class SpecialCaseSpecies { e, HM, DM, Mg, Al, S, Fe };

/// @brief A builder for configuring the model of actively evolved species
///
/// Essentially you use the methods to specify the different kinds of species
/// that are **ACTIVELY EVOLVED** (and the relative order of species of a given
/// kind).
///
/// The following code snippet shows a simple snippet corresponding to the
/// network corresponding to the classic primordial_chemistry == 0
///
/// @code{C++}
/// auto rslt = SpeciesBuilder()
///     .chemical<SpKind::PRIMORDIAL>(SpecialCaseSpecies::e)
///     .chemical<SpKind::PRIMORDIAL>("HI")
///     .chemical<SpKind::PRIMORDIAL>("HII")
///     .chemical<SpKind::PRIMORDIAL>("HeI")
///     .chemical<SpKind::PRIMORDIAL>("HeII")
///     .chemical<SpKind::PRIMORDIAL>("HeIII")
///     .build();
/// @endcode
///
/// Why make this a builder?
/// ------------------------
/// Longer term, it would be great if we could more easily allow easy extension
/// of a dynamic chemical network. Dynamically specifying the species (and
/// their properties) is a very important part of doing that. We would also
/// need to be able to dynamically specify the relvant reactions. To achieve
/// that we might extend this builder for that purpose (or create a new builder
/// composed out of this one if things get too complex).
///
/// To be clear, this is currently just intended for internal use (making it
/// easier to customise the relevant species would be useful for extending
/// networks or testing purposes). We could eventually consider exposing an
/// iterface to end users, but that would require a lot of consideration about
/// API-design
class SpeciesBuilder {
  /// @brief a type-safe union encoding information about a chemical Species
  ///
  /// This type-safe union holds either
  /// - a string encoding a "regular" chemical formula (i.e. we can parse it and
  ///   use `NuclideModel` to infer the chemical property)
  /// - a enum encoding a particular special case
  ///
  /// @note In the near future, we may track a custom type holding a reference
  ///   to a string-literal and the result of a compile-time attempt at parsing
  ///   that literal
  using ChemEntryVariant_ = std::variant<std::string, SpecialCaseSpecies>;

  // internal data members
  std::vector<ChemEntryVariant_> primordial_vec_;
  std::vector<ChemEntryVariant_> metal_vec_;
  std::vector<std::string> dust_vec_;

  // a helper method returning a reference to primordial_vec_ or metal_vec_
  template <SpKind kind>
  std::vector<ChemEntryVariant_>& get_chem_vec_() {
    static_assert(kind == SpKind::PRIMORDIAL || kind == SpKind::METAL,
                  "kind parameter doesn't correspond to a chemical species");
    return (kind == SpKind::PRIMORDIAL) ? primordial_vec_ : metal_vec_;
  }

public:
  SpeciesBuilder() {
    primordial_vec_.reserve(PrimordialSpLUT::NUM_ENTRIES);
    metal_vec_.reserve(MetalSpLUT::NUM_ENTRIES);
    dust_vec_.reserve(DustSpLUT::NUM_ENTRIES);
  }

  /// @brief specify a standard primordial or chemical species
  template <SpKind kind>
  SpeciesBuilder& chemical(std::string name) {
    get_chem_vec_<kind>().emplace_back(std::move(name));
    return *this;
  }

  /// @brief specify a special-case primordial or chemical species
  template <SpKind kind>
  SpeciesBuilder& chemical(SpecialCaseSpecies entry) {
    get_chem_vec_<kind>().emplace_back(entry);
    return *this;
  }

  /// @brief specify an actively evolved dust species
  SpeciesBuilder& dust(std::string name) {
    dust_vec_.emplace_back(name);
    return *this;
  }

  /// @brief build the resulting species initialization information
  ///
  /// Near Future Plan
  /// ----------------
  /// In the near future, the goal is to pass in a reference to a
  /// @ref NuclideModel as an argument and use it to construct a
  /// collection holding @ref ChemSpeciesProps objects that will be stored as
  /// part of the output object
  Expected<SpInitializeInfo, Error> build() const;
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // CHEM_MODEL_SPECIES_BUILDER_HPP