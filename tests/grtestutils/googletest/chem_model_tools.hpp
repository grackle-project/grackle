//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define some custom tools for testing chem_model machinery
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_GOOGLETEST_CHEM_MODEL_TOOLS_HPP
#define GRTESTUTILS_GOOGLETEST_CHEM_MODEL_TOOLS_HPP

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <cstdint>
#include <ostream>
#include <string>

#include "chem_model/chemsp_props.hpp"
#include "chem_model/parse.hpp"
#include "chem_model/rxn_eqn_props.hpp"
#include "chem_model/species.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "support/inplace_vector.hpp"
#include "support/error.hpp"
#include "support/expected.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// teach GoogleTest how to print InplaceVector
template <typename T, std::size_t N>
void PrintTo(const InplaceVector<T, N>& v, std::ostream* os) {
  *os << '{';
  for (std::size_t i = 0; i < v.size(); i++) {
    if (i > 0) {
      *os << ", ";
    }
    *os << testing::PrintToString(v[i]);
  }
  *os << '}';
}

inline void PrintTo(const IdCoefPair& p, std::ostream* os) {
  *os << "{id=" << p.id << ",coef=" << int64_t(p.coef) << '}';
}

// teach GoogleTest how to print our SpKind enum for more informative errors
inline void PrintTo(const SpKind& kind, std::ostream* os) {
  switch (kind) {
    case SpKind::PRIMORDIAL:
      *os << "SpKind::PRIMORDIAL";
      return;
    case SpKind::METAL:
      *os << "SpKind::METAL";
      return;
    case SpKind::DUST:
      *os << "SpKind::DUST";
      return;
    case SpKind::FILLER_A:
      *os << "SpKind::FILLER_A";
      return;
    case SpKind::FILLER_B:
      *os << "SpKind::FILLER_B";
      return;
    default:
      *os << "SpKind::<UNKNOWN>";
      return;
  }
}

// teach GoogleTest how to print our SpKind enum for more informative errors
inline void PrintTo(const EqnSide& side, std::ostream* os) {
  *os << ((side == EqnSide::LEFT) ? "EqnSide::LEFT" : "EqnSide::RIGHT");
}

}  // namespace GRIMPL_NAMESPACE_DECL

namespace grtest {

struct ParameterSpec {
  int primordial_chemistry;
  int metal_chemistry;
  int dust_species;
};

// teach GoogleTest how to print our sample partition descriptors for more
// informative errors
inline void PrintTo(ParameterSpec spec, std::ostream* os) {
  *os << "{prim=" << spec.primordial_chemistry
      << ",metal=" << spec.metal_chemistry << ",dustsp=" << spec.dust_species
      << '}';
}

// a helper function to help with constrrucing the species maps for the tests
inline GRIMPL_NS::Expected<GRIMPL_NS::SpInitializeInfo, GRIMPL_NS::Error>
infer_species_maps(ParameterSpec spec) {
  using GRIMPL_NS::GrainSpeciesInfo;
  std::optional<GrainSpeciesInfo> maybe_grainsp_info;
  const GrainSpeciesInfo* grain_info = nullptr;
  if (spec.dust_species > 0) {
    using GRIMPL_NS::GrainSpeciesInfo;
    GRIMPL_NS::Expected<GrainSpeciesInfo, GRIMPL_NS::Error> rslt =
        GrainSpeciesInfo::create(spec.dust_species);
    if (!rslt.has_value()) {
      return GRIMPL_NS::Unexpected(
          rslt.error().context_literal("unable to init GrainSpeciesInfo"));
    }
    maybe_grainsp_info = std::move(rslt).value();
    grain_info = &maybe_grainsp_info.value();
  }

  return GRIMPL_NS::infer_species_maps(spec.primordial_chemistry,
                                       spec.metal_chemistry, grain_info);
}

struct NameCoefPair {
  std::string name;
  int64_t coef;

  bool operator==(const NameCoefPair&) const = default;
};

inline void PrintTo(const NameCoefPair& p, std::ostream* os) {
  *os << "{name=" << p.name << ",coef=" << int64_t(p.coef) << '}';
}

/// compares a NameCoefPair objects against a of GRIMPL_NS::IdCoefPair object
///
/// Without this, checking validity of a GRIMPL_NS::IdCoefPair is messy:
/// - when specifying the expected value, you need to figure out what the
///   expected id is (usually involves a function call or 2)
/// - when there is an error, it can be hard to judge what went wrong from id
///   values, alone
///
/// This works around both of these values by mapping between names and ids
class NameCoefMatcher {
  GRIMPL_NS::FrozenKeyIdxBiMap name_map_;
  NameCoefPair ref_;
  uint16_t ref_id_;

public:
  using is_gtest_matcher = void;

  NameCoefMatcher(GRIMPL_NS::FrozenKeyIdxBiMap name_map, NameCoefPair ref)
      : name_map_(name_map), ref_(ref), ref_id_() {
    std::optional<uint16_t> maybe_ref_id = name_map.find(ref.name);
    if (!maybe_ref_id.has_value()) {
      GR_INTERNAL_ERROR("name_map doesn't have the \"%s\" entry",
                        ref_.name.c_str());
    }
    ref_id_ = maybe_ref_id.value();
  }

  bool MatchAndExplain(const GRIMPL_NS::IdCoefPair& pair,
                       std::ostream* os) const {
    if (pair.id != ref_id_) {
      if (os != nullptr) {
        const char* name = name_map_.inverse_find(pair.id);
        if (name == nullptr) {
          *os << "id doesn't correspond to known entry";
        } else {
          *os << "id corresponds to \"" << name << '"';
        }
        *os << ": id should be " << ref_id_;
      }
      return false;
    } else if (static_cast<int64_t>(pair.coef) != ref_.coef) {
      if (os != nullptr) {
        *os << "coef should be " << ref_.coef;
      }
      return false;
    } else {
      return true;
    }
  }

  // Describes the property of a value matching this matcher.
  void DescribeTo(std::ostream* os) const {
    *os << "consistent with {name=\"" << ref_.name << "\" (id=" << ref_id_
        << "), coef=" << ref_.coef << '}';
  }

  // Describes the property of a value NOT matching this matcher.
  void DescribeNegationTo(std::ostream* os) const {
    *os << "not ";
    DescribeTo(os);
  }
};

// create a matcher to checks if an IdCoefPair is consistent with a specified
// NameCoefPair (the name_map is used to translate between a name and an id)
inline ::testing::Matcher<const GRIMPL_NS::IdCoefPair&> IsNameCoefPair(
    const GRIMPL_NS::FrozenKeyIdxBiMap& name_map, const NameCoefPair& pair) {
  return NameCoefMatcher(name_map, pair);
}

// create a matcher to checks if a sequence of IdCoefPair is consistent with
// the contents of a vector of NameCoefPair (the name_map is used to translate
// between a name and an id)
//
// the return type is auto here because the gmock doesn't establish a stable
// return-type for `::testing::ElementsAreArray` as part of the API
inline auto NameCoefSeq(const GRIMPL_NS::FrozenKeyIdxBiMap& name_map,
                        const std::vector<NameCoefPair>& v) {
  std::vector<::testing::Matcher<const GRIMPL_NS::IdCoefPair&>> tmp;
  for (const NameCoefPair& element : v) {
    tmp.push_back(IsNameCoefPair(name_map, element));
  }
  return ::testing::ElementsAreArray(std::move(tmp));
}

// define a nice matcher for standard cases of GRIMPL_NS::ChemSpeciesProps
inline ::testing::Matcher<GRIMPL_NS::ChemSpeciesProps> HasChemProps(
    const GRIMPL_NS::NuclideModel& nuclide_model, int elementary_charge,
    const std::vector<NameCoefPair>& nuclide_subscript_pairs) {
  using ::GRIMPL_NS::ChemSpeciesProps;
  using ::testing::Field;

  double mass_factor = 0.0;
  bool composed_of_nuclides = !nuclide_subscript_pairs.empty();
  for (const NameCoefPair& pair : nuclide_subscript_pairs) {
    const GRIMPL_NS::NuclideProp* props =
        nuclide_model.try_get(pair.name.c_str());
    if (props == nullptr) {
      GR_INTERNAL_ERROR("nuclide model doesn't have entry for: %s",
                        pair.name.c_str());
    }
    mass_factor += pair.coef * props->mass_factor;
  }

  auto name_coef_seq_matcher =
      NameCoefSeq(nuclide_model.symbol_map(), nuclide_subscript_pairs);

  return ::testing::AllOf(
      Field("elementary_charge", &ChemSpeciesProps::elementary_charge,
            ::testing::Eq(elementary_charge)),
      Field("mass_factor", &ChemSpeciesProps::mass_factor,
            ::testing::DoubleEq(mass_factor)),
      Field("nuclide_id_coef_pairs", &ChemSpeciesProps::nuclide_id_coef_pairs,
            name_coef_seq_matcher),
      ::testing::Property("composed_of_nuclides",
                          &ChemSpeciesProps::composed_of_nuclides,
                          ::testing::Eq(composed_of_nuclides)));
}

// define a nice matcher for standard cases of GRIMPL_NS::RxnEqnProps
inline ::testing::Matcher<GRIMPL_NS::RxnEqnProps> HasRxnEqnProps(
    const GRIMPL_NS::FrozenKeyIdxBiMap& species_name_map,
    const std::vector<NameCoefPair>& lhs_entries,
    const std::vector<NameCoefPair>& rhs_entries,
    std::optional<GRIMPL_NS::EqnSide> has_photon = std::nullopt) {
  using ::GRIMPL_NS::RxnEqnProps;
  using ::testing::Field;

  auto lhs_matcher = NameCoefSeq(species_name_map, lhs_entries);
  auto rhs_matcher = NameCoefSeq(species_name_map, rhs_entries);

  return ::testing::AllOf(
      Field("lhs_entries", &RxnEqnProps::lhs_entries, lhs_matcher),
      Field("rhs_entries", &RxnEqnProps::rhs_entries, rhs_matcher),
      Field("has_photon", &RxnEqnProps::has_photon, ::testing::Eq(has_photon)));
}

}  // namespace grtest

#endif  // GRTESTUTILS_GOOGLETEST_CHEM_MODEL_TOOLS_HPP