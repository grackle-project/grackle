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

#include <ostream>
#include <string>

#include "chem_model/chemsp_props.hpp"
#include "support/inplace_vector.hpp"

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

}  // namespace GRIMPL_NAMESPACE_DECL

namespace grtest {

struct NameCoefPair {
  std::string name;
  int64_t coef;
};

// a polymorphic matcher for comparing a vector of NameCoefPair objects
// against a collection of IdCoefPair objects
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

  // for the sake of convenience, the plan is to introduce a second overload of
  // bool MatchAndExplain to handle an argument that directly specifies a
  // name-map pair

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

// ok, so the better plan is to replace NameCoefSeqMatcher with NameCoefMatcher
// (i.e. that just compares a single NameCoefPair) and then define

inline auto NameCoefSeq(const GRIMPL_NS::FrozenKeyIdxBiMap& name_map,
                        const std::vector<NameCoefPair>& v) {
  std::vector<NameCoefMatcher> tmp;
  for (const NameCoefPair& element : v) {
    tmp.emplace_back(name_map, element);
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

}  // namespace grtest

#endif  // GRTESTUTILS_GOOGLETEST_CHEM_MODEL_TOOLS_HPP