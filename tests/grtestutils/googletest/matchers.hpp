//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file Define custom matchers
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_GOOGLETEST_MATCHERS_HPP
#define GRTESTUTILS_GOOGLETEST_MATCHERS_HPP

#include <gmock/gmock.h>
#include <set>

namespace grtest {

/// Implements the polymorphic IsUnique matcher
///
/// This can be used as a Matcher<T> as long as T is a container (like
/// std::vector or std::string) that defines:
/// - the T::begin() and T::end() member functions
/// - the T::const_iterator member type
/// - the T::value_type member type
///
/// @note
/// The implementation is modelled on the way that gmock implements matchers
/// for C++ standard library containers
///
/// @todo
/// Currently, there is an additional implicit requirement that the elements
/// are copy-constructible. If we want a more general implementation, we should
/// consider using std::reference_wrapper
class HoldsUniqueMatcher {
public:
  template <typename MatcheeContainerType>
  bool MatchAndExplain(const MatcheeContainerType& c,
                       ::testing::MatchResultListener* listener) const {
    typename MatcheeContainerType::const_iterator it = c.begin();
    typename MatcheeContainerType::const_iterator last = c.end();

    using val_type = typename MatcheeContainerType::value_type;
    std::set<val_type> set;

    for (; it != last; ++it) {
      auto insertion_result = set.emplace(*it);
      if (!insertion_result.second) {
        *listener << "who contains at least 2 occurrences of the element, "
                  << *it;
        return false;
      }
    }
    return true;
  }

  // Describes what this matcher matches.
  void DescribeTo(std::ostream* os) const { *os << "contains unique elements"; }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "contains more than one occurrence of an element";
  }
};

/// Matches an STL-style container that holds unique values
///
/// @par Example
/// In the following example, each call to `EXPECT_THAT` should pass
/// @code{.cpp}
///   std::vector<int> v;
///   EXPECT_THAT(v, grtest::HoldsUnique());
///   v.push_back(1);
///   v.push_back(3);
///   EXPECT_THAT(v, grtest::HoldsUnique());
///   v.push_back(1);
///   EXPECT_THAT(v, ::testing::Not(grtest::HoldsUnique()));
/// @endcode
inline testing::PolymorphicMatcher<HoldsUniqueMatcher> HoldsUnique() {
  return testing::MakePolymorphicMatcher(HoldsUniqueMatcher());
}

}  // namespace grtest

#endif  // GRTESTUTILS_GOOGLETEST_MATCHERS_HPP
