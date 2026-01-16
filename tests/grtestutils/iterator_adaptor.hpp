//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare and implement the IteratorAdaptor
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_ITERATOR_ADAPTOR_HPP
#define GRTESTUTILS_ITERATOR_ADAPTOR_HPP

#include <iterator>
#include <ostream>
#include <string>

#include "grackle.h"
#include "preset.hpp"

namespace grtest {

/// the standard value-type that an IteratorAdaptor instantiation refers to
struct NameIdPair {
  std::string name;
  long long id;
};

/// teach std::ostream how to format NameIdPair
///
/// The motivation is to make it easier write detailed error messages
inline std::ostream& operator<<(std::ostream& os, const NameIdPair& pair) {
  os << "{name=\"" << pair.name << "\", id=" << pair.id << "}";
  return os;
}

/// implements a C++ style InputIterator by adapting a simple Plugin type
/// that wraps a set of Grackle functions
///
/// This is useful for making use of C++ standard library algorithms and
/// (arguably more importantly) making use of range-based for-loops
template <class Plugin>
class IteratorAdaptor {
  unsigned long long counter_;
  unsigned long long n_rates_;
  Plugin plugin_;
  NameIdPair current_pair_;

  /// Updates current_pair_ and returns `*this`
  IteratorAdaptor& update_pair_and_ret_(unsigned long long current_count) {
    if (current_count < this->n_rates_) {
      this->current_pair_ = this->plugin_(current_count);
    }
    return *this;
  }

public:
  using iterator_category = std::input_iterator_tag;
  using value_type = NameIdPair;
  using difference_type = std::ptrdiff_t;
  using pointer = const NameIdPair*;
  using reference = const NameIdPair;

  /// construct a new instance
  IteratorAdaptor(unsigned long long counter, unsigned long long n_rates,
                  Plugin plugin)
      : counter_(counter), n_rates_(n_rates), plugin_(plugin) {
    update_pair_and_ret_(counter);
  }

  /// implements the equality operation
  bool operator==(const IteratorAdaptor& other) const {
    return (counter_ == other.counter_) && (plugin_ == other.plugin_);
  }

  /// implements the inequality operation
  bool operator!=(const IteratorAdaptor& other) const {
    return !(*this == other);
  }

  /// implements the dereference operation
  reference operator*() const { return current_pair_; }

  /// implements the prefix increment operation
  ///
  /// This effectively implements `++x`, which increments the value of `x`
  /// before determining the returned value. In other words, `++x` returns the
  /// value of `x` from **after** after the increment
  IteratorAdaptor& operator++() { return update_pair_and_ret_(++counter_); }

  /// implements the prefix increment operation
  ///
  /// This effectively implements `x++`, which increments the value of `x`
  /// after determining the returned value. In other words, `x++` returns the
  /// value of `x` from **before** the increment
  IteratorAdaptor operator++(int) {
    IteratorAdaptor ret = *this;
    ++(*this);
    return ret;
  }
};

// Now lets use this machinery to implement logic iterating over the names
// accessible through the ratequery api

struct RateQueryPlugin {
  chemistry_data_storage* my_rates;

  NameIdPair operator()(unsigned long long i) const {
    grunstable_rateid_type tmp;
    const char* name = grunstable_ith_rate(my_rates, i, &tmp);
    return NameIdPair{name, tmp};
  }

  bool operator==(const RateQueryPlugin& other) const {
    return my_rates == other.my_rates;
  }
};

/// used for creating the iterator and within range-based for-loops
class RateQueryRange {
  RateQueryPlugin plugin_;
  using iterator = IteratorAdaptor<RateQueryPlugin>;
  long long n_rates_;

public:
  explicit RateQueryRange(grtest::GrackleCtxPack& pack)
      : plugin_(RateQueryPlugin{pack.my_rates()}),
        n_rates_(grunstable_ratequery_nrates(pack.my_rates())) {}

  iterator begin() { return iterator(0, n_rates_, plugin_); }
  iterator end() { return iterator(n_rates_, n_rates_, plugin_); }
};

}  // namespace grtest

#endif  // GRTESTUTILS_ITERATOR_ADAPTOR_HPP
