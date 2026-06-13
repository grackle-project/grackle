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
/// a common value-type that an @ref IteratorAdaptor instantiation refers to
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
/// that wraps random-access functionality.
///
/// A plugin type should look something like the following
/// @code{cpp}
/// struct MyPlugin {
///   // define the "ValueType". The is the type returned by the function call
///   // method. (for an instance of ``MyPlugin`` called  ``obj``,
///   // ``obj(index)`` will return an instance of the ValueType)
///   using ValueType = /* my_value_type */;
///
///   // This overloads the function call operator (its directly analogous to
///   // pythons __call__ method). It returns the value corresponding to the
///   // provided index ``i``
///   ValueType operator()(unsigned long long i) const;
///
///   // overloads the equality check operation (it's invoked when you write
///   // `objA == objB`)
///   bool operator==(const RateQueryPlugin& other) const;
/// };
/// @endcode
///
/// This is useful for making use of C++ standard library algorithms and
/// (arguably more importantly) making use of range-based for-loops
template <class Plugin>
class IteratorAdaptor {
public:
  // define the type aliases that make up the common interface expected by C++
  // function templates in the standard library for conveying properties about
  // the iterator.
  using iterator_category = std::input_iterator_tag;
  using value_type = typename Plugin::ValueType;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type*;
  using reference = const value_type;

private:
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
  using ValueType = NameIdPair;
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

// Now lets use this machinery to implement logic for iterating over the names
// of chemistry parameters
extern "C" {
typedef const char* param_name_fn(unsigned int);
}

namespace grtest {

struct ParamItrPlugin {
  using ValueType = NameIdPair;
  param_name_fn* fn;

  NameIdPair operator()(unsigned long long i) const {
    return {fn(static_cast<unsigned>(i)), static_cast<long long>(i)};
  }

  bool operator==(const ParamItrPlugin& other) const { return fn == other.fn; }
};

enum class ChemParamType { INT, DOUBLE, STRING };

/// used for creating the iterator and within range-based for-loops
class ChemParamRange {
  using iterator = IteratorAdaptor<ParamItrPlugin>;
  ParamItrPlugin plugin_;
  long long n_rates_;

public:
  explicit ChemParamRange(grtest::ChemParamType param) {
    switch (param) {
      case ChemParamType::INT:
        n_rates_ = grackle_num_params("int");
        plugin_ = ParamItrPlugin{&param_name_int};
        break;
      case ChemParamType::DOUBLE:
        n_rates_ = grackle_num_params("double");
        plugin_ = ParamItrPlugin{&param_name_double};
        break;
      case ChemParamType::STRING:
        n_rates_ = grackle_num_params("string");
        plugin_ = ParamItrPlugin{&param_name_string};
        break;
      default:  // should NEVER come up (but the behavior is well-defined)
        n_rates_ = 0;
        plugin_ = ParamItrPlugin{nullptr};
    }
  }

  iterator begin() { return iterator(0, n_rates_, plugin_); }
  iterator end() { return iterator(n_rates_, n_rates_, plugin_); }
};

}  // namespace grtest

#endif  // GRTESTUTILS_ITERATOR_ADAPTOR_HPP
