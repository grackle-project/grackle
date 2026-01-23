//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// <add a short description>
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_PARAM_HPP
#define GRTESTUTILS_PARAM_HPP

#include "grackle.h"
#include <initializer_list>
#include <iosfwd>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace grtest {

namespace param_detail {
// a forward declaration
class StrAllocTracker;
}  // namespace param_detail

/// Used to hold a chemistry_data parameter value (of arbitrary type)
///
/// This is primarily intended to help construct a chemistry_data instance
class ParamVal {
  friend bool set_param(chemistry_data& my_chem, const std::string& name,
                        const ParamVal& par_val,
                        param_detail::StrAllocTracker* str_allocs);

  using val_type = std::variant<std::string, nullptr_t, int, double>;

  // private attribute
  val_type val_;

  static val_type coerce_c_string_(const char* s) {
    return (s == nullptr) ? val_type(std::in_place_type<std::nullptr_t>)
                          : val_type(std::in_place_type<std::string>, s);
  }

public:
  // it isn't possible to have an empty parameter value
  ParamVal() = delete;

  explicit ParamVal(int val) : val_(val) {}
  explicit ParamVal(double val) : val_(val) {}
  explicit ParamVal(std::string val) : val_(val) {}
  explicit ParamVal(std::string_view val) : val_(std::string(val)) {}
  explicit ParamVal(const char* val) : val_(ParamVal::coerce_c_string_(val)) {}

  /// return a string-representation of `this`
  std::string to_string() const;

  /// teach googletest how to print the value
  friend void PrintTo(const ParamVal& p, std::ostream* os);

  bool operator==(const ParamVal& other) const { return val_ == other.val_; }
  bool operator!=(const ParamVal& other) const { return val_ != other.val_; }

  // all the following methods are primarily intended for debugging purposes
  // (they are not essential to the core functionality)

  /// Try to access a pointer the contained value (intended for debugging)
  ///
  /// For example, `param_val.try_get<int>()`, tries to access the value as an
  /// integer
  template <class T>
  std::optional<T> try_get() const {
    T* ptr = std::get_if<T>(&val_);
    return (ptr == nullptr) ? std::nullopt : std::optional<T>{*ptr};
  }

  /// Returns whether `this` holds an empty string
  bool is_empty_string() const {
    return std::holds_alternative<nullptr_t>(val_);
  }

  /*!{*/
  /// Checks equivalence with the specified value
  bool is_equal(int val) const;
  bool is_equal(double val) const;
  bool is_equal(std::string_view val) const;
  bool is_equal(const char* val) const;
  /*!}*/
};

/// Use this whenever you want to make a sequence of key-value pairs
using ParamPair = std::pair<std::string, ParamVal>;

/// Nicely format a @ref ParamPair as a string
inline std::string to_string(const ParamPair& pair) {
  std::string str_val = pair.second.to_string();
  std::string out;
  out.reserve(pair.first.size() + str_val.size() + 6);
  out += "{\"";
  out += pair.first;
  out += "\", ";
  out += str_val;
  out += '}';
  return out;
}

namespace param_detail {

/// this **ONLY** exists to help implement @ref make_ParamPair_vec
struct ParamPairInit {
  ParamPair pair;

  template <class V>
  ParamPairInit(const char* k, V&& v) : pair(k, ParamVal(std::forward<V>(v))) {}
};

}  // namespace param_detail

/// A convenience function for easier construction of a vector of Param pairs
///
/// This function lets you write code like
/// @code{C++}
/// std::vector<ParamPair> v = make_ParamPair_vec({
///   {"use_grackle", 1},
///   {"with_radiative_cooling", 1},
///   {"metal_cooling", 1},
///   {"UVbackground", 1},
///   {"primordial_chemistry", 2}
/// });
/// @endcode
///
/// We could definitely get rid of this function, but then we would need to
/// rewrite the above snippet as:
/// @code{C++}
/// std::vector<ParamPair> v = {
///   {"use_grackle", ParamVal(1)},
///   {"with_radiative_cooling", ParamVal(1)},
///   {"metal_cooling", ParamVal(1)},
///   {"UVbackground", ParamVal(1)},
///   {"primordial_chemistry", ParamVal(2)}
/// };
inline std::vector<ParamPair> make_ParamPair_vec(
    std::initializer_list<param_detail::ParamPairInit> l) {
  // if we were a little more clever, we might be able to get rid of an
  // intermediate allocation when constructing ParamPairInit...
  std::vector<ParamPair> v;
  v.reserve(l.size());
  for (const param_detail::ParamPairInit& p : l) {
    v.push_back(p.pair);
  }
  return v;
}

/// Tries to update @p my_chem to hold the value specified by @p par_val
///
/// @param[in,out] my_chem Tracks various Grackle parameters
/// @param[in] name The name of the parameter getting updated
/// @param[in] par_val The value of the parameter
/// @param[in] str_allocs Tracks allocations of the strings held by @p my_chem
/// @returns true if successful and `false` if there was an error (e.g. @p name
///     isn't a known parameter or is a parameter that doesn't expect a string)
///
/// @note
/// When @p par_val holds a non-null string and @p name is a known parameter
/// holding a string:
/// - this function ordinarily updates @p str_allocs to hold a copy of
///   @p par_val and @p my_chem will be updated to store a pointer to that copy.
/// - When @p str_allocs is a nullptr, @p my_chem is directly updated to track
///   the data contained within @p par_val.
bool set_param(chemistry_data& my_chem, const std::string& name,
               const ParamVal& par_val,
               param_detail::StrAllocTracker* str_allocs);

/// Tries to update @p my_chem to hold the value specified by @p par_val
///
/// This is designed work with various standard library container types holding
/// key-value parameter name pairs. For concreteness, suppose we have a variable
/// called `my_params` with one of the following types:
/// - `std::vector<std::pair<std::string, ParamVal>>` (this is equivalent to
///   `std::vector<ParamPair>`)
/// - `std::map<std::string, ParamVal>`
/// - `std::unordered_map<std::string, ParamVal>`
///
/// For each of the cases, you would invoke
/// @code{C++}
/// set_params(my_params.cbegin(), my_params.cend(), my_chem, str_allocs);
/// @endcode
///
/// @param[in] first, last the pair of iterators defining the range of parameter
///     name-value pairs that will be used to update my_chem.
/// @param[in,out] my_chem Tracks various Grackle parameters
/// @param[in] str_allocs Tracks allocations of the strings held by @p my_chem
/// @returns An empty optional if successful. Otherwise, it returns the name
///     of the parameter that couldn't be updated.
///
/// @tparam It The type of the iterator. When str_allocs is a `nullptr` this
///     must be a ForwardIterator (i.e. after incrementing the iterator,
///     previous values must remain valid). Otherwise, this can also be an
///     InputIterator
template <class It>
std::optional<std::string> set_params(
    It first, It last, chemistry_data& my_chem,
    param_detail::StrAllocTracker* str_allocs) {
  for (It it = first; it != last; ++it) {
    const std::string& name = it->first;
    const ParamVal& val = it->second;
    if (!set_param(my_chem, name, val, str_allocs)) {
      return std::optional<std::string>(name);
    }
  }
  return std::nullopt;  // we succeeded
}

}  // namespace grtest

#endif  // GRTESTUTILS_PARAM_HPP
