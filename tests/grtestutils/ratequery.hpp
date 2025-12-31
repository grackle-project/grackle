//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define assorted machinery to help with testing the ratequery machinery
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_RATEQUERY_UTILS_HPP
#define GRTESTUTILS_RATEQUERY_UTILS_HPP

#include "grackle.h"
#include "status_reporting.h"

#include <algorithm>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace grtest {

/// returns the rateid used to denote invalid rate names
inline grunstable_rateid_type get_invalid_rateid(
    const chemistry_data_storage* my_rates) {
  return grunstable_ratequery_id(my_rates, nullptr);
}

inline grunstable_rateid_type get_invalid_rateid(
    const grtest::GrackleCtxPack& pack) {
  return get_invalid_rateid(pack.my_rates());
}

inline std::string stringify_prop_kind(
    enum grunstable_ratequery_prop_kind kind) {
  switch (kind) {
    case GRUNSTABLE_QPROP_NDIM:
      return "GRUNSTABLE_QPROP_NDIM";
    case GRUNSTABLE_QPROP_SHAPE:
      return "GRUNSTABLE_QPROP_SHAPE";
    case GRUNSTABLE_QPROP_MAXITEMSIZE:
      return "GRUNSTABLE_QPROP_MAXITEMSIZE";
    case GRUNSTABLE_QPROP_WRITABLE:
      return "GRUNSTABLE_QPROP_WRITABLE";
    case GRUNSTABLE_QPROP_DTYPE:
      return "GRUNSTABLE_QPROP_DTYPE";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

inline std::string stringify_type(enum grunstable_types kind) {
  switch (kind) {
    case GRUNSTABLE_TYPE_F64:
      return "GRUNSTABLE_TYPE_F64";
    case GRUNSTABLE_TYPE_STR:
      return "GRUNSTABLE_TYPE_STR";
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

inline std::optional<enum grunstable_types> safe_type_enum_cast(long long val) {
  if (static_cast<long long>(GRUNSTABLE_TYPE_F64) == val) {
    return std::make_optional(GRUNSTABLE_TYPE_F64);
  } else if (static_cast<long long>(GRUNSTABLE_TYPE_STR) == val) {
    return std::make_optional(GRUNSTABLE_TYPE_STR);
  } else {
    return std::nullopt;
  }
}

/// A helper function that is used to help implement PrintTo for
/// ExpectedRateProperties and RateProperties
template <typename T>
void print_standard_props_(const T& props, std::ostream* os) {
  *os << "shape={";
  for (std::size_t i = 0; i < props.shape.size(); i++) {
    if (i != 0) {
      *os << ", ";
    }
    *os << props.shape[i];
  }
  *os << "}, dtype=" << stringify_type(props.dtype)
      << ", writable=" << props.writable;
}

/// summarizes details about rate properties
struct RateProperties {
  std::vector<long long> shape;
  std::size_t maxitemsize;
  enum grunstable_types dtype;
  bool writable;

  long long n_items() const {
    long long n_items = 1LL;
    for (std::size_t i = 0; i < shape.size(); i++) {
      n_items *= shape[i];
    }
    return n_items;  // this is correct even for scalars
  }

  /// teach googletest how to print this type
  friend void PrintTo(const RateProperties& props, std::ostream* os) {
    *os << "RateProperties{";
    print_standard_props_(props, os);
    *os << ", maxitemsize=" << props.maxitemsize << '}';
  }
};

/// analogous to RateProperties, but it doesn't have the maxitemsize member
///
/// The premise is that you construct this to express expectations since you
/// can't know the maxitemsize for an arbitrary array of strings
///
/// @todo
/// This is an ugly kludge. We may want to replace this with some kind of
/// custom "matcher"...
struct ExpectedRateProperties {
  std::vector<long long> shape;
  enum grunstable_types dtype;
  bool writable;

  /// teach googletest how to print this type
  friend void PrintTo(const ExpectedRateProperties& props, std::ostream* os) {
    *os << "ExpectedRateProperties{";
    print_standard_props_(props, os);
    *os << '}';
  }
};

inline bool operator==(const RateProperties& a, const RateProperties& b) {
  return a.maxitemsize == b.maxitemsize && a.shape == b.shape &&
         a.dtype == b.dtype && a.writable == b.writable;
}

static bool has_appropriate_maxitemsize_(const RateProperties& props) {
  switch (props.dtype) {
    case GRUNSTABLE_TYPE_F64:
      return props.maxitemsize == sizeof(double);
    case GRUNSTABLE_TYPE_STR:
      return true;
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

inline bool operator==(const RateProperties& a,
                       const ExpectedRateProperties& b) {
  return has_appropriate_maxitemsize_(a) && a.shape == b.shape &&
         a.dtype == b.dtype && a.writable == b.writable;
}

inline bool operator==(const ExpectedRateProperties& a,
                       const RateProperties& b) {
  return b == a;
}

/// construct a RateProperties instance for the specified rate
///
/// returns an empty optional if any there are any issues. to be clear, this
/// function shouldn't actually be used to test whether the property-querying
/// machinery works (it just tries to gracefully handle cases issues to avoid
/// distracting error-messages)
///
/// @note
/// After we query each property, we encode an extra sanity-check. We **only**
/// encode this in case there are other underlying issues (so that we don't end
/// up with an extremely crazy set of errors). The Property tests should
/// generally provide more details about these tests. To be clear, user-code
/// should never need to include these sanity check!
inline std::optional<RateProperties> try_query_RateProperties(
    chemistry_data_storage* my_rates, grunstable_rateid_type rateid) {
  long long ndim = -1LL;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_NDIM,
                                &ndim) != GR_SUCCESS) {
    return std::nullopt;
  }
  if (ndim < 0LL) {
    return std::nullopt;  // sanity-check failed!
  }

  std::vector<long long> shape;
  if (ndim > 0LL) {
    shape.assign(ndim, 0LL);
    if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_SHAPE,
                                  shape.data()) != GR_SUCCESS) {
      return std::nullopt;
    }
  }
  if (std::count_if(shape.begin(), shape.end(),
                    [](long long x) { return x <= 0; }) > 0) {
    return std::nullopt;  // sanity check failed!
  }

  long long maxitemsize = -1LL;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_MAXITEMSIZE,
                                &maxitemsize) != GR_SUCCESS) {
    return std::nullopt;
  }
  if (maxitemsize <= 0LL) {
    return std::nullopt;  // sanity check failed!
  }

  long long dtype_tmp;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_DTYPE,
                                &dtype_tmp) != GR_SUCCESS) {
    return std::nullopt;
  }
  std::optional<enum grunstable_types> dtype = safe_type_enum_cast(dtype_tmp);
  if (!dtype.has_value()) {
    return std::nullopt;  // sanity check failed!
  }

  long long writable;
  if (grunstable_ratequery_prop(my_rates, rateid, GRUNSTABLE_QPROP_WRITABLE,
                                &writable) != GR_SUCCESS) {
    return std::nullopt;
  }
  if ((writable != 0LL) && (writable != 1LL)) {
    return std::nullopt;
  }

  return {RateProperties{shape, static_cast<std::size_t>(maxitemsize),
                         dtype.value(), static_cast<bool>(writable)}};
}

/// A convenience function to try to query a 1d array of strings
///
/// returns an empty optional if any there are any issues. to be clear, this
/// function shouldn't actually be used to test the ratequiery API machinery,
/// itself. It's intended to be used to query that queried strings hold
/// expected values or as a building-block in more sophisticated calculations.
inline std::optional<std::vector<std::string>> try_query_str_arr1d(
    chemistry_data_storage* my_rates, std::string name) {
  grunstable_rateid_type id = grunstable_ratequery_id(my_rates, name.c_str());
  if (id == get_invalid_rateid(my_rates)) {
    return std::nullopt;
  }

  std::optional<RateProperties> props_maybe =
      try_query_RateProperties(my_rates, id);
  if (!props_maybe.has_value()) {
    return std::nullopt;
  }
  RateProperties props = props_maybe.value();

  if ((props.shape.size() != 1) || props.dtype != GRUNSTABLE_TYPE_STR) {
    return std::nullopt;
  }

  // allocate output storage
  std::size_t bytes_per_string = props.maxitemsize;
  std::size_t max_chars_per_string = bytes_per_string - 1;
  std::size_t n_strings = static_cast<std::size_t>(props.n_items());
  // NOTE: C++ standards older than C++17 can't directly write to strings
  //       without triggering undefined behavior
  std::vector<std::string> out_arr;
  out_arr.reserve(n_strings);
  for (std::size_t i = 0; i < n_strings; i++) {
    out_arr.emplace_back(max_chars_per_string, '\0');
  }

  // now, we are actually going to make the query
  // WARNING: do NOT touch `out_arr` until after we are done with this.
  // Otherwise, undefined behavior may arise since std::string may (must?) be
  // implemented using the "small-string optimization" (SSO).
  // -> for context, most standard library C++ containers (e.g. std::vector)
  //    are implemented as a class/struct with 3 data-members:
  //    1. a pointer to heap-memory that is the storage for contained items
  //    2. the container's current size
  //    3. the container's current capacity
  // -> SSO is the idea that sufficiently small strings may directly embed the
  //    stored elements as part of the struct/class as long as the total size
  //    of the struct/class doesn't exceed the size of the container's required
  //    members (outlined above)
  // -> this is helpful since it avoids needless heap-allocations
  // -> issues arise for us when `bytes_per_string` is small enough to take
  //    advantage of the SSO. In this scenario, the characters of the strings
  //    are being tracked as part of the storage allocating by `out_arr`. If we
  //    do something (e.g. append an element) that causes `out_arr` to
  //    reallocate its storage pointer, then the pointers in tmp_buf would be
  //    invalidated
  {
    // set up the buffer to actually pass to Grackle
    std::vector<char*> tmp_buf;
    tmp_buf.reserve(n_strings);
    for (std::size_t i = 0; i < n_strings; i++) {
      tmp_buf.push_back(out_arr[i].data());
    }

    if (grunstable_ratequery_get_str(my_rates, id, tmp_buf.data()) !=
        GR_SUCCESS) {
      return std::nullopt;
    }
  }

  // since we allocated more memory than was necessary, let's resize each
  // element in out_arr so it has the appropriate size
  for (std::string& element : out_arr) {
    element.resize(std::strlen(element.c_str()));
  }

  return {out_arr};
}

/// A convenience function to try to query a 1d array of `double`s
///
/// returns an empty optional if any there are any issues. to be clear, this
/// function shouldn't actually be used to test the ratequiery API machinery,
/// itself. It's intended to be used to query that queried strings hold
/// expected values or as a building-block in more sophisticated calculations.
inline std::optional<std::vector<double>> try_query_f64_arr1d(
    chemistry_data_storage* my_rates, std::string name) {
  grunstable_rateid_type id = grunstable_ratequery_id(my_rates, name.c_str());
  if (id == get_invalid_rateid(my_rates)) {
    return std::nullopt;
  }

  std::optional<RateProperties> props_maybe =
      try_query_RateProperties(my_rates, id);
  if (!props_maybe.has_value()) {
    return std::nullopt;
  }
  RateProperties props = props_maybe.value();

  if ((props.shape.size() != 1) || props.dtype != GRUNSTABLE_TYPE_F64) {
    return std::nullopt;
  }
  std::vector<double> out_arr;
  out_arr.assign(props.n_items(), -1.0);
  if (grunstable_ratequery_get_f64(my_rates, id, out_arr.data()) !=
      GR_SUCCESS) {
    return std::nullopt;
  }

  return {out_arr};
}

}  // namespace grtest

#endif  // GRTESTUTILS_RATEQUERY_UTILS_HPP
