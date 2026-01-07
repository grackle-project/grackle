// This defines a series of utilities to aid with testing grackle

#ifndef GRTEST_UTILS_HPP
#define GRTEST_UTILS_HPP

#include <random>
#include <grackle.h>
#include <string>
#include <string_view>

namespace grtest {

/// Returns whether the `s` begins with the provided prefix
///
/// @note
/// When we start using C++ 20, this can be replaced with the `starts_with`
/// member functions of `std::string` and `std::string_view`
inline bool starts_with(std::string_view s, std::string_view prefix) {
#if defined(__cpp_lib_starts_ends_with) && __cpp_lib_starts_ends_with >= 201711L
  return s.starts_with(prefix);
#else
  return s.substr(0, prefix.size()) == prefix;
#endif
}

inline bool starts_with(std::string_view s, const char* prefix) {
  return starts_with(s, std::string_view(prefix));
}

inline bool starts_with(const std::string& s, std::string_view prefix) {
  return starts_with(std::string_view(s), prefix);
}

inline bool starts_with(const std::string& s, const char* prefix) {
  return starts_with(std::string_view(s), prefix);
}

/// this function records the desired standard datafile within the
/// chemistry_data struct. It deals with the minutia of making sure that
/// grackle can find the standardized data-file
///
/// @returns true if successful or false if unsuccessful
///
/// @note
/// For the sake of forward compatability (we will probably change the
/// implementation if we merge PRs 235, 237, and 246) you should:
/// - only pass string literals (or the addresses of string-literals) as the
///   this function's datafile arg
/// - AND never deallocate my_chemistry.datafile after calling this function
///   (there won't be a memory leak)
bool set_standard_datafile(chemistry_data& my_chemistry, const char* datafile);

}

/// In this namespace, we define tools to help with drawing random numbers
/// in a reproducible manner (for a given seed) across platforms
///
/// In this regard, there issues with common built-in choices:
/// -> functions like drand48 apparently aren't particularly portable
///    (https://charm.readthedocs.io/en/latest/faq/manual.html#how-much-does-getting-a-random-number-generator-right-matter)
/// -> machinery in the C++ standard library for sampling distributions (e.g.
///    std::normal_distribution or std::uniform_real_distribution) are allowed
///    to produce differente sequences of output values (for a given generator
///    and seed) in different implementations of the standard library
///    -> https://stackoverflow.com/a/24554535
///
/// The solution:
/// -> the pseudo-random number generators themselves are very portable in the
///    C++ standard library
/// -> consequently, this namespace provides functions sample distributions
///    from a given generator so we can be sure that the mapping is portable
namespace grtest::random {

/// return randomly drawn double-precision value drawn from the uniform 
/// distribution over the interval [0.0, 1.0).
inline double uniform_dist_transform(std::minstd_rand &generator) {

  // this static_assert breaks things on macOS's apple-clang compiler
  //static_assert((generator.max() <= UINT32_MAX) && (generator.min() == 1), 
  //              "Unexpected PRNG Property"); // sanity-check!

  // cast to double since they can perfectly represent all values of uint32_t
  double raw = static_cast<double>(generator()) - 1.0;
  double range = static_cast<double>(generator.max());
  return raw / range;
}

} // namespace grtest::random

#endif /* GRTEST_UTILS_HPP */
