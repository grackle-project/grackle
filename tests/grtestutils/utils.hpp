// This defines a series of utilities to aid with testing grackle

#ifndef GRTEST_UTILS_HPP
#define GRTEST_UTILS_HPP

#include <grackle.h>

#include <random>
#include <optional>
#include <string>
#include <type_traits>  // std::is_floating_point_v

namespace grtest {

/// this function returns the desired standard datafile. It deals with the minutia of making sure that
/// grackle can find the standardized data-file
///
/// @returns An empty optional if unsuccessful
std::optional<std::string> get_standard_datafile(const char* datafile);

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


namespace grtest {

/// equivalent of %g (the value should be fully specified)
/// 
/// The current implementation is extremely crude!
std::string to_pretty_string(double val);

inline std::string to_pretty_string(float val) {
  // this is crude!
  return to_pretty_string(static_cast<double>(val));
}

/// formats a std::vector as a string
///
/// @note
/// This is highly inefficient, partially because it consists of code written
/// from before we adopted googletest
template <typename T>
std::string ptr_to_string(const T* ptr, std::size_t len) {
  static_assert(std::is_floating_point_v<T>);

  std::string out = "{";

  std::size_t pause_start;
  std::size_t pause_stop;

  if (len > 30){
    pause_start = 3;
    pause_stop = len - 3;
  } else {
    pause_start = len *2;
    pause_stop = pause_start;
  }

  for (std::size_t i = 0; i < len; i++) {
    if ((i > pause_start) && (i < pause_stop)) { continue; }

    if (i == pause_stop) {
      out += ", ... ";
    } else if (i != 0) {
      out += ", ";
    }

    const int BUF_SIZE = 30;
    char buf[BUF_SIZE];
    snprintf(buf, BUF_SIZE, "%g", ptr[i]);
    out += buf;
  }
  return out + "}";
}

}  // namespace grtest

#endif /* GRTEST_UTILS_HPP */
