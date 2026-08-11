//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares various C++ utilities used in routines transcribed from Fortran
///
//===----------------------------------------------------------------------===//

#ifndef UTILS_HPP
#define UTILS_HPP

#ifndef __cplusplus
#error "This file must be used by a c++ compiler"
#endif

#include "fortran_func_decls.h" // gr_mask_type
#include "support/status_reporting.hpp"

// for historical reasons, some files that include the current header are
// looking for macros defined in the following pair of headers
#include "support/config.hpp"
#include "support/View.hpp"

#include <cmath>
#include <cstdio> // printf
#include <type_traits> // std::is_floating_point_v, std::is_same_v


// ---------------------------------------------
// define some functions used in the translation
// ---------------------------------------------

/// convenience function that acts like printf, but prints to stderr
///
/// @note
/// It may make sense to move this type to the C layer (this is the primary
/// reason that we don't put it inside a namespace)
int eprintf(const char* format, ...);

namespace grackle::impl {

/// crude implementation of some logic for printing arrays
///
/// @note
/// The choice to have a switch statements that pass varying numbers of
/// arguments to printf, rather than formatting with sprintf or something
/// similar was motivated by our desire to write a simple function that will
/// probably work when compiled as CUDA (to my knowledge, CUDA doesn't have
/// other flavors of printf).
///
/// @par
/// If we deem this function to be useful, we could come up with something
/// better and more efficient (in that case, we may end up reimplemnting some
/// formatting options by hand to support it on GPUs). It seems more likely
/// that we would just avoid this function on GPUs.
template<typename T>
void print_contiguous_row_(const T* ptr, int start_idx, int stop_idx) {
  // for gr_mask_type, we cast to `int` before passing to printf
  using castT = std::conditional_t<std::is_same_v<T, gr_mask_type>, int, T>;

  const int max_elem = 7;
  const char* fmtline; // string with max_elem occurences of a formatter
                       // (delimited by ' ') & a trailing \n
  int fmtstep; // width of a single fmtspec plus 1

  if constexpr (std::is_same_v<T, gr_mask_type>) {
    fmtline = "%d %d %d %d %d %d %d\n";
    fmtstep = 2 + 1;
  } else if constexpr (std::is_floating_point_v<T>) {
    // total fixed width is 8 larger than # of decimal digits to accound for:
    //   minus-sign, 1st digit, decimal-point,
    //   'e' (for exponent), exponent-sign, 3 exponent digits
    fmtline = "%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n";
    fmtstep = 6 + 1;
  } else {
    printf("can't print specified type\n");
    return;
  }

  if (start_idx >= stop_idx) { printf("\n"); return; }

  const T* stop = ptr + stop_idx;
  ptr += start_idx;
  while (ptr < stop) {
    const int step = (ptr+max_elem < stop) ? max_elem : (int)(stop - ptr);
    const int off = (max_elem-step) * fmtstep;
    switch (step) {
      case 7:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]), castT(ptr[2]),
               castT(ptr[3]), castT(ptr[4]), castT(ptr[5]), castT(ptr[6]));
        break;
      case 6:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]),
               castT(ptr[2]), castT(ptr[3]), castT(ptr[4]), castT(ptr[5]));
        break;
      case 5:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]),
               castT(ptr[2]), castT(ptr[3]), castT(ptr[4]));
        break;
      case 4:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]),
               castT(ptr[2]), castT(ptr[3]));
        break;
      case 3:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]),
               castT(ptr[2]));
        break;
      case 2:
        printf(fmtline + off, castT(ptr[0]), castT(ptr[1]));
        break;
      case 1:
        printf(fmtline + off, castT(ptr[0]));
    }
    ptr += step;
  }
}

///@{
/// Implements alternatives to fmax that take more args
///
/// @note
/// These functions **only** exists to support transcription
inline float fmax(float a, float b, float c) {
  return std::fmax(a, std::fmax(b, c));
}
inline double fmax(double a, double b, double c) {
  return std::fmax(a, std::fmax(b, c));
}
///@}


///@{
/// Implements alternatives to fmin that take more args
///
/// @note
/// These functions **only** exists to support transcription
inline float fmin(float a, float b, float c) {
  return std::fmin(a, std::fmin(b, c));
}
inline float fmin(float a, float b, float c, float d) {
  return std::fmin(std::fmin(a, b), std::fmin(c, d));
}
inline double fmin(double a, double b, double c) {
  return std::fmin(a, std::fmin(b, c));
}
inline double fmin(double a, double b, double c, double d) {
  return std::fmin(std::fmin(a, b), std::fmin(c, d));
}
///@}


/// Crude implementation of Fortran's
/// [MOD function](https://gcc.gnu.org/onlinedocs/gfortran/MOD.html)
///
/// @note
/// We currently only implement behavior for positive lengths
inline int mod(int a, int p) {
  GRIMPL_REQUIRE((a >= 0) && (p>0),
      "a must be non-negative and p must be positive");
  return a % p;
}

/// Implementation of Fortran's
/// [DAbs intrinsic](https://gcc.gnu.org/onlinedocs/gcc-3.4.6/g77/DAbs-Intrinsic.html)
inline double dabs(double a) { return std::fabs(a); }

} // namespace grackle::internal

// ----------------------
// other useful machinery
// ----------------------

namespace grackle::impl {

///@{
/// Implements alternatives to std::clamp that always works with floating point
/// values
///
/// For context, std::clamp, produces undefined behavior when a NaN is passed
/// into the function
inline int clamp(int v, int lo, int hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
inline long clamp(long v, long lo, long hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
inline long long clamp(long long v, long long lo, long long hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
inline float clamp(float v, float lo, float hi) {
  return std::fmax(lo, std::fmin(v, hi));
}
inline double clamp(double v, double lo, double hi) {
  return std::fmax(lo, std::fmin(v, hi));
}
///@}

} // namespace grackle::impl

#endif /* UTILS_HPP */
