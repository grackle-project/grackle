//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define functionality for working with/testing multidimensional "views"
///
//===----------------------------------------------------------------------===//

#include "view.hpp"

#include <cstddef>  // std::size_t
#include <sstream>
#include <string>
#include <type_traits>  // std::is_floating_point_v<T>

namespace grtest {

std::string to_pretty_string(double val) {
  // the current implementation is EXTREMELY crude
  char buf[30];
  snprintf(buf, 30, "%g", val);
  return std::string(buf);
}

std::string to_pretty_string(float val) {
  // this is crude!
  return to_pretty_string(static_cast<double>(val));
}

template <typename T>
static std::string ptr_to_string_(const T* ptr, std::size_t len) {
  static_assert(std::is_floating_point_v<T>);

  std::size_t pause_start;
  std::size_t pause_stop;

  if (len > 30) {
    pause_start = 3;
    pause_stop = len - 3;
  } else {
    pause_start = len * 2;
    pause_stop = pause_start;
  }

  std::stringstream s;
  s << '{';

  for (std::size_t i = 0; i < len; i++) {
    if ((i > pause_start) && (i < pause_stop)) {
      continue;
    }

    if (i == pause_stop) {
      s << ", ... ";
    } else if (i != 0) {
      s << ", ";
    }

    const int BUF_SIZE = 30;
    char buf[BUF_SIZE];
    snprintf(buf, BUF_SIZE, "%g", ptr[i]);
    s << buf;
  }
  s << '}';
  return s.str();
}

std::string ptr_to_string(const float* ptr, std::size_t len) {
  return ptr_to_string_(ptr, len);
}
std::string ptr_to_string(const double* ptr, std::size_t len) {
  return ptr_to_string_(ptr, len);
}

}  // namespace grtest