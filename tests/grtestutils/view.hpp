//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares functionality for working with/testing multidimensional "views"
///
/// For additional context, we use the word "view" in the same sense as
/// std::string_view or a std::span (it is a view onto externally owned data)
///
//===----------------------------------------------------------------------===//

#ifndef GRTESTUTILS_VIEW_HPP
#define GRTESTUTILS_VIEW_HPP

#include <cstddef>  // std::size_t
#include <string>

namespace grtest {

/// equivalent of converting output of printf("%g", val) to std::string
std::string to_pretty_string(float val);
std::string to_pretty_string(double val);

/// formats a pointer as a string
std::string ptr_to_string(const float* ptr, std::size_t len);
std::string ptr_to_string(const double* ptr, std::size_t len);

}  // namespace grtest

#endif  // GRTESTUTILS_VIEW_HPP