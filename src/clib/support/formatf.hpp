//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare machinery related to printf-style string formatting
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_FORMATF_HPP
#define SUPPORT_FORMATF_HPP

#include <cstdarg>
#include <string>

#include "./config.hpp"

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Construct a std::string from a printf-style message
///
/// @note
/// For less experienced c++ developers: ``[[gnu::format(printf, 1, 2)]]``
/// is a compiler-specific attribute that instructs gcc (or clang++) to check
/// at compile-time that the arguments are consistent with the printf style
/// format string ``s``. It **SHOULD** also check that ``s`` is a
/// string-literal. If a compiler doesn't recognize the attribute, it simply
/// ignores (that's mandated by the C++ standard)
[[gnu::format(printf, 1, 2)]] std::string str_formatf(const char* s, ...);

/// @brief Construct a std::string from a printf-style message using the data
///     contained by @p vlist
std::string vstr_formatf(const char* s, std::va_list vlist);

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_FORMATF_HPP