//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// This is a lightweight header file that defines a few macros that are used
/// throughout the codebase
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_CONFIG_HPP
#define SUPPORT_CONFIG_HPP

/// Expands to Grackle's global internal namespace enclosing the implementation
///
/// This is used when referring to the qualified name of a PREVIOUSLY DECLARED
/// entity (a function, class/struct, etc.) that is a member of this namespace.
///
/// For example, suppose that a function <tt> int foo(const char* s) </tt> was
/// as part of this namespace. From global scope, we could use this macro to
/// call the function with <tt> GRIMPL_NS::foo(s) </tt>.
///
/// @important
/// This macro should NOT be used to open namespace blocks. You should use
/// GRIMPL_NAMESPACE_DECL for that purpose.
///
/// @par Motivation
/// The existence of this macro has 3 motivations:
/// 1. Once we use it everywhere, the name will be easier to change. The current
///    choice, <tt> grackle::impl </tt> was picked because it was relatively
///    short (there are probably slightly more descriptive alternatives).
/// 2. Relatedly, it might be convenient to encode the current version number
///    in the namespace. This comes up in the context of header-only libraries.
///    (It would also open the door to linking and testing 2 Grackle versions
///    in the same program).
/// 3. It also goes hand-in-hand with the GRIMPL_NAMESPACE_DECL macro. See its
///    docstring for the benefits that it brings.
#define GRIMPL_NS grackle::impl

/// Used for opening the namespace block of Grackle's global internal namespace
/// that is used to enclose Grackle's implementation
///
/// To declare a function <tt> int foo(const char* s) </tt> in this namespace,
/// you might write something like:
/// @code{.cpp}
/// #include "support/config.hpp"
///
/// namespace GRIMPL_NAMESPACE_DECL {
/// int foo(const char* s)
/// }  //  namespace GRIMPL_NAMESPACE_DECL
/// @endcode
///
/// @par Motivation
/// This is primarily inspired by a similar macro, \c LIBC_NAMESPACE_DECL from
/// LLVM-libc (LLVM's implementation of the C standard library), which is
/// described at https://libc.llvm.org/dev/code_style.html#libc-namespace-decl
///
/// @par
/// Once we adopt this everywhere, we would change this macro's definition to
/// <tt> [[gnu::visibility("hidden")]] GRIMPL_NS </tt> when compiling Grackle
/// as a shared library, which would declare all enclosed symbols as having
/// hidden visibility. A brief primer on symbol visibility is provided by
/// https://cs.dartmouth.edu/~sergey/cs258/ABI/UlrichDrepper-How-To-Write-Shared-Libraries.pdf
/// This would effectively reduce the size of the Grackle shared library,
/// improve startup time of programs linked against a Grackle shared library
/// and enforce that private functions shouldn't be used outside of Grackle
#define GRIMPL_NAMESPACE_DECL GRIMPL_NS

#endif  // SUPPORT_CONFIG_HPP
