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

// we may want to find a better home for some of these macros in the future
// - the central theme for all of these macros is that they are totally
//   unrelated to physics

/// @def OMP_PRAGMA
/// Macro used to wrap OpenMP's pragma directives.
///
/// When the program:
///  * is compiled with OpenMP, the pragma contents are honored.
///  * is NOT compiled with OpenMP, the pragma contents are ignored.
///
/// @note
/// This macro is implemented using the ``_Pragma`` operator, described
/// [here](https://en.cppreference.com/w/cpp/preprocessor/impl). More details
/// can be found [here](https://gcc.gnu.org/onlinedocs/cpp/Pragmas.html).
#ifdef _OPENMP
#define OMP_PRAGMA(x) _Pragma(#x)
#else
#define OMP_PRAGMA(x) /* ... */
#endif

#ifdef _OPENMP
#define OMP_PRAGMA_CRITICAL _Pragma("omp critical")
#else
#define OMP_PRAGMA_CRITICAL /* ... */
#endif

/// @macro GRIMPL_FORCE_INLINE
/// @brief replacement for ``inline`` that forces inlining (on some compilers)
///
/// This macro should be used sparingly (if you force inlining of too much, you
/// will slow down the code. On unsupported macros, this just becomes the
/// ``inline`` keyword
///
/// While the initial purpose of ``inline`` keyword was to encourage optimizers
/// to inline annotated functions, they are free to ignore this hint.
/// - In C++, the ``inline`` keyword has come to mean that a single function
///   can be defined in separate translation units (the definitions must be the
///   same unless it is also a ``static`` function).
///   - This behavior is what allows you to put the definition of an ordinary,
///     ``inline`` function in a header and safely include that file in
///     multiple source files (function template behave kinda like an
///     ``inline`` function, even without the keyword).
///   - Unless you are using link-time-optimization, you **need** to follow
///     this pattern for defining a given function if you want to the optimizer
///     to be able to inline calls to that function in multiple source files.
///     (But, the optimizer is still not to inline).
/// - Be aware, the semantics of the ``inline`` keyword are slightly different
///   in C (C actually adopted the keyword from C++).
///   - The easiest way to achieve roughly the same behavior as C++ (with
///     putting the definitions of ``inline`` functions in headers that are
///     included in multiple source files) is to make the function
///     ``static inline``
///   - in C, ``inline`` functions, that aren't ``static inline``, also have
///     special rules/restrictions about using ``static`` local variables
///     inside of the definition.
#if defined(__GNUC__)
#define GRIMPL_FORCE_INLINE __attribute__((always_inline)) inline
#else
#define GRIMPL_FORCE_INLINE inline
#endif

/// @macro GRIMPL_RESTRICT
/// @brief Equivalent to C99's ``restrict`` qualifier (on supported compilers)
///
/// Background
/// ----------
/// Fortran's reputation for producing faster numerical code involving arrays
/// than C/C++ comes from the fact that Fortran forbids (by default) overlap
/// in the memory regions for arrays that are represented by distinct
/// variables. that use to store array elements. A C/C++ compiler typically
/// must be much more conservative about optimizations because pointers are
/// allowed to freely overlap.
/// - consider a function that performs loops over pairs of elements in arrays
///   tracked by pointers ``a`` and ``b`` and write the result to an array
///   tracked by pointer ``c``.
/// - Outside of special cases, the compiler commonly has to assume that the
///   memory locations accessed through pointer ``c`` may overlap with the
///   elements accessed through ``a`` or ``b``. In that case, the precise order
///   of the loop evaluation (and the timing of when values are accessed from
///   ``a`` or ``b`` or are written to ``c``) matters.
/// - This significantly inhibits auto-vectorization.
/// - (There is also overhead with pointers in general)
///
/// C99's ``restrict`` qualifier
/// ----------------------------
/// To facillitate faster code, the C99 standard introduced the ``restrict``
/// qualifier that can be attached to a pointer to object-type (e.g. a pointer
/// to an ``int`` or ``float``)
/// - When developers add this qualifier to a pointer, ``ptr``, they are making
///   making a promise to the compiler for a given scope. They essentially
///   promise that if an object accessible through ``ptr`` is modified in that
///   scope, all access/modifications to that object are performed with ``ptr``
/// - This facillitates Fortran-like optimizations (the compiler can also
///   choose to ignore the information)
/// - Obviously, violating this promise can obviously produces all manner of
///   undefined behavior.
///
/// C++ does **NOT** provide a ``restrict`` qualifier, but many compilers
/// support it as an extension (usually called ``__restrict__``). This macro is
/// intended to expand to ``__restrict__`` where applicable)
///
/// > [!note]
/// > The ``restrict`` qualifier is a lesser known and lesser used feature of
/// > C. Consequently, machinery in backends of C/C++ compilers that use this
/// > information to make optimization may not be rigorously tested unless they
/// > are also used as backends for other languages that more regularly perform
/// > this kind of optimization.
/// > - gcc's backend is probably pretty robust since it is also used by
/// >   gfortran
/// > - in contrast, bugs were found (they've been fixed) over the last several
/// >   years in LLVM (clang's backend) by developers of the Rust lanugague.
/// >
/// > With that in mind, we currently just use this macro as a placeholder that
/// > we can always try later.
/// > - I think it's useful to use this in transcribed code so that we properly
/// >   retain as much contextual information as possible (it's easy to remove
/// >   later!). The Fortran code is written in such a way that we can apply
/// >   this very liberally.
/// > - To be clear, I have a lot more faith in enabling this feature than
/// >   passing the ``-ffast-math`` flag to gcc (__restrict__ semantics are
/// >   well defined and its opt-in)
#if defined(__GNUC__)
// #define GRIMPL_RESTRICT __restrict__
#define GRIMPL_RESTRICT /* ... */
#else
#define GRIMPL_RESTRICT /* ... */
#endif

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
