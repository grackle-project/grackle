/// @file utils.hpp
/// @brief Header file defining various C++ utilities used in the routines
///        transcribed from Fortran

#ifndef UTILS_HPP
#define UTILS_HPP

#ifndef __cplusplus
#error "This file must be used by a c++ compiler"
#endif

#include "fortran_func_decls.h" // gr_mask_type

#include <cmath>
#include <cstdio> // printf
#include <type_traits> // std::is_floating_point_v, std::is_same_v

// ---------------------------------------------
// first, we define some generally useful macros
// ---------------------------------------------
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
#if defined (__GNUC__)
//#define GRIMPL_RESTRICT __restrict__
#define GRIMPL_RESTRICT /* ... */
#else
#define GRIMPL_RESTRICT /* ... */
#endif

// -----------------------------------------------------------------------
// define functions/function-like-macros for internal unrecoverable errors
// -----------------------------------------------------------------------
// -> this section could be removed in the near future

namespace grackle::impl {

/// helper function that prints an error message & aborts the program.
///
/// This is called by function-like macros. Instead macros
///
/// @note
/// On the next line ``[[noreturn]]`` is syntax used to introduce a function
/// attribute (the double-bracket syntax was supported by C++11 and C23).
/// In this case, ``noreturn`` informs the compiler that the function will not
/// return (i.e. the program exits). Without this information the compiler may
/// produce spurious compiler warnings.
[[noreturn]] void Abort_With_Err_(const char* func_name, const char* file_name, int line_num, const char* msg, ...);

} // namespace grackle::internal

/// @def      __GRIMPL_FUNC__
/// @brief    a magic contant like __LINE__ or __FILE__ used to specify the name
///           of the current function
///
/// In more detail:
/// - The C++11 standard ensures __func__ is provided on all platforms, but it
///   only provides limited information (just the name of the function).
/// - note that __func__ is technically not a macro. It's a static constant
///   string implicitly defined by the compiler within each function definition
/// - Where available, we prefer to use compiler-specific features that provide
///   more information about the function (like the scope of the function, the
///   the function signature, any template specialization, etc.).
#ifdef __GNUG__
  #define __GRIMPL_PRETTY_FUNC__ __PRETTY_FUNCTION__
#else
  #define __GRIMPL_PRETTY_FUNC__ __func__
#endif


/// @def GRIMPL_ERROR
/// @brief function-like macro that handles a (lethal) error message
///
/// This macro should be treated as a function with the signature:
///
///   [[noreturn]] void GRIMPL_ERROR(const char* fmt, ...);
///
/// The ``msg`` arg is printf-style format argument specifying the error
/// message. The remaining args arguments are used to format error message
///
/// @note
/// the ``msg`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``msg`` doesn't format
/// any arguments). There is no portable way around this until C++ 20.
///
/// @note
/// the ``msg`` string is part of the variadic args so that there is always
/// at least 1 variadic argument (even in cases when ``msg`` doesn't format
/// any arguments). There is no portable way around this until C++ 20.
#define GRIMPL_ERROR(...)                                                      \
  { grackle::impl::Abort_With_Err_                                             \
      (__GRIMPL_PRETTY_FUNC__, __FILE__, __LINE__, __VA_ARGS__); }

/// @def GRIMPL_REQUIRE
/// @brief implements functionality analogous to the assert() macro
///
/// if the condition is false, print an error-message (with printf
/// formatting) & abort the program.
///
/// This macro should be treated as a function with the signature:
///
///   void GRIMPL_REQUIRE(bool cond, const char* fmt, ...);
///
/// - The 1st arg is a boolean condition. When true, this does nothing
/// - The 2nd arg is printf-style format argument specifying the error message
/// - The remaining args arguments are used to format error message
///
/// @note
/// The behavior is independent of the ``NDEBUG`` macro
#define GRIMPL_REQUIRE(cond, ...)                                              \
  {  if (!(cond))                                                              \
      { grackle::impl::Abort_With_Err_                                         \
         (__GRIMPL_PRETTY_FUNC__, __FILE__, __LINE__, __VA_ARGS__); } }

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


// ---------------------------------------------------
// define machinery related to the View class template
// ---------------------------------------------------

namespace grackle::impl {

/// this simply helps us implement the View class template
template <typename T>
struct MDPtrProps_ {
  // try to strip off 3-pointer layers (it's ok to have fewer layers)
  using strip0_ = std::remove_pointer_t<T>;
  using strip1_ = std::remove_pointer_t<strip0_>;
  using type = std::remove_pointer_t<strip1_>;

  static constexpr int rank = (
    1 + int(std::is_pointer_v<strip0_>) + int(std::is_pointer_v<strip1_>)
    + int(std::is_pointer_v<type>)
  );

  static_assert(0 < rank && rank <4, "template arg must be 1D, 2D, or 3D");
};

/// Implements the View class template
///
/// @tparam T defines the datatype and dimensionality of the span. This is
///     easiest to explain with examples.
///       - ``float**`` specifies a 2D View of ``float`` values
///       - ``const double***`` specifies a 3D View of ``const double`` values
///     Be aware that the use of multiple of multiple pointer indirection is
///     purely a symbolic shorthand. Under the hood, ``int*`` is used whether
///     this parameter is ``int*``, ``int**``, or ``int***``.
///
/// Overview
/// --------
/// This file holds the declaration/definition of the grackle::Imple::View
/// class template. You should think of instances of the class template as a
/// special kind of pointer:
///  - instances can be empty (i.e. they encode a nullptr) or store the address
///    to memory that is used to store a multidimensional array
///  - instances also track the shape of the multidimensional array in order 
///    to support. For example, an to access the value stored at index 
///    ``i``, ``j``, ``k`` in data represented by an instance ``view`` using
///    ``view(i,j,k)`` where ``i`` is the index along the fast access (we could
///    add support for customizing data layout in the future if we deem it
///    useful)
///  - it has pointer semantics (more on this below)
///
/// The idea of a ``View`` is common in various C++ HPC libraries (e.g. see
/// Kokkos or Raja). Enzo-E makes use of a similar ``CelloView``. In modern C++
/// lingo this is a kind of span. If we were using C++23, we might use
/// std::mdspan instead of defining a custom type.
///
/// Motivation
/// ----------
/// This class template primarily exists to ease the process of transcribing
/// fortran logic multidimensional arrays.
///
/// More detailed description
/// -------------------------
/// It is important to understand that this has all of the semantics of a
/// pointer and not the "value semantics" of a container implemented by
/// C++'s standard library.
///
/// Let's consider a few key scenarios:
/// 1. Const-semantics:
///    - When you have a View of ``const`` values, ``View<const int*>``, you
///      can't modify the values, similar to ``const int*``. There is no
///      container equivalent (e.g. ``std::vector<const int>`` doesn't exist)
///    - Like a variable holding a ``const`` pointer to an integer,
///      ``int const *``, a variable holding a ``const View<int*>`` can be used
///      to freely modify referenced values. In both cases ``const`` means that
///      the properties of the array (shape/memory-address) can't change. (In
///      contrast, you can't mutate elements of a ``const std::vector<int>``)
/// 2. Copying/assignment:
///    - Like with copying a pointer, copying a View just copies the properites
///      of the underlying data (shape/memory-address). If you store a copy of
///      ``view_a`` and store it in a variable ``view_b``, then ``view_a`` and
///      ``view_b`` can be used to access/modify data at the same memory
///      location. If ``view_b`` previously held information about a different
///      view, the act of copying has no impact on the values in the old view.
///      (Both behaviors contrast with ``std::vector`` where copy operations
///      always involve making a deepcopy).
///
/// At this time, a View can not allocate its own data. For example of how to
/// do this, see Enzo-E's CelloView class template.
///
/// Considerations
/// --------------
/// In the long-term, it may make sense to use a ``View`` template class that
/// wraps Kokkos::View. It is a relatively elegant way to attach information
/// about where memory is allocated.
/// - We might want to remove all use of this as a 2D View
/// - We might also want to remove all use of this as a 3D View and just use it
///   as a 1D View (it depends on our thoughts about self-shielding)
///
template<typename T>
struct View {

  // first, we define useful types used by instances of the class template
private:
  using ptrprops_ = MDPtrProps_<T>;
public:
  /// the element type
  using element_type = typename ptrprops_::type;
  using reference_type = element_type &;
  using size_type = int; // maybe revisit this?
  static constexpr int rank = ptrprops_::rank;

private: // attributes
  element_type* data_;
  size_type extent_[rank];
  size_type strides_[rank];

private: // helper methods

  void check_invariants_() const {
    if ((data_ != nullptr) && extent_[0] <= 0) {
      GRIMPL_ERROR("ilen can't be 0 for non nullptr");
    }
    for (int i = 0; i < rank; i++) {
      GRIMPL_REQUIRE(extent_[i] >= 0, "extent can't be negative");
    }
  }

public:
  
  /// Default constructor
  ///
  /// @note
  /// The syntax ensures that the contents of extent_ are all initialized to zero
  /// since extent_ is an array of non-class types
  View() : data_{nullptr}, extent_{}, strides_{} {}

  ///@{
  /// Construct a view from an existing pointer `ptr`. Every arg after the pointer
  /// specifies the extent of an access (from the fastest axis to the slowest axis)
  View(element_type* ptr, int ilen)
    : data_(ptr),
      extent_{ilen},
      strides_{1}
  {
    // we may need to enforce this check in a different way
    static_assert(rank==1, "constructor only works with 1D views");
    check_invariants_();
  }

  View(element_type* ptr, int ilen, int jlen)
    : data_(ptr),
      extent_{ilen, jlen},
      strides_{1, ilen}
  {
    static_assert(rank==2, "constructor only works with 2D views");
    check_invariants_();
  }

  View(element_type* ptr, int ilen, int jlen, int klen)
    : data_(ptr),
      extent_{ilen, jlen, klen},
      strides_{1, ilen, ilen*jlen}
  {
    static_assert(rank==3, "constructor only works with 3D views");
    check_invariants_();
  }
  ///@}

  // explicitly use defaults for a handful of cases
  ~View()=default;
  View(const View&)=default; // copy constructor
  View(View&&)=default; // move constructor
  View& operator=(const View&)=default; // copy assignment
  View& operator=(View&&)=default; // move assignment

  element_type* data() const noexcept { return data_; }
  size_type extent(int i) const {
    GRIMPL_REQUIRE(i >= 0 && i <= rank,
                   "i must be non-negative and can't exceed %d", rank-1);
    return extent_[i];
  }

  ///@{
  /// Implements multi-dimensional indexing. The first argument corresponds
  /// to the contiguous axis
  GRIMPL_FORCE_INLINE element_type& operator()(int i) const {
    static_assert(rank==1, "1 index should only be specified for 1D views");
    return data_[i]; // strides_[0] == 1
  }

  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j) const {
    static_assert(rank==2, "2 indices should only be specified for 2D views");
    return data_[i + j*strides_[1]]; // strides_[0] == 1
  }

  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j, int k) const {
    static_assert(rank==3, "3 indices should only be specified for 3D views");
    return data_[i + j*strides_[1] + k*strides_[2]]; // strides_[0] == 1
  }
  ///@}

};

} // namespace grackle::internal

#endif /* UTILS_HPP */
