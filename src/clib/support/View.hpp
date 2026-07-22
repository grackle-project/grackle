//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the GeneralView construct
///
//===----------------------------------------------------------------------===//

#ifndef SUPPORT_VIEW_HPP
#define SUPPORT_VIEW_HPP

#include "support/config.hpp"
#include "support/status_reporting.hpp"

#include <type_traits>  // std::remove_pointer_t, std::is_pointer_v

// A brief discussion on Views and Data Layout
// ===========================================
//
// Background
// ----------
// For some background, C++23 introduced `std::mdspan` to describe
// multi-dimensional views. A `std::mdspan` is parameterized by
// - the data's extents (aka the shape)
// - the data's layout, which dictates how a multidimensional index is mapped
//   to a 1D pointer offset
//
// For views of contiguous data there are 2 obvious layouts:
// 1. layout-right: where the stride is `1` along the rightmost extent.
//    - for extents `{a,b,c}`, an optimal nested for-loop will iterates from
//      `0` up to `a` in the outermost loop and from `0` up to `c`
//      in the innermost loop
//    - this is the "natural layout" for a multidimensional c-style array
//      `arr[a][b][c]`
// 2. layout-left: where the stride is `1` along the leftmost extent
//    - for extents `{a,b,c}`, an optimal nested for-loop will iterates from
//      `0` up to `c` in the outermost loop and from `0` up to `a`
//      in the innermost loop
//    - this is the "natural layout" for a multidimensional fortran array
//      `arr(a, b, c)`
//
// Views In Grackle
// ----------------
// Since Grackle doesn't yet use C++23, we define our own custom view type.
//
// Because Grackle was originally written in Fortran, all 3d arrays
// (they represent spatial grids have left layout). However, there are also
// some arrays with right layout (e.g. and gaussj after transcription and
// interpolation grids). It would nice to be as consistent as possible.
// Since Grackle is now written in C++, we will be using Right Layout
//
// The Plan
// --------
// We are starting to transition away from treating fields as full 3D spatial
// grids and moving towards treating them as 1D spatial grids. While we do
// this, this is a natural place for us to start transitioning the standard
// layout type
//
// Once we finish transitioning the layout type, we should delete the
// `DataLayout` type and simplify GeneralView down to View
namespace GRIMPL_NAMESPACE_DECL {

/// @brief describes the data layout
///
/// For contiguous multidimensional view with:
/// - LEFT layout means that the leftmost component of a multidimensional index
///   is the contiguous axis. This corresponds the "natural" layout of Fortran
///   arrays
/// - RIGHT layout means that the rightmost component of a multidimensional
///   index is the contiguous axis. This corresponds the "natural" layout of
///   arrays in C/C++
enum struct DataLayout {
  LEFT,  ///< the leftmost dimension has a stride 1
  RIGHT  ///< the rightmost dimension has a stride 1
};

namespace view_detail {

// define template and some partial specializations to implement MkPtr_t
template <typename value_type, int Rank>
struct MkPtr_ {
  static_assert(0 < Rank && Rank < 4, "rank isn't 1, 2 or 3");
  using type = value_type*;
};

template <typename value_type>
struct MkPtr_<value_type, 2> {
  using type = value_type**;
};

template <typename value_type>
struct MkPtr_<value_type, 3> {
  using type = value_type***;
};

/// machinery that provides the appropriate pointer to ``value_type``
///
/// Examples: ``MkPtr_t<double, 3>::type`` is ``double***``
///           ``MkPtr_t<const int, 2>::type`` is ``const int**``
///           ``MkPtr_t<const int, 1>::type`` is ``const int*``
template <typename value_type, int Rank>
using MkPtr_t = typename MkPtr_<value_type, Rank>::type;

/// this simply helps us implement the GeneralView class template
template <typename T>
struct MDPtrProps_ {
  // try to strip off 3-pointer layers (it's ok to have fewer layers)
  using strip0_ = std::remove_pointer_t<T>;
  using strip1_ = std::remove_pointer_t<strip0_>;
  using type = std::remove_pointer_t<strip1_>;

  static constexpr int rank =
      (1 + int(std::is_pointer_v<strip0_>) + int(std::is_pointer_v<strip1_>) +
       int(std::is_pointer_v<type>));

  static_assert(0 < rank && rank < 4, "template arg must be 1D, 2D, or 3D");
};

}  // namespace view_detail

/// A general purpose multidimensional view
///
/// @tparam T defines the datatype and dimensionality of the span. This is
///     easiest to explain with examples.
///       - ``float**`` specifies a 2D GeneralView of ``float`` values
///       - ``const double***`` specifies a 3D GeneralView of ``const double``
///         values
///     Be aware that the use of multiple of multiple pointer indirection is
///     purely a symbolic shorthand. Under the hood, ``int*`` is used whether
///     this parameter is ``int*``, ``int**``, or ``int***``.
///
/// Overview
/// --------
/// You should think of instances of the class template as a special kind of
/// pointer:
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
/// The idea of a ``GeneralView`` is common in various C++ HPC libraries (e.g.
/// see Kokkos or Raja). Enzo-E makes use of a similar ``CelloGeneralView``. In
/// modern C++ lingo this is a kind of span. If we were using C++23, we might
/// use std::mdspan instead of defining a custom type.
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
///    - When you have a GeneralView of ``const`` values,
///      ``GeneralView<const int*>``, you can't modify the values, similar to
///      ``const int*``. There is no container equivalent (e.g.
///      ``std::vector<const int>`` doesn't exist)
///    - Like a variable holding a ``const`` pointer to an integer,
///      ``int const *``, a variable holding a ``const GeneralView<int*>`` can
///      be used to freely modify referenced values. In both cases ``const``
///      means that the properties of the array (shape/memory-address) can't
///      change. (In contrast, you can't mutate elements of a
///      ``const std::vector<int>``)
/// 2. Copying/assignment:
///    - Like with copying a pointer, copying a GeneralView just copies the
///      properites of the underlying data (shape/memory-address). If you store
///      a copy of ``view_a`` and store it in a variable ``view_b``, then
///      ``view_a`` and ``view_b`` can be used to access/modify data at the same
///      memory location. If ``view_b`` previously held information about a
///      different view, the act of copying has no impact on the values in the
///      old view. (Both behaviors contrast with ``std::vector`` where copy
///      operations
///      always involve making a deepcopy).
///
/// At this time, a GeneralView can not allocate its own data. For example of
/// how to do this, see Enzo-E's CelloView class template.
///
/// Considerations
/// --------------
/// In the long-term, it may make sense to use a ``GeneralView`` template class
/// that wraps Kokkos::View. It is a relatively elegant way to attach
/// information about where memory is allocated.
///
/// We may want to remove all use of the 3D GeneralView.
template <typename T, DataLayout data_layout>
class GeneralView {
  // first, we define useful types used by instances of the class template
  using ptrprops_ = view_detail::MDPtrProps_<T>;

public:
  /// the element type
  using element_type = typename ptrprops_::type;
  using reference_type = element_type&;
  using size_type = int;  // maybe revisit this?
  static constexpr int rank = ptrprops_::rank;
  static constexpr DataLayout layout = data_layout;

private:
  // these entries exist to help us support implicit casts of views of mutable
  // data to the corresponding view of const data
  using non_const_ptr_ =
      view_detail::MkPtr_t<std::remove_const_t<element_type>, rank>;
  using const_ptr_ = view_detail::MkPtr_t<std::add_const_t<element_type>, rank>;
  friend class GeneralView<const_ptr_, data_layout>;

  // attributes:
  element_type* data_;
  size_type extent_[rank];
  size_type strides_[rank];

private:  // helper methods
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
  /// The syntax ensures that the contents of extent_ are all initialized to
  /// zero since extent_ is an array of non-class types
  GeneralView() : data_{nullptr}, extent_{}, strides_{} {}

  ///@{
  /// Construct a view from an existing pointer `ptr`. Every arg after the
  /// pointer specifies the extent of an access (from the fastest axis to the
  /// slowest axis)
  GeneralView(element_type* ptr, int ilen)
      : data_(ptr), extent_{ilen}, strides_{1} {
    // we may need to enforce this check in a different way
    static_assert(rank == 1, "constructor only works with 1D views");
    check_invariants_();
  }

  GeneralView(element_type* ptr, int ilen, int jlen)
      : data_(ptr), extent_{ilen, jlen} {
    static_assert(rank == 2, "constructor only works with 2D views");
    if constexpr (layout == DataLayout::LEFT) {
      strides_[0] = 1;
      strides_[1] = ilen;
    } else {
      strides_[0] = jlen;
      strides_[1] = 1;
    }
    check_invariants_();
  }

  GeneralView(element_type* ptr, int ilen, int jlen, int klen)
      : data_(ptr), extent_{ilen, jlen, klen} {
    static_assert(rank == 3, "constructor only works with 3D views");
    if constexpr (layout == DataLayout::LEFT) {
      strides_[0] = 1;
      strides_[1] = ilen;
      strides_[2] = ilen * jlen;
    } else {
      strides_[0] = jlen * klen;
      strides_[1] = klen;
      strides_[2] = 1;
    }
    check_invariants_();
  }
  ///@}

  /// conversion constructor that facilitates implicit casts from views of
  /// non-constant values to views of constant values
  ///
  /// For example, this allows implicit creation of ``GeneralView<const
  /// double**>`` from ``GeneralView<double**>``
  ///
  /// @note
  /// This is only defined for instances of GeneralView for which T is a pointer
  /// to a a const (e.g. `const int*`, `const double**`). If it were defined
  /// when T is not a pointer-to-const, then it would duplicate the
  /// copy-constructor.
  template <class = std::enable_if<std::is_same<T, const_ptr_>::value>>
  GeneralView(const GeneralView<non_const_ptr_, layout>& other) {
    data_ = other.data_;
    for (int i = 0; i < rank; i++) {
      extent_[i] = other.extent_[i];
      strides_[i] = other.strides_[i];
    }
  }

  // explicitly use defaults for a handful of cases
  ~GeneralView() = default;
  GeneralView(const GeneralView&) = default;             // copy constructor
  GeneralView(GeneralView&&) = default;                  // move constructor
  GeneralView& operator=(const GeneralView&) = default;  // copy assignment
  GeneralView& operator=(GeneralView&&) = default;       // move assignment

  element_type* data() const noexcept { return data_; }
  size_type extent(int i) const {
    GRIMPL_REQUIRE(i >= 0 && i <= rank,
                   "i must be non-negative and can't exceed %d", rank - 1);
    return extent_[i];
  }

  ///@{
  /// Implements multi-dimensional indexing. The first argument corresponds
  /// to the contiguous axis
  GRIMPL_FORCE_INLINE element_type& operator()(int i) const {
    static_assert(rank == 1, "1 index should only be specified for 1D views");
    return data_[i];  // strides_[0] == 1
  }

  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j) const {
    static_assert(rank == 2, "2 indices should only be specified for 2D views");
    if constexpr (layout == DataLayout::LEFT) {
      return data_[i + j * strides_[1]];  // strides_[0] == 1
    } else {
      return data_[i * strides_[0] + j];  // strides_[1] == 1
    }
  }

  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j, int k) const {
    static_assert(rank == 3, "3 indices should only be specified for 3D views");
    if constexpr (layout == DataLayout::LEFT) {
      return data_[i + j * strides_[1] + k * strides_[2]];  // strides_[0] == 1
    } else {
      return data_[i * strides_[0] + j * strides_[1] + k];  // strides_[2] == 1
    }
  }
  ///@}
};

/// @brief a multidimensional view with a data layout equivalent to fortran
///
/// REMINDER: we use 0-indexing
template <typename T>
using FortranView = GeneralView<T, DataLayout::LEFT>;

/// @brief a multidimensional view with C++'s "natural" data layout
template <typename T>
using View = GeneralView<T, DataLayout::RIGHT>;

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_VIEW_HPP