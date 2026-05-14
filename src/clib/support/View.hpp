//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the View construct
///
//===----------------------------------------------------------------------===//

#ifndef SUPPORT_VIEW_HPP
#define SUPPORT_VIEW_HPP

#include "support/config.hpp"
#include "support/status_reporting.hpp"

#include <type_traits>  // std::remove_pointer_t, std::is_pointer_v

namespace GRIMPL_NAMESPACE_DECL {

namespace view_detail {
/// this simply helps us implement the View class template
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
template <typename T>
struct View {
  // first, we define useful types used by instances of the class template
private:
  using ptrprops_ = view_detail::MDPtrProps_<T>;

public:
  /// the element type
  using element_type = typename ptrprops_::type;
  using reference_type = element_type&;
  using size_type = int;  // maybe revisit this?
  static constexpr int rank = ptrprops_::rank;

private:  // attributes
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
  View() : data_{nullptr}, extent_{}, strides_{} {}

  ///@{
  /// Construct a view from an existing pointer `ptr`. Every arg after the
  /// pointer specifies the extent of an access (from the fastest axis to the
  /// slowest axis)
  View(element_type* ptr, int ilen) : data_(ptr), extent_{ilen}, strides_{1} {
    // we may need to enforce this check in a different way
    static_assert(rank == 1, "constructor only works with 1D views");
    check_invariants_();
  }

  View(element_type* ptr, int ilen, int jlen)
      : data_(ptr), extent_{ilen, jlen}, strides_{1, ilen} {
    static_assert(rank == 2, "constructor only works with 2D views");
    check_invariants_();
  }

  View(element_type* ptr, int ilen, int jlen, int klen)
      : data_(ptr), extent_{ilen, jlen, klen}, strides_{1, ilen, ilen * jlen} {
    static_assert(rank == 3, "constructor only works with 3D views");
    check_invariants_();
  }
  ///@}

  // explicitly use defaults for a handful of cases
  ~View() = default;
  View(const View&) = default;             // copy constructor
  View(View&&) = default;                  // move constructor
  View& operator=(const View&) = default;  // copy assignment
  View& operator=(View&&) = default;       // move assignment

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
    return data_[i + j * strides_[1]];  // strides_[0] == 1
  }

  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j, int k) const {
    static_assert(rank == 3, "3 indices should only be specified for 3D views");
    return data_[i + j * strides_[1] + k * strides_[2]];  // strides_[0] == 1
  }
  ///@}
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_VIEW_HPP