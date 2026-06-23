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

  /// returns whether `*this` has associated data (it doesn't wrap a nullptr)
  explicit operator bool() const noexcept { return data_ != nullptr; }

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

  /// Return the pointer to the start of the ``j``th contiguous 1d span.
  ///
  /// This method exists purely to provide feature parity with @ref Multi1DView.
  /// This is useful if you want to change the underlying type of the data
  /// that is being passed to a function. For that reason, this method only
  /// supports view with a rank of 2.
  ///
  /// Performance
  /// -----------
  /// Using "normal" multidimensional indexing (i.e. not this method) to access
  /// elements generally produces instructions with comparable or better
  /// performance than using this method. This directly constrasts the advise
  /// given for the methods of @ref Multi1DView.
  GRIMPL_FORCE_INLINE element_type* contig1d_ptr(int j) const {
    static_assert(rank == 2, "This method is only provided for 2D views");
    return data_ + (j * strides_[1]);
  }
};

/// Represents a sequence of pointers.
///
/// At this time, instances never own the underlying memory (those allocations
/// are managed outside of this type). The primary feature of this type is that
/// it implements the interface of a `View<T**>`. If you have an instance
/// called `mview`, then you can invoke `mview(i,j)` to access an element.
///
/// - (Just like for `View<T**>`) in the expression `mview(i,j)`:
///   - `i` specifies the index along the contiguous axis
///   - `j` specifies the index along the non-contiguous axis
/// - in practice, you may want to avoid using `mview(i,j)` within a tight
///   loop for the aliasing-related reasons discussed down below. With that
///   said, this interface is quite useuf in certain contexts
/// - the initial impetus for providing this interface will be for setting up
///   the vector used in the Newton-Rapshon solver
///
/// Data Layout Considerations
/// ==========================
/// While it may be obvious, it's worth emphasizing that under this type
/// effectively remaps the "natural data layout" of an array of pointers in
/// C/C++ in order to be consistent with the data layout of `View<T**>`.
///
/// Let's elaborate and consider a variable ``a`` that represents a 2D matrix
/// with shape ``(m,n)``. It's useful to recall that the mathematical notation
/// for specifying a matrix element is ``aᵢⱼ`` (where ``i`` increases as you go
/// down a column and ``j`` increases as you move along a row).
///
/// Suppose ``a`` is represented by a type ``T**``:
/// - To access an element, you would write an ``a[i][j]``.
/// - Because this expression is equivalent to ``(a[i])[j]``, it's easy to see
///   that incrementing the ``j`` index points to an adjacent memory location.
/// - This has row-major ordering or right-layout (i.e. data is contiguous
///   along the rightmost extent).
///
/// Now, instead suppose ``a`` is represented by ``Multi1DView<T>``:
/// - To access an element, you would write ``a(i,j)``.
/// - The type is designed such that incrementing the ``i`` index will point
///   to an adjacent memory location.
/// - This has column-major ordering or left-layout (i.e. data is contiguous
///   along the leftmost extent)
///
/// All of this is to say that you should think of individual 1d views in this
/// type as if they are filling different columns of data.
///
/// @note
/// If we were willing to make View<T**> use the natural right-layout,
///
/// Array of Pointers vs Pointer of Pointers
/// ========================================
/// An important implementation question is whether we should declare this
/// type's `data_` member as a fixed size c-style array of pointers or a
/// pointer to pointers
///
/// 1. a fixed size c-style array of pointers.
///    - the declaration of the `data_` data-member would look like
///      `T* data_[BIG_NUM]`, where `BIG_NUM` is a compile-time constant (we
///      *could* instead make it a template-parameter of the type)
///    - currently, ``BIG_NUM`` needs to be ~50 to be able to track pointers for
///      all known species densities
///    - in this case, the `data_` data-member would take up approximately
///      ``BIG_NUM * sizeof(void*)`` bytes (the amount of memory associated
///      with each stored pointer is irrelevant to this discussion)
/// 2. OR a pointer of pointers
///    - the declaration of the `data_` data-member would look like `T** data_`
///    - in this case, the `data_` data-member would take up approximately
///      ``sizeof(void*)`` bytes (the amount of memory associated
///      with each stored pointer is irrelevant to this discussion)
///
/// The main tradeoff:
/// - For a stack-allocated ``Multi1DView<T>`` instance called ``mview``,
///   accessing a value with ``mview(i, j)`` under the first approach involves
///   1 fewer pointer dereferences
/// - the size of the ``data_`` member is significantly larger for the first
///   approach. Thus, a stack-allocated ``Multi1DView`` instance requires much
///   more stack-space. This *might* be important on GPUs (per-thread stack
///   usage can limit the number of usable GPU threads). This may encourage
///   passing views around by reference (which effectively adds another pointer
///   dereference).
///
/// The big question pertains to a for-loop when we we call `mview(i,j)`. We
/// have taken steps to ensure that the underlying `data_[j][idx]` access is
/// always inlined. Generally, a C/C++ compiler would insert extra instructions
/// to check whether `data_[j]` is changes between accesses. I suspect, that
/// first approach may make it easier to conclude at compile time that the
/// value won't change (but I don't actually know -- it almost certainly
/// depends on whether there are other function calls inside the loop and what
/// assumptions the compiler can make about them).
///
/// We defer this point to the future since
/// 1. this construct will be at least as performant as the code it replaces
/// 2. If we really want to treat the data as 2D view in a tight loop, we
///    should really be copying the data into a contiguous array so we can
///    track it directly in a `View<T**>`
template <typename T>
class Multi1DView {
public:
  // define type-aliases consistent with View<T**>
  using element_type = T;
  using reference_type = element_type&;
  using size_type = int;
  static constexpr int rank = 2;

private:
  // the choice to declare `data_` as `element_type*const*`, rather than
  // `element_type**` is intended to make sure that we never accidently mutate
  // a pointer address
  element_type* const* data_;
  /// number of pointers tracked by data_
  size_type n_ptr_;
  /// a constant offset applied to each pointer in data_
  size_type ptr_offset_;
  /// number of elements accessible through each pointer
  ///
  /// For a pointer in `data_` called `ptr`, a `Multi1DView` instance provides
  /// access to values at `ptr[ptr_offset_ + i]` for `0 <= i < elem_in_ptr_`.
  size_type elem_in_ptr_;

public:
  /// Default constructor
  ///
  /// @note
  /// The syntax ensures that the contents of extent_ are all initialized to
  /// zero since extent_ is an array of non-class types
  Multi1DView() : data_{nullptr}, n_ptr_{0}, ptr_offset_{0}, elem_in_ptr_{0} {}

  /// Construct an instance a from an array of pointers, @p arr_of_ptr. Every
  /// subsequent arg specifies the extent of an axis (from the fastest axis to
  /// the slowest axis)
  ///
  /// The provided @p arr_of_ptr should not be a nullptr or hold a nullptr.
  /// Furthermore, undefined behavior may arise if the constructed object
  /// outlives the memory associated associated with the argument (it is of
  /// course ok for the arg's memory to be freed before the object's destructor
  /// is called)
  ///
  /// @note
  /// Be aware that the number of elements in @p arr_of_ptr is specified by the
  /// the final argument passed to the constructor. This is counter-intuitive
  /// given the way you would directly access memory from @p arr_of_ptr (but
  /// it is consistent with argument order for the constructor of `View<T**>`)
  Multi1DView(element_type* const* arr_of_ptr, int ilen, int jlen)
      : data_(arr_of_ptr), n_ptr_{jlen}, ptr_offset_{0}, elem_in_ptr_{ilen} {
    GRIMPL_REQUIRE(ilen > 0, "ilen must be positive, received: %d", ilen);
    GRIMPL_REQUIRE(jlen > 0, "jlen must be positive, received: %d", jlen);
    GRIMPL_REQUIRE(arr_of_ptr != nullptr, "arr_of_ptr is a nullptr");

    // we are not ready to enforce the following invariant yet (it will be a
    // useful safety check once we are ready)

    // for (int j = 0; j < jlen; j++) {
    //   GRIMPL_REQUIRE(arr_of_ptr[j] != nullptr, "arr_of_ptr holds nullptr");
    // }
  }

  // This is a quick and dirty solution. Ideally, we would set this in the
  // constructor (I'm just not clear on the best way to differentiate the
  // ptr_offset from the extents)
  void override_ptr_offset_and_ilen(int ptr_offset, int ilen) {
    GRIMPL_REQUIRE(ptr_offset >= 0, "ptr_offset can't be negative");
    GRIMPL_REQUIRE(ilen > 0, "ilen must be positive");
    ptr_offset_ = ptr_offset;
    elem_in_ptr_ = ilen;
  }

  // explicitly use defaults for a handful of cases
  ~Multi1DView() = default;
  Multi1DView(const Multi1DView&) = default;             // copy constructor
  Multi1DView(Multi1DView&&) = default;                  // move constructor
  Multi1DView& operator=(const Multi1DView&) = default;  // copy assignment
  Multi1DView& operator=(Multi1DView&&) = default;       // move assignment

  /// returns whether `*this` has associated data (it doesn't wrap a nullptr)
  explicit operator bool() const noexcept { return data_ != nullptr; }

  size_type extent(int i) const {
    switch (i) {
      case 0:
        return elem_in_ptr_;  // <- contiguous axis
      case 1:
        return n_ptr_;
      default:
        GRIMPL_ERROR("i must be 0 or 1");
    }
  }

  /// Access the element corresponding to the provided multidimensional index
  ///
  /// @note
  /// For consistency with the element-access method provided by View<T**>, the
  /// the first argument corresponds to the contiguous axis.
  GRIMPL_FORCE_INLINE element_type& operator()(int i, int j) const {
    // reminder: the order of indices is reversed from the order of that they
    //           are passed to this member-function in order to follow the
    //           indexing convention used in View (which inheritted the
    //           convention from Fortran during transcription of code to C++)
    return data_[j][i + ptr_offset_];
  }

  /// Return the pointer to the start of the ``j``th contiguous 1d span
  ///
  /// @note
  /// Ideally, this method wouldn't exist since it's a departure from the
  /// basic interface of a 2d array (i.e. its a leaky abstraction). But, for
  /// performance purposes, it's a necessary evil
  ///
  /// Purpose
  /// -------
  /// This method serves as an alternative mechanism for accessing values when
  /// iterating over a multiview, as opposed to direct multidimensional
  /// indexing.
  ///
  /// Consider the following code snippets how you might loop over 2 contiguous
  /// slices at `jA` and `jB` from an instance called `mview`. For the sake of
  /// this discussion, we use `OP` to denote the operations using these values
  /// (it could represent a function call or multiple operations directly
  /// inscribed into the loop that may or may not involve function calls).
  /// Suppose we store the results in an output buffer, `buf`.
  ///
  /// Using normal multidimensional-array access:
  /// @code{cpp}
  /// for (int i = i_start; i < i_stop; i++) {
  ///   buf[i] = OP(mview(i, jA), mview(i, jB) /* , <other args...> */);
  /// }
  /// @endcode
  ///
  /// Using this method:
  /// @code{cpp}
  /// double* ptr_A = mview.contig1d_ptr(jA);
  /// double* ptr_B = mview.contig1d_ptr(jB);
  /// for (int i = i_start; i < i_stop; i++) {
  ///   buf[i] = OP(ptr_A[i], ptr_B[j] /* , <other args...> */);
  /// }
  /// @endcode
  ///
  /// The second approach will generally compile to instructions with
  /// perfomance that is comparable or superior to the former.
  /// - This is fundamentally a consequence of the fact that each occurence of
  ///   `mview(i,j)` corresponds to 2 levels of pointer-dereferences. In the
  ///   second approach, we are manually caching one-level of dereferences
  /// - While a compiler *may* be able to automatically cache one level of the
  ///   pointer dereferences, but that depends on the precise details of the
  ///   operations performed in the loop. For example, the presence of a
  ///   non-inlined function-call may prevent such an optimization.
  GRIMPL_FORCE_INLINE element_type* contig1d_ptr(int j) const {
    // I don't love the name of this method (I'm open to alternatives).
    // Relevant requirements:
    // -> we want to be able to add it to the View<T**> method.
    // -> while we could call it something like get_column, I think referring to
    //    rows/columns may be confusing since this type remaps mapping compared
    //    to the natural data mapping
    return data_[j] + ptr_offset_;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_VIEW_HPP