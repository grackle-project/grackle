// See LICENSE file for license and copyright information

/// @file common.hpp
/// @brief Common utilities used with ScratchBuf Data Structures
///
/// Just below, we first include a discussion of "ScratchBuf Data Structures,"
/// and subsequently discuss the visitor machinery.
///
/// ScratchBuf Data Structures
/// ==========================
/// At the time of writing this, Grackle defines a number of data structures
/// with special characteristics to aid the process of converting Fortran to
/// C/C++. These data structures are effectively just structs that are used
/// to group pointers to memory buffers together that were previously passed
/// through Fortran subroutines  as separate arguments.
/// - The older Fortran approach lead to exceptionally long argument lists
/// - The current groupings of buffers in each struct are very experimental (we
///   expect the group. We expect the groups to change over time as we port
///   more and more code.
/// - Each of these data structures will support the visitor design pattern.
///   More about that down below.
///
/// Pointer Semantics
/// -----------------
/// Because
///  * our goal is to support performance-portable logic (that should work
///    on CPUs and GPUs)
///  * we have not settled on how we will implement express the logic in
///    performance-portable manner (a Kokkos-style system vs something that
///    looks more like C)
///  * and we're not ready to fully embrace C++ (code mostly resembles
///    a C-like subset with a few C++ features)
/// it therefore makes sense for these data-structures to have "pointer
/// semantics."
///
/// Containers/collections with "pointer semantics" act just like a pointer that
/// represents an array of values. The key point is that when you "copy the
/// container" (or copy a pointer address), you aren't actually copying the
/// values stored within the collection; the underlying values are shared and
/// accessible through both copies of the container (i.e. it is a "shallow
/// copy"). A numpy array and Kokkos view are example of containers with
/// pointer semantics.
///
/// For context, containers/collections without "pointer semantics" generally
/// have "value semantics" instead. When copying a container with "value
/// semantics," the newly created copy holds stored within the container (a
/// deepcopy). Examples of a container/collection "value semantics" include
/// C++'s ``std::vector`` or something like the following struct
/// ```{.c}
/// struct MyContainer{ double data[4] };
/// ```
///
/// > [!note]
/// > Currently, the idea is to require manual initialization and cleanup of
/// > these data structures. This means that care is required to avoid
/// > unitialized memory, memory leaks and dangling pointers
///
/// Avoiding Constructors and Destructors
/// -------------------------------------
/// Since we are currently writing C++ code, there is a temptation to leverage
/// Constructors and Destructors (that we could easily convert back to C code).
/// However, this is not a good use of time right now (and will probably waste
/// more time in the future as our implementation strategy changes) unless we
/// are willing to fully embrace C++.
///
/// In more detail, a container with "pointer semantics" that implements
/// (con|de)structors generally must implement semantics like C++'s
/// ``std::shared_ptr`` or ``std::unique_ptr``. The former is a little tedious
/// to implement without embracing C++ semantics, and we probably want to avoid
/// the latter (at least for now).
///
/// > [note]
/// > It is definitely possible to make a container work with semantics similar
/// > to a ``std::unique_ptr``, but you would need to make some extensions to
/// > the logic to get it working on GPUs (especially since passing an object
/// > to a GPU involves an operation like memcpy). Furthermore, to support
/// > code that runs on many platforms we would probably want to pass this
/// > around by reference. We might ultimately want to pursue this route, but
/// > we should be very deliberate and consistent about doing this and I think
/// > we should defer this choice until we have a coherent strategy)
///
/// Other Thoughts
/// --------------
/// It is probably a bad idea to publicly expose any of these data structures
/// as a part of the API. There probably will be some (maybe a lot of) value to
/// using these structures to organize internal data, but these structures
/// should not be visible (they can be opaque types or held within an opaque
/// type). If we must provide access to the internal values, we probably want
/// to use some kind of dynamic API.
///
/// Summary
/// -------
/// There are a lot of things we can improve about these data structures if we
/// fully embrace C++, but there is probably a lot of value to waiting until we
/// have transcribed most code from Fortran and have a coherent plan for
/// refactoring Grackle (we don't want to waste time refactoring and then
/// undoing the refactoring)
///
/// Overview of the Visitor Design Pattern
/// ======================================
/// Because all of these internal data structures are effectively structs of
/// pointers there is value to supporting a variation on the
/// [visitor-design pattern](https://en.wikipedia.org/wiki/Visitor_pattern).
/// Historically, this was described in terms of object-oriented programming,
/// but it is more general than that.
///
/// For the uninitiated, there are 4 parts to the visitor pattern:
///   1. an object-structure composed of 1 or more elements of various kinds
///      (the set of element-kinds must be well-described)
///   2. a visitor entity that provides a set of logic to perform a particular
///      operation for each element-kind
///   3. when applying a visitor to an element, a mechanism is required to
///      invoke the logic based on the element-kind (this can be implicit,
///      e.g. with function overloads)
///   4. logic to apply the visitor to all elements of the object-structure
///
/// How we'll use it
/// ----------------
/// In our case:
/// - the different object-structures correspond to different ScratchBuf Data
///   Structures and the elements are the data-members.
/// - given our general reticence to fully embrace C++:
///   - each visitor entity is specified by a function pointer (with optional
///     context-data)
///   - the visitor is internally responsible for dispatching the appropriate
///     logic based on provided information about the member.
/// - each ScratchBuf Data Structure will have an associated function,
///   ```{.cpp}
///     template<typename BinaryVisitor>
///     visit_member(ScratchBufType* obj, UnaryVisitor fn);
///   ```
///   that apply the visitor to each member
///
/// An earlier version of this machinery involved function pointers. However,
/// I decided to fully embrace templates for this particular aspect of the code
/// since you can't portably use function pointers on GPUs. (And there are some
/// special cases, at least for now, where need to use the visitor machinery in
/// compute-logic). We can revisit this choice at a later date.
///
/// If we are willing to more fully embrace C++, we can have visit_member
/// accept a reference.
///
/// Assorted Thoughts
/// -----------------
/// Right now, this machinery is mostly used to simplify logic for allocation
/// and deallocation. The machinery is also used to create Printer, which will
/// print the contents of a struct (it prints the values held by pointers
/// rather than pointer addresses) -- this is intended for debugging. There are
/// 1 or 2 special cases where we use the visitor machinery for other stuff.
///
/// In the future, this could be used to:
///   - refactoring functions to support GPUs in a piecemeal fashion.
///     - we would make new allocator/deallocator visitors that would
///       allocate/deallocate memory on the GPU
///     - we would implement a visitor to perform GPU-memcpy (e.g. hipMemcpy or
///       cudaMemcpy). We currently have a proof-of-concept CopyMember visitor
///       that invokes the standard memcpy function (it won't be hard to use it
///       as a template for implementing this other idea)
///   - aggregating all allocations and deallocations (way down the road, this
///     could plausibly be important for attaining performance with GPUs -- but
///     that's a way off and the code structure could change a bunch by then)
///
/// For reasons related to transcribing step_rate_newton_raphson, we are going
/// to implement the visit_member_<...> function in terms of a
/// template-function that visits pairs of members. We should be able to get
/// rid of it in the future (essentially, we need to maintain some consistent
/// behavior during transcription, but I think we should change that behavior
/// over the long-term)
///
/// If we are ever willing to more fully embrace C++:
/// - attaching the vist_member functions to the corresponding structs would
///   make a lot of sense.
///
/// Future Thoughts
/// ---------------
/// it may make sense to adjust the visit_member function to accept info about
/// grackle's current configuration and specify whether or not a visited member
/// is actually active... Alternatively, that may be too granular (and maybe
/// the structs should be organized so that either all a struct's members are
/// active or none are)

#ifndef VISITOR_COMMON_HPP
#define VISITOR_COMMON_HPP

#include <cstddef>  // std::size_t

/// VIS_MEMBER_NAME is intended to wrap the name of a member. Essentially, the
/// idea is to make it possible to provide visitors (like the printer) with more
/// more information in some-compilation-configurations (e.g. DEBUG) and to
/// avoid storing string-literals in the compiled library in other cases
#if 1
#define VIS_MEMBER_NAME(name) name
#else
#define VIS_MEMBER_NAME(name) nullptr
#endif

namespace grackle::impl::visitor {

/// This is basic context that almost all visitors should have
struct VisitorCtx {
  unsigned int idx_range_len;
};

/// This specifies the buffer-length of a given ScratchBuf member. Instances of
/// this struct are passed as an argument to the visitor.
///
/// An instance ``spec`` specifies that the number of elements in a ScratchBuf
/// member is ``spec.idxrange_len_multiple * idxrange_len + spec.constant``,
/// where ``idxrange_len`` denotes the length of the calculation's IdxRange.
///
/// @note
/// In almost all cases, ``idxrange_len_multiple`` is 1 and ``constant`` is 0.
struct BufLenSpec {
  unsigned int idxrange_len_multiple;
  unsigned int constant;
};

/// Returns a BufLenSpec expressing that the a ScratchBuf member is just
inline BufLenSpec idx_range_len_multiple(unsigned int idxrange_len_multiple) {
  return BufLenSpec{idxrange_len_multiple, 0};
}

/// Utility function used by visitors. It computes the length of a buffer.
inline std::size_t get_buf_len(const BufLenSpec& spec, const VisitorCtx& ctx) {
  return ctx.idx_range_len * spec.idxrange_len_multiple + spec.constant;
}

/// informs the visitor that we are beginning the visit with the given type
///
/// This can be specialized for different Visitors, but the default
/// implementation does nothing.
template <typename Visitor>
inline void begin_visit(const char* type_name, Visitor& visitor) { /* NO-OP */ }

/// informs the visitor that we are ending the visit with a given type.
template <typename Visitor>
inline void end_visit(Visitor& visitor) { /* NO-OP */ }

/// wraps a visitor that expects a single argument in a way that it can be
/// used in contexts where we expect a single member
///
/// This is only intended to be used within the GRIMPL_IMPL_VISIT_MEMBER macro
/// (it shouldn't be used directly).
template <typename T>
struct UnaryAdaptor {
  T& adaptee;

  explicit UnaryAdaptor(T& adaptee) : adaptee(adaptee) {}

  void operator()(const char* name, long long*& val, long long*& dummy,
                  const BufLenSpec& spec) {
    adaptee(name, val, spec);
  }

  void operator()(const char* name, double*& val, double*& dummy,
                  const BufLenSpec& spec) {
    adaptee(name, val, spec);
  }
};

template <typename T>
inline void begin_visit(const char* type_name, UnaryAdaptor<T>& visitor) {
  begin_visit(type_name, visitor.adaptee);
}

template <typename T>
inline void end_visit(UnaryAdaptor<T>& visitor) {
  end_visit(visitor.adaptee);
}

/// this is a macro used to help implement the visit_member function
///
/// The basic idea is that we fully implement an overload for visit_member_pair
/// for a given ScratchBuf struct. Then, the internal implementation of the
/// overload for the visit_member ScratchBuf struct can be replaced with this
/// macro
#define GRIMPL_IMPL_VISIT_MEMBER(visit_mempair_fn, Tobj, objptr, visitor_obj)  \
  {                                                                            \
    /* create a stack allocated instance of Tobj, which holds garbage data */  \
    Tobj dummy;                                                                \
    visit_mempair_fn(*objptr, dummy,                                           \
                     grackle::impl::visitor::UnaryAdaptor(visitor_obj));       \
  }

}  // namespace grackle::impl::visitor

#endif  // VISITOR_COMMON_HPP