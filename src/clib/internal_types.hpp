// See LICENSE file for license and copyright information

/// @file internal_types.hpp
/// @brief Declares some types used internally by Grackle
///
/// ScratchBuf Data Structures
/// ==========================
/// At the time of writing this, this file contains a number of data structures
/// used while converting Fortran to C/C++. A number of these structures are
/// effectively just structs that are being used to group pointers to memory
/// buffers together that were previously passed through Fortran subroutines
/// as separate arguments.
/// - The older Fortran approach lead to exceptionally long argument lists
/// - The current groupings of buffers in each struct are very experimental (we
///   expect the group. We expect the groups to change over time as we port
///   more and more code.
/// - Each of these data structures will support the visitor design pattern.
///   More about that down below.
///
/// Pointer Semantics
/// -----------------
/// For the following reasons:
///  * our goal is to support performance-portable logic (that should work
///    on CPUs and GPUs)
///  * we have not settled on how we will implement express the logic in
///    performance-portable manner (a kokkos-style system vs something that
///    looks more like C)
///  * and because we are not ready to fully embrace C++ (code mostly resembles
///    a C-like subset with a few C++ features)
/// it therefore makes sense for these data-structures to have "pointer
/// semantics."
///
/// Containers/collections with "pointer semantics" act just like a pointer that
/// represents an array of values. The key point is that when you "copy the
/// container" (or copy a pointer address), you aren't actually copying the
/// values stored within the collection; the underlying values are shared and
/// accessible through both copies of the container (i.e. it is a "shallow
/// copy"). A numpy array and kokkos view are example of containers with
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
/// > these data structures. This means that the care is required to avoid
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
/// to implement without embracing C++ semantics and we probably want to avoid
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
/// as a part of the API. There probably will be some (maybe a lot) value to
/// using these structures to organize internal data, but these structures
/// should not be visible (they can be opaque types or held within an opaque
/// type). If we must provide access to the internal values, we probably want
/// to use some kind of dynamic API.
///
/// In the future, we probably want to define these datatypes in close
/// proximity to where they are used.
///
/// Summary
/// -------
/// There are a lot of things we can improve about these data structures if we
/// fully embrace C++, but there is probably a lot of value to waiting until we
/// have transcribed most code from Fortran and have a coherent plan for
/// refactoring Grackle (we don't want to waste time refactoring and then
/// undoing the refactoring)

#ifndef INTERNAL_TYPES_HPP
#define INTERNAL_TYPES_HPP

#ifndef __cplusplus
#error "This file must be used by a c++ compiler"
#endif

#include "LUT.hpp"

namespace grackle::impl {

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
///   ```{.c}
///     visit_member_ScratchBufType(
///       ScratchBufType* obj, visitor_callback* fn, void* visitor_ctx
///     );
///   ```
///   that apply the visitor to each member
///
/// The branching and use of function pointers will introduce runtime overhead,
/// but that generally shouldn't be an issue.
///
/// If we were willing to embrace C++ and have visit_member_ScratchBufType,
/// accept a template argument, we could avoid totally any runtime overhead.
///
/// Assorted Thoughts
/// -----------------
/// Right now, this is mostly just used to simplify logic for allocation and
/// deallocation.
///
/// In the future, this could be used to:
///   - help with pretty-printing values during debugging
///   - refactoring functions to support GPUs in a piecemeal fashion
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
/// - we it would make more sense for visit_member to accept a template arg
///   rather than a function pointer. Doing so would:
///   - reduce the runtime overhead of applying a visitor
///   - avoid casting the struct members to and from void*. (This is actually
///     quite significant! It's currently very easy to mess this up without
///     the compiler stopping us)
///
/// Future Thoughts
/// ---------------
/// it may make sense to adjust the visit_member function to accept info about
/// grackle's current configuration and specify whether or not a visited member
/// is actually active... Alternatively, that may be too granular (and maybe
/// the structs should be organized so that either all of the members are
/// or none are)

/// this specifies the kind of member being visited
///
/// @note
/// This is a scoped enum, which are more type-safe than regular enums. If we
/// choose to use a regular enum, we should stick the definition inside of its
/// own namespace
enum struct MemberKind { i64_buffer, f64_buffer };

/// an instance of this struct is passed to the visitor
/// -> for a subset of ScratchBuf Data Structures, `name` may be a nullptr
///
/// @note
/// In the future, it may be nice to add another member to this struct in order
/// for the visitor to know whether the member is actively being used under
/// grackle's current configuration (there are a couple of ways to accomplish
/// this)
struct MemberInfo {
  const char* name;
  MemberKind kind;
};

/// the visitor function signature
typedef void visitor_callback(
  MemberInfo member_info, void* member_ptr, void* visitor_ctx
);

/// this is a macro used to help implement the visit_member_... function
#define GRIMPL_IMPL_VISIT_MEMBER(visit_mempair_fn, Tobj, objptr, fn, ctxptr)  \
  {                                                                           \
    /* create a stack allocated instance of Tobj, which holds garbage data */ \
    Tobj dummy;                                                               \
                                                                              \
    /* we create a lambda function to wrap fn and ctx_ptr                     \
     * - it captures both the function-pointer and context pointer by value   \
     * - it recieves references to corresponding members of objptr & dummy    \
     */                                                                       \
    auto wrapper = [ctxptr, fn](MemberInfo member_info,                       \
                                auto& reference_to_obj_member,                \
                                auto& reference_to_dummy_member) {            \
      /* reference_to_obj_member is always a reference to the member.         \
       * I need to make it a pointer to the member type                       \
       *  -> if `reference_to_obj_member.kind == MemberKind.f64_buffer`,      \
       *     then we effectively are converting from double*& to double**     \
       *     before casting to void*                                          \
       */                                                                     \
      void* mem_ptr = (void*)(&reference_to_obj_member);                      \
      fn(member_info, mem_ptr, ctxptr);                                       \
    };                                                                        \
                                                                              \
    visit_mempair_fn(*objptr, dummy, wrapper);                                \
  }

// =================================================================
// Start implementing typical ScratchBuf data structures
// =================================================================

/// Holds 1D arrays used for cooling and heating
///
/// @note
/// If we adopt a strategy with callbacks, I imagine that most (all?) of these
/// pointers will be removed from here.
struct CoolHeatScratchBuf {

  double* ceHI = nullptr;
  double* ceHeI = nullptr;
  double* ceHeII = nullptr;
  double* ciHI = nullptr;
  double* ciHeI = nullptr;
  double* ciHeIS = nullptr;
  double* ciHeII = nullptr;
  double* reHII = nullptr;
  double* reHeII1 = nullptr;
  double* reHeII2 = nullptr;
  double* reHeIII = nullptr;
  double* brem = nullptr;
  double* cieco = nullptr;
  double* hyd01k = nullptr;
  double* h2k01 = nullptr;
  double* vibh = nullptr;
  double* roth = nullptr;
  double* rotl = nullptr;
  double* gpldl = nullptr;
  double* gphdl = nullptr;
  double* hdlte = nullptr;
  double* hdlow = nullptr;

};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  CoolHeatScratchBuf& obj0, CoolHeatScratchBuf& obj1, BinaryFn f
) {
  f(MemberInfo{"ceHI", MemberKind::f64_buffer}, obj0.ceHI, obj1.ceHI);
  f(MemberInfo{"ceHeI", MemberKind::f64_buffer}, obj0.ceHeI, obj1.ceHeI);
  f(MemberInfo{"ceHeII", MemberKind::f64_buffer}, obj0.ceHeII, obj1.ceHeII);
  f(MemberInfo{"ciHI", MemberKind::f64_buffer}, obj0.ciHI, obj1.ciHI);
  f(MemberInfo{"ciHeI", MemberKind::f64_buffer}, obj0.ciHeI, obj1.ciHeI);
  f(MemberInfo{"ciHeIS", MemberKind::f64_buffer}, obj0.ciHeIS, obj1.ciHeIS);
  f(MemberInfo{"ciHeII", MemberKind::f64_buffer}, obj0.ciHeII, obj1.ciHeII);
  f(MemberInfo{"reHII", MemberKind::f64_buffer}, obj0.reHII, obj1.reHII);
  f(MemberInfo{"reHeII1", MemberKind::f64_buffer}, obj0.reHeII1, obj1.reHeII1);
  f(MemberInfo{"reHeII2", MemberKind::f64_buffer}, obj0.reHeII2, obj1.reHeII2);
  f(MemberInfo{"reHeIII", MemberKind::f64_buffer}, obj0.reHeIII, obj1.reHeIII);
  f(MemberInfo{"brem", MemberKind::f64_buffer}, obj0.brem, obj1.brem);
  f(MemberInfo{"cieco", MemberKind::f64_buffer}, obj0.cieco, obj1.cieco);
  f(MemberInfo{"hyd01k", MemberKind::f64_buffer}, obj0.hyd01k, obj1.hyd01k);
  f(MemberInfo{"h2k01", MemberKind::f64_buffer}, obj0.h2k01, obj1.h2k01);
  f(MemberInfo{"vibh", MemberKind::f64_buffer}, obj0.vibh, obj1.vibh);
  f(MemberInfo{"roth", MemberKind::f64_buffer}, obj0.roth, obj1.roth);
  f(MemberInfo{"rotl", MemberKind::f64_buffer}, obj0.rotl, obj1.rotl);
  f(MemberInfo{"gpldl", MemberKind::f64_buffer}, obj0.gpldl, obj1.gpldl);
  f(MemberInfo{"gphdl", MemberKind::f64_buffer}, obj0.gphdl, obj1.gphdl);
  f(MemberInfo{"hdlte", MemberKind::f64_buffer}, obj0.hdlte, obj1.hdlte);
  f(MemberInfo{"hdlow", MemberKind::f64_buffer}, obj0.hdlow, obj1.hdlow);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_CoolHeatScratchBuf(
  CoolHeatScratchBuf* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, CoolHeatScratchBuf, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new CoolHeatScratchBuf
///
/// @param nelem The number of elements in each buffer
CoolHeatScratchBuf new_CoolHeatScratchBuf(int nelem);

/// performs cleanup of the contents of CoolHeatScratchBuf
///
/// This effectively invokes the destructor
void drop_CoolHeatScratchBuf(CoolHeatScratchBuf*);

// -----------------------------------------------------------------

/// holds other assorted scratch buffers used within cool1d_multi_g
///
/// @note
/// Prior to transcription, there was a distinction between the buffers that
/// are now inside of this struct and the buffers in CoolHeatScratchBuf. The
/// distinction has been preserved during transcription, but it is not clear
/// how real the distinction truly is.
struct Cool1DMultiScratchBuf {
  /// unlike the othe members in this struct, tgasold is retained between
  /// iterations. Thus, it may be better to remove it from this struct
  double* tgasold = nullptr;
  double* mynh = nullptr;
  double* myde = nullptr;
  double* gammaha_eff = nullptr;
  double* gasgr_tdust = nullptr;
  double* regr = nullptr;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  Cool1DMultiScratchBuf& obj0, Cool1DMultiScratchBuf& obj1, BinaryFn f
) {
  f(MemberInfo{"tgasold", MemberKind::f64_buffer}, obj0.tgasold, obj1.tgasold);
  f(MemberInfo{"mynh", MemberKind::f64_buffer}, obj0.mynh, obj1.mynh);
  f(MemberInfo{"myde", MemberKind::f64_buffer}, obj0.myde, obj1.myde);
  f(MemberInfo{"gammaha_eff", MemberKind::f64_buffer},
    obj0.gammaha_eff, obj1.gammaha_eff);
  f(MemberInfo{"gasgr_tdust", MemberKind::f64_buffer},
    obj0.gasgr_tdust, obj1.gasgr_tdust);
  f(MemberInfo{"regr", MemberKind::f64_buffer}, obj0.regr, obj1.regr);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_Cool1DMultiScratchBuf(
  Cool1DMultiScratchBuf* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, Cool1DMultiScratchBuf, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new Cool1DMultiScratchBuf
///
/// @param nelem The number of elements in each buffer
Cool1DMultiScratchBuf new_Cool1DMultiScratchBuf(int nelem);

/// performs cleanup of the contents of Cool1DMultiScratchBuf
///
/// This effectively invokes the destructor
void drop_Cool1DMultiScratchBuf(Cool1DMultiScratchBuf*);

// -----------------------------------------------------------------

/// Holds 1D arrays used for linear interpolation
///
/// A common idiom within Grackle is to construct 1D tables of different
/// quantities (I think its mostly different types of rates) sampled at various
/// log temperature values during setup.
/// - when performing calculations on simulation data, Grackle will compute
///   values from these tables.
/// - Since these tables are all sampled at the same log-temperature values, we
///   can reuse information about the current location in the table between
///   different interpolations to speed up the calculation
/// - the values of the buffers in this data structure at a given location are
///   used encode this information about logT table location.
///
/// @note
/// Logic related to this struct is a prime candidate for logic that we probably
/// want to refactor after we complete transcription from Fortran
///
/// @todo
/// Once we finish transcribing, we may want to make the naming a little more
/// generic since it can be used for more than just temperature
struct LogTLinInterpScratchBuf{
  long long* indixe = nullptr;
  double* t1 = nullptr;
  double* t2 = nullptr;
  double* logtem = nullptr;
  double* tdef = nullptr;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  LogTLinInterpScratchBuf& obj0, LogTLinInterpScratchBuf& obj1, BinaryFn f
) {
  f(MemberInfo{"indixe", MemberKind::i64_buffer}, obj0.indixe, obj1.indixe);
  f(MemberInfo{"t1", MemberKind::f64_buffer}, obj0.t1, obj1.t1);
  f(MemberInfo{"t2", MemberKind::f64_buffer}, obj0.t2, obj1.t2);
  f(MemberInfo{"logtem", MemberKind::f64_buffer}, obj0.logtem, obj1.logtem);
  f(MemberInfo{"tdef", MemberKind::f64_buffer}, obj0.tdef, obj1.tdef);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_LogTLinInterpScratchBuf(
  LogTLinInterpScratchBuf* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, LogTLinInterpScratchBuf, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new LogTLinInterpScratchBuf
///
/// @param nelem The number of elements in each buffer
LogTLinInterpScratchBuf new_LogTLinInterpScratchBuf(int nelem);

/// performs cleanup of the contents of LogTLinInterpScratchBuf
///
/// This effectively invokes the destructor
void drop_LogTLinInterpScratchBuf(LogTLinInterpScratchBuf*);

// -----------------------------------------------------------------

/// holds radiative reaction rate buffers
///
/// @note
/// In the future, it would be nice to use this within the general rate
/// storage struct (and maybe also within the photo_rate_storage struct).
/// Since the storage struct only needs to store these rates as scalars we
/// would need to adopt the LUT strategy (or templates)
struct PhotoRxnRateCollection {
  /* Radiative rates for 6-species. */
  double* k24;
  double* k25;
  double* k26;

  /* Radiative rates for 6-species. */
  double* k27;
  double* k28;
  double* k29;
  double* k30;
  double* k31;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  PhotoRxnRateCollection& obj0, PhotoRxnRateCollection& obj1, BinaryFn f
) {
  f(MemberInfo{"k24", MemberKind::f64_buffer}, obj0.k24, obj1.k24);
  f(MemberInfo{"k25", MemberKind::f64_buffer}, obj0.k25, obj1.k25);
  f(MemberInfo{"k26", MemberKind::f64_buffer}, obj0.k26, obj1.k26);
  f(MemberInfo{"k27", MemberKind::f64_buffer}, obj0.k27, obj1.k27);
  f(MemberInfo{"k28", MemberKind::f64_buffer}, obj0.k28, obj1.k28);
  f(MemberInfo{"k29", MemberKind::f64_buffer}, obj0.k29, obj1.k29);
  f(MemberInfo{"k30", MemberKind::f64_buffer}, obj0.k30, obj1.k30);
  f(MemberInfo{"k31", MemberKind::f64_buffer}, obj0.k31, obj1.k31);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_PhotoRxnRateCollection(
  PhotoRxnRateCollection* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, PhotoRxnRateCollection, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new PhotoRxnRateCollection
///
/// @param nelem The number of elements in each buffer
PhotoRxnRateCollection new_PhotoRxnRateCollection(int nelem);

/// performs cleanup of the contents of PhotoRxnRateCollection
///
/// This effectively invokes a destructor
void drop_PhotoRxnRateCollection(PhotoRxnRateCollection*);

// -----------------------------------------------------------------

/// holds reaction rates chemical heating reaction rates
///
/// @note
/// In the future, it would be nice to use this within the general rate
/// storage struct
///
/// @note
/// For now, the modelled rates are fairly limited, but the idea is this could
/// be extended in the future
struct ChemHeatingRates{
  // Chemical heating from H2 formation.
  // numerator and denominator of Eq 23 of Omukai ea. 2000.
  double *n_cr_n;
  double *n_cr_d1;
  double *n_cr_d2;
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  ChemHeatingRates& obj0, ChemHeatingRates& obj1, BinaryFn f
) {
  f(MemberInfo{"n_cr_n", MemberKind::f64_buffer}, obj0.n_cr_n, obj1.n_cr_n);
  f(MemberInfo{"n_cr_d1", MemberKind::f64_buffer}, obj0.n_cr_d1, obj1.n_cr_d1);
  f(MemberInfo{"n_cr_d2", MemberKind::f64_buffer}, obj0.n_cr_d2, obj1.n_cr_d2);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_ChemHeatingRates(
  ChemHeatingRates* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, ChemHeatingRates, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new ChemHeatingRates
///
/// @param nelem The number of elements in each buffer
ChemHeatingRates new_ChemHeatingRates(int nelem);

/// performs cleanup of the contents of ColRecRxnRateCollection
///
/// This effectively invokes a destructor
void drop_ChemHeatingRates(ChemHeatingRates*);

// =================================================================
// The next set of data structures don't really follow all of the
// standard conventions of the other ScratchBuf Data Structures.
// - Instead, they are each implemented as a struct that contains a
//   member called `data`, which is a statically sized array of
//   double* pointers
// - The idea is to use a compile-time lookup table to access members
//   of the struct (the compile-time lookup table is implemented in
//   terms of enumerators)
// =================================================================

/// holds properties about each kind of species
///
/// The basic premise is that we can use `SpLUT::<entry>` to lookup values
/// for the desired species
///
/// @note
/// The way this is currently a lot like a struct of arrays (i.e. there is an
/// extra level of indirection).
/// - For concreteness, to access index `idx` of the data for HeI, you
///   would write `obj.data[SpLUT::HeI][idx]`.
/// - in terms of performance, this is very similar to what came before (and
///   is good enough for now)
/// - in the long term, it would be even better to replace this with a
///   `View<double**>`, but I want to wait until after we have transcribed the
///   bulk of grackle before we do that.
///
/// @note
/// At the time of writing, both this data structure and the
/// GrainSpeciesCollection data structure both reserve space for each of the
/// dust species. Maybe this data structure shouldn't do this? (If we remove
/// the grain species, maybe we rename this so that it is called
/// ChemSpeciesCollection?)
struct SpeciesCollection {
  double* data[SpLUT::NUM_ENTRIES];
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  SpeciesCollection& obj0, SpeciesCollection& obj1, BinaryFn f
) {
  // it's okay that we aren't specifing the name of the visited member.
  // -> if we need this info, we could easily write a function mapping the
  //    LUT index to the species name
  // -> frankly, this function is totally unnecessary. We're only implementing
  //    it for the sake of consistency/convenience
  for (int i = 0; i < SpLUT::NUM_ENTRIES; i++) {
    f(MemberInfo{nullptr, MemberKind::f64_buffer}, obj0.data[i], obj1.data[i]);
  }
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_SpeciesCollection(
  SpeciesCollection* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, SpeciesCollection, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new SpeciesCollection
///
/// @param nelem The number of elements in each buffer
SpeciesCollection new_SpeciesCollection(int nelem);

/// performs cleanup of the contents of SpeciesCollection
///
/// This effectively invokes a destructor
void drop_SpeciesCollection(SpeciesCollection*);

// -----------------------------------------------------------------

/// holds properties about each kind of dust grain
///
/// @note
/// This operates in a similar manner to SpeciesCollection (i.e. we use
/// `OnlyGrainSpLUT::<entry>` to lookup values for the desired rate).
///
/// @note
/// This is something we may want to reuse. If we are willing to embrace C++,
/// then we may want to use templates
struct GrainSpeciesCollection {
  double* data[OnlyGrainSpLUT::NUM_ENTRIES];
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  GrainSpeciesCollection& obj0, GrainSpeciesCollection& obj1, BinaryFn f
) {
  for (int i = 0; i < OnlyGrainSpLUT::NUM_ENTRIES; i++) {
    f(MemberInfo{nullptr, MemberKind::f64_buffer}, obj0.data[i], obj1.data[i]);
  }
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_GrainSpeciesCollection(
  GrainSpeciesCollection* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, GrainSpeciesCollection, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new GrainSpeciesCollection
///
/// @param nelem The number of elements in each buffer
GrainSpeciesCollection new_GrainSpeciesCollection(int nelem);

/// performs cleanup of the contents of GrainSpeciesCollection
///
/// This effectively invokes a destructor
void drop_GrainSpeciesCollection(GrainSpeciesCollection*);

// -----------------------------------------------------------------

/// holds properties about Collisional and Recombination Reaction Rates that
/// behave in the "standard" way (i.e. we interpolate in 1D with respect to
/// log T)
///
/// This operates in a similar manner to SpeciesCollection (i.e. we use
/// `ColRecRxnLUT::<entry>` to lookup values for the desired rate).
struct ColRecRxnRateCollection {
  double* data[ColRecRxnLUT::NUM_ENTRIES];
};

/// used to help implement the visitor design pattern
///
/// (avoid using this unless you really have to)
template<class BinaryFn>
void visit_member_pair(
  ColRecRxnRateCollection& obj0, ColRecRxnRateCollection& obj1, BinaryFn f
) {
  for (int i = 0; i < ColRecRxnLUT::NUM_ENTRIES; i++) {
    f(MemberInfo{nullptr, MemberKind::f64_buffer}, obj0.data[i], obj1.data[i]);
  }
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
/// @param visitor_ctx[in,out] User-defined callback function context
inline void visit_member_ColRecRxnRateCollection(
  ColRecRxnRateCollection* ptr, visitor_callback* fn, void* visitor_ctx
) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, ColRecRxnRateCollection, ptr, fn,
                           visitor_ctx)
}

/// allocates the contents of a new ColRecRxnRateCollection
///
/// @param nelem The number of elements in each buffer
ColRecRxnRateCollection new_ColRecRxnRateCollection(int nelem);

/// performs cleanup of the contents of ColRecRxnRateCollection
///
/// This effectively invokes a destructor
void drop_ColRecRxnRateCollection(ColRecRxnRateCollection*);

} // namespace grackle::impl

#undef GRIMPL_IMPL_VISIT_MEMBER

#endif /* INTERNAL_TYPES_HPP */
