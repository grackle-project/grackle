// See LICENSE file for license and copyright information

/// @file internal_types.hpp
/// @brief Declares some types used internally by Grackle
///
/// This files defines a bunch of datatypes that all behave in a very similar
/// manner. See the discussion at the top of visitor/common.hpp for an extended
/// description.
///
/// Going forward, we plan to distribute the definitions of a lot of these data
/// types among other header/implementation files so that the struct is defined
/// in proximity to where it is used.

#ifndef INTERNAL_TYPES_HPP
#define INTERNAL_TYPES_HPP

#ifndef __cplusplus
#error "This file must be used by a c++ compiler"
#endif

#include "LUT.hpp"

namespace grackle::impl {

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
  /// unlike the other members in this struct, tgasold is retained between
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
