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
#include "visitor/common.hpp"
#include "utils-cpp.hpp"

namespace grackle::impl {

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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("CoolHeatScratchBuf", f);
  f(VIS_MEMBER_NAME("ceHI"), obj0.ceHI, obj1.ceHI, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ceHeI"), obj0.ceHeI, obj1.ceHeI, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ceHeII"), obj0.ceHeII, obj1.ceHeII, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ciHI"), obj0.ciHI, obj1.ciHI, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ciHeI"), obj0.ciHeI, obj1.ciHeI, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ciHeIS"), obj0.ciHeIS, obj1.ciHeIS, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("ciHeII"), obj0.ciHeII, obj1.ciHeII, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("reHII"), obj0.reHII, obj1.reHII, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("reHeII1"), obj0.reHeII1, obj1.reHeII1, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("reHeII2"), obj0.reHeII2, obj1.reHeII2, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("reHeIII"), obj0.reHeIII, obj1.reHeIII, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("brem"), obj0.brem, obj1.brem, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("cieco"), obj0.cieco, obj1.cieco, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("hyd01k"), obj0.hyd01k, obj1.hyd01k, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("h2k01"), obj0.h2k01, obj1.h2k01, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("vibh"), obj0.vibh, obj1.vibh, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("roth"), obj0.roth, obj1.roth, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("rotl"), obj0.rotl, obj1.rotl, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("gpldl"), obj0.gpldl, obj1.gpldl, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("gphdl"), obj0.gphdl, obj1.gphdl, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("hdlte"), obj0.hdlte, obj1.hdlte, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("hdlow"), obj0.hdlow, obj1.hdlow, vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(CoolHeatScratchBuf* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, CoolHeatScratchBuf, ptr, fn);
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("Cool1DMultiScratchBuf", f);
  f(VIS_MEMBER_NAME("tgasold"), obj0.tgasold, obj1.tgasold, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("mynh"), obj0.mynh, obj1.mynh, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("myde"), obj0.myde, obj1.myde, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("gammaha_eff"), obj0.gammaha_eff, obj1.gammaha_eff, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("gasgr_tdust"), obj0.gasgr_tdust, obj1.gasgr_tdust, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("regr"), obj0.regr, obj1.regr, vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(Cool1DMultiScratchBuf* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, Cool1DMultiScratchBuf, ptr, fn);
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("LogTLinInterpScratchBuf", f);
  f(VIS_MEMBER_NAME("indixe"), obj0.indixe, obj1.indixe, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("t1"), obj0.t1, obj1.t1, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("t2"), obj0.t2, obj1.t2, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("logtem"), obj0.logtem, obj1.logtem, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("tdef"), obj0.tdef, obj1.tdef, vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(LogTLinInterpScratchBuf* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, LogTLinInterpScratchBuf, ptr, fn);
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("PhotoRxnRateCollection", f);
  f(VIS_MEMBER_NAME("k24"), obj0.k24, obj1.k24, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k25"), obj0.k25, obj1.k25, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k26"), obj0.k26, obj1.k26, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k27"), obj0.k27, obj1.k27, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k28"), obj0.k28, obj1.k28, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k29"), obj0.k29, obj1.k29, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k30"), obj0.k30, obj1.k30, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("k31"), obj0.k31, obj1.k31, vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(PhotoRxnRateCollection* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, PhotoRxnRateCollection, ptr, fn)
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("ChemHeatingRates", f);
  f(VIS_MEMBER_NAME("n_cr_n"), obj0.n_cr_n, obj1.n_cr_n, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("n_cr_d1"), obj0.n_cr_d1, obj1.n_cr_d1, vis::idx_range_len_multiple(1));
  f(VIS_MEMBER_NAME("n_cr_d2"), obj0.n_cr_d2, obj1.n_cr_d2, vis::idx_range_len_multiple(1));
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(ChemHeatingRates* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, ChemHeatingRates, ptr, fn)
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
  // it's okay that we aren't specifying the name of the visited member.
  // -> if we need this info, we could easily write a function mapping the
  //    LUT index to the species name
  // -> frankly, this function is totally unnecessary. We're only implementing
  //    it for the sake of consistency/convenience
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("SpeciesCollection", f);
  for (int i = 0; i < SpLUT::NUM_ENTRIES; i++) {
    f(VIS_MEMBER_NAME("data[...]"), obj0.data[i], obj1.data[i], vis::idx_range_len_multiple(1));
  }
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(SpeciesCollection* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, SpeciesCollection, ptr, fn);
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("GrainSpeciesCollection", f);
  for (int i = 0; i < OnlyGrainSpLUT::NUM_ENTRIES; i++) {
    f(VIS_MEMBER_NAME("data[...]"), obj0.data[i], obj1.data[i], vis::idx_range_len_multiple(1));
  }
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(GrainSpeciesCollection* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, GrainSpeciesCollection, ptr, fn)
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
  namespace vis = ::grackle::impl::visitor;

  vis::begin_visit("ColRecRxnRateCollection", f);
  for (int i = 0; i < ColRecRxnLUT::NUM_ENTRIES; i++) {
    f(VIS_MEMBER_NAME("data[...]"), obj0.data[i], obj1.data[i], vis::idx_range_len_multiple(1));
  }
  vis::end_visit(f);
}

/// implements the visitor design pattern
///
/// @param ptr[in,out] Members of the specified object will be visited
/// @param fn[in] Calls function that will be applied to each function
template <class UnaryVisitor>
inline void visit_member(ColRecRxnRateCollection* ptr, UnaryVisitor fn) {
  GRIMPL_IMPL_VISIT_MEMBER(visit_member_pair, ColRecRxnRateCollection, ptr, fn)
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

#endif /* INTERNAL_TYPES_HPP */
