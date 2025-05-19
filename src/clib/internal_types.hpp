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

/// allocates the contents of a new CoolHeatScratchBuf
///
/// @param nelem The number of elements in each buffer
CoolHeatScratchBuf new_CoolHeatScratchBuf(int nelem);

/// performs cleanup of the contents of CoolHeatScratchBuf
///
/// This effectively invokes the destructor
void drop_CoolHeatScratchBuf(CoolHeatScratchBuf*);

/// holds other assorted scratch buffers used within cool1d_multi_g
///
/// @note
/// Prior to transcription, there was a distinction between the buffers that
/// are now inside of this struct and the buffers in CoolHeatScratchBuf. The
/// distinction has been preserved during transcription, but it is not clear
/// how real the distinction truly is.
struct Cool1DMultiScratchBuf {
  double* tgasold = nullptr;
  double* mynh = nullptr;
  double* myde = nullptr;
  double* gammaha_eff = nullptr;
  double* gasgr_tdust = nullptr;
  double* regr = nullptr;
};

/// allocates the contents of a new Cool1DMultiScratchBuf
///
/// @param nelem The number of elements in each buffer
Cool1DMultiScratchBuf new_Cool1DMultiScratchBuf(int nelem);

/// performs cleanup of the contents of Cool1DMultiScratchBuf
///
/// This effectively invokes the destructor
void drop_Cool1DMultiScratchBuf(Cool1DMultiScratchBuf*);

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

/// allocates the contents of a new LogTLinInterpScratchBuf
///
/// @param nelem The number of elements in each buffer
LogTLinInterpScratchBuf new_LogTLinInterpScratchBuf(int nelem);

/// performs cleanup of the contents of LogTLinInterpScratchBuf
///
/// This effectively invokes the destructor
void drop_LogTLinInterpScratchBuf(LogTLinInterpScratchBuf*);


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

/// allocates the contents of a new GrainSpeciesCollection
///
/// @param nelem The number of elements in each buffer
GrainSpeciesCollection new_GrainSpeciesCollection(int nelem);

/// performs cleanup of the contents of GrainSpeciesCollection
///
/// This effectively invokes a destructor
void drop_GrainSpeciesCollection(GrainSpeciesCollection*);


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


/// allocates the contents of a new SpeciesCollection
///
/// @param nelem The number of elements in each buffer
SpeciesCollection new_SpeciesCollection(int nelem);

/// performs cleanup of the contents of SpeciesCollection
///
/// This effectively invokes a destructor
void drop_SpeciesCollection(SpeciesCollection*);


/// holds properties about Collisional and Recombination Reaction Rates that
/// behave in the "standard" way (i.e. we interpolate in 1D with respect to
/// log T)
///
/// This operates in a similar manner to SpeciesCollection (i.e. we use
/// `ColRecRxnLUT::<entry>` to lookup values for the desired rate).
struct ColRecRxnRateCollection {
  double* data[ColRecRxnLUT::NUM_ENTRIES];
};

/// allocates the contents of a new ColRecRxnRateCollection
///
/// @param nelem The number of elements in each buffer
ColRecRxnRateCollection new_ColRecRxnRateCollection(int nelem);

/// performs cleanup of the contents of ColRecRxnRateCollection
///
/// This effectively invokes a destructor
void drop_ColRecRxnRateCollection(ColRecRxnRateCollection*);

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

/// allocates the contents of a new PhotoRxnRateCollection
///
/// @param nelem The number of elements in each buffer
PhotoRxnRateCollection new_PhotoRxnRateCollection(int nelem);

/// performs cleanup of the contents of PhotoRxnRateCollection
///
/// This effectively invokes a destructor
void drop_PhotoRxnRateCollection(PhotoRxnRateCollection*);

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

/// allocates the contents of a new ChemHeatingRates
///
/// @param nelem The number of elements in each buffer
ChemHeatingRates new_ChemHeatingRates(int nelem);

/// performs cleanup of the contents of ColRecRxnRateCollection
///
/// This effectively invokes a destructor
void drop_ChemHeatingRates(ChemHeatingRates*);



} // namespace grackle::impl

#endif /* INTERNAL_TYPES_HPP */
