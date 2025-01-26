// See LICENSE file for license and copyright information

/// @file LUT.hpp
/// @brief Declares some lookup-tables (LUTs) used internally by Grackle

#ifndef LUT_HPP
#define LUT_HPP

#ifndef __cplusplus
#error "This file can only be read by a c++ compiler"
#endif

// once we have transcribed more code, we should really put this header's
// contents inside of the grackle::impl namespace
// - if a lookup-table is called {name}, then its fully qualified name will be
//   `grackle::impl::{name}`.
// - when using the LUT in a function, we can shorten its name, within the
//   function to just {name} by inserting `using grackle::impl::{name};` near
//   the start of the function
// - we should hold off on doing this until more code is transcribed (since it
//   will be hard for the transcription tools to automatically handle the
//   shortenning of the fully qualified name)



/// This is collection of enumerators (localized to the `SpLUT::` scope), with
/// an enumerator named for EVERY species (primordial-species, metal-species,
/// grain-species, etc). The enumerator values are intended to be used in a
/// compile-time lookup table (LUT), that incur no overhead at runtime.
///
/// To understand this entity's purpose, imagine that we want to temporarily
/// track data for all N species known to Grackle (the data for each species is
/// stored in a buffer of common size determined at runtime). Consider the
/// following 3 ways we could accomplish this:
///   1. define a separate variable holding a pointer for every single species
///      - this is what the fortran routines historically did
///      - this won't scale well with N (e.g. a function that requires all of
///        this info needs to take N arguments)
///   2. define a single struct with a separate data member holding a pointer
///      for every single species.
///      - still doesn't scale well if for applying a common operation to the
///        data of each species (e.g. allocating/deallocating/zeroing-out)
///   3. track each data-member in a statically sized array that holds pointers
///      for each species
///      - the array length is `SpLUT::NUM_ENTRIES` (known at compile time)
///      - it is now trivial to apply a common operation to each pointer
///      - we use the enumerator values in this LUT to lookup the value
///        corresponding to a given species (for example, we could use
///        `SpLUT::HI` to look up values related to HI)
///
/// > [!important]
/// > This is a private implementation detail. We should never unintentionally
/// > expose this information directly in the public API (i.e. we want to
/// > maintain our flexibility). With that said, there are hypothetical
/// > scenarios where we might deliberately expose it (for performance purposes)
/// > -- but we should think long and hard about it before we do
///
/// Implementation
/// --------------
/// We implement choose to implement this entity by defining the enumerators
/// within a "stateless struct" (i.e. has a size of 0) that won't be
/// instantiated. The enclosing struct acts exactly like a pseudo-namespace
///  - doing this lets us use short/simple names for each enumerator without
///    fear of collision
///  - there are 2 advantages of using a struct over a namespace:
///    1. it forces us to localize the definition in 1 location (in comparison,
///       C++ lets us add identifiers to a namespace in other files)
///    2. it provides the flexibility to experiment with passing the LUT as a
///       template parameter
///       - a function that take accepts a specialized LUT (e.g. picked based
///         on primordial_chemistry) could reduce branching within loops and
///         and plausibly produce faster code.
///       - **WE DON'T CURRENTLY HAVE PLANS TO PURSUE THIS.** There are
///         much more immediately useful things to pursue first
///
/// Other Thoughts
/// --------------
/// - maybe we should store the dust species in a separate lookup table?
/// - adopting data-structures that leverage this LUT, is probably an important
///   stepping stone to implementing a system where reaction rates are
///   dynamically specified. After adopting such a system, we may not even need
///   a lookup table any more (this all assumes, of course, that such a system
///   is adequately performant)
struct SpLUT {

  // in the future, we may want to reimplement the following in terms of the
  // XMacros provided in grackle_field_data_fdatamembers.def (or we may need to
  // slightly revise the system?)
  enum {
    e,
    HI,
    HII,
    HeI,
    HeII,
    HeIII,
 
    HM,
    H2I,
    H2II,
 
    DI,
    DII,
    HDI,
 
    DM,
    HDII,
    HeHII,
 
    CI,
    CII,
    CO,
    CO2,
    OI,
    OH,
    H2O,
    O2,
    SiI,
    SiOI,
    SiO2I,
    CH,
    CH2,
    COII,
    OII,
    OHII,
    H2OII,
    H3OII,
    O2II,
 
    Mg,
    Al,
    S,
    Fe,
 
    SiM,
    FeM,
    Mg2SiO4,
    MgSiO3,
    Fe3O4,
    AC,
    SiO2D,
    MgO,
    FeS,
    Al2O3,
    reforg,
    volorg,
    H2Oice,
 
    NUM_ENTRIES // <- always last (so it specifies the number of species)
  }; // enum

}; // SpLUT struct


/// Defines the LUT for Standard Collisional and Recombination reaction rates
struct ColRecRxnLUT {

  enum {
    #define ENTRY(NAME) NAME,
    #include "col_rec_rxn_rate_members.def"
    #undef ENTRY

    NUM_ENTRIES // <- always last (so it specifies the number of species)
  }; // enum

}; // ColRecRxnLUT struct


#endif /* LUT_HPP */
