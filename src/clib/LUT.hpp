// See LICENSE file for license and copyright information

/// @file LUT.hpp
/// @brief Declares some lookup-tables used internally by Grackle

#ifndef LUT_HPP
#define LUT_HPP

#ifndef __cplusplus
#error "This file can only be read by a c++ compiler"
#endif

// This namespace holds enumerators defined for every species. The enumerator
// values can be used in a lookup table.
//  - the idea here is that we can create arrays to hold data for every species
//    (rather than structs).
//  - in this picture, we can use the enumerator values to lookup the value
//    corresponding to a given species (for example, we could use `SpLUT::HI`
//    to look up values related to HI). This essentially gives us all the
//    benefits of a struct (instead of an array).
//  - if we are using an array to organize data instead of a struct, we can
//    then write certain operations that are identical for all species (e.g.
//    allocating buffers, deallocating buffers, initializing buffers) in a much
//    more concise (more maintainable) way.
//
// IMPORTANT: this is a private implementation detail. We should never
// unintentionally expose this information directly in the public API (i.e. we
// want to maintain our flexibility). With that said, there are hypothetical
// scenarios where we might deliberately expose it (for performance purposes)
// -- but we should think long and hard about it before we do
//
// Why use a namespace?
// --------------------
//  - it lets us use short/simple names for each enumerator without feat of
//    collision
//  - we later have the option to convert SpLUT from a namespace to a struct
//    (it still contains the enumerators).
//    - This would be useful if we wanted to write a template function that
//      accepts different choices of LUTs as a template argument (without
//      modifying any existing code). This would be hypothetically useful if we
//      wanted to create specializations of the code that only have LUTs for
//      a subset of species (e.g. based on primordial_chemistry == 1)
//    - we currently have no plans to do this.
//
// Other Thoughts
// --------------
// - maybe we should store the dust species in a separate lookup table?
// - adopting data-structures that leverage this LUT, is probably an important
//   stepping stone to implementing a system where reaction rates are
//   dynamically specified. After adopting such a system, we may not even need
//   a lookup table any more (this all assumes, of course, that such a system
//   is adequately performant)
namespace SpLUT{

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

  NUM_ENTRIES // <- this is always last (so it specifies the number of species)
}; // enum

} // SpLUT namespace


// Defines the LUT for Standard Collisional and Recombination reaction rates
namespace ColRecRxnLUT {

enum {
  #define ENTRY(NAME) NAME,
  #include "col_rec_rxn_rate_members.def"
  #undef ENTRY

  NUM_ENTRIES // <- this is always last (so it specifies the number of species)
}; // enum

} // ColRecRxnLUT namespace


#endif /* LUT_HPP */
