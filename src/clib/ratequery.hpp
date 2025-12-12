//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares functionality for querying rate data
///
//===----------------------------------------------------------------------===//
#ifndef RATEQUERY_HPP
#define RATEQUERY_HPP

#include "grackle.h"

namespace grackle::impl::ratequery {

/// @defgroup Dynamic Rate Query Machinery
///
/// This group of entities provides the machinery for implementing the API for
/// dynamically accessing rate data.
///
/// The design of this machinery is based on a simple model:
/// - the query machinery queries a Registry
/// - you can think of a Registry as a container. Each item is an Entry that
///   describes a unique queryable rate.
///
/// The actualy implementation is a little more sophisticated. In practice is
/// divided into subsets of `Entry` instances. Moreover, subsets are free to
/// treat `Entry` instances as ephemeral objects (i.e. a subset can lazily
/// construct a new `Entry` instance and is allowed to forget about the
/// instance).
///
/// Design Considerations
/// =====================
/// Ideally, this machinery should balance 3 design considerations: (1) memory
/// usage, (2) impact on grackle solver initialization, (3) Query Performance.
///
/// The relative of importance of these considerations should be informed by
/// the fact that dynamic rate API logic is not central to querying logic is
/// not of central important to Grackle as a library:
/// - currently, none of the rates ever need to queried for regular Grackle
///   usage (i.e. they are only queried for debugging/experimentation)
/// - moreover, if we do need to query some information in certain
///   configurations (e.g. assumed grain species yields for different injection
///   pathways), the data probably isn't in the most useful format for
///   everyone. Even if queries were ultra-fast, a performance-oriented user
///   would only query that information once, repack the data so that it's
///   structured in a more useful way, and cache the repacked data.
///
/// In principle, if we had an idea for an delivered significantly better
/// performance, at the cost of increased memory usage and/or longer
/// initialization, we could consider making construction of the registry (or
/// non-essential registry-entries) a runtime parameter.
///
/// ASIDE: The current design is very much prioritizes doing something easy
/// over runtime performance.
/** @{ */

/// A queryable entity
struct Entry {
  double* data;
  const char* name;
};

/// Constructs an entry
inline Entry new_Entry(double* rate, const char* name) {
  Entry out;
  out.data = rate;
  out.name = name;
  return out;
}

/// a recipe for querying 1 or more entries from a chemistry_data_storage
/// pointer given an index. A recipe generally lazily creates an Entry as
/// the underlying data is fetched.
///
/// If `N` denotes the number of entries that can be queried by a given recipe,
/// then this function should produce unique entries for each unique index that
/// satisfies `0 <= index <= (N-1)`
typedef Entry fetch_Entry_recipe_fn(chemistry_data_storage*, int);

/// Describes the set of entries that can be accessed through a given recipe
struct RecipeEntrySet {
  int id_offset;
  int len;
  fetch_Entry_recipe_fn* fn;
};

#define ENTRY_SET_COUNT 2

/// Describes a registry of queryable entries
struct EntryRegistry {
  int len;
  RecipeEntrySet sets[ENTRY_SET_COUNT];
};

/** @}*/  // end of group

}  // namespace grackle::impl::ratequery

#endif  // RATEQUERY_HPP
