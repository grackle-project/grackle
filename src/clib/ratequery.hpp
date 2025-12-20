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

/// Describes properties about the data in an entry
struct EntryProps {
  int ndim;
  int shape[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
};

inline EntryProps mk_invalid_EntryProps() {
  EntryProps out;
  out.ndim = -1;
  return out;
}

inline bool EntryProps_is_valid(EntryProps obj) { return obj.ndim >= 0; }

/// A queryable entity
struct Entry {
  double* data;
  const char* name;
  EntryProps props;
};

inline Entry mk_invalid_Entry() {
  return Entry{nullptr, nullptr, mk_invalid_EntryProps()};
}

/// Constructs an Entry
inline Entry new_Entry(double* rate, const char* name) {
  return Entry{rate, name, mk_invalid_EntryProps()};
}

/// a recipe for querying 1 or more entries from a chemistry_data_storage
/// pointer given an index. A recipe generally lazily creates an Entry as
/// the underlying data is fetched.
///
/// If `N` denotes the number of entries that can be queried by a given recipe,
/// then this function should produce unique entries for each unique index that
/// satisfies `0 <= index <= (N-1)`
typedef Entry fetch_Entry_recipe_fn(chemistry_data_storage*, int);

// temporary forward declaration
struct EntrySet;

/// Describes a registry of queryable entries
struct Registry {
  /// number of entries
  int n_entries;
  /// number of contained EntrySets
  int n_sets;
  /// stores the minimum rate_id for each EntrySet
  int* id_offsets;
  /// stores sets of entries
  EntrySet* sets;
};

/// construct a new registry
Registry new_Registry(const chemistry_data&);

/// deallocate the contents of a registry
void drop_Registry(Registry* ptr);

/// An interface for gradually configuring a Registry
///
/// @par Context within the codebase
/// This is intended to be an ephemeral object that only lives during within
/// the function that initializes a Grackle solver from user-specified
/// parameters
/// - the instance is created at the start of this function
/// - a pointer to the instance is passed to various initialization functions.
///   These functions can use the instance to register entries that will be
///   accessible in the resulting Registry
/// - if all configuration has gone well, this instance will be used to
///   initialize the Registry
/// - the instance is **always** destroyed when it is time to exit the solver
///   initialization function
///
/// @par Motivation
/// There are 2 main sources of motivation
/// 1. It lets us locate a "recipe-routine" for accessing data next to the
///    routines that initialize the same data. This makes sense from a
///    code-organization perspective. Moreover, if we conditionally initialize
///    data, it will be easier to ensure that recipies won't try to provide
///    access to the uninitialized data
/// 2. this will make it easier for us to handle data that Grackle only
///    initializes for the sake of supporting user queries (at the moment,
///    this is hypothetical)
///
/// @important
/// Other parts of grackle should refrain from directly accessing the internals
/// of this function (i.e. they should only use the associated methods)
struct RegBuilder {
  int capacity;
  int len;
  EntrySet* sets;
};

/// initialize a new instance
inline RegBuilder new_RegBuilder() { return {0, 0, nullptr}; }

/// deallocates all storage within a RegBuilder instance
void drop_RegBuilder(RegBuilder* ptr);

/// register a recipe for accessing scalar values
///
/// @param[inout] ptr The RegBuilder that will be updated
/// @param[in] n_entries The number of entries accessible through the recipe
/// @param[in] recipe_fn The recipe being registered
/// @param[in] common_props The properties shared by each Entry in the recipe
///
/// @returns GR_SUCCESS if successful, otherwise returns a different value
int RegBuilder_recipe_scalar(RegBuilder* ptr, int n_entries,
                             fetch_Entry_recipe_fn* recipe_fn);

/// register a recipe for accessing 1D arrays
///
/// @param[inout] ptr The RegBuilder that will be updated
/// @param[in] n_entries The number of entries accessible through the recipe
/// @param[in] recipe_fn The recipe being registered
/// @param[in] common_len The length shared by each 1D array accessible
///     through this recipe.
///
/// @returns GR_SUCCESS if successful, otherwise returns a different value
int RegBuilder_recipe_1d(RegBuilder* ptr, int n_entries,
                         fetch_Entry_recipe_fn* recipe_fn, int common_len);

/// build a new Registry.
///
/// In the process, the current Registry is consumed; it's effectively reset to
/// the state immediately after it was initialized. (This lets us avoid
/// reallocating lots of memory)
///
/// @note
/// For safety, the caller should still plan to call drop_RegBuilder
Registry RegBuilder_consume_and_build(RegBuilder* ptr);

/// Describes the set of entries that can be access through a given recipe
struct EntrySet {
  /// number of entries in the current set
  int len;

  /// a function pointer that can be used to access entries through a recipe
  fetch_Entry_recipe_fn* recipe_fn;

  /// properties used by all entries accessed through a recipe
  ///
  /// In more detail, an entry returned by `recipe_fn` has its `props` member
  /// overwritten by this value
  EntryProps common_props;
};

/** @}*/  // end of group

}  // namespace grackle::impl::ratequery

#endif  // RATEQUERY_HPP
