//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define/declare @ref PartMap
///
//===----------------------------------------------------------------------===//
#ifndef SUPPORT_PARTMAP_HPP
#define SUPPORT_PARTMAP_HPP
#include "support/config.hpp"
#include "support/FrozenKeyIdxBiMap.hpp"
#include "status_reporting.h"

namespace GRIMPL_NAMESPACE_DECL {
namespace partmap_detail {
/// we are starting with an arbitrarily low number
inline constexpr int MAX_LEN = 4;
}  // namespace partmap_detail

using partition_descr_type = int;

/// This type encodes a table of partitions
///
/// The premise of this type is extremely simple:
/// - we may work with sequences of data that we need to access by index.
///   (the indices may be described by a FrozenKeyIdxBiMap)
/// - the sequences are commonly subdivided into different partitions that
///   have special semantic meaning (or explicitly don't have a meaning).
///   For our purposes:
///   - each index must lie in exactly 1 partition
///   - each partition spans 0 or more contiguous indices
/// - Instances of this type exist to provide information about the indices
///   bounding a partition **AND** to find the partition containing an index
///
/// @par Basic Vocabulary
/// To avoid confusion (especially with abbreviations):
/// - an index is an index in array that is partitioned
/// - the identifier of a partition is always a *partition descriptor* (and is
///   abbreviated as `pd`)
///
/// @par Example
/// Let's consider a simple example. This example breaks up a 9 element array
/// into 3 partitions, PartitionName::A, PartitionName::B, & PartitionName::C,
/// in a manner sketched below:
///
/// @code{unparsed}
///        ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐
/// data:  │  0  │  1  │  2  │  3  │  4  │  5  │  6  │  7  │  8  │
///        └─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘
///        ┗━━━━━━━━━━━┳━━━━━━━━━━━┻━━━━━┳━━━━━┻━━━━━━━━┳━━━━━━━━┛
///                    ┃                 ┃              ┃
/// PartitionName::A  ━┛                 ┃              ┃
/// PartitionName::C  ━━━━━━━━━━━━━━━━━━━┛              ┃
/// PartitionName::B  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
/// @endcode
///
/// The names of the partitions are typically be declared as an enum (ideally,
/// they would be "scoped enums," but that is tricky for reasons explained
/// later). The declaration might look like:
///
/// @code{C++}
/// namespace PartitionName { enum {A, B, C}; }
/// @endcode
///
/// The code that constructs the enum would look something like:
///
/// @code{C++}
///   /* Construct the partition map */
///   int pds[3] = {PartitionName::A, PartitionMap::C, PartitionMap::B};
///   int size[3] = {4, 2, 3};
///   PartMap m = new_PartitionMap(pds, sizes, 3);
///   if (!PartitionMap_is_ok(&m)) { /* <error-propagation ...> */ }
///
///   /* query the bounds associated with PartitionMap::C */
///   IdxInterval bounds = PartMap_part_bounds(&m, PartitionMap::C);
///   assert(bounds.start == 4);
///   assert(bounds.stop == 6); /
///
///   /* query the partition associated with index 7 */
///   partmap::IdxSearch search_rslt = PartMap_search_idx(&m, 7);
///   assert(search_rslt.has_val);
///   assert(search_rslt.pd == PartitionName::B);
/// @endcode
///
/// @par Motivation
/// The concept modeled by this datatype models is **EXTREMELY** common in
/// scientific software written in languages like C or Fortran, and it's
/// almost always handled very implicitly.
///
/// In contrast, this type primarily exists for the sake of being explicit.
/// In fact, we are explicitly trading off a tiny amount of performance (not
/// in any performance-critical loops) for the benefit of explicitness. While
/// we continue refactoring Grackle, this tradeoff is warranted. (We can
/// revisit this tradeoff once we finish refactoring)
///
/// @par Relevant Applications
/// For added context, I can think of at least 3 contexts where this construct
/// is relevant:
/// 1. For tracking which indices in a lookup table refer to chemical species
///    vs dust species
/// 2. I hope to consolidate the various arrays used to track the computed
///    arrays into a single row.
/// 2. (Hypothetical) when we adopt a new interface type to replace
///    @ref grackle_field_data (similar in spirit to the one in PR #271), the
///    idea is to use have users use string keys to access key-value pairs.
///    When we do this, it would be really useful to have our API support
///    queries of field category or ranges of keys in a given category (in
///    particular: chemical-species-fields, dust-species fields, injection
///    path metal density, etc.)
///
/// @par Ideas for improvement
/// There are 2 ideas:
/// 1. Make this machinery compatible with "scoped enums."
///    - For context, regular C-style enums implicitly converts to and from
///      integer datatypes, whereas scoped enums (aka "class enums") require
///      explicit casts. This produces nice behaviors:
///      1. if a function argument has type `E`, where `E` is a scoped enum,
///         the compiler reports an error if you try to pass anything other
///         than an enumerator declared within the declaration of `E`.
///      2. if you try to pass an enumerator that was declared as part of a
///         scoped enum, `E`, to a function argument with a type other than
///         `E`, the compiler will report an error.
///    - To make this machinery compatible with scoped enums, we would need to
///      make PartMap into a template struct where `partition_descr_type` is a
///      template parameter
/// 2. Better Performance: We can almost always assume a particular ordering
///    of the partition descriptors.
///    - There are a few ways we can take advantage of this.
///    - It's probably wise to hold off on this until we start using this type
///      in a bunch of places and performance is a demonstrated issue (I'm a
///      little skeptical, since this probably won't be used deep within any
///      nested loops)
///
/// @note
/// The contents of this struct should be considered an implementation
/// detail! Always prefer the associated functions (they are defined in such
/// a way that they should be inlined)
struct PartMap {
  /// number of partitions
  int n_parts;
  /// the list of partition descriptors associated with each partition
  partition_descr_type pd_array[partmap_detail::MAX_LEN];
  /// the upper bounds on each partition
  int right_idx_bounds[partmap_detail::MAX_LEN];
};

namespace partmap_detail {
// helper function (that I wish didn't exist)
inline PartMap mk_invalid() {
  PartMap out;
  out.n_parts = -1;
  for (int i = 0; i < partmap_detail::MAX_LEN; i++) {
    out.pd_array[i] = 0;
    out.right_idx_bounds[i] = 0;
  }
  return out;
}
}  // namespace partmap_detail

/// Construct a PartMap from the sizes of each partition.
///
/// @param[in] pds Array of unique partition descriptors
/// @param[in] sizes Holds the number of indices for each partition.
/// @param[in] n_parts The number of partitions
///
/// @note
/// Callers should pass the returned value to @ref PartMap_isok
/// to check whether there was an error during creation. This is pretty
/// ugly/clunky, but it's the only practical way to achieve comparable behavior
/// to other internal data types. The best alternatives involve things like
/// std::optional or converting this type to a simple C++ class.
inline PartMap new_PartMap(const partition_descr_type* pds, const int* sizes,
                           int n_parts) {
  if (n_parts == 0) {
    PartMap out;
    out.n_parts = n_parts;
    out.pd_array[0] = -1;
    out.right_idx_bounds[0] = 0;
    return out;
  }

  // (in reality, any error here points to an internal logic-error)
  if (pds == nullptr || sizes == nullptr) {
    GrPrintErrMsg("pds and sizes can only be a nullptr when n_parts is 0");
    return partmap_detail::mk_invalid();
  } else if (n_parts < 0 || n_parts > partmap_detail::MAX_LEN) {
    GrPrintErrMsg("n_parts doesn't satisfy 0 <= n_parts <= %d",
                  partmap_detail::MAX_LEN);
  }

  PartMap out;
  out.n_parts = n_parts;
  int running_sum = 0;
  for (int i = 0; i < n_parts; i++) {
    // error checks:
    for (int j = 0; j < i; j++) {
      if (pds[i] == pds[j]) {
        GrPrintErrMsg("pds[%d] and pds[%d] hold the same descriptor", i, j);
        return partmap_detail::mk_invalid();
      }
    }
    if (sizes[i] < 0) {
      GrPrintErrMsg("sizes[%d] is negative", i);
      return partmap_detail::mk_invalid();
    }

    out.pd_array[i] = pds[i];
    running_sum += sizes[i];
    out.right_idx_bounds[i] = running_sum;
  }
  return out;
}

/// checks whether a creational function produced a valid partition map
///
/// @param[in] ptr Points to the partition map being checked
/// @return true if the value is ok or false if the value is invalid
///
/// @important
/// The interface of @ref PartMap sets values in a very particular way to
/// signal that an instance is in an invalid state. This function @b ONLY
/// checks for that particular signature.
inline bool PartMap_is_ok(const PartMap* ptr) { return (ptr->n_parts >= 0); }

/// Destroys the internal data tracked by an instance
///
/// @param[in] ptr A non-null pointer to a valid partition map
///
/// @note
/// The @b ONLY exists for consistency with other data types. (It's
/// unnecessary to call this function)
inline void drop_PartMap(PartMap* ptr) {
  return;  // NO-OP
}

/// number of partitions in the partition map
inline int PartMap_n_partitions(const PartMap* m) { return m->n_parts; }

/// number of indices bounded by the partition map
inline int PartMap_n_idx(const PartMap* m) {
  return (m->n_parts == 0) ? 0 : m->right_idx_bounds[m->n_parts - 1];
}

/// @todo Perhaps we should reconcile with field_flat_index_range?
struct IdxInterval {
  int start;
  int stop;
};

/// Query the interval of indices that bound a partition
///
/// @param[in] m Valid pointer to a partition map
/// @param[in] pd The partition descriptor to query
inline IdxInterval PartMap_part_bounds(const PartMap* m,
                                       partition_descr_type pd) {
  // simple, stupid, linear search
  for (int i = 0; i < m->n_parts; i++) {
    if (pd == m->pd_array[i]) {
      return IdxInterval{/*start=*/(i == 0) ? 0 : m->right_idx_bounds[i - 1],
                         /*stop=*/m->right_idx_bounds[i]};
    }
  }
  return IdxInterval{-1, -1};
}

namespace partmap {
/// holds info pertaining to the partition holding an index
struct IdxSearch {
  /// indicates whether the index was found
  bool has_val;

  /// The index being searched
  ///
  /// @note
  /// This is primarily tracked in order to make the key_partition_search
  /// convenience function provide more useful results
  int index;

  /// the partition descriptor
  partition_descr_type pd;

  /// offset of the index relative to the start of the partition
  int start_offset;
};
}  // namespace partmap

/// search for the partition containing an index
///
/// @param[in] m Valid pointer to a partition map
/// @param[in] idx The index to search for
inline partmap::IdxSearch PartMap_search_idx(const PartMap* m, int idx) {
  if (idx >= 0) {
    // simple, stupid, linear search
    for (int i = 0; i < m->n_parts; i++) {
      if (idx < m->right_idx_bounds[i]) {
        int part_start = (i == 0) ? 0 : m->right_idx_bounds[i - 1];
        return partmap::IdxSearch{true, idx, m->pd_array[i], idx - part_start};
      }
    }
  }
  return {false, idx, m->pd_array[0], -1};
}

/// search for the partition containing @p key
///
/// This is a convenience function that can be used when the keys of @p bimap
/// are partitioned according to @p part_map
///
/// @param[in] part_map Describes the partitioning of indices
/// @param[in] bimap Describes bi-directional mapping between indices & keys
/// @param[in] key The key to that is being searched
inline partmap::IdxSearch key_partition_search(const PartMap* part_map,
                                               const FrozenKeyIdxBiMap* bimap,
                                               const char* key) {
  bimap::AccessRslt maybe_idx = FrozenKeyIdxBiMap_find(bimap, key);
  return PartMap_search_idx(part_map,
                            maybe_idx.has_value ? -1 : maybe_idx.value);
}
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_PARTMAP_HPP
