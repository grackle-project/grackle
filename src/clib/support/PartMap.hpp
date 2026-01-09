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
#ifndef SUPPORT_PARTSEQ_HPP
#define SUPPORT_PARTSEQ_HPP
#include "support/config.hpp"
#include "status_reporting.h"

namespace GRIMPL_NAMESPACE_DECL {
/// we are starting with an arbitrarily low number
inline constexpr int MAX_PART_SEQ_LEN = 4;

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
/// For context, the following diagram sketches an example where an array of 9
/// elements has been partitioned into 3 unequally sized partitions
///
/// @code{unparsed}
///        ┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐
/// data:  │  0  │  1  │  2  │  3  │  4  │  5  │  6  │  7  │  8  │
///        └─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘
///        ┗━━━━━━━━━━━┳━━━━━━━━━━━┻━━━━━┳━━━━━┻━━━━━━━━┳━━━━━━━━┛
///                    ┃                 ┃              ┃
///               partition 0       partition 1    partition 2
///                start: 0          start 4        start 6
///                stop: 4           stop 6         stop 9
/// @endcode
///
/// @par Basic Vocabulary
/// To avoid confusion (especially with abbreviations):
/// - an index is an index in array that is partitioned
/// - the identifier of a partition is always a *partition descriptor* (and is
///   abbreviated as `pd`)
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
struct PartMap {
  /// number of partitions
  int n_parts;
  /// the upper bounds on each partition
  int right_idx_bounds[MAX_PART_SEQ_LEN];
};

/// Construct a PartMap from the sizes of each partition.
///
/// @param[in] part_sizes Holds the number of indices for each partition.
/// @param[in] n_parts The number of partitions
inline PartMap PartMap_from_part_sizes(const int* part_sizes, int n_parts) {
  if (part_sizes == nullptr && n_parts == 0) {
    PartMap out;
    out.n_parts = n_parts;
    out.right_idx_bounds[0] = 0;
    return out;
  }

  // todo: gracefully handle errors
  // (in reality, any error here points to an internal logic-error)
  if (part_sizes == nullptr) {
    GR_INTERNAL_ERROR("part_sizes is nullptr when n_parts isn't 0");
  } else if (n_parts == 0) {
    GR_INTERNAL_ERROR("n_parts is 0 when part_sizes isn't nullptr");
  } else if (n_parts < 0 || n_parts > MAX_PART_SEQ_LEN) {
    GR_INTERNAL_ERROR("n_parts doesn't satisfy 0 <= n_parts <= %d",
                      MAX_PART_SEQ_LEN);
  }

  PartMap out;
  out.n_parts = n_parts;
  int running_sum = 0;
  for (int part_id = 0; part_id < n_parts; part_id++) {
    GR_INTERNAL_REQUIRE(part_sizes[part_id] >= 0,
                        "partition size can't be negative");
    running_sum += part_sizes[part_id];
    out.right_idx_bounds[part_id] = running_sum;
  }
  return out;
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
inline IdxInterval PartMap_part_bounds(const PartMap* m, int pd) {
  if (pd < 0 || pd >= m->n_parts) {
    return IdxInterval{-1, -1};
  }
  return IdxInterval{
      /*start=*/(pd == 0) ? 0 : m->right_idx_bounds[pd - 1],
      /*stop=*/m->right_idx_bounds[pd]};
}

/// holds info pertaining to the partition holding an index
struct IdxPartSearch {
  /// indicates whether the index was found
  bool is_valid;
  /// the partition descriptor
  int pd;
  /// offset of the index relative to the start of the partition
  int start_offset;
};

/// search for the partition containing an index
///
/// @param[in] m Valid pointer to a partition map
/// @param[in] idx The index to search for
inline IdxPartSearch PartMap_search_idx(const PartMap* m, int idx) {
  if (idx >= 0) {  // simple, stupid, linear search
    for (int pd = 0; pd < m->n_parts; pd++) {
      if (idx < m->right_idx_bounds[pd]) {
        int part_start = (pd == 0) ? 0 : m->right_idx_bounds[pd - 1];
        return IdxPartSearch{true, pd, idx - part_start};
      }
    }
  }
  return {false, -1, -1};
}
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_PARTSEQ_HPP