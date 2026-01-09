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

/// This type encodes a sequence of partitions
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
/// Important invariants:
/// - every index lies in exactly 1 partition
/// - each partition has a non-zero start
/// - each partitions has a stop value that is greater or equal it the start
///
/// @note
/// While this seems like overkill, I'm convinced this will come up a bunch
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

inline int PartMap_n_partitions(const PartMap* s) { return s->n_parts; }

inline int PartMap_n_idx(const PartMap* s) {
  return (s->n_parts == 0) ? 0 : s->right_idx_bounds[s->n_parts - 1];
}

/// @todo Perhaps we should reconcile with field_flat_index_range?
struct IdxInterval {
  int start;
  int stop;
};

inline IdxInterval PartMap_part_prop(const PartMap* s, int part_id) {
  if (part_id < 0 || part_id >= s->n_parts) {
    return IdxInterval{-1, -1};
  }
  return IdxInterval{
      /*start=*/(part_id == 0) ? 0 : s->right_idx_bounds[part_id - 1],
      /*stop=*/s->right_idx_bounds[part_id]};
}

/// holds info pertaining to the partition holding an index
struct RelativePartitionInfo {
  /// id of the partition
  int part_id;
  /// offset of the index relative to the start of the partition
  int start_offset;
};

inline RelativePartitionInfo PartMap_find_part_info(const PartMap* s, int idx) {
  if (idx >= 0) {  // simple, stupid, linear search
    for (int part_id = 0; part_id < s->n_parts; part_id++) {
      if (idx < s->right_idx_bounds[part_id]) {
        int part_start = (part_id == 0) ? 0 : s->right_idx_bounds[part_id - 1];
        return RelativePartitionInfo{part_id, idx - part_start};
      }
    }
  }
  return {-1, -1};
}
}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_PARTSEQ_HPP