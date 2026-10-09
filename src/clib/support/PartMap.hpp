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
#include "support/error.hpp"
#include "support/expected.hpp"
#include "support/index_helper.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace partmap_detail {
/// we are starting with an arbitrarily low number
inline constexpr int MAX_LEN = 4;
}  // namespace partmap_detail

/// @todo Perhaps we should reconcile with FieldFlatIndexRange?
struct IdxInterval {
  int start;
  int stop;
};

namespace partmap {
/// @brief encodes the result of a search for an index
///
/// See @ref PartMap::search_idx for more details
template <class PartitionDescrT>
struct IdxSearch {
  /// The index being searched
  ///
  /// @note
  /// This is primarily tracked in order to make the key_partition_search
  /// convenience function provide more useful results
  int index;

  /// the partition descriptor
  PartitionDescrT pd;

  /// offset of the index relative to the start of the partition
  int start_offset;

  /// @brief overloads equality operation (for testing purposes)
  bool operator==(const IdxSearch&) const = default;
};
}  // namespace partmap

/// @brief This class template encodes a table of partitions
///
/// The premise of this type is extremely simple:
/// - we may work with sequences of data that we need to access by index.
///   (the indices may be described by a FrozenKeyIdxBiMap, see PR #492)
/// - the sequences are commonly subdivided into different partitions that
///   have special semantic meaning (or explicitly don't have a meaning).
///   For our purposes:
///   - each index must lie in exactly 1 partition
///   - each partition spans 0 or more contiguous indices
/// - Instances of this type exist to provide information about the indices
///   bounding a partition **AND** to find the partition containing an index
///
/// @tparam PartitionDescrT the type of the partition descriptior. This should
///      generally be a type like an `int` or a scoped enum. Using a scoped
///      enum is generally preferable as it reduces the likelihood of
///      programming errors.
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
/// The names of the partitions typically be declared as an enum (ideally,
/// they would be "scoped enums," but that is tricky for reasons explained
/// later). The declaration might look like:
///
/// @code{C++}
/// enum class PartitionName {A, B, C};
/// @endcode
///
/// The code that constructs the enum would look something like:
///
/// @code{C++}
///   /* Construct the partition map */
///   PartitionName pds[3] = {
///       PartitionName::A, PartitionName::C, PartitionName::B};
///   int size[3] = {4, 2, 3};
///   Expected<PartMap<PartitionName>, Error> rslt
///       = PartMap<PartitionName>::create(pds, sizes, 3);
///   if (!rslt.has_value()) { /* <error-propagation ...> */ }
///
///   /* move PartMap out of rslt (deepcopies of PartMap are also cheap) */
///   PartMap<PartitionName> m = std::move(rslt).value();
///
///   /* query the bounds associated with PartitionName::C */
///   IdxInterval bounds = m.part_bounds(PartitionName::C);
///   assert(bounds.start == 4);
///   assert(bounds.stop == 6);
///
///   /* query the partition associated with index 7 */
///   std::optional<partmap::IdxSearch<PartitionName>> search_rslt =
///       m.search_idx(&m, 7);
///   assert(search_rslt.has_value());
///   assert(search_rslt.value().pd == PartitionName::B);
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
/// Better Performance: We can almost always assume a particular ordering
/// of the partition descriptors.
/// - There are a few ways we can take advantage of this.
/// - It's probably wise to hold off on this until we start using this type
///   in a bunch of places and performance is a demonstrated issue (I'm a
///   little skeptical, since this probably won't be used deep within any
///   nested loops)
template <class PartitionDescrT>
class PartMap {
  /// number of partitions
  int n_parts_;
  /// the list of partition descriptors associated with each partition
  PartitionDescrT pd_array_[partmap_detail::MAX_LEN];
  /// the upper bounds on each partition
  int right_idx_bounds_[partmap_detail::MAX_LEN];

public:
  /// @brief default constructor (constructs a map with 0 partitions)
  PartMap() {
    n_parts_ = 0;
    for (int i = 0; i < partmap_detail::MAX_LEN; i++) {
      right_idx_bounds_[i] = 0;
    }
  }

  /// @brief Construct a PartMap from the sizes of each partition.
  ///
  /// @param[in] pds Array of unique partition descriptors
  /// @param[in] sizes Holds the number of indices for each partition.
  /// @param[in] n_parts The number of partitions
  static Expected<PartMap, Error> create(const PartitionDescrT* pds,
                                         const int* sizes, int n_parts) {
    // (in reality, any error here points to an internal logic-error)
    if (n_parts != 0 && (pds == nullptr || sizes == nullptr)) {
      return Unexpected(Error::msg_literal(
          "pds and sizes can only be a nullptr when n_parts is 0"));
    } else if (n_parts < 0 || n_parts > partmap_detail::MAX_LEN) {
      return Unexpected(
          Error::msgf("n_parts doesn't satisfy 0 <= n_parts <= %d",
                      partmap_detail::MAX_LEN));
    }

    PartMap out;
    if (n_parts == 0) {
      return out;
    }

    int running_sum = 0;
    out.n_parts_ = n_parts;
    for (int i = 0; i < n_parts; i++) {
      // error checks:
      for (int j = 0; j < i; j++) {
        if (pds[i] == pds[j]) {
          return Unexpected(Error::msgf(
              "pds[%d] and pds[%d] hold the same descriptor", i, j));
        }
      }
      if (sizes[i] < 0) {
        return Unexpected(Error::msgf("sizes[%d] is negative", i));
      }

      out.pd_array_[i] = pds[i];
      running_sum += sizes[i];
      out.right_idx_bounds_[i] = running_sum;
    }
    return out;
  }

  // use default move/copy constructors & assignment operations
  PartMap(const PartMap&) = default;
  PartMap(PartMap&&) = default;
  PartMap& operator=(const PartMap&) = default;
  PartMap& operator=(PartMap&&) = default;

  /// number of partitions in the partition map
  int n_partitions() const { return n_parts_; }

  /// number of indices bounded by the partition map
  int n_idx() const {
    return (n_parts_ == 0) ? 0 : right_idx_bounds_[n_parts_ - 1];
  }

  /// @brief Query the interval of indices that bound a partition
  ///
  /// This function is explicitly written such that the start and stop values
  /// of the returned interval are equal to each other if @p pd
  /// - corresponds to a contained partition with a length of 0
  /// - isn't a contained partition
  /// Often times, the caller won't care about this distinction
  ///
  /// @param pd The partition descriptor to query
  /// @return An interval of indices that bound a partition. If @p pd is not
  ///     contained, the start and stop values of the range are set to a
  ///     negative value.
  IdxInterval part_bounds(PartitionDescrT pd) const {
    // simple, stupid, linear search
    for (int i = 0; i < n_parts_; i++) {
      if (pd == pd_array_[i]) {
        return IdxInterval{/*start=*/(i == 0) ? 0 : right_idx_bounds_[i - 1],
                           /*stop=*/right_idx_bounds_[i]};
      }
    }
    return IdxInterval{-1, -1};
  }

  /// @brief search for the partition containing an index
  ///
  /// @param[in] idx The index to search for
  ///
  /// @note
  /// Originally I was a little hesitant to return a std::optional from this
  /// function, since it could theoretically give GPUs some trouble (in
  /// practice, it would probably work on most platforms, but might be slow).
  /// After giving it some thought, I was reminded that unlike many other
  /// methods of this type, you don't want to execute this function on GPUs
  /// (you really just want to use it at startup while initializing Grackle)
  std::optional<partmap::IdxSearch<PartitionDescrT>> search_idx(int idx) const {
    if (idx >= 0) {
      // simple, stupid, linear search
      for (int i = 0; i < n_parts_; i++) {
        if (idx < right_idx_bounds_[i]) {
          int part_start = (i == 0) ? 0 : right_idx_bounds_[i - 1];
          return partmap::IdxSearch<PartitionDescrT>{idx, pd_array_[i],
                                                     idx - part_start};
        }
      }
    }
    return std::nullopt;
  }
};

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // SUPPORT_PARTMAP_HPP
