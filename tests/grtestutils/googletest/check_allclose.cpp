//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements assertions for comparing arrays
///
//===----------------------------------------------------------------------===//

#include "./check_allclose.hpp"
#include "../view.hpp"
#include "gtest/gtest.h"
#include "status_reporting.h"

#include <cmath>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include <gtest/gtest.h>

namespace grtest::arraycmp_detail {

namespace {  // stuff inside an anonymous namespace is local to this file

/// calls the unary function once for each element
///
/// if selected is provided, we skip locations where selected[i] is false
template <typename Fn>
void flat_for_each(Fn fn, int n_elements, const bool* selected) {
  if (selected == nullptr) {
    for (int flat_idx = 0; flat_idx < n_elements; flat_idx++) {
      fn(flat_idx);
    }
  } else {
    for (int flat_idx = 0; flat_idx < n_elements; flat_idx++) {
      if (selected[flat_idx]) {
        fn(flat_idx);
      }
    }
  }
}

/// Specifies a noteworthy detail about a mismatching pairs of pointers
///
/// There are 2 kinds of details:
/// - a detail describing a particular pair of mismatching elements (e.g. the
///   very first mismatch, the place where the size of the mismatch is most
///   significant). In these cases, the flattened index is tracked so that
///   the location of the mismatch and the values of the elements can be
///   printed
/// - a generic string that doesn't have any single associated location
struct MismatchDetail {
  std::string description;
  std::optional<int> flat_idx;

  // the following constructor is defined in order to make this work with
  // std::vector::emplace_back. Delete it, once we require C++20 or newer
#if __cpp_aggregate_paren_init < 201902L
  MismatchDetail(const std::string& description,
                 const std::optional<int>& flat_idx)
      : description(description), flat_idx(flat_idx) {}
#endif
};

/// collects interesting details about mismatched elements in a pair of pointers
///
/// This is called by @ref compare_ptrs_, if we determine that the compared
/// pointers contain at least one pair of mismatching elements. This loops back
/// over all pairs of elements and collects a vector of noteworthy mismatches.
///
/// @note
/// It's ok if this is a little slow, as long as it provides useful messages.
/// (After all, this logic only gets invoked when comparisons fail).
template <typename T, DataLayout Layout, typename Cmp>
std::vector<MismatchDetail> collect_details_(
    const T* actual, const T* desired, const IdxMapping<Layout>& idx_mapping,
    const bool* selection_mask, Cmp cmp_fn) {
  // define some variables that we will fill as we loop over the array
  int first_mismatch_idx = -1;

  int first_nan_mismatch_idx = -1;
  int nan_mismatch_count = 0;
  bool any_nan = false;

  double max_absDiff = 0.0;
  int max_absDiff_idx = -1;

  double max_relDiff = 0.0;
  int max_relDiff_idx = -1;

  auto fn = [&](int flat_idx) {
    bool either_nan =
        std::isnan(actual[flat_idx]) || std::isnan(desired[flat_idx]);
    any_nan = any_nan || either_nan;  // <- record whether we have seen a NaN

    // record properties if there is a mismatch
    if (!cmp_fn(actual[flat_idx], desired[flat_idx])) {
      if (first_mismatch_idx == -1) {
        first_mismatch_idx = flat_idx;
      }
      nan_mismatch_count += either_nan;
      if (either_nan && first_nan_mismatch_idx == -1) {
        first_nan_mismatch_idx = flat_idx;
      } else if (!either_nan) {
        double cur_absDiff = std::fabs(actual[flat_idx] - desired[flat_idx]);
        if (cur_absDiff > max_absDiff) {
          max_absDiff = cur_absDiff;
          max_absDiff_idx = flat_idx;
        }

        if (cur_absDiff > (max_relDiff * std::fabs(desired[flat_idx]))) {
          max_relDiff = cur_absDiff / std::fabs(desired[flat_idx]);
          max_relDiff_idx = flat_idx;
        }
      }
    }
  };
  flat_for_each(fn, idx_mapping.n_elements(), selection_mask);

  // now, let's construct the vector of details
  std::vector<MismatchDetail> details;

  if (first_mismatch_idx == -1) {
    return details;  // <- this is probably indicative of an error
  } else {
    details.emplace_back("first mismatch",
                         std::optional<int>{first_mismatch_idx});
  }

  if (max_absDiff_idx == -1) {
    details.emplace_back("Max abs diff: NaN (i.e. each mismatch involves NaN)",
                         std::nullopt);
  } else {
    details.emplace_back("Max abs diff: " + to_pretty_string(max_absDiff),
                         std::optional<int>{max_absDiff_idx});
  }

  if (max_relDiff_idx == -1) {
    details.emplace_back(
        "Max rel diff: NaN (i.e. each mismatch involves NaN or has actual=0.0)",
        std::nullopt);
  } else {
    details.emplace_back("Max rel diff: " + to_pretty_string(max_relDiff),
                         std::optional<int>{max_relDiff_idx});
  }

  if (first_nan_mismatch_idx == -1) {
    details.emplace_back(any_nan ? "all NaNs match" : "there are no NaNs",
                         std::nullopt);
  } else {
    details.emplace_back(
        "First (of " + std::to_string(nan_mismatch_count) + ") NaN mismatch",
        first_nan_mismatch_idx);
  }

  return details;
}

/// Returns a `testing::AssertionResult` instance specifying whether all pairs
/// of values from @p actual and @p desired pointers satisfy the comparison
/// operation specified by @p cmp_fn
///
/// @tparam T is either `float` or `double`
/// @tparam Layout specifies the data-layout
/// @tparam Cmp Function-like type that does the underlying comparison. See the
///     description of the @p cmp_fn function for more details
///
/// @param actual,desired The pointers being compared
/// @param idx_mapping Specifies information for treating the pointers as
///     contiguous multi-dimensional arrays. It maps between multi-dimensional
///     indices & pointer 1d offsets, and specifies all relevant information
///     for this mapping (i.e. extents and data layout)
/// @param selection_mask When specified, only the locations holding `true`
///     values are compared
/// @param cmp_fn "Callable" object that implements a function signature
///     equivalent to `bool fun(T actual, T desired)`. This signature is called
///     by passing pairs of values from the @p actual and @p desired pointers.
///     This should implement a member function called `describe_false` that
///     returns a `std::string`
template <typename T, DataLayout Layout, typename Cmp>
testing::AssertionResult compare_ptrs_(const T* actual, const T* desired,
                                       const IdxMapping<Layout>& idx_mapping,
                                       const bool* selection_mask, Cmp cmp_fn) {
  GR_INTERNAL_REQUIRE(actual != nullptr && desired != nullptr,
                      "it's illegal to compare nullptr");
  // Part 1: perform the comparison (this is as fast as possible)
  const int n_elements = idx_mapping.n_elements();
  int mismatch_num = 0;
  int n_comparisons = 0;
  auto loop_callback = [=, &mismatch_num, &n_comparisons](int flat_idx) {
    n_comparisons++;
    mismatch_num += !cmp_fn(actual[flat_idx], desired[flat_idx]);
  };
  flat_for_each(loop_callback, n_elements, selection_mask);

  if (mismatch_num == 0) {
    return testing::AssertionSuccess();
  }

  // Part 2: build the failure result and construct the detailed error message
  // -> it's ok if this isn't extremely optimized. This logic shouldn't come up
  //    very frequently
  testing::AssertionResult out = testing::AssertionFailure();

  out << '\n'
      << "arrays are " << cmp_fn.describe_false() << '\n'
      << "index mapping: " << testing::PrintToString(idx_mapping) << '\n';
  out << "Mismatched elements: " << mismatch_num << " (" << n_comparisons
      << " were compared";
  if (n_comparisons != n_elements) {
    out << ", " << n_elements - n_comparisons << "ignored from masking";
  }
  out << ")\n";

  std::vector<MismatchDetail> detail_vec =
      collect_details_(actual, desired, idx_mapping, selection_mask, cmp_fn);
  if (detail_vec.empty()) {
    GR_INTERNAL_ERROR("something went wrong with finding mismatch details");
  }

  // now let's append the interesting mismatch details
  for (const MismatchDetail& detail : detail_vec) {
    if (!detail.flat_idx.has_value()) {
      out << detail.description << '\n';
      continue;
    }
    int flat_idx = detail.flat_idx.value();

    out << detail.description << ", ";
    // write the index
    int idx_components[IdxMapping<Layout>::MAX_RANK];
    idx_mapping.offset_to_md_idx(flat_idx, idx_components);
    out << "idx: {";
    for (int i = 0; i < idx_mapping.rank(); i++) {
      out << idx_components[i];
      out << ((i + 1) < idx_mapping.rank() ? ',' : '}');
    }
    out << ", ";

    // write the actual and description value
    out << "actual = " << to_pretty_string(actual[flat_idx]) << ", "
        << "desired = " << to_pretty_string(desired[flat_idx]) << '\n';
  }

  // print out final summary details
  bool has_mask = selection_mask != nullptr;
  out << "Flattened Ptr Details"
      << (has_mask ? " (selection mask is ignored):\n" : ":\n")
      << "  actual:  " << ptr_to_string(actual, idx_mapping) << '\n'
      << "  desired: " << ptr_to_string(desired, idx_mapping);
  return out;
}

}  // anonymous namespace

testing::AssertionResult compare_(CmpPack pack) {
  // this function launches the appropriate specialization of compare_ptrs_
  // -> there are 3 template parameters to consider
  // -> (see the docstring in the header file for a little more context)

  // load either (f32_actual, f32_desired) OR (f64_actual, f64_desired)
  const float *f32_actual, *f32_desired;
  const double *f64_actual, *f64_desired;
  bool use_f32 =
      std::holds_alternative<PtrPair<float>>(pack.actual_desired_pair);
  if (use_f32) {
    f32_actual = std::get<PtrPair<float>>(pack.actual_desired_pair).first;
    f32_desired = std::get<PtrPair<float>>(pack.actual_desired_pair).second;
  } else {
    f64_actual = std::get<PtrPair<double>>(pack.actual_desired_pair).first;
    f64_desired = std::get<PtrPair<double>>(pack.actual_desired_pair).second;
  }

  // Either idx_map_L OR idx_map_R will not be a nullptr
  const IdxMapping<DataLayout::LEFT>* idx_map_L =
      std::get_if<IdxMapping<DataLayout::LEFT>>(&pack.idx_mapping);
  const IdxMapping<DataLayout::RIGHT>* idx_map_R =
      std::get_if<IdxMapping<DataLayout::RIGHT>>(&pack.idx_mapping);

  // dispatcher_ is a "generic lambda"
  // -> it acts as if the type of cmp_fn is a template parameter.
  // -> when we pass it to std::visit, the cmp_fn argument is a copy of the
  //    alternative currently held by the `CmpPack::cmp_fn` variant
  auto dispatcher_ = [&](auto cmp_fn) -> testing::AssertionResult {
    const bool* smask = pack.selection_mask;
    if (use_f32 && idx_map_L != nullptr) {
      return compare_ptrs_(f32_actual, f32_desired, *idx_map_L, smask, cmp_fn);
    } else if (use_f32 && idx_map_R != nullptr) {
      return compare_ptrs_(f32_actual, f32_desired, *idx_map_R, smask, cmp_fn);
    } else if (idx_map_L != nullptr) {
      return compare_ptrs_(f64_actual, f64_desired, *idx_map_L, smask, cmp_fn);
    } else if (idx_map_R != nullptr) {
      return compare_ptrs_(f64_actual, f64_desired, *idx_map_R, smask, cmp_fn);
    } else {
      GR_INTERNAL_ERROR("should be unreachable");
    }
  };
  return std::visit(dispatcher_, pack.cmp_fn);
}

}  // namespace grtest::arraycmp_detail