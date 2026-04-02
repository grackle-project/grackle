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

#include <cmath>
#include <string>
#include <vector>
#include <gtest/gtest.h>

testing::AssertionResult check_allclose(const std::vector<double>& actual,
                                        const std::vector<double>& desired,
                                        double rtol, double atol) {
  if (actual.size() != desired.size()) {
    return testing::AssertionFailure()
           << "the compared arrays have different lengths";
  }

  std::size_t num_mismatches = 0;
  double max_absDiff = 0.0;
  std::size_t max_absDiff_ind = 0;
  double max_relDiff = 0.0;
  std::size_t max_relDiff_ind = 0;
  bool has_nan_mismatch = false;

  for (std::size_t i = 0; i < actual.size(); i++) {
    double cur_absDiff = std::fabs(actual[i] - desired[i]);

    bool isnan_actual = std::isnan(actual[i]);
    bool isnan_desired = std::isnan(desired[i]);

    if ((cur_absDiff > (atol + rtol * std::fabs(desired[i]))) ||
        (isnan_actual != isnan_desired)) {
      num_mismatches++;
      if (isnan_actual != isnan_desired) {
        has_nan_mismatch = true;
        max_absDiff = NAN;
        max_absDiff_ind = i;
        max_relDiff = NAN;
        max_relDiff_ind = i;
      } else if (!has_nan_mismatch) {
        if (cur_absDiff > max_absDiff) {
          max_absDiff = cur_absDiff;
          max_absDiff_ind = i;
        }

        if (cur_absDiff > (max_relDiff * std::fabs(desired[i]))) {
          max_relDiff = cur_absDiff / std::fabs(desired[i]);
          max_relDiff_ind = i;
        }
      }
    }
  }

  if (num_mismatches == 0) {
    return testing::AssertionSuccess();
  }

  std::string actual_vec_str =
      grtest::ptr_to_string(actual.data(), actual.size());
  std::string ref_vec_str =
      grtest::ptr_to_string(desired.data(), desired.size());

  using grtest::to_pretty_string;

  return testing::AssertionFailure()
         << "\narrays are unequal for the tolerance: "
         << "rtol = " << to_pretty_string(rtol) << ", "
         << "atol = " << to_pretty_string(atol) << '\n'
         << "Mismatched elements: " << num_mismatches << " / " << actual.size()
         << '\n'
         << "Max absolute difference: " << to_pretty_string(max_absDiff) << ", "
         << "ind = " << max_absDiff_ind << ", "
         << "actual = " << to_pretty_string(actual[max_absDiff_ind]) << ", "
         << "reference = " << to_pretty_string(desired[max_absDiff_ind]) << '\n'
         << "Max relative difference: " << to_pretty_string(max_relDiff) << ", "
         << "ind = " << max_absDiff_ind << ", "
         << "actual = " << to_pretty_string(actual[max_relDiff_ind]) << ", "
         << "desired = " << to_pretty_string(desired[max_relDiff_ind]) << '\n'
         << "actual:  " << actual_vec_str << '\n'
         << "desired: " << ref_vec_str << '\n';
}