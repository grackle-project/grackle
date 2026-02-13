#ifndef UTEST_HELPERS_HPP
#define UTEST_HELPERS_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "grtestutils/utils.hpp"

/// this compares 2 std::vectors
///
/// This draws a lot of inspiration from numpy.testing.assert_allclose
///
/// Parts of this are fairly inefficient, partially because it is adapted from
/// code written from before we adopted googletest
inline testing::AssertionResult check_allclose(
  const std::vector<double>& actual,
  const std::vector<double>& desired,
  double rtol = 0.0, double atol = 0.0,
  std::string err_msg = "")
{
  if (actual.size() != desired.size()){
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
    double cur_absDiff = fabs(actual[i]-desired[i]);

    bool isnan_actual = isnan(actual[i]);
    bool isnan_desired = isnan(desired[i]);

    if ( (cur_absDiff > (atol + rtol * fabs(desired[i]))) ||
         (isnan_actual != isnan_desired) ) {

      num_mismatches++;
      if (isnan_actual != isnan_desired){
        has_nan_mismatch = true;
        max_absDiff = NAN;
        max_absDiff_ind = i;
        max_relDiff = NAN;
        max_relDiff_ind = i;
      } else if (!has_nan_mismatch) {
        if (cur_absDiff > max_absDiff){
          max_absDiff = cur_absDiff;
          max_absDiff_ind = i;
        }

        if ( cur_absDiff > (max_relDiff * fabs(desired[i])) ) {
          max_relDiff = cur_absDiff / fabs(desired[i]);
          max_relDiff_ind = i;
        }
      }
    }
  }

  if (num_mismatches == 0) { return testing::AssertionSuccess(); }

  std::string actual_vec_str = grtest::ptr_to_string(actual.data(),
                                                     actual.size());
  std::string ref_vec_str = grtest::ptr_to_string(desired.data(),
                                                  desired.size());

  using grtest::to_pretty_string;

  return testing::AssertionFailure()
    << "\narrays are unequal for the tolerance: "
       << "rtol = " << to_pretty_string(rtol) << ", "
       << "atol = " << to_pretty_string(atol) << '\n'
    << err_msg << '\n' // custom error message
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

#endif /* UTEST_HELPERS */
