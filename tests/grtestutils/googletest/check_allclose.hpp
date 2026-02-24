#ifndef UTEST_HELPERS_HPP
#define UTEST_HELPERS_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <gtest/gtest.h>

/// equivalent of %g (this is extremely crude!)
std::string pretty_format_(double val) {
  char buf[30];
  snprintf(buf, 30, "%g", val);
  return std::string(buf);
}

/// formats a std::vector as a string
///
/// @note
/// This is highly inefficient, partially because it consists of code written
/// from before we adopted googletest
inline std::string vec_to_string(const std::vector<double>& vec) {
  std::string out = "{";

  std::size_t len = vec.size();

  std::size_t pause_start;
  std::size_t pause_stop;

  if (len > 30){
    pause_start = 3;
    pause_stop = len - 3;
  } else {
    pause_start = len *2;
    pause_stop = pause_start;
  }

  for (std::size_t i = 0; i < len; i++) {
    if ((i > pause_start) && (i < pause_stop)) { continue; }

    if (i == pause_stop) {
      out += ", ... ";
    } else if (i != 0) {
      out += ", ";
    }

    const int BUF_SIZE = 30;
    char buf[BUF_SIZE];
    snprintf(buf, BUF_SIZE, "%g", vec[i]);
    out += buf;
  }
  return out + "}";
}

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

  std::string actual_vec_str = vec_to_string(actual);
  std::string ref_vec_str = vec_to_string(desired);

  return testing::AssertionFailure()
    << "\narrays are unequal for the tolerance: "
       << "rtol = " << pretty_format_(rtol) << ", "
       << "atol = " << pretty_format_(atol) << '\n'
    << err_msg << '\n' // custom error message
    << "Mismatched elements: " << num_mismatches << " / " << actual.size()
       << '\n'
    << "Max absolute difference: " << pretty_format_(max_absDiff) << ", "
       << "ind = " << max_absDiff_ind << ", "
       << "actual = " << pretty_format_(actual[max_absDiff_ind]) << ", "
       << "reference = " << pretty_format_(desired[max_absDiff_ind]) << '\n'
    << "Max relative difference: " << pretty_format_(max_relDiff) << ", "
       << "ind = " << max_absDiff_ind << ", "
       << "actual = " << pretty_format_(actual[max_relDiff_ind]) << ", "
       << "desired = " << pretty_format_(desired[max_relDiff_ind]) << '\n'
    << "actual:  " << actual_vec_str << '\n'
    << "desired: " << ref_vec_str << '\n';
}

#endif /* UTEST_HELPERS */
