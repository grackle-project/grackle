// See LICENSE file for license and copyright information

/// @file interpolate.hpp
/// @brief Defines C++ versions of the interpolate functions
///
/// There is a huge potential for optimization in these functions. For example:
///  - we could take advantage of the regular spacing between parameter values
///    and precompute 1/dgridPar1, 1/dgridPar2, 1/dgridPar3, ... ahead of time.
///    This would allow us to avoid a lot of expensive division operations
///  - we could restructure the ordering of the tables (maybe when we read them
///    in?) so that we can further reduce the number of calls to the log
///    function, which is generally very slow.

#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include <cmath> // log
#include "fortran_func_decls.h" // gr_i64
#include "utils-cpp.hpp" // GRIMPL_RESTRICT

namespace grackle::impl{

/// helper function that determines the 1-indexed interpolation index
///
/// This assumes that parameter is evenly spaced on the grid
static inline gr_i64 get_index_(double input, gr_i64 parLen,
                                  const double *gridPar, double dgridPar)
{
  gr_i64 index = (gr_i64)((input-gridPar[0])/dgridPar)+1;
  gr_i64 index_with_floor = (index > 1) ? index : 1;
  return ((parLen - 1) < index_with_floor) ? (parLen - 1) : index_with_floor;
};

/// helper function to interpolate along a single dimension
///
/// @note
/// we could probably make this a 3 argument function where the second & third
/// arguments each expect a pointer to a pair of values. The main reason I have
/// avoided this right now is that it may interfere with vectorization
///
/// @note
/// We may want to use some compiler specific extensions to force inlining of
/// this function (OR define it as a macro)
static inline double interp_(double x,
                             double xref0, double xref1,
                             double yref0, double yref1)
{
  double slope = (yref1 - yref0) / (xref1 - xref0);
  return (x - xref0) * slope + yref0;
}

void interpolate_1d_g(double input1,
                      const gr_i64 * GRIMPL_RESTRICT gridDim, // 1 elements
                      const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                      gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                      double * GRIMPL_RESTRICT value)
{
  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);

  // interpolate over parameter 1
  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   dataField[index1-1], dataField[index1]);
}

void interpolate_2d_g(double input1, double input2,
                      const gr_i64 * GRIMPL_RESTRICT gridDim, // 2 elements
                      const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                      const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
                      gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                      double * GRIMPL_RESTRICT value)
{
  double value2[2];

  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index2 = get_index_(input2, gridDim[1], gridPar2, dgridPar2);

  for (gr_i64 q=0; q < 2; q++) {

    // interpolate over parameter 2
    gr_i64 int_index = (q+index1-1) * gridDim[1] + index2;

    value2[q] = interp_(input2, gridPar2[index2-1], gridPar2[index2],
                        dataField[int_index-1], dataField[int_index]);
  }

  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   value2[0], value2[1]);
}

void interpolate_3d_g(double input1, double input2, double input3,
                      const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
                      const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                      const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
                      const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
                      gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                      double * GRIMPL_RESTRICT value)
{
  double value3[2], value2[2];

  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index2 = get_index_(input2, gridDim[1], gridPar2, dgridPar2);
  const gr_i64 index3 = get_index_(input3, gridDim[2], gridPar3, dgridPar3);

  for (gr_i64 q = 0; q < 2; q++) {
    for (gr_i64 w = 0; w < 2; w++) {

      // interpolate over parameter 3
      gr_i64 int_index = ((q+index1-1) * gridDim[1] +
                            (w+index2-1)) * gridDim[2] + index3;

      value3[w] = interp_(input3, gridPar3[index3-1], gridPar3[index3],
                          dataField[int_index-1], dataField[int_index]);
    }

    // interpolate over parameter 2
    value2[q] = interp_(input2, gridPar2[index2-1], gridPar2[index2],
                        value3[0], value3[1]);
  }

  // interpolate over parameter 1
  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   value2[0], value2[1]);
}

void interpolate_4d_g(double input1, double input2, double input3,
                      double input4,
                      const gr_i64 * GRIMPL_RESTRICT gridDim, // 4 elements
                      const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                      const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
                      const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
                      const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
                      gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                      double * GRIMPL_RESTRICT value)
{

  double value4[2], value3[2], value2[2];

  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index2 = get_index_(input2, gridDim[1], gridPar2, dgridPar2);
  const gr_i64 index3 = get_index_(input3, gridDim[2], gridPar3, dgridPar3);
  const gr_i64 index4 = get_index_(input4, gridDim[3], gridPar4, dgridPar4);

  for (gr_i64 q = 0; q < 2; q++) {
    for (gr_i64 w = 0; w < 2; w++) {
      for (gr_i64 e = 0; e < 2; e++) {

        // interpolate over parameter 4
        gr_i64 int_index = (((q+index1-1) * gridDim[1] +
                               (w+index2-1)) * gridDim[2] +
                              (e+index3-1)) * gridDim[3] + index4;

        value4[e] = interp_(input4, gridPar4[index4-1], gridPar4[index4],
                            dataField[int_index-1], dataField[int_index]);
      }

      // interpolate over parameter 3
      value3[w] = interp_(input3, gridPar3[index3-1], gridPar3[index3],
                          value4[0], value4[1]);
    }

    // interpolate over parameter 2
    value2[q] = interp_(input2, gridPar2[index2-1], gridPar2[index2],
                        value3[0], value3[1]);
  }

  // interpolate over parameter 1
  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   value2[0], value2[1]);
}

void interpolate_5d_g(double input1, double input2, double input3,
                      double input4, double input5,
                      const gr_i64 * GRIMPL_RESTRICT gridDim, // 5 elements
                      const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                      const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
                      const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
                      const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
                      const double * GRIMPL_RESTRICT gridPar5, double dgridPar5,
                      gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                      double * GRIMPL_RESTRICT value)
{
  double value5[2], value4[2], value3[2], value2[2];

  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index2 = get_index_(input2, gridDim[1], gridPar2, dgridPar2);
  const gr_i64 index3 = get_index_(input3, gridDim[2], gridPar3, dgridPar3);
#define INDEX_4_BISECTION
#ifdef INDEX_4_BISECTION
  // get index 4 with bisection, since not evenly spaced
  gr_i64 index4;
  if (input4 <= gridPar4[0]) {
    index4 = 1;
  } else if (input4 >= gridPar4[gridDim[3]-2]) { // -2 isn't a typo
    index4 = gridDim[3] - 1;
  } else {
    index4 = 1;
    gr_i64 highPt = gridDim[3];
    while ((highPt - index4) > 1) {
      gr_i64 midPt = (gr_i64)((highPt + index4) / 2);
      if (input4 >= gridPar4[midPt-1]) {
        index4 = midPt;
      } else {
        highPt = midPt;
      }
    }
  }
#else
  gr_i64 index4 = get_index_(input4, gridDim[3], gridPar4, dgridPar4);
#endif /* INDEX_4_BISECTION */
  const gr_i64 index5 = get_index_(input5, gridDim[4], gridPar5, dgridPar5);

  for (gr_i64 q = 0; q < 2; q++) {
    for (gr_i64 w = 0; w < 2; w++) {
      for (gr_i64 e = 0; e < 2; e++) {
        for (gr_i64 r = 0; r < 2; r++) {

          // interpolate over parameter 5
          gr_i64 int_index = ((((q+index1-1) * gridDim[1] +
                                  (w+index2-1)) * gridDim[2] + (e+index3-1)) *
                                gridDim[3] + (r+index4-1)) * gridDim[4] +
                                index5;

          value5[r] = interp_(input5, gridPar5[index5-1], gridPar5[index5],
                              dataField[int_index-1], dataField[int_index]);
        }

        // interpolate over parameter 4
        value4[e] = interp_(input4, gridPar4[index4-1], gridPar4[index4],
                            value5[0], value5[1]);
      }

      // interpolate over parameter 3
      value3[w] = interp_(input3, gridPar3[index3-1], gridPar3[index3],
                          value4[0], value4[1]);
    }

    // interpolate over parameter 2
    value2[q] = interp_(input2, gridPar2[index2-1], gridPar2[index2],
                        value3[0], value3[1]);
  }

  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   value2[0], value2[1]);
}


// OTHER FUNCTIONS
// Interpolation in 2 dimensions but with a 3D grid.
// This is used for interpolating from just the last
// slice in the datacube before the redshift where
// the UV background turns on.
static inline double interpolate_2Df3D_g(double input1, double input3,
                                  const gr_i64* GRIMPL_RESTRICT gridDim, // 3 elements
                                  const double* GRIMPL_RESTRICT gridPar1, double dgridPar1,
                                  gr_i64 index2,
                                  const double* GRIMPL_RESTRICT gridPar3, double dgridPar3,
                                  gr_i64 dataSize,
                                  const double* dataField)
{
  double value3[2];

  // Calculate interpolation indices
  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index3 = get_index_(input3, gridDim[2], gridPar3, dgridPar3);

  for (gr_i64 q = 0; q < 2; q++) { // interpolate over parameter 3

    gr_i64 int_index = ((q+index1-1) * gridDim[1] + (index2-1)) *
      gridDim[2] + index3;

    value3[q] = interp_(input3, gridPar3[index3-1], gridPar3[index1],
                        dataField[int_index-1], dataField[int_index]);
  }

  // interpolate over parameter 1
  return interp_(input1, gridPar1[index1-1], gridPar1[index1],
                 value3[0], value3[1]);
}


void interpolate_3dz_g(double input1, double input2, double input3,
                       const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
                       const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
                       const double * GRIMPL_RESTRICT gridPar2, gr_i64 index2,
                       const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
                       gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
                       gr_i64 end_int, double * GRIMPL_RESTRICT value)
{
  if (end_int == 1) {
    *value = interpolate_2Df3D_g(input1, input3,
                                 gridDim,
                                 gridPar1, dgridPar1,
                                 index2,
                                 gridPar3, dgridPar3,
                                 dataSize, dataField);
    return;
  }

  double value3[2], value2[2];

  // Calculate interpolation indices
  const gr_i64 index1 = get_index_(input1, gridDim[0], gridPar1, dgridPar1);
  const gr_i64 index3 = get_index_(input3, gridDim[2], gridPar3, dgridPar3);


  // it turns out that precomputing the following 2 variables reduces runtime
  // appreciably (because the C compiler can't automatically hoist these
  // calculations out of the loop)
  const double par2_slope_denom = log((1+gridPar2[index2]) /
                                      (1+gridPar2[index2-1]));
  const double par2_offset_from_grid = log((1+input2)/(1+gridPar2[index2-1]));


  // preliminary testing on gcc 9.4 suggests that unrolling the outer loop
  // could speed this function up by ~10%
  for (gr_i64 q = 0; q < 2; q++) {
    for (gr_i64 w = 0; w < 2; w++) {

      // interpolate over parameter 3
      gr_i64 int_index = ((q+index1-1) * gridDim[1] +
                            (w+index2-1)) * gridDim[2] + index3;

      value3[w] = interp_(input3, gridPar3[index3-1], gridPar3[index3],
                          dataField[int_index-1], dataField[int_index]);
    }

    // interpolate over parameter 2
    double slope = (value3[1] - value3[0]) / par2_slope_denom;
    value2[q] = par2_offset_from_grid * slope + value3[0];
  }

  // interpolate over parameter 1
  *value = interp_(input1, gridPar1[index1-1], gridPar1[index1],
                   value2[0], value2[1]);
}

} // namespace grackle::impl

#endif /* INTERPOLATE_HPP */
