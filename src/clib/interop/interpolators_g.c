#include <stdint.h>

#include <math.h> // log

//#include "../grackle_macros"
// we steal the following from  "../grackle_macros.h"
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

typedef int64_t gr_int64;
typedef int64_t gr_dint; // equivalent of int(val, DIKIND)

void interpolate_1d_g(double input1,
                      const gr_int64* gridDim, // 1 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p;

  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));

  // interpolate over parameter 1
  double slope = (dataField[index1] - dataField[index1-1]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + dataField[index1-1];
}

void interpolate_2d_g(double input1, double input2,
                      const gr_int64* gridDim, // 2 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p; const double dgridPar2 = *dgridPar2_p;

  double value2[2];

  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index2 = min(gridDim[1]-1,
                        max(1, (gr_dint)((input2-gridPar2[0])/dgridPar2)+1));

  for (gr_int64 q=1; q <= 2; q++) {

    // interpolate over parameter 2
    gr_int64 int_index = (q+index1-2) * gridDim[1] + index2;

    double slope = (dataField[int_index] - dataField[int_index-1]) /
      (gridPar2[index2] - gridPar2[index2-1]);

    value2[q-1] = (input2 - gridPar2[index2-1]) * slope +
      dataField[int_index-1];
  }

  double slope = (value2[1] - value2[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + value2[0];
}

void interpolate_3d_g(double input1, double input2, double input3,
                      const gr_int64* gridDim, // 3 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p; const double dgridPar2 = *dgridPar2_p;
  const double dgridPar3 = *dgridPar3_p;

  double value3[2], value2[2];

  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index2 = min(gridDim[1]-1,
                        max(1, (gr_dint)((input2-gridPar2[0])/dgridPar2)+1));
  gr_int64 index3 = min(gridDim[2]-1,
                        max(1, (gr_dint)((input3-gridPar3[0])/dgridPar3)+1));

  for (gr_int64 q=1; q <= 2; q++) {
    for (gr_int64 w=1; w <= 2; w++) {

      // interpolate over parameter 3
      gr_int64 int_index = ((q+index1-2) * gridDim[1] +
                            (w+index2-2)) * gridDim[2] + index3;

      double slope = (dataField[int_index] - dataField[int_index-1]) /
        (gridPar3[index3] - gridPar3[index3-1]);

      value3[w-1] = (input3 - gridPar3[index3-1]) * slope +
          dataField[int_index-1];
    }

    // interpolate over parameter 2
    double slope = (value3[1] - value3[0]) /
      (gridPar2[index2] - gridPar2[index2-1]);

    value2[q-1] = (input2 - gridPar2[index2-1]) * slope + value3[0];
  }

  // interpolate over parameter 1
  double slope = (value2[1] - value2[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + value2[0];
}

void interpolate_4d_g(double input1, double input2, double input3,
                      double input4,
                      const gr_int64* gridDim, // 4 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const double* gridPar4, const double* dgridPar4_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p; const double dgridPar2 = *dgridPar2_p;
  const double dgridPar3 = *dgridPar3_p; const double dgridPar4 = *dgridPar4_p;

  double value4[2], value3[2], value2[2];

  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index2 = min(gridDim[1]-1,
                        max(1, (gr_dint)((input2-gridPar2[0])/dgridPar2)+1));
  gr_int64 index3 = min(gridDim[2]-1,
                        max(1, (gr_dint)((input3-gridPar3[0])/dgridPar3)+1));
  gr_int64 index4 = min(gridDim[3]-1,
                        max(1, (gr_dint)((input4-gridPar4[0])/dgridPar4)+1));

  for (gr_int64 q=1; q <= 2; q++) {
    for (gr_int64 w=1; w <= 2; w++) {
      for (gr_int64 e=1; e <= 2; e++) {

        // interpolate over parameter 4

        gr_int64 int_index = (((q+index1-2) * gridDim[1] +
                               (w+index2-2)) * gridDim[2] +
                              (e+index3-2)) * gridDim[3] + index4;

        double slope = (dataField[int_index] - dataField[int_index-1]) /
          (gridPar4[index4] - gridPar4[index4-1]);

        value4[e-1] = (input4 - gridPar4[index4-1]) * slope +
          dataField[int_index-1];
      }

      // interpolate over parameter 3
      double slope = (value4[1] - value4[0]) /
        (gridPar3[index3] - gridPar3[index3-1]);

      value3[w-1] = (input3 - gridPar3[index3-1]) * slope + value4[0];
    }

    // interpolate over parameter 2
    double slope = (value3[1] - value3[0]) /
      (gridPar2[index2] - gridPar2[index2-1]);

    value2[q-1] = (input2 - gridPar2[index2-1]) * slope + value3[0];
  }

  // interpolate over parameter 1
  double slope = (value2[1] - value2[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + value2[0];
}

void interpolate_5d_g(double input1, double input2, double input3,
                      double input4, double input5,
                      const gr_int64* gridDim, // 5 elements
                      const double* gridPar1, const double* dgridPar1_p,
                      const double* gridPar2, const double* dgridPar2_p,
                      const double* gridPar3, const double* dgridPar3_p,
                      const double* gridPar4, const double* dgridPar4_p,
                      const double* gridPar5, const double* dgridPar5_p,
                      const gr_int64* dataSize, const double* dataField,
                      double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p; const double dgridPar2 = *dgridPar2_p;
  const double dgridPar3 = *dgridPar3_p; const double dgridPar4 = *dgridPar4_p;
  const double dgridPar5 = *dgridPar5_p;
  
  double value5[2], value4[2], value3[2], value2[2];

  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index2 = min(gridDim[1]-1,
                        max(1, (gr_dint)((input2-gridPar2[0])/dgridPar2)+1));
  gr_int64 index3 = min(gridDim[2]-1,
                        max(1, (gr_dint)((input3-gridPar3[0])/dgridPar3)+1));
#define INDEX_4_BISECTION
#ifdef INDEX_4_BISECTION
  // get index 4 with bisection, since not evenly spaced
  gr_int64 index4;
  if (input4 <= gridPar4[0]) {
    index4 = 1;
  } else if (input4 >= gridPar4[gridDim[3]-2]) { // -2 isn't a typo
    index4 = gridDim[3] - 1;
  } else {
    index4 = 1;
    gr_int64 highPt = gridDim[3];
    while ((highPt - index4) > 1) {
      gr_int64 midPt = (gr_dint)((highPt + index4) / 2);
      if (input4 >= gridPar4[midPt-1]) {
        index4 = midPt;
      } else {
        highPt = midPt;
      }
    }
  }
#else
  gr_int64 index4 = min(gridDim[3]-1,
                        max(1, (gr_dint)((input4-gridPar4[0])/dgridPar4)+1));
#endif /* INDEX_4_BISECTION */
  gr_int64 index5 = min(gridDim[4]-1,
                        max(1, (gr_dint)((input5-gridPar5[0])/dgridPar5)+1));


  for (gr_int64 q=1; q <= 2; q++) {
    for (gr_int64 w=1; w <= 2; w++) {
      for (gr_int64 e=1; e <= 2; e++) {
        for (gr_int64 r=1; r <= 2; r++) {
          // interpolate over parameter 5
          gr_int64 int_index = ((((q+index1-2) * gridDim[1] + 
                                  (w+index2-2)) * gridDim[2] + (e+index3-2)) * 
                                gridDim[3] + (r+index4-2)) * gridDim[4] +
                                index5;
          double slope = (dataField[int_index] - dataField[int_index-1]) /
                         (gridPar5[index5] - gridPar5[index5-1]);

          value5[r-1] = (input5 - gridPar5[index5-1]) * slope +
            dataField[int_index-1];
        }

        // interpolate over parameter 4
        double slope = (value5[1] - value5[0]) /
          (gridPar4[index4] - gridPar4[index4-1]);

        value4[e-1] = (input4 - gridPar4[index4-1]) * slope + value5[0];
      }

      // interpolate over parameter 3
      double slope = (value4[1] - value4[0]) /
        (gridPar3[index3] - gridPar3[index3-1]);

      value3[w-1] = (input3 - gridPar3[index3-1]) * slope + value4[0];
    }

    // interpolate over parameter 2
    double slope = (value3[1] - value3[0]) /
      (gridPar2[index2] - gridPar2[index2-1]);

    value2[q-1] = (input2 - gridPar2[index2-1]) * slope + value3[0];
  }

  double slope = (value2[1] - value2[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + value2[0];
}


// OTHER FUNCTIONS
// Interpolation in 2 dimensions but with a 3D grid.
// This is used for interpolating from just the last
// slice in the datacube before the redshift where
// the UV background turns on.
static double interpolate_2Df3D_g(double input1, double input3,
                                  const gr_int64* gridDim, // 2 elements
                                  const double* gridPar1, double dgridPar1,
                                  gr_int64 index2,
                                  const double* gridPar3, double dgridPar3,
                                  const gr_int64 dataSize,
                                  const double* dataField)
{
  double value3[2];

  // Calculate interpolation indices
  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index3 = min(gridDim[2]-1,
                        max(1, (gr_dint)((input3-gridPar3[0])/dgridPar3)+1));

  for (gr_int64 q=1; q <= 2; q++) { // interpolate over parameter 3

    // NOTE: if there's a bug, it's probably here!
    gr_int64 int_index = ((q+index1-2) * gridDim[1] + (index2-1)) *
      gridDim[2] + index3;

    double slope = (dataField[int_index] - dataField[int_index-1]) /
      (gridPar3[index3] - gridPar3[index3-1]);

    value3[q-1] = (input3 - gridPar3[index3-1]) * slope +
      dataField[int_index-1];
  }

  // interpolate over parameter 1
  double slope = (value3[1] - value3[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  return (input1 - gridPar1[index1-1]) * slope + value3[0];
}


void interpolate_3dz_g(double input1, double input2, double input3,
                       const gr_int64* gridDim, // 3 elements
                       const double* gridPar1, const double* dgridPar1_p,
                       const double* gridPar2, const gr_int64* index2_p,
                       const double* gridPar3, const double* dgridPar3_p,
                       const gr_int64* dataSize, const double* dataField,
                       const gr_int64* end_int_p, double* value)
{
  // store some arguments that were passed by pointer
  const double dgridPar1 = *dgridPar1_p;
  const gr_int64 index2 = *index2_p;
  const double dgridPar3 = *dgridPar3_p;

  const gr_int64 end_int = *end_int_p;

  if (end_int == 1) {
    *value = interpolate_2Df3D_g(input1, input3,
                                 gridDim,
                                 gridPar1, dgridPar1,
                                 index2,
                                 gridPar3, dgridPar3,
                                 *dataSize, dataField);
    return;
  }

  double value3[2];
  double value2[2];

  // Calculate interpolation indices
  gr_int64 index1 = min(gridDim[0]-1,
                        max(1, (gr_dint)((input1-gridPar1[0])/dgridPar1)+1));
  gr_int64 index3 = min(gridDim[2]-1,
                        max(1, (gr_dint)((input3-gridPar3[0])/dgridPar3)+1));

  for (gr_int64 q=1; q <= 2; q++) {
    for (gr_int64 w=1; w <= 2; w++) {

      // interpolate over parameter 3
      gr_int64 int_index = ((q+index1-2) * gridDim[1] +
                            (w+index2-2)) * gridDim[2] + index3;

      double slope = (dataField[int_index] - dataField[int_index-1]) /
        (gridPar3[index3] - gridPar3[index3-1]);

      value3[w-1] = (input3 - gridPar3[index3-1]) * slope +
          dataField[int_index-1];
    }

    // interpolate over parameter 2
    double slope = (value3[1] - value3[0]) /
      log((1+gridPar2[index2]) / (1+gridPar2[index2-1]));

    value2[q-1] = log((1+input2)/(1+gridPar2[index2-1])) * slope + value3[0];
  }

  // interpolate over parameter 1
  double slope = (value2[1] - value2[0]) /
    (gridPar1[index1] - gridPar1[index1-1]);

  *value = (input1 - gridPar1[index1-1]) * slope + value2[0];
}


// NOTE: we could probably abstract out the internals here
// static inline calc_slope(const double* y_pair, const double* x_pair)
// { return (y_pair[1] - y_pair[0]) / (x_pair[1] - x_pair[0]) };
/*
inline double interpolate_along_dim_(const double* y_val_pair,
                                     const double* x_grid_vals,
                                     gr_int64 lo_x_index,
                                     double x_val,
                                     double offset)
{
  double slope = (y_val_pair[1] - y_val_pair[0]) /
    (x_grid_vals[lo_x_index+1] - x_grid_vals[lo_x_index]);
  return (x_val - x_grid_vals[lo_x_index]) * slope + offset;
}
*/
