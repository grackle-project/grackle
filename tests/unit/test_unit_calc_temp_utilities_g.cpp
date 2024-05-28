#include "gtest/gtest.h"

extern "C" {
    #include "grackle_macros.h"
}

// TODO: FORTRAN_NAME for unit tests not needed
extern "C" void FORTRAN_NAME(calc_parameter_value_slopes)(double* dclPar, long long* clGridRank, 
                                                          long long* clGridDim, double* clPar1, 
                                                          double* clPar2, double* clPar3);

TEST(CalcTempUtilitiesTest, CalcParameterValueSlopes) {

    long long clGridRank = 3;
    long long clGridDim[clGridRank] = {3, 2, 3}; 
    double clPar1[clGridDim[0]] = {1., 2., 3.};
    double clPar2[clGridDim[1]] = {4., 5.};
    double clPar3[clGridDim[2]] = {6., 7., 8.};
    double dclPar[clGridRank]; 

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim, clPar1, clPar2, clPar3);

    ASSERT_EQ(dclPar[0], 1.);
    ASSERT_EQ(dclPar[1], 1.);
    ASSERT_EQ(dclPar[2], 1.);

}
