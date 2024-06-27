#include "gtest/gtest.h"

extern "C" {
    #include "grackle_macros.h"
}

// TODO: FORTRAN_NAME for unit tests not needed
extern "C" void FORTRAN_NAME(calc_parameter_value_slopes)(double* dclPar, long long* clGridRank, long long* clGridDim, 
                                                          double* clPar1, double* clPar2, double* clPar3, 
                                                          double* clPar4, double* clPar5);


TEST(CalcTempUtilitiesTest, CalcParameterValueSlopes) {

    // clGridRank = 1
    long long clGridRank = 1;
    double* dclPar = new double[clGridRank];
    long long* clGridDim = new long long[clGridRank]{3};
    double* clPar1 = new double[clGridDim[0]]{1., 2., 3.};

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim,
                                              clPar1, nullptr, nullptr,
                                              nullptr, nullptr);

    ASSERT_NEAR(dclPar[0], 1., 1e-15);

    delete [] dclPar;
    delete [] clGridDim;

    // clGridRank = 2
    clGridRank = 2;
    dclPar = new double[clGridRank];
    clGridDim = new long long[clGridRank]{3, 2};
    double* clPar2 = new double[clGridDim[1]]{4., 4.5};

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim,
                                              clPar1, clPar2, nullptr,
                                              nullptr, nullptr);

    ASSERT_NEAR(dclPar[0], 1., 1e-15);
    ASSERT_NEAR(dclPar[1], 0.5, 1e-15);


    delete [] dclPar;
    delete [] clGridDim;


    // clGridRank = 3
    clGridRank = 3;
    dclPar = new double[clGridRank];
    clGridDim = new long long[clGridRank]{3, 2, 4};
    double* clPar3 = new double[clGridDim[2]]{1., 1.25, 1.5, 1.75};

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim,
                                              clPar1, clPar2, clPar3,
                                              nullptr, nullptr);

    ASSERT_NEAR(dclPar[0], 1., 1e-15);
    ASSERT_NEAR(dclPar[1], 0.5, 1e-15);
    ASSERT_NEAR(dclPar[2], 0.25, 1e-15);


    delete [] dclPar;
    delete [] clGridDim;


    // clGridRank = 4
    clGridRank = 4;
    dclPar = new double[clGridRank];
    clGridDim = new long long[clGridRank]{3, 2, 4, 3};
    double* clPar4 = new double[clGridDim[3]]{10., 15., 20.};

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim,
                                              clPar1, clPar2, clPar3,
                                              clPar4, nullptr);

    ASSERT_NEAR(dclPar[0], 1., 1e-15);
    ASSERT_NEAR(dclPar[1], 0.5, 1e-15);
    ASSERT_NEAR(dclPar[2], 0.25, 1e-15);
    ASSERT_NEAR(dclPar[3], 5., 1e-15);


    delete [] dclPar;
    delete [] clGridDim;

    // clGridRank = 5
    clGridRank = 5;
    dclPar = new double[clGridRank];
    clGridDim = new long long[clGridRank]{3, 2, 4, 3, 3};
    double* clPar5 = new double[clGridDim[4]]{1.4, 1.8, 2.2};

    FORTRAN_NAME(calc_parameter_value_slopes)(dclPar, &clGridRank, clGridDim,
                                              clPar1, clPar2, clPar3,
                                              clPar4, clPar5);

    ASSERT_NEAR(dclPar[0], 1., 1e-15);
    ASSERT_NEAR(dclPar[1], 0.5, 1e-15);
    ASSERT_NEAR(dclPar[2], 0.25, 1e-15);
    ASSERT_NEAR(dclPar[3], 5., 1e-15);
    ASSERT_NEAR(dclPar[4], 0.4, 1e-15);


    delete [] dclPar;
    delete [] clGridDim;

    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;

}
