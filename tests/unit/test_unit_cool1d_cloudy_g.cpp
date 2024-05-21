#include "gtest/gtest.h"
#include <cmath>

extern "C" {
    #include "grackle_macros.h"
}

extern "C" void FORTRAN_NAME(cool_rank1_interpolation)(int* get_heat, int* icmbTfloor,
                                                       double* log10tem_i, double* log10_tCMB,
                                                       long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                       double* clPar1, double* clCooling, double* clHeating,
                                                       double* dclPar_1,
                                                       double*  edot_met_i);

extern "C" void FORTRAN_NAME(cool_rank2_interpolation)(int* get_heat, int* icmbTfloor,
                                                       double* log_n_h_i, double* log10tem_i, double* log10_tCMB,
                                                       long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                       double* clPar1, double* clPar2, double* clCooling, double* clHeating,
                                                       double* dclPar_1, double* dclPar_2,
                                                       double* edot_met_i);

extern "C" void FORTRAN_NAME(cool_rank3_interpolation)(int* get_heat, int* icmbTfloor,
                                                       double* log10tem_i, double* edot_met_i,
                                                       double* log_cool_i, double* log_cool_cmb_i,
                                                       double* log_heat_i, double* log_n_h_i,
                                                       long long* zindex, long long* clGridRank,
                                                       long long* clDataSize, long long* clGridDim,
                                                       double* clPar1, double* clPar2, double* clPar3,
                                                       double* clCooling, double* clHeating,
                                                       double* dclPar_1, double* dclPar_3,
                                                       double* log10_tCMB);

TEST(UnitCool1DCloudy, CoolRank1InterpolationTest) {

    int get_heat = 0;
    int icmbTfloor = 0;
    double log10tem_i = 3.4;
    double log10_tCMB = log10tem_i - 1.3;
    long long clGridRank = 1;
    long long clDataSize = 10;
    long long clGridDim[clGridRank] = {clDataSize};
    double clPar1[clDataSize] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
    double clCooling[clDataSize] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
    double clHeating[clDataSize] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
    double dclPar_1 = 1.0;
    double edot_met_i = 0.;

    FORTRAN_NAME(cool_rank1_interpolation)(&get_heat, &icmbTfloor,
                                           &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim, 
                                           clPar1, clCooling,
                                           clHeating, &dclPar_1,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -2511.8864315095798);

    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank1_interpolation)(&get_heat, &icmbTfloor,
                                           &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim, 
                                           clPar1, clCooling,
                                           clHeating, &dclPar_1,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -2385.9938903301631);

    get_heat = 1;
    icmbTfloor = 0;
    FORTRAN_NAME(cool_rank1_interpolation)(&get_heat, &icmbTfloor,
                                           &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim, 
                                           clPar1, clCooling,
                                           clHeating, &dclPar_1,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, 5431.3959157332338);

    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank1_interpolation)(&get_heat, &icmbTfloor,
                                           &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim, 
                                           clPar1, clCooling,
                                           clHeating, &dclPar_1,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, 5557.2884569126509);
}

TEST(UnitCool1DCloudy, CoolRank2InterpolationTest) {

    int get_heat = 0;
    int icmbTfloor = 0;
    double log_n_h_i = 2.5;
    double log10tem_i = 0.75;
    double log10_tCMB = log10tem_i - 0.5;
    long long clGridRank = 2;
    long long clDataSize1 = 5;
    long long clDataSize2 = 3;
    long long clDataSize = clDataSize1*clDataSize2;
    long long clGridDim[clGridRank] = {clDataSize1, clDataSize2};
    double clPar1[clDataSize1] = {0., 1., 2., 3., 4.,};
    double clPar2[clDataSize2] = {0.5, 1., 1.5};
    double clCooling[clDataSize] = {0., 0.3, 0.6, 1., 1.3, 1.6, 2., 2.3, 2.6,
                                    3., 3.3, 3.6, 4., 4.3, 4.6};
    double clHeating[clDataSize] = {0.5, 0.6, 0.7, 0.55, 0.65, 0.75, 0.6, 0.65, 0.7,
                                    0.65, 0.7, 0.75, 0.7, 0.75, 0.8};
    double dclPar_1 = 1.0;
    double dclPar_2 = 0.5;
    double edot_met_i = 0.;


    FORTRAN_NAME(cool_rank2_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2,
                                           clCooling, clHeating,
                                           &dclPar_1, &dclPar_2,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -446.68359215096302);

    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank2_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2,
                                           clCooling, clHeating,
                                           &dclPar_1, &dclPar_2,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -222.81147829412902);

    get_heat = 1;
    icmbTfloor = 0;
    FORTRAN_NAME(cool_rank2_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2,
                                           clCooling, clHeating,
                                           &dclPar_1, &dclPar_2,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -442.21675622945338);

    get_heat = 1;
    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank2_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2,
                                           clCooling, clHeating,
                                           &dclPar_1, &dclPar_2,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -218.34464237261938);

}

TEST(UnitCool1DCloudy, CoolRank3InterpolationTest) {
/*
    int get_heat = 2;
    int icmbTfloor = 1;
    double log10tem_i = 3.4066747435375704;
    double edot_met_i = 1.8116928740079170E-316;
    double log_cool_i = 1.8116897119877836E-316;
    double  log_cool_cmb_i = 1.8116897119877836E-316;
    double  log_heat_i = 1.8116897119877836E-316;
    double log_n_h_i = -0.12552852219224844;
    long long clGridRank = 3;
    long long clDataSize = 121394;
    long long* clGridDim = [29, 26, 161];
    long long zindex = 1;
    double* clPar1 = [4.0000000000000000, 14.849000000000000, 9.0000000000000000];
    double clCooling;
    double clHeating;
    double dclPar_1 = 0.50000000000000000;
    double dclPar_3 = 5.0000000000000003E-002;
    double log10_tCMB = 0.43616264704075602;

    FORTRAN_NAME(cool_rank3_interpolation)(get_heat, icmbTfloor,
                                           log10tem_i, edot_met_i, log_cool_i,
                                           log_cool_cmb_i, log_heat_i, log_n_h_i,
                                           zindex, clGridRank, clDataSize, clGridDim,
                                           clPar1, clPar2, clPar3,
                                           clCooling, clHeating,
                                           dclPar_1, dclPar_3,
                                           log10_tCMB);
                                           */

}
