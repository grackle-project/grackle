#include "gtest/gtest.h"
#include <cmath>

extern "C" {
    #include "grackle_macros.h"
}

extern "C" void FORTRAN_NAME(compute_redshift_dimension)(double* zr, long long* clGridDim_i, double* clPar_i,
                                                         long long* zindex, int* end_int, int* get_heat);


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
                                                       double* log_n_h_i, double* zr, double* log10tem_i, double* log10_tCMB,
                                                       long long* clGridRank, long long* clDataSize, long long* clGridDim,
                                                       double* clPar1, double* clPar2, double* clPar3,
                                                       double* clCooling, double* clHeating,
                                                       double* dclPar_1, long long* zindex, double* dclPar_3,
                                                       double* edot_met_i);

TEST(UnitCool1DCloudy, ComputeRedshiftDimensionsTest) {

    double zr = 0.5;
    long long clGridDim_i = 5;
    double clPar_i[clGridDim_i] = {1., 2., 3., 4., 5.};
    long long zindex = 0;
    int end_int = 0;
    int get_heat = 0;

    FORTRAN_NAME(compute_redshift_dimension)(&zr, &clGridDim_i, clPar_i,
                                             &zindex, &end_int, &get_heat);
    ASSERT_EQ(zindex, 1);
    ASSERT_EQ(end_int, 0);
    ASSERT_EQ(get_heat, 0);

    zr = clPar_i[clGridDim_i-2] + 0.5;
    FORTRAN_NAME(compute_redshift_dimension)(&zr, &clGridDim_i, clPar_i,
                                             &zindex, &end_int, &get_heat);
    ASSERT_EQ(zindex, clGridDim_i);
    ASSERT_EQ(end_int, 1);
    ASSERT_EQ(get_heat, 0);

    zr = clPar_i[clGridDim_i-3] + 0.5;
    FORTRAN_NAME(compute_redshift_dimension)(&zr, &clGridDim_i, clPar_i,
                                             &zindex, &end_int, &get_heat);
    ASSERT_EQ(zindex, clGridDim_i-2);

    zr = clPar_i[1] + 0.25;
    FORTRAN_NAME(compute_redshift_dimension)(&zr, &clGridDim_i, clPar_i,
                                             &zindex, &end_int, &get_heat);
    ASSERT_EQ(zindex, 2);

}

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

    int get_heat = 0;
    int icmbTfloor = 0;
    double log_n_h_i = 2.5;
    double zr = 0.5;
    double log10tem_i = 1.5;
    double log10_tCMB = log10tem_i - 0.5;
    long long clGridRank = 3;
    long long clDataSize1 = 5;
    long long clDataSize2 = 3;
    long long clDataSize3 = 2;
    long long clDataSize = clDataSize1*clDataSize2*clDataSize3;
    long long clGridDim[clGridRank] = {clDataSize1, clDataSize2, clDataSize3};
    double clPar1[clDataSize1] = {0., 1., 2., 3., 4.,};
    double clPar2[clDataSize2] = {0.5, 1., 1.5};
    double clPar3[clDataSize3] = {1., 2.};
    double clCooling[clDataSize] = {0., 0.3, 0.6, 1., 1.3, 1.6, 2., 2.3, 2.6,
                                    3., 3.3, 3.6, 4., 4.3, 4.6, 4., 4.3, 4.6,
                                    3., 3.3, 3.6, 2., 2.3, 2.6, 1., 1.3, 1.6,
                                    0., 0.3, 0.6};
    double clHeating[clDataSize] = {0.5, 0.6, 0.7, 0.55, 0.65, 0.75, 0.6, 0.65, 0.7,
                                    0.65, 0.7, 0.75, 0.7, 0.75, 0.8, 0.7, 0.75, 0.8,
                                    0.65, 0.7, 0.75, 0.6, 0.65, 0.7, 0.55, 0.65, 0.75,
                                    0.5, 0.6, 0.7};
    double dclPar_1 = 1.0;
    long long zindex = 1;
    double dclPar_3 = 0.5;
    double edot_met_i = 0.;



    FORTRAN_NAME(cool_rank3_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &zr, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2, clPar3,
                                           clCooling, clHeating,
                                           &dclPar_1, &zindex, &dclPar_3,
                                           &edot_met_i);
    // TODO: use near with a tolerance
    ASSERT_EQ(edot_met_i, -4466.8359215096352);

    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank3_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &zr, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2, clPar3,
                                           clCooling, clHeating,
                                           &dclPar_1, &zindex, &dclPar_3,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -1304.5582613412557);

    get_heat = 1;
    icmbTfloor = 0;
    FORTRAN_NAME(cool_rank3_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &zr, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2, clPar3,
                                           clCooling, clHeating,
                                           &dclPar_1, &zindex, &dclPar_3,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -4461.8240491733623);

    icmbTfloor = 1;
    FORTRAN_NAME(cool_rank3_interpolation)(&get_heat, &icmbTfloor,
                                           &log_n_h_i, &zr, &log10tem_i, &log10_tCMB,
                                           &clGridRank, &clDataSize, clGridDim,
                                           clPar1, clPar2, clPar3,
                                           clCooling, clHeating,
                                           &dclPar_1, &zindex, &dclPar_3,
                                           &edot_met_i);
    ASSERT_EQ(edot_met_i, -1299.546389004983);
}
