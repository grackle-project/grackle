#include "gtest/gtest.h"

extern "C" {
    #include "grackle_macros.h"
}


extern "C" void FORTRAN_NAME(calc_tdust_1d_g)(double* tdust, double* tgas, double* nh, double* gasgr,
                                              double* gamma_isrfa, double* isrf, int* itmask, double* trad,
                                              int* in, int* is, int* ie, int* j, int* k);

TEST(FuncCalcTDust1DTest, CalcTDust1DTest) {

    int in = 10;
    int is = 0;
    int ie = in - 1;
    double tdust[in];
    double tgas[in] = {11.0, 11.2, 11.3, 11.4, 11.5, 11.4, 11.35, 11.25, 11.15, 11.05};
    double nh[in] = {1.1, 1.2, 2.1, 3.2, 2.3, 4.3, 3.4, 5.4, 4.5, 6.0};
    double gasgr[in] = {0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59};
    double gamma_isrfa = 1.0;
    double isrf[in] = {0.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    int itmask[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double trad = 2.73;
    int j = 0;
    int k = 0;

    double tdust_reference[in] = {10.90543148329828, 12.18628896658262, 11.808040273194752,
                                  11.677495935870839, 11.78768899170481, 11.521945069945298,
                                  11.452894618067857, 11.286738480633945, 11.15972409842734,
                                  11.032215177747794};


    FORTRAN_NAME(calc_tdust_1d_g)(tdust, tgas, nh, gasgr, &gamma_isrfa,
                                  isrf, itmask, &trad, &in, &is, &ie, &j, &k);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(tdust[i], tdust_reference[i]);
        //ASSERT_NEAR(tdust[i], tdust_reference[i], 1e-8);
    }
}
