#include "gtest/gtest.h"
#include <cmath>

extern "C" {
    #include "grackle_macros.h"
}

// TODO: FORTRAN_NAME not needed for unit tests
extern "C" void FORTRAN_NAME(set_iteration_mask_and_initial_guess)(double *gamma_isrf, double *isrf, double *gamma_isrfa,
                                                                   int *itmask, int *nm_itmask, int *bi_itmask,
                                                                   double *tdustnow, double *tgas, double *trad,
                                                                   double *t_subl, double *radf, double *kgr1,
                                                                   double *pert_i, double *pert, int *c_done, int *nm_done,
                                                                   int *is, int *ie);

extern "C" void FORTRAN_NAME(calc_tdplus)(int* in, int *is, int *ie, int *nm_itmask,
                                          double *pert, double *tdustnow, double *tdplus);

extern "C" void FORTRAN_NAME(run_newton_method_iteration)(int* in, int* is, int* ie, int* nm_itmask, int* bi_itmask,
                                                          double* solplus, double* sol, double* slope, double* tdustnow,
                                                          double* tdustold, double* pert, double* minpert, double* trad,
                                                          double* tol, int* nm_done, int* c_done);

extern "C" void FORTRAN_NAME(run_bisection_method_iteration)(int* in, int* is, int* ie, int* bi_itmask, double* tdustnow,
                                                             double* bi_t_high, int* iter, double* t_subl,
                                                             double* kgr, double* tgas, double* trad4, double* gasgr,
                                                             double* gamma_isrf, double* nh, double* sol, double* bi_tol,
                                                             int* c_done, int* c_total);

extern "C" void FORTRAN_NAME(calc_gr_balance_g)(double *tdust, double *tgas, double *kgr, double *trad4,
                                                double *gasgr, double *gamma_isrf, double *nh,
                                                int *itmask, double *sol, int *in, int *is, int *ie);

extern "C" void FORTRAN_NAME(calc_kappa_gr_g)(double *tdust, double *kgr, int *itmask,
                                              int *in, int *is, int *ie, double *t_subl);

TEST(UnitCalcTDust1DTest, CalcGRBalanceTest) {

    int in = 10;
    double tdust[in] = {3.0, 3.2, 3.3, 3.4, 3.5, 3.6, 3.35, 3.25, 3.15, 3.05};
    double tgas[in] = {11.0, 11.2, 11.3, 11.4, 11.5, 11.4, 11.35, 11.25, 11.15, 11.05};
    double kgr[in] = {4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9};
    double trad4 = std::pow(2.73, 4.);
    double gasgr[in] = {0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59};
    double nh[in] = {1.1, 1.2, 2.1, 3.2, 2.3, 4.3, 3.4, 5.4, 4.5, 6.0};
    double gamma_isrf[in] = {0.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    int itmask[in] = {1, 1, 0, 1, 1, 0, 1, 1, 1, 1};
    int is = 0;
    int ie = in - 1;
    double sol[in];
    double sol_reference[in] = {4.4769063566300273,
                                5.7501428510406694,
                                6.9532564830293191e-310,
                                14.191840564552038,
                                10.441673608478029,
                                5.9287877500949585e-323,
                                15.558549521012679,
                                24.864280259195358,
                                21.033283071182687,
                                28.385557051983135};

    FORTRAN_NAME(calc_gr_balance_g)(tdust, tgas, kgr, &trad4,
                                    gasgr, gamma_isrf, nh,
                                    itmask, sol, &in, &is, &ie);


    for (int i=0;i<in;++i) {
        ASSERT_NEAR(sol[i], sol_reference[i], 1e-300);
    }
}

TEST(UnitCalcTDust1DTest, CalcKappaGRBalanceTest) {

    int in = 5;
    int is = 0;
    int ie = in - 1;
    double t_subl = 200.0;
    double tdust[in] = {100.0, 102.0, 103.0, 104.0, 105.0};
    int itmask[in] = {1, 1, 1, 1, 1};
    double kgr[in];
    double kgr_reference[in] = {4.0, 4.1616, 4.2436, 4.3264000000000005, 4.41};

    FORTRAN_NAME(calc_kappa_gr_g)(tdust, kgr, itmask, &in, &is, &ie, &t_subl);

    for (int i=0;i<in;++i) {
        ASSERT_EQ(kgr[i], kgr_reference[i]);
    }
}

TEST(UnitCalcTDust1DTest, SetIterationMaskAndInitialGuessTest) {

    const int size = 10;
    double gamma_isrf[size], isrf[size], gamma_isrfa = 1.0;
    int itmask[size], nm_itmask[size], bi_itmask[size];
    double tdustnow[size], tgas[size], pert[size], trad = 100.0;
    double t_subl = 1500.0, radf = 4.0, kgr1 = 4.0e-4, pert_i = 0.001;
    int c_done = 0, nm_done = 0, is = 0, ie = size - 1;

    for (int i = 0; i < size; ++i) {
        gamma_isrf[i] = 0.0;
        isrf[i] = 0.0;
        itmask[i] = 1;
        nm_itmask[i] = 1;
        bi_itmask[i] = 1;
        tdustnow[i] = 0.0;
        tgas[i] = 0.0;
    }

    FORTRAN_NAME(set_iteration_mask_and_initial_guess)(
        gamma_isrf, isrf, &gamma_isrfa, itmask, nm_itmask, bi_itmask,
        tdustnow, tgas, &trad, &t_subl, &radf, &kgr1,
        &pert_i, pert, &c_done, &nm_done, &is, &ie);

    for (int i = 0; i < size; ++i) {
        ASSERT_EQ(tdustnow[i], trad);
        ASSERT_FALSE(nm_itmask[i]);
        ASSERT_FALSE(bi_itmask[i]);
        ASSERT_EQ(c_done, size);
        ASSERT_EQ(nm_done, size);
    }
}

TEST(UnitCalcTDust1DTest, CalcTDPlusTest) {

    int in = 10;
    int is = 0;
    int ie = in - 1;
    int nm_itmask[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double pert[in] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1};
    double tdustnow[in] = {273.0, 275.0, 277.0, 279.0, 0.0, 283.0, 285.0, 287.0, 289.0, 291.0};
    double tdplus[in];
    double tdplus_reference[in] = {300.3, 330.0, 360.1, 390.6, 0.001,
                                   452.8, 484.5, 516.6, 549.1, 582.0};

    FORTRAN_NAME(calc_tdplus)(&in, &is, &ie, nm_itmask, pert, tdustnow, tdplus);

    for (int i=0;i<in; ++i) {
        EXPECT_DOUBLE_EQ(tdplus[i], tdplus_reference[i]);
    }
}


TEST(UnitCalcTDust1DTest, RunNewtonMethodIterationTest) {

    int in = 11;
    int is = 0;
    int ie = in - 1;
    int nm_itmask[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int bi_itmask[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double solplus[in] = {0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1};
    double sol[in] = {10.5, 9.4, 8.3, 7.2, 6.1, 5.0, 6.1, 7.2, 8.3, 9.4, 10.5};
    double slope[in];
    double tdustnow[in] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
    double tdustold[in];
    double pert[in] = {1.01, 1.02, 1.01, 1.03, 1.01, 1.04, 1.01, 1.01, 1.02, 1.01, 1.01};
    double minpert = 1.e-10;
    double trad = 1.0;
    double tol = 1.0e-5;
    int nm_done = 0;
    int c_done = 0;

    double tdustnow_reference[in] = {2.019712, 2.386195, 2.876600, 3.837053,
                                     6.790875, -11.500000, -14.829333,
                                     -18.904000, -23.598000, 62.028667,
                                     0.000000};
    int nm_itmask_reference[in] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0};
    int bi_itmask_reference[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


    FORTRAN_NAME(run_newton_method_iteration)(&in, &is, &ie, nm_itmask, bi_itmask,
                                              solplus, sol, slope, tdustnow, tdustold,
                                              pert, &minpert, &trad, &tol, &nm_done, &c_done);

    EXPECT_EQ(nm_done, 5);
    EXPECT_EQ(c_done, 0);
    for(int i=0; i<in; i++){
        ASSERT_NEAR(tdustnow[i], tdustnow_reference[i], 1.e-6);
    }
    for(int i=0; i<in; i++){
        EXPECT_EQ(nm_itmask[i], nm_itmask_reference[i]);
    }
    for(int i=0; i<in; i++){
        EXPECT_EQ(bi_itmask[i], bi_itmask_reference[i]);
    }
}

TEST(UnitCalcTDust1DTest, RunBisectionMethodTest) {


    int in = 11;
    int is = 0;
    int ie = in - 1;
    int bi_itmask[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double tdustnow[in] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
    double tgas[in] = {11.0, 11.2, 11.3, 11.4, 11.5, 11.4, 11.35, 11.25, 11.15, 11.05, 11.35};
    double bi_t_high[in];
    int iter;
    double t_subl = 1.5e3;
    double kgr[in];
    double trad4 = std::pow(2.73, 4.);
    double gasgr[in] = {0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60};
    double gamma_isrf[in] = {0.1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1};
    double nh[in] = {1.1, 1.2, 2.1, 3.2, 2.3, 4.3, 3.4, 5.4, 4.5, 6.5, 6.6};
    double sol[in] = {10.5, 9.4, 8.3, 7.2, 6.1, 5.0, 6.1, 7.2, 8.3, 9.4, 10.5};
    double bi_tol = 1.0e-3;
    int c_done = 0;
    int c_total = ie - is + 1;
    for(int i=0;i<in;i++){
        bi_t_high[i] = tgas[i];
    }

    double tdustnow_reference[in] = {6.000000, 6.150000, 6.250000, 6.350000,
                                     6.450000, 6.450000, 6.475000, 6.475000,
                                     6.475000, 6.475000, 5.675000};

    int bi_itmask_reference[in] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    iter = 1;
    FORTRAN_NAME(run_bisection_method_iteration)(&in, &is, &ie, bi_itmask, tdustnow, bi_t_high, &iter, &t_subl,
                                                 kgr, tgas, &trad4, gasgr, gamma_isrf, nh, sol, &bi_tol, &c_done, &c_total);


    EXPECT_EQ(c_done, 0);
    for(int i=0; i<in; i++){
        ASSERT_NEAR(tdustnow[i], tdustnow_reference[i], 1.e-6);
    }
    for(int i=0; i<in; i++){
        EXPECT_EQ(bi_itmask[i], bi_itmask_reference[i]);
    }

}
