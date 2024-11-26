#include <cmath>
#include "gtest/gtest.h"

extern "C" {
    #include "grackle_macros.h"
}
// TODO: ispecies is not used in the subroutine, remove?
extern "C" void FORTRAN_NAME(cool1d_cloudy_old_tables_g)(double* d, double* de, double* rhoH, double* metallicity,
                                                         int* in, int* jn, int* kn, int* is, int* ie, int* j, int* k,
                                                         double* logtem, double* edot, double* comp2, int* ispecies,
                                                         double* dom, double* zr, int* icmbTfloor, int* iClHeat, double* clEleFra,
                                                         long long* clGridRank, long long*  clGridDim,
                                                         double* clPar1, double* clPar2, double* clPar3, double* clPar4, double* clPar5,
                                                         long long* clDataSize, double* clCooling, double* clHeating, int* itmask);


TEST(FuncCool1DCloudyOldTablesTest, Cool1DCloudyOldTablesTest) {

    int in = 2;
    int jn = 2;
    int kn = 2;
    int is = 0;
    int ie = in - 1;
    int j = 1; // TODO: not used in the subroutine, remove?
    int k = 1; // TODO: not used in the subroutine, remove?

    int itmask[in];
    double d[in*jn*kn];
    double de[in*jn*kn];
    double logtem[in], rhoH[in], zr, dom; // Input parameters used as independent variables ("x"-like) for interpolation
    double comp2 = 2.73;
    double clEleFra = 0.5;
    int iZscale = 1;
    int ispecies = 0; // TODO: ispecies is not used in the subroutine, remove?
    double metallicity[in]; // Used to scale cooling by metallicity
    double edot[in]; // Output parameter to be tested
    double edot_reference[in];

    // Grid descriptions
    long long clGridRank;
    long long* clGridDim;
    long long clDataSize;
    double* clPar1;
    double* clPar2;// = new double[clGridDim[1]] = {0.0, 0.12202000000000000, 0.25892999999999999,
                   //                0.41254000000000002, 0.58489000000000002};
    double* clPar3;// = new double[clGridDim[2]] = {1.0000000000000000, 1.0499998241554094, 1.0999998579423673,
                                                //1.1500001400088755 ,1.1999999472615537, 1.2499998998595714};
    double* clPar4;//[clGridDim[3]] = new double[clGridDim[3]] = {0.0, 1.0};
    double* clPar5;//[clGridDim[4]] = new double[clGridDim[4]]{0.0, 1.0};
    double* clCooling;
    double* clHeating;

    // ??
    int icmbTfloor;
    int iClHeat;

    for (int i=0; i<in; ++i) {
        metallicity[i] = 1.0;
        itmask[i] = 1;
        for (int j=0; j<jn; ++j) {
            for (int k=0; k<kn; ++k) {
                d[i*jn*kn+j*kn+k] = 1.0;
                de[i*jn*kn+j*kn+k] = 3.0;
            }
        }
    }

    ////////// clGridDim = 1 \\\\\\\\\\\

    clGridRank = 1;
    clGridDim = new long long[5]{3,1,1,1,1};// This must always be full size due to Fortran interface
    clDataSize = clGridDim[0]*clGridDim[1]*clGridDim[2]*clGridDim[3]*clGridDim[4];
    clPar1 = new double[clGridDim[0]]{0.6, 0.7, 0.8};//C++11
    clPar2 = new double[clGridDim[1]];
    clPar3 = new double[clGridDim[2]];
    clPar4 = new double[clGridDim[3]];
    clPar5 = new double[clGridDim[4]];
    clCooling = new double[clDataSize]{33.435745212863388, 33.461567803779815, 33.487283203685813};
    clHeating = new double[clDataSize]{39.224844181639753, 39.224844181639753, 39.202567762461001};

    // Independent variables for interpolation
    for (int i=0; i<in; ++i) {
        rhoH[i] = 1.0;
        logtem[i] = log(comp2 + 2. + (1.*i)/in);
    }
    zr = 0.;
    dom = 0.99843257445223532;


    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -2.8515191924225000E+033;
    edot_reference[1] = -2.9263454131841624E+033;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }


    //// icmbTfloor = 1 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -3.7730004524134467E+032;
    edot_reference[1] = -4.5212626600300711E+032;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 0;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 1.6781989448879827E+039;
    edot_reference[1] = 1.6623478396918585E+039;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }


    //// icmbTfloor = 1 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 1;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 1.6782014191071299E+039;
    edot_reference[1] = 1.6623503139110057E+039;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    delete [] clGridDim;
    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;
    delete [] clCooling;
    delete [] clHeating;

    ////////// clGridDim = 2 \\\\\\\\\\\

    clGridRank = 2;
    clGridDim = new long long[5]{5,3,1,1,1};// This must always be full size due to Fortran interface
    clDataSize = clGridDim[0]*clGridDim[1]*clGridDim[2]*clGridDim[3]*clGridDim[4];
    clPar1 = new double[clGridDim[0]]{0.6, 0.65, 0.7, 0.75, 0.8};//C++11
    clPar2 = new double[clGridDim[1]]{0.6, 0.7, 0.8};
    clPar3 = new double[clGridDim[2]];
    clPar4 = new double[clGridDim[3]];
    clPar5 = new double[clGridDim[4]];
    clCooling = new double[clDataSize]{33.435745212863388, 33.461567803779815, 33.487283203685813,
                                      33.535745212863388, 33.561567803779815, 33.587283203685813,
                                      33.635745212863388, 33.661567803779815, 33.687283203685813,
                                      33.735745212863388, 33.761567803779815, 33.787283203685813,
                                      33.835745212863388, 33.861567803779815, 33.887283203685813};
    clHeating = new double[clDataSize]{39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001};

    // Independent variables for interpolation
    for (int i=0; i<in; ++i) {
        rhoH[i] = 5.0 + (1.*i)/in;
        logtem[i] = log(comp2 + 2. + (1.*i)/in);
    }
    zr = 0.;
    dom = 0.99843257445223532;

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -2.2419390180022109e+34;
    edot_reference[1] = -3.0623240063001747e+34;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -2.966431701278334e+33;
    edot_reference[1] = -4.7313523277941653e+33;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 0;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 1.3194439316812583E+040;
    edot_reference[1] = 1.7395922140204456E+040;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 1;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 1.3194458769771063E+040;
    edot_reference[1] = 1.7395948032092192E+040;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }


    delete [] clGridDim;
    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;
    delete [] clCooling;
    delete [] clHeating;

    ////////// clGridDim = 3 \\\\\\\\\\\

    clGridRank = 3;
    clGridDim = new long long[5]{5,3,2,1,1};// This must always be full size due to Fortran interface
    clDataSize = clGridDim[0]*clGridDim[1]*clGridDim[2]*clGridDim[3]*clGridDim[4];
    clPar1 = new double[clGridDim[0]]{0.6, 0.65, 0.7, 0.75, 0.8};//C++11
    clPar2 = new double[clGridDim[1]]{0.6, 0.7, 0.8};
    clPar3 = new double[clGridDim[2]]{0.6, 0.8};
    clPar4 = new double[clGridDim[3]];
    clPar5 = new double[clGridDim[4]];
    clCooling = new double[clDataSize]{33.435745212863388, 33.461567803779815, 33.487283203685813,
                                      33.535745212863388, 33.561567803779815, 33.587283203685813,
                                      33.635745212863388, 33.661567803779815, 33.687283203685813,
                                      33.735745212863388, 33.761567803779815, 33.787283203685813,
                                      33.835745212863388, 33.861567803779815, 33.887283203685813,
                                      33.935745212863388, 33.961567803779815, 33.987283203685813,
                                      33.135745212863388, 33.161567803779815, 33.187283203685813,
                                      33.235745212863388, 33.261567803779815, 33.287283203685813,
                                      33.335745212863388, 33.361567803779815, 33.387283203685813,
                                      33.435745212863388, 33.461567803779815, 33.487283203685813};
    clHeating = new double[clDataSize]{39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001,
                                       39.724844181639753, 39.724844181639753, 39.702567762461001,
                                       39.824844181639753, 39.824844181639753, 39.802567762461001,
                                       39.924844181639753, 39.924844181639753, 39.902567762461001,
                                       39.124844181639753, 39.124844181639753, 39.102567762461001,
                                       39.6224844181639753, 39.224844181639753, 39.202567762461001};

    // Independent variables for interpolation
    for (int i=0; i<in; ++i) {
        rhoH[i] = 5.0 + (1.*i)/in;
        metallicity[i] = 5.0 + (1.*i)/in;
        logtem[i] = log(comp2 + 2. + (1.*i)/in);
    }
    zr = 0.;
    dom = 0.99843257445223532;

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -3.9530260562451938e+34;
    edot_reference[1] = -1.33463762488875e+34;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -4.9069223101035244e+33;
    edot_reference[1] = -1.5998148327232414e+33;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 2.1888638944122671e+40;
    edot_reference[1] = 3.8805183328955476e+40;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 1;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 2.1888673567460923e+40;
    edot_reference[1] = 3.8805195075516887e+40;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);


    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    delete [] clGridDim;
    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;
    delete [] clCooling;
    delete [] clHeating;

    ////////// clGridDim = 4 \\\\\\\\\\\

    clGridRank = 4;
    clGridDim = new long long[5]{5,3,2,2,1};// This must always be full size due to Fortran interface
    clDataSize = clGridDim[0]*clGridDim[1]*clGridDim[2]*clGridDim[3]*clGridDim[4];
    clPar1 = new double[clGridDim[0]]{0.6, 0.65, 0.7, 0.75, 0.8};//C++11
    clPar2 = new double[clGridDim[1]]{0.6, 0.7, 0.8};
    clPar3 = new double[clGridDim[2]]{0.6, 0.8};
    clPar4 = new double[clGridDim[3]]{0.6, 0.8};
    clPar5 = new double[clGridDim[4]];
    clCooling = new double[clDataSize]{33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.535745212863388, 33.561567803779815, 33.587283203685813,
                                       33.635745212863388, 33.661567803779815, 33.687283203685813,
                                       33.735745212863388, 33.761567803779815, 33.787283203685813,
                                       33.835745212863388, 33.861567803779815, 33.887283203685813,
                                       33.935745212863388, 33.961567803779815, 33.987283203685813,
                                       33.135745212863388, 33.161567803779815, 33.187283203685813,
                                       33.235745212863388, 33.261567803779815, 33.287283203685813,
                                       33.335745212863388, 33.361567803779815, 33.387283203685813,
                                       33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.335745212863388, 33.361567803779815, 33.387283203685813,
                                       33.235745212863388, 33.261567803779815, 33.287283203685813,
                                       33.135745212863388, 33.161567803779815, 33.187283203685813,
                                       33.035745212863388, 33.061567803779815, 33.087283203685813,
                                       33.935745212863388, 33.961567803779815, 33.987283203685813,
                                       33.835745212863388, 33.861567803779815, 33.887283203685813,
                                       33.735745212863388, 33.761567803779815, 33.787283203685813,
                                       33.635745212863388, 33.661567803779815, 33.687283203685813,
                                       33.535745212863388, 33.561567803779815, 33.587283203685813};
    clHeating = new double[clDataSize]{39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001,
                                       39.724844181639753, 39.724844181639753, 39.702567762461001,
                                       39.824844181639753, 39.824844181639753, 39.802567762461001,
                                       39.924844181639753, 39.924844181639753, 39.902567762461001,
                                       39.124844181639753, 39.124844181639753, 39.102567762461001,
                                       39.9224844181639753, 39.924844181639753, 39.902567762461001,
                                       39.824844181639753, 39.824844181639753, 39.802567762461001,
                                       39.724844181639753, 39.724844181639753, 39.702567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.124844181639753, 39.124844181639753, 39.102567762461001,
                                       39.024844181639753, 39.024844181639753, 39.002567762461001,
                                       39.1224844181639753, 39.124844181639753, 39.102567762461001};

    // Independent variables for interpolation
    for (int i=0; i<in; ++i) {
        rhoH[i] = 5.0 + (1.*i)/in;
        metallicity[i] = 5.0 + (1.*i)/in;
        logtem[i] = log(comp2 + 2. + (1.*i)/in);
    }
    zr = 0.;
    dom = 0.99843257445223532;

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -1.1354844090587073e+35;
    edot_reference[1] = -4.3014010232448388e+34;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 1;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -7.4513497032104963e+33;
    edot_reference[1] = -4.1659718396574788e+34;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 0;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 3.6078199934769295e+41;
    edot_reference[1] = 1.6266705682839908e+41;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }


    //// icmbTfloor = 1 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 1;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 3.6078210544478416e+41;
    edot_reference[1] = 1.6266705818269093e+41;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    delete [] clGridDim;
    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;
    delete [] clCooling;
    delete [] clHeating;

    ////////// clGridDim = 5 \\\\\\\\\\\

    clGridRank = 5;
    clGridDim = new long long[5]{5,3,2,2,2};// This must always be full size due to Fortran interface
    clDataSize = clGridDim[0]*clGridDim[1]*clGridDim[2]*clGridDim[3]*clGridDim[4];
    clPar1 = new double[clGridDim[0]]{0.6, 0.65, 0.7, 0.75, 0.8};//C++11
    clPar2 = new double[clGridDim[1]]{0.6, 0.7, 0.8};
    clPar3 = new double[clGridDim[2]]{0.6, 0.8};
    clPar4 = new double[clGridDim[3]]{0.6, 0.8};
    clPar5 = new double[clGridDim[4]]{0.6, 0.8};
    clCooling = new double[clDataSize]{33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.535745212863388, 33.561567803779815, 33.587283203685813,
                                       33.635745212863388, 33.661567803779815, 33.687283203685813,
                                       33.735745212863388, 33.761567803779815, 33.787283203685813,
                                       33.835745212863388, 33.861567803779815, 33.887283203685813,
                                       33.935745212863388, 33.961567803779815, 33.987283203685813,
                                       33.135745212863388, 33.161567803779815, 33.187283203685813,
                                       33.235745212863388, 33.261567803779815, 33.287283203685813,
                                       33.335745212863388, 33.361567803779815, 33.387283203685813,
                                       33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.435745212863388, 33.461567803779815, 33.487283203685813,
                                       33.335745212863388, 33.361567803779815, 33.387283203685813,
                                       33.235745212863388, 33.261567803779815, 33.287283203685813,
                                       33.135745212863388, 33.161567803779815, 33.187283203685813,
                                       33.035745212863388, 33.061567803779815, 33.087283203685813,
                                       33.935745212863388, 33.961567803779815, 33.987283203685813,
                                       33.835745212863388, 33.861567803779815, 33.887283203685813,
                                       33.735745212863388, 33.761567803779815, 33.787283203685813,
                                       33.635745212863388, 33.661567803779815, 33.687283203685813,
                                       33.535745212863388, 33.561567803779815, 33.587283203685813,
                                       32.435745212863388, 32.461567803779815, 32.487283203685813,
                                       32.535745212863388, 32.561567803779815, 32.587283203685813,
                                       32.635745212863388, 32.661567803779815, 32.687283203685813,
                                       32.735745212863388, 32.761567803779815, 32.787283203685813,
                                       32.835745212863388, 32.861567803779815, 32.887283203685813,
                                       32.935745212863388, 32.961567803779815, 32.987283203685813,
                                       32.135745212863388, 32.161567803779815, 32.187283203685813,
                                       32.235745212863388, 32.261567803779815, 32.287283203685813,
                                       32.335745212863388, 32.361567803779815, 32.387283203685813,
                                       32.435745212863388, 32.461567803779815, 32.487283203685813,
                                       32.435745212863388, 32.461567803779815, 32.487283203685813,
                                       32.335745212863388, 32.361567803779815, 32.387283203685813,
                                       32.235745212863388, 32.261567803779815, 32.287283203685813,
                                       32.135745212863388, 32.161567803779815, 32.187283203685813,
                                       32.035745212863388, 32.061567803779815, 32.087283203685813,
                                       32.935745212863388, 32.961567803779815, 32.987283203685813,
                                       32.835745212863388, 32.861567803779815, 32.887283203685813,
                                       32.735745212863388, 32.761567803779815, 32.787283203685813,
                                       32.635745212863388, 32.661567803779815, 32.687283203685813,
                                       32.535745212863388, 32.561567803779815, 32.587283203685813};
    clHeating = new double[clDataSize]{39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001,
                                       39.724844181639753, 39.724844181639753, 39.702567762461001,
                                       39.824844181639753, 39.824844181639753, 39.802567762461001,
                                       39.924844181639753, 39.924844181639753, 39.902567762461001,
                                       39.124844181639753, 39.124844181639753, 39.102567762461001,
                                       39.9224844181639753, 39.924844181639753, 39.902567762461001,
                                       39.824844181639753, 39.824844181639753, 39.802567762461001,
                                       39.724844181639753, 39.724844181639753, 39.702567762461001,
                                       39.624844181639753, 39.624844181639753, 39.602567762461001,
                                       39.524844181639753, 39.524844181639753, 39.502567762461001,
                                       39.424844181639753, 39.424844181639753, 39.402567762461001,
                                       39.324844181639753, 39.324844181639753, 39.302567762461001,
                                       39.224844181639753, 39.224844181639753, 39.202567762461001,
                                       39.124844181639753, 39.124844181639753, 39.102567762461001,
                                       39.024844181639753, 39.024844181639753, 39.002567762461001,
                                       39.1224844181639753, 39.124844181639753, 39.102567762461001,
                                       38.224844181639753, 38.224844181639753, 38.202567762461001,
                                       38.324844181639753, 38.324844181639753, 38.302567762461001,
                                       38.424844181639753, 38.424844181639753, 38.402567762461001,
                                       38.524844181639753, 38.524844181639753, 38.502567762461001,
                                       38.624844181639753, 38.624844181639753, 38.602567762461001,
                                       38.724844181639753, 38.724844181639753, 38.702567762461001,
                                       38.824844181639753, 38.824844181639753, 38.802567762461001,
                                       38.924844181639753, 38.924844181639753, 38.902567762461001,
                                       38.124844181639753, 38.124844181639753, 38.102567762461001,
                                       38.9224844181639753, 38.924844181639753, 38.902567762461001,
                                       38.824844181639753, 38.824844181639753, 38.802567762461001,
                                       38.724844181639753, 38.724844181639753, 38.702567762461001,
                                       38.624844181639753, 38.624844181639753, 38.602567762461001,
                                       38.524844181639753, 38.524844181639753, 38.502567762461001,
                                       38.424844181639753, 38.424844181639753, 38.402567762461001,
                                       38.324844181639753, 38.324844181639753, 38.302567762461001,
                                       38.224844181639753, 38.224844181639753, 38.202567762461001,
                                       38.124844181639753, 38.124844181639753, 38.102567762461001,
                                       38.024844181639753, 38.024844181639753, 38.002567762461001,
                                       38.1224844181639753, 38.124844181639753, 38.102567762461001};

    // Independent variables for interpolation
    for (int i=0; i<in; ++i) {
        rhoH[i] = 5.0 + (1.*i)/in;
        metallicity[i] = 5.0 + (1.*i)/in;
        logtem[i] = log(comp2 + 2. + (1.*i)/in);
    }
    zr = 0.7;
    dom = 0.99843257445223532;

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 0;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = -1.8035835890887779e+38;
    edot_reference[1] = -1.7290148698993495e+34;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 0 \\\\

    icmbTfloor = 1;
    iClHeat = 0;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 2.9469425985396465e+38;
    edot_reference[1] = -3.4385534804228648e+33;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 0 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 0;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 9.1199865036544127e+42;
    edot_reference[1] = 3.6540571532577576e+41;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    //// icmbTfloor = 1 \\\\
    //// iClHeat = 1 \\\\

    icmbTfloor = 1;
    iClHeat = 1;
    for (int i=0; i<in; ++i) {
        edot[i] = 0.0;
    }
    edot_reference[0] = 9.1204615562731772e+42;
    edot_reference[1] = 3.6540572917737097e+41;
    FORTRAN_NAME(cool1d_cloudy_old_tables_g)(d, de, rhoH, metallicity,
                                             &in, &jn, &kn, &is, &ie, &j, &k,
                                             logtem, edot, &comp2, &ispecies, &dom, &zr,
                                             &icmbTfloor, &iClHeat, &clEleFra,
                                             &clGridRank, clGridDim,
                                             clPar1, clPar2, clPar3, clPar4, clPar5,
                                             &clDataSize, clCooling, clHeating, itmask);

    for (int i=0; i<in; ++i) {
        EXPECT_DOUBLE_EQ(edot[i], edot_reference[i]);
        //ASSERT_NEAR(edot[i], edot_reference[i], 1e-8);
    }

    delete [] clGridDim;
    delete [] clPar1;
    delete [] clPar2;
    delete [] clPar3;
    delete [] clPar4;
    delete [] clPar5;
    delete [] clCooling;
    delete [] clHeating;

}
