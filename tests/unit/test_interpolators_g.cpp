#include "gtest/gtest.h"

extern "C" {
    #include "grackle_macros.h"
}


// TODO: FORTRAN_NAME for unit tests not needed
extern "C" void FORTRAN_NAME(interpolate_1d_g)(double *input1,
                                               long long *gridDim,
                                               double *gridPar1, double *dgridPar1,
                                               long long *dataSize, double *dataField,
                                               double *value);

extern "C" void FORTRAN_NAME(interpolate_2d_g)(double *input1, double *input2,
                                               long long *gridDim,
                                               double *gridPar1, double *dgridPar1,
                                               double *gridPar2, double *dgridPar2,
                                               long long *dataSize, double *dataField,
                                               double *value);

extern "C" void FORTRAN_NAME(interpolate_3d_g)(double *input1, double *input2, double *input3,
                                               long long *gridDim,
                                               double *gridPar1, double *dgridPar1,
                                               double *gridPar2, double *dgridPar2,
                                               double *gridPar3, double *dgridPar3,
                                               long long *dataSize, double *dataField,
                                               double *value);

extern "C" void FORTRAN_NAME(interpolate_3dz_g)(double *input1, double *input2, double *input3,
                                                long long *gridDim,
                                                double *gridPar1, double *dgridPar1,
                                                double *gridPar2, long long *index2,
                                                double *gridPar3, double *dgridPar3,
                                                long long *dataSize, double *dataField,
                                                int *end_int,
                                                double *value);

extern "C" void FORTRAN_NAME(interpolate_2df3d_g)(double *input1, double *input3,
                                                  long long *gridDim,
                                                  double *gridPar1, double *dgridPar1,
                                                  long long *index2,
                                                  double *gridPar3, double *dgridPar3,
                                                  long long *dataSize, double *dataField,
                                                  double *value);

extern "C" void FORTRAN_NAME(interpolate_4d_g)(double *input1, double *input2, double *input3, double *input4,
                                               long long *gridDim,
                                               double *gridPar1, double *dgridPar1,
                                               double *gridPar2, double *dgridPar2,
                                               double *gridPar3, double *dgridPar3,
                                               double *gridPar4, double *dgridPar4,
                                               long long *dataSize, double *dataField,
                                               double *value);

extern "C" void FORTRAN_NAME(interpolate_5d_g)(double *input1, double *input2, double *input3, double *input4, double *input5,
                                               long long *gridDim,
                                               double *gridPar1, double *dgridPar1,
                                               double *gridPar2, double *dgridPar2,
                                               double *gridPar3, double *dgridPar3,
                                               double *gridPar4, double *dgridPar4,
                                               double *gridPar5, double *dgridPar5,
                                               long long *dataSize, double *dataField,
                                               double *value);

TEST(InterpolationTest, Interpolate1D) {
    double input1 = 5.0;
    std::vector<long long> gridDim = {10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_1d_g)(&input1,
                                   gridDim.data(),
                                   gridPar1.data(), &dgridPar1,
                                   &dataSize, dataField.data(),
                                   &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate2D) {
    double input1 = 5.0;
    double input2 = 7.0;
    std::vector<long long> gridDim = {10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar2(10, 0.0);
    double dgridPar2 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_2d_g)(&input1, &input2,
                                   gridDim.data(),
                                   gridPar1.data(), &dgridPar1,
                                   gridPar2.data(), &dgridPar2,
                                   &dataSize, dataField.data(),
                                   &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate3D) {
    double input1 = 5.0;
    double input2 = 7.0;
    double input3 = 9.0;
    std::vector<long long> gridDim = {10, 10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar2(10, 0.0);
    double dgridPar2 = 1.0;
    std::vector<double> gridPar3(10, 0.0);
    double dgridPar3 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_3d_g)(&input1, &input2, &input3,
                                   gridDim.data(),
                                   gridPar1.data(), &dgridPar1,
                                   gridPar2.data(), &dgridPar2,
                                   gridPar3.data(), &dgridPar3,
                                   &dataSize, dataField.data(),
                                   &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate3Dz) {
    double input1 = 5.0;
    double input2 = 7.0;
    double input3 = 9.0;
    std::vector<long long> gridDim = {10, 10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar2(10, 0.0);
    double dgridPar2 = 1.0;
    long long index2 = 5;
    std::vector<double> gridPar3(10, 0.0);
    double dgridPar3 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    int end_int;
    double value;

    FORTRAN_NAME(interpolate_3dz_g)(&input1, &input2, &input3,
                                    gridDim.data(),
                                    gridPar1.data(), &dgridPar1,
                                    gridPar2.data(), &index2,
                                    gridPar3.data(), &dgridPar3,
                                    &dataSize, dataField.data(),
                                    &end_int,
                                    &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate2Df3D) {
    double input1 = 5.0;
    double input3 = 7.0;
    std::vector<long long> gridDim = {10, 10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar3(10, 0.0);
    double dgridPar3 = 1.0;
    long long index2 = 5;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_2df3d_g)(&input1, &input3,
                                      gridDim.data(),
                                      gridPar1.data(), &dgridPar1,
                                      &index2,
                                      gridPar3.data(), &dgridPar3,
                                      &dataSize, dataField.data(),
                                      &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate4D) {
    double input1 = 5.0;
    double input2 = 7.0;
    double input3 = 10.0;
    double input4 = 15.0;
    std::vector<long long> gridDim = {10, 10, 10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar2(10, 0.0);
    double dgridPar2 = 1.0;
    std::vector<double> gridPar3(10, 0.0);
    double dgridPar3 = 1.0;
    std::vector<double> gridPar4(10, 0.0);
    double dgridPar4 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_4d_g)(&input1, &input2, &input3, &input4,
                                   gridDim.data(),
                                   gridPar1.data(), &dgridPar1,
                                   gridPar2.data(), &dgridPar2,
                                   gridPar3.data(), &dgridPar3,
                                   gridPar4.data(), &dgridPar4,
                                   &dataSize, dataField.data(),
                                   &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}

TEST(InterpolationTest, Interpolate5D) {
    double input1 = 5.0;
    double input2 = 7.0;
    double input3 = 10.0;
    double input4 = 15.0;
    double input5 = 20.0;
    std::vector<long long> gridDim = {10, 10, 10, 10, 10};
    std::vector<double> gridPar1(10, 0.0);
    double dgridPar1 = 1.0;
    std::vector<double> gridPar2(10, 0.0);
    double dgridPar2 = 1.0;
    std::vector<double> gridPar3(10, 0.0);
    double dgridPar3 = 1.0;
    std::vector<double> gridPar4(10, 0.0);
    double dgridPar4 = 1.0;
    std::vector<double> gridPar5(10, 0.0);
    double dgridPar5 = 1.0;
    long long dataSize = 10;
    std::vector<double> dataField(dataSize, 0.0);
    double value;

    FORTRAN_NAME(interpolate_5d_g)(&input1, &input2, &input3, &input4, &input5,
                                   gridDim.data(),
                                   gridPar1.data(), &dgridPar1,
                                   gridPar2.data(), &dgridPar2,
                                   gridPar3.data(), &dgridPar3,
                                   gridPar4.data(), &dgridPar4,
                                   gridPar5.data(), &dgridPar5,
                                   &dataSize, dataField.data(),
                                   &value);

    EXPECT_DOUBLE_EQ(value, 0.);
}
