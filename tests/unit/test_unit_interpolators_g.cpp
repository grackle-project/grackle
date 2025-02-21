#include "gtest/gtest.h"

#include "fortran_func_wrappers.hpp"

extern "C" {
    #include "grackle_macros.h"
}



TEST(InterpolationTest, Interpolate1D) {

    long long dataSize = 10;
    double input1 = 5.0;
    double dgridPar1 = 1.3;
    std::vector<long long> gridDim(1,dataSize);
    std::vector<double> gridPar1(dataSize);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize;i++){
        gridPar1[i]=1.0 + i*dgridPar1;
        dataField[i]=1.0*i;
    }

    
    double value = grackle::impl::fortran_wrapper::interpolate_1d_g(
        input1,
        gridDim.data(),
        gridPar1.data(), dgridPar1,
        dataSize, dataField.data()
    );

    EXPECT_DOUBLE_EQ(value, 3.0769230769230766);
}

TEST(InterpolationTest, Interpolate2D) {

    long long dataSize1 = 10;
    long long dataSize2 = 7;
    long long dataSize = dataSize1*dataSize2;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    double input1 = 2.3;
    double input2 = 1.7;
    std::vector<long long> gridDim = {dataSize1, dataSize2};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=1.0+i*dgridPar2;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            dataField[i*dataSize2+j]=1.0*i+2.0*j;
        }
    }


    double value = grackle::impl::fortran_wrapper::interpolate_2d_g(
        input1, input2,
        gridDim.data(),
        gridPar1.data(), dgridPar1,
        gridPar2.data(), dgridPar2,
        dataSize, dataField.data()
    );

    EXPECT_DOUBLE_EQ(value, 1.2333333333333329);
}

TEST(InterpolationTest, Interpolate3D) {

    long long dataSize1 = 10;
    long long dataSize2 = 7;
    long long dataSize3 = 5;
    long long dataSize = dataSize1*dataSize2*dataSize3;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    double dgridPar3 = 2.5;
    double input1 = 2.4;
    double input2 = 5.3;
    double input3 = 1.5;
    std::vector<long long> gridDim = {dataSize1, dataSize2, dataSize3};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> gridPar3(dataSize3);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=0.5+i*dgridPar2;
    }
    for(auto i=0;i<dataSize3;i++){
        gridPar3[i]=1.0+i*dgridPar3;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            for(auto k=0;k<dataSize3;k++){
                dataField[i*dataSize2*dataSize3+j*dataSize3+k]=1.*i+2.*j+3.*k;
            }
        }
    }

    double value = grackle::impl::fortran_wrapper::interpolate_3d_g(
        input1, input2, input3,
        gridDim.data(),
        gridPar1.data(), dgridPar1,
        gridPar2.data(), dgridPar2,
        gridPar3.data(), dgridPar3,
        dataSize, dataField.data()
    );

    EXPECT_DOUBLE_EQ(value, 7.3999999999999986);

}

TEST(InterpolationTest, Interpolate3Dz) {

    long long dataSize1 = 10;
    long long dataSize2 = 7;
    long long dataSize3 = 5;
    long long dataSize = dataSize1*dataSize2*dataSize3;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    long long index2 = 3;
    double dgridPar3 = 2.5;
    int end_int = 0;
    double input1 = 2.4;
    double input2 = 5.3;
    double input3 = 1.5;
    std::vector<long long> gridDim = {dataSize1, dataSize2, dataSize3};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> gridPar3(dataSize3);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=0.5+i*dgridPar2;
    }
    for(auto i=0;i<dataSize3;i++){
        gridPar3[i]=1.0+i*dgridPar3;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            for(auto k=0;k<dataSize3;k++){
                dataField[i*dataSize2*dataSize3+j*dataSize3+k]=1.*i+2.*j+3.*k;
            }
        }
    }

    double value_end_int_0 = grackle::impl::fortran_wrapper::interpolate_3dz_g(
      input1, input2, input3,
      gridDim.data(),
      gridPar1.data(), dgridPar1,
      gridPar2.data(), index2,
      gridPar3.data(), dgridPar3,
      dataSize, dataField.data(),
      end_int
    );

    end_int = 1;
    double value_end_int_1 = grackle::impl::fortran_wrapper::interpolate_3dz_g(
      input1, input2, input3,
      gridDim.data(),
      gridPar1.data(), dgridPar1,
      gridPar2.data(), index2,
      gridPar3.data(), dgridPar3,
      dataSize, dataField.data(),
      end_int
    );

    EXPECT_DOUBLE_EQ(value_end_int_0, 7.3391950270214341);
    EXPECT_DOUBLE_EQ(value_end_int_1, 5.0);


}

TEST(InterpolationTest, Interpolate2Df3D) {

    long long dataSize1 = 10;
    long long dataSize2 = 7;
    long long dataSize3 = 5;
    long long dataSize = dataSize1*dataSize2*dataSize3;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    long long index2 = 1;
    double dgridPar3 = 2.5;
    double input1 = 2.4;
    double input3 = 1.5;
    std::vector<long long> gridDim = {dataSize1, dataSize2, dataSize3};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> gridPar3(dataSize3);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=0.5+i*dgridPar2;
    }
    for(auto i=0;i<dataSize3;i++){
        gridPar3[i]=1.0+i*dgridPar3;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            for(auto k=0;k<dataSize3;k++){
                dataField[i*dataSize2*dataSize3+j*dataSize3+k]=1.*i+2.*j+3.*k;
            }
        }
    }

    double value = grackle::impl::fortran_wrapper::interpolate_2df3d_g(
      input1, input3,
      gridDim.data(),
      gridPar1.data(), dgridPar1,
      index2,
      gridPar3.data(), dgridPar3,
      dataSize, dataField.data());

    EXPECT_DOUBLE_EQ(value, 0.99999999999999989);
}

TEST(InterpolationTest, Interpolate4D) {

    long long dataSize1 = 3;
    long long dataSize2 = 6;
    long long dataSize3 = 5;
    long long dataSize4 = 4;
    long long dataSize = dataSize1*dataSize2*dataSize3*dataSize4;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    double dgridPar3 = 2.5;
    double dgridPar4 = 3.5;
    double input1 = 2.4;
    double input2 = 5.3;
    double input3 = 1.5;
    double input4 = 5.7;
    std::vector<long long> gridDim = {dataSize1, dataSize2, dataSize3, dataSize4};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> gridPar3(dataSize3);
    std::vector<double> gridPar4(dataSize4);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=0.5+i*dgridPar2;
    }
    for(auto i=0;i<dataSize3;i++){
        gridPar3[i]=1.0+i*dgridPar3;
    }
    for(auto i=0;i<dataSize4;i++){
        gridPar4[i]=3.5+i*dgridPar4;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            for(auto k=0;k<dataSize3;k++){
                for(auto l=0;l<dataSize4;l++){
                    dataField[i*dataSize2*dataSize3*dataSize4+
                              j*dataSize3*dataSize4+k*dataSize4+l]=1.*i+2.*j+3.*k+4.*l;
                }
            }
        }
    }

    double value = grackle::impl::fortran_wrapper::interpolate_4d_g(
        input1, input2, input3, input4,
        gridDim.data(),
        gridPar1.data(), dgridPar1,
        gridPar2.data(), dgridPar2,
        gridPar3.data(), dgridPar3,
        gridPar4.data(), dgridPar4,
        dataSize, dataField.data()
    );

    EXPECT_DOUBLE_EQ(value, 9.9142857142857146);
}

TEST(InterpolationTest, Interpolate5D) {

    long long dataSize1 = 3;
    long long dataSize2 = 6;
    long long dataSize3 = 5;
    long long dataSize4 = 4;
    long long dataSize5 = 3;
    long long dataSize = dataSize1*dataSize2*dataSize3*dataSize4*dataSize5;
    double dgridPar1 = 1.0;
    double dgridPar2 = 1.5;
    double dgridPar3 = 2.5;
    double dgridPar4 = 3.5;
    double dgridPar5 = 3.0;
    double input1 = 2.4;
    double input2 = 5.3;
    double input3 = 1.5;
    double input4 = 5.7;
    double input5 = 2.7;
    std::vector<long long> gridDim = {dataSize1, dataSize2, dataSize3, dataSize4, dataSize5};
    std::vector<double> gridPar1(dataSize1);
    std::vector<double> gridPar2(dataSize2);
    std::vector<double> gridPar3(dataSize3);
    std::vector<double> gridPar4(dataSize4);
    std::vector<double> gridPar5(dataSize5);
    std::vector<double> dataField(dataSize);

    for(auto i=0;i<dataSize1;i++){
        gridPar1[i]=2.0+i*dgridPar1;
    }
    for(auto i=0;i<dataSize2;i++){
        gridPar2[i]=0.5+i*dgridPar2;
    }
    for(auto i=0;i<dataSize3;i++){
        gridPar3[i]=1.0+i*dgridPar3;
    }
    for(auto i=0;i<dataSize4;i++){
        gridPar4[i]=3.5+i*dgridPar4;
    }
    for(auto i=0;i<dataSize5;i++){
        gridPar5[i]=2.5+i*dgridPar5;
    }
    for(auto i=0;i<dataSize1;i++){
        for(auto j=0;j<dataSize2;j++){
            for(auto k=0;k<dataSize3;k++){
                for(auto l=0;l<dataSize4;l++){
                    for(auto m=0;m<dataSize5;m++){
                        dataField[i*dataSize2*dataSize3*dataSize4*dataSize5+
                                  j*dataSize3*dataSize4*dataSize5+
                                  k*dataSize4*dataSize5+l*dataSize5+m]=1.*i+2.*j+3.*k+4.*l+5.*m;
                    }
                }
            }
        }
    }

    double value = grackle::impl::fortran_wrapper::interpolate_5d_g(
        input1, input2, input3, input4, input5,
        gridDim.data(),
        gridPar1.data(), dgridPar1,
        gridPar2.data(), dgridPar2,
        gridPar3.data(), dgridPar3,
        gridPar4.data(), dgridPar4,
        gridPar5.data(), dgridPar5,
        dataSize, dataField.data()
    );

    EXPECT_DOUBLE_EQ(value, 10.247619047619049);
}
