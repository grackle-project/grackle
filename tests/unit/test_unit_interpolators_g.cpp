#include "gtest/gtest.h"
#include <cstdint>

#include <initializer_list>
#include <utility>  // std::pair

#include "grtestutils/coord_grid.hpp"
#include "interp_grid.hpp"
#include "interpolate.hpp"
#include "utils-cpp.hpp"  // GRIMPL_NS::clamp


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

    
    double value = GRIMPL_NS::interpolate_1d(
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


    double value = GRIMPL_NS::interpolate_2d(
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

    double value = GRIMPL_NS::interpolate_3d(
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
    std::int64_t end_int = 0;
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

    double value_end_int_0 = GRIMPL_NS::interpolate_3dz(
      input1, input2, input3,
      gridDim.data(),
      gridPar1.data(), dgridPar1,
      gridPar2.data(), index2,
      gridPar3.data(), dgridPar3,
      dataSize, dataField.data(),
      end_int
    );

    end_int = 1;
    double value_end_int_1 = GRIMPL_NS::interpolate_3dz(
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

    double value = GRIMPL_NS::interpolate_2Df3D(
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

    double value = GRIMPL_NS::interpolate_4d(
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

    double value = GRIMPL_NS::interpolate_5d(
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

// =============================================================================
// Down below we move onto making more exhaustive property tests
// =============================================================================

static double perform_interp(const grtest::Coordinate& c,
                             const GRIMPL_NS::InterpGrid& interp_grid) {
  const GRIMPL_NS::InterpGridProps& props = interp_grid.props;
  switch (props.rank) {
    case 1:
      return GRIMPL_NS::interpolate_1d(
          c.components[0], props.dimension, props.parameters[0],
          props.parameter_spacing[0], props.data_size, interp_grid.data);
    case 2:
      return GRIMPL_NS::interpolate_2d(
          c.components[0], c.components[1], props.dimension,
          props.parameters[0], props.parameter_spacing[0], props.parameters[1],
          props.parameter_spacing[1], props.data_size, interp_grid.data);
    case 3:
      return GRIMPL_NS::interpolate_3d(
          c.components[0], c.components[1], c.components[2], props.dimension,
          props.parameters[0], props.parameter_spacing[0], props.parameters[1],
          props.parameter_spacing[1], props.parameters[2],
          props.parameter_spacing[2], props.data_size, interp_grid.data);
    case 4:
      return GRIMPL_NS::interpolate_4d(
          c.components[0], c.components[1], c.components[2], c.components[3],
          props.dimension, props.parameters[0], props.parameter_spacing[0],
          props.parameters[1], props.parameter_spacing[1], props.parameters[2],
          props.parameter_spacing[2], props.parameters[3],
          props.parameter_spacing[3], props.data_size, interp_grid.data);
    case 5:
      return GRIMPL_NS::interpolate_5d(
          c.components[0], c.components[1], c.components[2], c.components[3],
          c.components[4], props.dimension, props.parameters[0],
          props.parameter_spacing[0], props.parameters[1],
          props.parameter_spacing[1], props.parameters[2],
          props.parameter_spacing[2], props.parameters[3],
          props.parameter_spacing[3], props.parameters[4],
          props.parameter_spacing[4], props.data_size, interp_grid.data);
    default:
      GRIMPL_ERROR("unexpected rank: %d", (int)props.rank);
  }
}

template <typename Fn>
static std::pair<GRIMPL_NS::InterpGrid, grtest::CoordGrid> setup(
    const std::vector<GRIMPL_NS::InterpDimScale>& dim_scales, Fn fn) {
  int rank = dim_scales.size();
  GRIMPL_NS::InterpGridProps grid_props(rank, dim_scales.data());
  std::vector<std::vector<double>> coords;
  for (int i = 0; i < rank; i++) {
    int n_elements = grid_props.dimension[i];
    std::vector<double>& cur_ax = coords.emplace_back(n_elements);
    for (int j = 0; j < n_elements; j++) {
      cur_ax[j] = grid_props.parameters[i][j];
    }
  }

  grtest::CoordGrid coord_grid(coords);

  // set up the data_field
  double* data_field = new double[coord_grid.size()];
  int i = 0;
  for (const grtest::Coordinate& coord : coord_grid) {
    data_field[i] = fn(coord);
    // std::string s = testing::PrintToString(coord);
    // printf("  fn(%s)=%g, \n", s.c_str(), data_field[i]);
    // fflush(stdout);
    i++;
  }

  GRIMPL_NS::InterpGrid interp_grid(std::move(grid_props), data_field);
  return {std::move(interp_grid), coord_grid};
}

/// This namespace aggregates different "policies" that are used as parameters
/// for the @ref InterpGridTest test fixture.
///
/// In slightly more detail, each policy is a distinct type that defines a
/// function and a set of grid dimension scales, from which the tests
/// associated with @ref InterpGridTest will construct an interpolation grid.
/// And then the tests verify that interpolations can reproduce the initial
/// values
///
/// See the docstring of @ref InterpGridTest for more details
namespace InterpTestPolicies {

struct Simple1D {
  // implements the function that will be interpolated by the test
  static double fn(const grtest::Coordinate& c) { return c[0] + 5.0; };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {GRIMPL_NS::InterpDimScale::Linear(2, 0, 1.0)};
  };
};

struct Piecewise1D {
  // implements the function that will be interpolated by the test
  static double fn(const grtest::Coordinate& c) {
    return (c[0] < 1.0) ? c[0] + 5.0 : (c[0] - 0.25) * 2 + 5.25;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    // important: make sure the interpolation grid includes the breakpoint of
    // the interpolation grid
    return {GRIMPL_NS::InterpDimScale::Linear(4, 0, 0.25)};
  };
};

struct Simple2D {
  // implements the function that will be interpolated in the test
  static double fn(const grtest::Coordinate& c) {
    return c[0] * 0.5 + c[1] * -1.0 + 5.0;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {
        GRIMPL_NS::InterpDimScale::Linear(3, 0, 1.0),
        GRIMPL_NS::InterpDimScale::Linear(3, -3.0, 0.5),
    };
  }
};

struct Alt2D {
  // implements the function that will be interpolated in the test
  static double fn(const grtest::Coordinate& c) {
    return 100.0 * c[0] * c[1] + c[0] * 0.5 + c[1] * -1.0 + 5.0;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {
        GRIMPL_NS::InterpDimScale::Linear(3, 0, 1.0),
        GRIMPL_NS::InterpDimScale::Linear(3, -3.0, 0.5),
    };
  }
};

struct Piecewise2D {
  // implements the function that will be interpolated in the test
  static double fn(const grtest::Coordinate& c) {
    double bp0 = 1.0;
    double bp1 = 0.5;

    double a = (c[0] < bp0) ? c[0] * 4 : (c[0] - bp0) + bp0 * 4;
    double b = (c[1] < bp1) ? c[1] * -0.5 : (c[1] - bp1) + bp1 * -0.5;

    return a + b;
  }

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    // important: make sure the interpolation grid includes the breakpoints of
    // the interpolation grid
    return {
        GRIMPL_NS::InterpDimScale::Linear(3, 0, 1.0),
        GRIMPL_NS::InterpDimScale::Linear(5, -1.0, 0.5),
    };
  }
};

struct Simple3D {
  // implements the function that will be interpolated during the test
  static double fn(const grtest::Coordinate& c) {
    return c[0] * 0.01 + c[1] * -1.0 + c[2] * 100.0 + 5.0;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {GRIMPL_NS::InterpDimScale::Linear(4, 0, 1.0),
            GRIMPL_NS::InterpDimScale::Linear(4, -3.0, 0.5),
            GRIMPL_NS::InterpDimScale::Linear(4, -1.0, 0.25)};
  }
};

struct Simple4D {
  // implements the function that will be interpolated during the test
  static double fn(const grtest::Coordinate& c) {
    return c[0] * 0.1 + c[1] * -1.0 + c[2] * 10.0 + c[3] * -100.0 + 5.0;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {
        GRIMPL_NS::InterpDimScale::Linear(2, 0, 1.0),
        GRIMPL_NS::InterpDimScale::Linear(2, -3.0, 0.5),
        GRIMPL_NS::InterpDimScale::Linear(2, -1.0, 0.25),
        GRIMPL_NS::InterpDimScale::Linear(2, 1.0, 0.5),
    };
  }
};

struct Simple5D {
  // implements the function that will be interpolated during the test
  static double fn(const grtest::Coordinate& c) {
    return c[0] * 0.1 + c[1] * -1.0 + c[2] * 10.0 + c[3] * -100.0 + 5.0;
  };

  // returns the dimension scales for producing the interpolation grid
  static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales() {
    return {
        GRIMPL_NS::InterpDimScale::Linear(2, 0, 1.0),
        GRIMPL_NS::InterpDimScale::Linear(2, -3.0, 0.5),
        GRIMPL_NS::InterpDimScale::Linear(2, -1.0, 0.25),
        GRIMPL_NS::InterpDimScale::Linear(2, 1.0, 0.5),
        GRIMPL_NS::InterpDimScale::Linear(2, -2.0, 2.0),
    };
  }
};

}  // namespace InterpTestPolicies

/// A fixture class template that defines logic for testing different
/// interpolation schemes on grids of values
///
/// We refer to the actual parameter as a "test-policy". A policy is a
/// unique type that each policy that defines a function and a set of grid
/// dimension scales. The tests associated with this fixture will use this
/// information to construct an interpolation grid for interpolating values of
/// that function. The tests will verify that interpolations can reproduce the
/// behavior of the function.
///
/// IMPORTANTLY, we must be able to exactly reproduce the value of the
/// function.
///
/// At this time, a policy should implement the following interface
/// @code{cpp}
///   struct MyPolicy {
///     // implements the function that is modelled by the tests
///     static double fn(const grtest::Coordinate& c);
///
///     // returns the dimension scales for producing the interpolation grid
///     static std::vector<GRIMPL_NS::InterpDimScale> grid_dim_scales();
///   };
/// @endcode
template <typename TestPolicy>
class InterpGridTest : public testing::Test {};

using MyPolicies = ::testing::Types<
    InterpTestPolicies::Simple1D, InterpTestPolicies::Piecewise1D,
    InterpTestPolicies::Simple2D, InterpTestPolicies::Alt2D,
    InterpTestPolicies::Piecewise2D, InterpTestPolicies::Simple3D,
    InterpTestPolicies::Simple4D, InterpTestPolicies::Simple5D>;
// the trailing comma on the next line can be removed once we start using c++20
TYPED_TEST_SUITE(InterpGridTest, MyPolicies, );

TYPED_TEST(InterpGridTest, AtGridPoints) {
  using Policy = TypeParam;
  auto my_fn = [](const grtest::Coordinate& c) -> double {
    return Policy::fn(c);
  };
  auto [interp_grid, coord_grid] = setup(Policy::grid_dim_scales(), my_fn);

  ASSERT_TRUE(interp_grid) << "error setting up coord grid";

  for (const grtest::Coordinate& c : coord_grid) {
    // allows a discrepancy of 4 ULPs
    EXPECT_DOUBLE_EQ(perform_interp(c, interp_grid), my_fn(c))
        << "for c = " << testing::PrintToString(c);
  }
}

TYPED_TEST(InterpGridTest, BetweenGridPoints) {
  using Policy = TypeParam;
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales = Policy::grid_dim_scales();
  auto my_fn = [](const grtest::Coordinate& c) -> double {
    return Policy::fn(c);
  };
  auto [interp_grid, _] = setup(Policy::grid_dim_scales(), my_fn);

  ASSERT_TRUE(interp_grid) << "error setting up coord grid";

  // construct a new coordinate grid where we query locations between each point
  std::vector<std::vector<double>> test_locs;
  for (const GRIMPL_NS::InterpDimScale& dim_scale : dim_scales) {
    std::vector<double>& v = test_locs.emplace_back();
    std::size_t dim_scale_count = dim_scale.count;
    v.reserve(2 * (dim_scale_count - 1));
    for (std::size_t i = 0; i < (dim_scale_count - 1); i++) {
      double left = dim_scale.access_value(i);
      double right = dim_scale.access_value(i + 1);
      v.push_back(0.75 * left + 0.25 * right);
      v.push_back(0.25 * left + 0.75 * right);
    }
  }

  grtest::CoordGrid check_coord_grid(test_locs);

  for (const grtest::Coordinate& c : check_coord_grid) {
    // allows a discrepancy of 4 ULPs
    EXPECT_DOUBLE_EQ(perform_interp(c, interp_grid), my_fn(c))
        << "for c = " << testing::PrintToString(c);
  }
}

// Tests what happens when we interpolate beyond the edge of the grid
// - at the time of writing, this simply extrapolates
// - in principle, we should probably default to snapping to the nearest
//   interpolated value within the domain of allowed points. We actually take
//   some steps to manually do this for certain quantities and in certain
//   locations of the codebase, but I'm a little skeptical that we do it
//   consistently enough
TYPED_TEST(InterpGridTest, BeyondGrid) {
  using Policy = TypeParam;
  std::vector<GRIMPL_NS::InterpDimScale> dim_scales = Policy::grid_dim_scales();
  auto my_fn = [](const grtest::Coordinate& c) -> double {
    return Policy::fn(c);
  };
  auto [interp_grid, _] = setup(Policy::grid_dim_scales(), my_fn);

  ASSERT_TRUE(interp_grid) << "error setting up coord grid";

  // construct a new coordinate grid where we query locations between each point
  // (we will also keep track of the most extreme values along each axis)
  constexpr int MAX_RANK = 5;
  double ax_extremes[MAX_RANK][2];
  std::vector<std::vector<double>> test_locs;

  int grid_rank = dim_scales.size();
  ASSERT_LE(grid_rank, MAX_RANK) << "sanity check failed";
  for (int i = 0; i < grid_rank; i++) {
    const GRIMPL_NS::InterpDimScale& dim_scale = dim_scales[i];
    std::size_t dim_scale_count = dim_scale.count;
    double leftmost = dim_scale.access_value(0);
    double leftmost_step = dim_scale.access_value(1) - leftmost;
    double rightmost = dim_scale.access_value(dim_scale_count - 1);
    double rightmost_step =
        rightmost - dim_scale.access_value(dim_scale_count - 2);

    std::vector<double> v{leftmost - 0.25 * leftmost_step, leftmost,
                          leftmost + 0.25 * leftmost_step, rightmost,
                          rightmost + 0.25 * rightmost_step};

    test_locs.emplace_back(std::move(v));
    ax_extremes[i][0] = leftmost;
    ax_extremes[i][1] = rightmost;
  }

  grtest::CoordGrid check_coord_grid(test_locs);

  for (const grtest::Coordinate& c : check_coord_grid) {
    grtest::Coordinate clamped_coord;
    clamped_coord.rank = grid_rank;
    for (int i = 0; i < grid_rank; i++) {
      clamped_coord.components[i] =
          GRIMPL_NS::clamp(c[i], ax_extremes[i][0], ax_extremes[i][1]);
    }

    // allows a discrepancy of 4 ULPs
    EXPECT_DOUBLE_EQ(perform_interp(c, interp_grid), my_fn(c))
        << "for c = " << testing::PrintToString(c)
        << ". Nearest location in domain: "
        << testing::PrintToString(clamped_coord);
  }
}