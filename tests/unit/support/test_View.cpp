//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// add some basic tests of the View type
///
//===----------------------------------------------------------------------===//

#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <utility>  // std::move

#include "gtest/gtest.h"
#include "support/View.hpp"
#include "support/config.hpp"

TEST(View, Empty1D) {
  GRIMPL_NS::View<double*> empty_view;
  EXPECT_EQ(empty_view.data(), nullptr);
}

TEST(View, Empty2D) {
  GRIMPL_NS::View<float**> empty_view;
  EXPECT_EQ(empty_view.data(), nullptr);
}

TEST(View, Empty3D) {
  GRIMPL_NS::View<const int***> empty_view;
  EXPECT_EQ(empty_view.data(), nullptr);
}

TEST(View, Simple1D) {
  int arr[3] = {8, -3, 6};
  GRIMPL_NS::View<int*> v(arr, 3);

  EXPECT_EQ(arr, v.data());
  EXPECT_EQ(v.extent(0), 3);
  EXPECT_EQ(v(0), 8);
  EXPECT_EQ(v(1), -3);
  EXPECT_EQ(v(2), 6);
}

TEST(View, Simple1DConstCast) {
  int arr[3] = {8, -3, 6};
  GRIMPL_NS::View<int*> v(arr, 3);

  GRIMPL_NS::View<const int*> v_const = v;

  EXPECT_EQ(arr, v_const.data());
  EXPECT_EQ(v_const.extent(0), 3);
  EXPECT_EQ(v_const(0), 8);
  EXPECT_EQ(v_const(1), -3);
  EXPECT_EQ(v_const(2), 6);

  arr[0] = 0;
  v(2) = -6;
  EXPECT_EQ(v_const(0), 0);
  EXPECT_EQ(v_const(1), -3);
  EXPECT_EQ(v_const(2), -6);
}

TEST(View, Simple2D) {
  // clang-format off:  formatter knows nothing about shape
  double arr[20] = {
    0, 0, 0, 8, 0,
    0, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 0, 0, 17
  };
  // clang-format on

  GRIMPL_NS::View<double**> v(arr, 5, 4);

  EXPECT_EQ(arr, v.data());
  EXPECT_EQ(v.extent(0), 5);
  EXPECT_EQ(v.extent(1), 4);
  EXPECT_EQ(v(3, 0), 8);
  EXPECT_EQ(v(1, 2), 1);
  EXPECT_EQ(v(4, 3), 17);
}

TEST(View, Simple2DConstCast) {
  // clang-format off:  formatter knows nothing about shape
  double arr[20] = {
    0, 0, 0, 8, 0,
    0, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 0, 0, 17
  };
  // clang-format on

  GRIMPL_NS::View<double**> v(arr, 5, 4);
  GRIMPL_NS::View<const double**> v_const = v;

  EXPECT_EQ(arr, v_const.data());
  EXPECT_EQ(v_const.extent(0), 5);
  EXPECT_EQ(v_const.extent(1), 4);
  EXPECT_EQ(v_const(3, 0), 8);
  EXPECT_EQ(v_const(1, 2), 1);
  EXPECT_EQ(v_const(4, 3), 17);
}

TEST(View, Simple3D) {
  // clang-format off:  formatter knows nothing about shape
  float arr[60] = {
    0, 0, 0, 0, 0,
    0, 0, 0, 9, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    -3, 0, 0, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1
  };
  // clang-format on

  GRIMPL_NS::View<float***> v(arr, 5, 4, 3);

  EXPECT_EQ(arr, v.data());
  EXPECT_EQ(v.extent(0), 5);
  EXPECT_EQ(v.extent(1), 4);
  EXPECT_EQ(v.extent(2), 3);
  EXPECT_EQ(v(3, 1, 0), 9);
  EXPECT_EQ(v(0, 2, 1), -3);
  EXPECT_EQ(v(4, 3, 2), 1);
}

TEST(View, Simple3DConstCast) {
  // clang-format off:  formatter knows nothing about shape
  float arr[60] = {
    0, 0, 0, 0, 0,
    0, 0, 0, 9, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    -3, 0, 0, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 1
  };
  // clang-format on

  GRIMPL_NS::View<float***> v(arr, 5, 4, 3);
  GRIMPL_NS::View<const float***> v_const = v;

  EXPECT_EQ(arr, v_const.data());
  EXPECT_EQ(v_const.extent(0), 5);
  EXPECT_EQ(v_const.extent(1), 4);
  EXPECT_EQ(v_const.extent(2), 3);
  EXPECT_EQ(v_const(3, 1, 0), 9);
  EXPECT_EQ(v_const(0, 2, 1), -3);
  EXPECT_EQ(v_const(4, 3, 2), 1);
}

//===----------------------------------------------------------------------===//

// this acts like a "plugin" for View2DInterfaceTest
struct View2DManager {
  using element_type = double;
  using ViewType = GRIMPL_NS::View<element_type**>;

private:
  std::vector<std::unique_ptr<element_type[]>> ptrs_;
  std::vector<std::unique_ptr<ViewType>> views_;

public:
  ViewType& make_zero_filled(int extent0, int extent1) {
    int size = extent0 * extent1;
    element_type* p =
        ptrs_.emplace_back(std::make_unique<element_type[]>(size)).get();
    // I think this loop is unnecessary
    for (int i = 0; i < size; i++) {
      p[i] = element_type{0};
    }
    return *views_.emplace_back(
        std::make_unique<ViewType>(p, extent0, extent1));
  }
};

// this is a generic test of the interface for a 2D view
// -> I decided to make this a type-parameterized test for reasons that will
//    become apparent in a followup PR
// -> ManagerT constructs an instance of the ViewType and ensures that all
//    associated memory allocations will be freed when we go out of scope
template <typename ManagerT>
class View2DInterfaceTest : public testing::Test {
protected:
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  // the manager constructs an instance of ViewType and retains ownership
  // of that type and makes sure all memory is freed in the destructor
  ManagerT manager_;

  // make a view that is filled with zeros
  ViewType& make_zero_filled(int extent0, int extent1) {
    return manager_.make_zero_filled(extent0, extent1);
  }

  // make a view with a well specified shape and well-specified values
  ViewType& make_sample() {
    ViewType& out = manager_.make_zero_filled(3, 2);

    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 3; i++) {
        out(i, j) = static_cast<element_type>(j * 10 + i);
      }
    }
    return out;
  }

  // confirm that v has all the properties you would expect for the sample view
  // (before it is mutated)
  testing::AssertionResult is_unchanged_sample_view(const ViewType& v) const {
    int extents[2] = {v.extent(0), v.extent(1)};
    if (extents[0] != 3 || extents[1] != 2) {
      testing::AssertionFailure()
          << "v.extent(0) and v.extent(1) are " << extents[0] << " & "
          << extents[1] << "; should be 3 & 2";
    }

    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 3; i++) {
        element_type expected = static_cast<element_type>(j * 10 + i);
        if (v(i, j) != expected) {
          return testing::AssertionFailure()
                 << "v(" << i << ',' << j << ") is " << v(i, j)
                 << "; should be " << expected;
        }
      }
    }

    return testing::AssertionSuccess();
  }
};

using MyManagerTypes = ::testing::Types<View2DManager>;
// until we adopt C++20, we need to insert a trailing comma in order to suppress
// warnings about legacy C++20 variadic macro behavior
// https://github.com/google/googletest/issues/2271#issuecomment-665742471
TYPED_TEST_SUITE(View2DInterfaceTest, MyManagerTypes, );

TYPED_TEST(View2DInterfaceTest, ElementAssign) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  ViewType& v = this->make_zero_filled(2, 5);
  EXPECT_EQ(v.extent(0), 2);  // <- sanity check!
  EXPECT_EQ(v.extent(1), 5);  // <- sanity check!

  for (int j = 0; j < 5; j++) {
    for (int i = 0; i < 2; i++) {
      EXPECT_EQ(v(i, j), element_type{0});
    }
  }

  v(0, 2) = element_type{17};
  v(1, 0) = element_type{-3};
  v(1, 1) = element_type{4};
  v(0, 4) = element_type{-100};

  EXPECT_EQ(v(0, 0), element_type{0});
  EXPECT_EQ(v(1, 0), element_type{-3});
  EXPECT_EQ(v(0, 1), element_type{0});
  EXPECT_EQ(v(1, 1), element_type{4});
  EXPECT_EQ(v(0, 2), element_type{17});
  EXPECT_EQ(v(1, 2), element_type{0});
  EXPECT_EQ(v(0, 3), element_type{0});
  EXPECT_EQ(v(1, 3), element_type{0});
  EXPECT_EQ(v(0, 4), element_type{-100});
  EXPECT_EQ(v(1, 4), element_type{0});
}

// this is a quick test to make sure our make_sample helper method does what
// we think it does
TYPED_TEST(View2DInterfaceTest, SampleSanityCheck) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;

  ViewType& v = this->make_sample();
  EXPECT_TRUE(this->is_unchanged_sample_view(v));
}

TYPED_TEST(View2DInterfaceTest, CopyConstruct) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  ViewType& v_orig = this->make_sample();
  ViewType v2(v_orig);

  EXPECT_TRUE(this->is_unchanged_sample_view(v2));

  v_orig(1, 1) = element_type{-29};
  EXPECT_EQ(v2(1, 1), element_type{-29});
  v2(0, 0) = -9;
  EXPECT_EQ(v_orig(0, 0), -9);
}

TYPED_TEST(View2DInterfaceTest, CopyConstructEmpty) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;

  ViewType empty_view;
  ViewType v2(empty_view);
  EXPECT_EQ(v2.data(), nullptr);
}

TYPED_TEST(View2DInterfaceTest, MoveConstruct) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  ViewType& v_orig = this->make_sample();
  ViewType v_tmp(v_orig);
  ViewType v2(std::move(v_tmp));
  // the current state of v_tmp is intentionally left unspecified (by us)

  EXPECT_TRUE(this->is_unchanged_sample_view(v2));

  v_orig(1, 1) = element_type{-29};
  EXPECT_EQ(v2(1, 1), element_type{-29});
  v2(0, 0) = -9;
  EXPECT_EQ(v_orig(0, 0), -9);
}

TYPED_TEST(View2DInterfaceTest, MoveConstructEmpty) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;

  ViewType empty_view;
  ViewType v2(std::move(empty_view));
  EXPECT_EQ(v2.data(), nullptr);
}

TYPED_TEST(View2DInterfaceTest, CopyAssign) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  ViewType& v = this->make_sample();
  ViewType v_copy_of_orig(v);
  ViewType& v_alt = this->make_zero_filled(3, 4);

  // a quick sanity check!
  EXPECT_EQ(v_alt.extent(0), 3);
  EXPECT_EQ(v_alt.extent(1), 4);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(v_alt(i, j), element_type{0});
    }
  }

  // perform the copy assignment
  v = v_alt;
  // check that v is now a copy of v_alt
  EXPECT_EQ(v.extent(0), 3);
  EXPECT_EQ(v.extent(1), 4);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(v(i, j), element_type{0});
    }
  }

  // demonstrate that v and v_alt reference the same memory
  v_alt(1, 1) = 91;
  EXPECT_EQ(v(1, 1), 91);
  v(0, 2) = 14;
  EXPECT_EQ(v_alt(0, 2), 14);

  // confirm that v_copy_of_orig was totally unaffected by all of this
  EXPECT_TRUE(this->is_unchanged_sample_view(v_copy_of_orig));
}

TYPED_TEST(View2DInterfaceTest, CopyAssignEmpty) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;

  ViewType& v = this->make_sample();
  ViewType v_copy_of_orig(v);
  ViewType empty_view;

  // perform copy assignment
  v = empty_view;

  // check that v is now like empty_view
  EXPECT_EQ(v.data(), nullptr);

  // confirm that v_copy_of_orig was totally unaffected by all of this
  EXPECT_TRUE(this->is_unchanged_sample_view(v_copy_of_orig));
}

TYPED_TEST(View2DInterfaceTest, MoveAssign) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;
  using element_type = typename ManagerT::element_type;

  ViewType& v = this->make_sample();
  ViewType v_copy_of_orig(v);
  ViewType& v_alt = this->make_zero_filled(3, 4);

  // a quick sanity check!
  EXPECT_EQ(v_alt.extent(0), 3);
  EXPECT_EQ(v_alt.extent(1), 4);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(v_alt(i, j), element_type{0});
    }
  }

  // perform the move assignment
  // -> to make the test slightly more robust, we move from a temporary copy of
  //    v_alt (since we haven't made a firm choice about the state of a
  //    variable on the RHS of move-assignment after the operation)
  {
    ViewType tmp(v_alt);
    v = std::move(tmp);
  }
  // check that v is now a copy of v_alt
  EXPECT_EQ(v.extent(0), 3);
  EXPECT_EQ(v.extent(1), 4);
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(v(i, j), element_type{0});
    }
  }

  // demonstrate that v and v_alt reference the same memory
  v_alt(1, 1) = 91;
  EXPECT_EQ(v(1, 1), 91);
  v(0, 2) = 14;
  EXPECT_EQ(v_alt(0, 2), 14);

  // confirm that v_copy_of_orig was totally unaffected by all of this
  EXPECT_TRUE(this->is_unchanged_sample_view(v_copy_of_orig));
}

TYPED_TEST(View2DInterfaceTest, MoveAssignEmpty) {
  using ManagerT = TypeParam;
  using ViewType = typename ManagerT::ViewType;

  ViewType& v = this->make_sample();
  ViewType v_copy_of_orig(v);
  ViewType empty_view;

  // perform move assignment
  v = std::move(empty_view);

  // check that v is now like empty_view
  EXPECT_EQ(v.data(), nullptr);

  // confirm that v_copy_of_orig was totally unaffected by all of this
  EXPECT_TRUE(this->is_unchanged_sample_view(v_copy_of_orig));
}