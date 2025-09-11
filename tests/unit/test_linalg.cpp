#include <vector>
#include <gtest/gtest.h>

#include "fortran_func_wrappers.hpp"
#include "utest_helpers.hpp"


/// Records the paramters for a linear algebra test-case
struct LinAlgCase {
  // attributes
  int n;
  std::vector<double> matrix_rowmajor;
  std::vector<double> rhs_vector;
  std::vector<double> solution_vector;

  // interface methods:
  LinAlgCase() = delete;

  /// Construct the test-case
  LinAlgCase(int n, std::vector<double> matrix_rowmajor,
             std::vector<double> rhs_vector,
             std::vector<double> solution_vector)
    : n(n), matrix_rowmajor(matrix_rowmajor), rhs_vector(rhs_vector),
      solution_vector(solution_vector)
  { }

};

static std::vector<double> transpose_matrix(
  int n_row, int n_col, std::vector<double> matrix
) {
  std::vector<double> out(n_row*n_col);

  for (int row_idx = 0; row_idx < n_row; row_idx++) {
    for (int col_idx = 0; col_idx < n_col; col_idx++) {
      int row_major_ind = row_idx + col_idx * n_row;
      int col_major_ind = col_idx + row_idx * n_row;
      // correctness is unaffected by whether matrix is row-major or col-major
      // (but performance is affected)
      out[col_major_ind] = matrix[row_major_ind];
    }
  }

  return out;
}

class LinAlgTestSolve : public testing::TestWithParam<LinAlgCase> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};

// Solve a basic linear algebra equation.
TEST_P(LinAlgTestSolve, CheckSuccessfulSolve) {
  LinAlgCase my_case = GetParam();
  std::vector<double> matrix_colmajor = transpose_matrix(
    my_case.n, my_case.n, my_case.matrix_rowmajor
  );
  std::vector<double> vec = my_case.rhs_vector;
  int rslt = grackle::impl::fortran_wrapper::gaussj_g(
    my_case.n, matrix_colmajor.data(), vec.data()
  );
  ASSERT_EQ(rslt, 0) << "expected a return-code of 0, which indicates "
                     << "that the linear equations were successfully solved";

  EXPECT_TRUE(check_allclose(/* actual: */ vec, my_case.solution_vector,
                             /* rtol: */1e-15, /* atol: */ 0.0));
}

INSTANTIATE_TEST_SUITE_P(
  LinAlgScenarios,
  LinAlgTestSolve,
  testing::Values(
    LinAlgCase(
      2,
      /* row-major matrix */ std::vector<double>{ 3.0, 4.0,
                                                 -6.0, 9.0},
      /* rhs-vector */       std::vector<double>{41.0, 3.0},
      /* solution-vector */  std::vector<double>{ 7.0, 5.0}
    ),
    LinAlgCase(
      3,
      /* row-major matrix */ std::vector<double>{2.0,  1.0, -1.0,
                                                -3.0, -1.0,  2.0,
                                                -2.0,  1.0,  2.0},
      /* rhs-vector */       std::vector<double>{8.0, -11.0, -3.0},
      /* solution-vector */  std::vector<double>{2.0, 3.0, -1.0}
    ),
    LinAlgCase(
      4,
      /* row-major matrix */ std::vector<double>{ 1.0,  4.0,  3.0, -2.0,
                                                  2.0,  3.0, -1.0,  4.0,
                                                 -1.0,  1.0, -1.0,  1.5,
                                                  2.0,  1.0,  0.0, -1.0},
      /* rhs-vector */       std::vector<double>{-8.0, 37.0, 12.5, -3.0},
      /* solution-vector */  std::vector<double>{ 1.0,  2.0, -1.0,  7.0}
    )
  )
);


TEST(LinAlgTestSolveSingular, SillyScenario) {
  int n = 4;
  std::vector<double> matrix_rowmajor{ 1.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0, 0.0};
  std::vector<double> matrix_colmajor = transpose_matrix(n, n, matrix_rowmajor);

  // it doesn't really matter what the values are in rhs_vec
  std::vector<double> rhs_vec = {0.0, 0.0, 0.0, 0.0};

  int rslt = grackle::impl::fortran_wrapper::gaussj_g(
    n, matrix_colmajor.data(), rhs_vec.data()
  );
  ASSERT_EQ(rslt, 1) << "expected a return-code of 1, to indicate that "
                     << "the matrix is singular";
}

TEST(LinAlgTestSolveSingular, AltScenario) {
  // we are picking a matrix without so many zeros (it is singular since the
  // determinant is 0)
  int n = 2;
  std::vector<double> matrix_rowmajor{ -2.0, 4.0,
                                        3.0,-6.0};
  std::vector<double> matrix_colmajor = transpose_matrix(n, n, matrix_rowmajor);

  // it doesn't really matter what the values are in rhs_vec
  std::vector<double> rhs_vec = {0.0, 0.0};

  int rslt = grackle::impl::fortran_wrapper::gaussj_g(
    n, matrix_colmajor.data(), rhs_vec.data()
  );
  ASSERT_EQ(rslt, 1) << "expected a return-code of 1, to indicate that "
                     << "the matrix is singular";
}



