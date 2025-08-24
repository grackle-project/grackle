// See LICENSE file for license and copyright information

/// @file gaussj_g.hpp
/// @brief Implementation of Gauss-Jordan elimination for solving linear systems


#ifndef GAUSSJ_G_HPP
#define GAUSSJ_G_HPP

namespace grackle::impl {

// Zero threshold for numerical stability
const double eps_zero = 1e-20;

/// Performs Gauss-Jordan elimination to solve the specified system of linear
/// equations and inverts the coefficient matrix.
///
/// In more detail, it solves the linear matrix `ax=b`, where `a` is the
/// square coefficient matrix, `b` is the right-hand side vector and `x`
/// is the solution vector
///
/// @param[in]     n The number of linear equations being solved
/// @param[in,out] coef_matrix An n by n array that initially specifies the
///    coefficient matrix. It's overwritten by the inverted matrix.
/// @param[in,out] vec An n element array that initially specifies the
///    right-hand side vector. It's overwritten by the solution vector.
///
/// @retval 0 indicates success
/// @retval 1 indicates that the matrix is singular
///
int gaussj_g(int n, double* coef_matrix, double* vector);

} // namespace grackle::impl

#endif // GAUSSJ_G_HPP