#include <cmath>
#include <vector>

#include "gaussj_g.hpp"

namespace grackle::impl {

int gaussj_g(int n, double* coef_matrix_fortran, double* vector) {

    // TODO: to be removed
    // Copy the matrix to a C-style layout (column-major order)
    std::vector<double> coef_matrix(n*n);
    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
        coef_matrix[j*n + i] = coef_matrix_fortran[i*n + j];
        }
    }

    // Loop over columns
    for(unsigned int col=0; col<n; col++) {

        // Find the largest (absolute value) element in the column
        double max_val = std::fabs(coef_matrix[col*n + col]);
        unsigned int max_row = col;
        for(unsigned int row = col+1; row < n; row++) {
            if(std::fabs(coef_matrix[row*n + col]) > max_val) {
                max_val = std::fabs(coef_matrix[row*n + col]);
                max_row = row;
            }
        }

        // Check for singularity
        if(max_val < eps_zero) {
            return 1;
        }

        // Do partial pivoting by swaping rows
        if(max_row != col) {
            for(unsigned int col_p=0; col_p<n; col_p++) {
                std::swap(coef_matrix[col*n + col_p], coef_matrix[max_row*n + col_p]);
            }
            std::swap(vector[col], vector[max_row]);
        }

        // Scale the pivot row with pivot value (optional, this should increase stability)
        const double pivot_val = coef_matrix[col*n + col];
        coef_matrix[col*n + col] = 1.0;
        for(unsigned int col_p = col+1; col_p<n; col_p++) {
            coef_matrix[col*n + col_p] /= pivot_val;
        }
        vector[col] /= pivot_val;

        // Eliminate elements of rows below the pivot
        for(unsigned int row = col+1; row<n; row++){
            const double scaling_factor = -coef_matrix[row*n+col]/coef_matrix[col*n + col];
            coef_matrix[row*n + col] = 0.0;
            for(unsigned int col_p = col+1; col_p<n; col_p++) {
                coef_matrix[row*n + col_p] += scaling_factor*coef_matrix[col*n + col_p];
            }
            vector[row] += scaling_factor*vector[col];
        }

    }

    // Compute back substitution
    for(unsigned int row=n-1; row<n; row--) {
        for(unsigned int col=row+1; col<n; col++) {
            vector[row] -= coef_matrix[row*n +col]*vector[col];
        }
        vector[row] /= coef_matrix[row*n + row];
    }

    return 0;
}

} // namespace grackle::impl
