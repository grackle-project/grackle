#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "gauss_j.hpp"

int gaussj_g_cpp(int n, double* coef_matrix, double* vector, double* solution) {

    // Find pivot position of each row
    std::vector<unsigned int> pivot(n,n-1);
    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            if(std::abs(coef_matrix[i*n+j]) > eps_zero) {
                pivot[i] = j;
                break;
            }
            if(pivot[i] == 0) {
                std::cerr << "Matrix is singular." << std::endl;
                return -1;
            }
        }
    }

#ifdef DEBUG
    // Debug: print pivot positions
    std::cout << "Pivot positions: ";
    for(unsigned int i = 0; i < n; ++i) {
        std::cout << pivot[i] << " ";
    }
    std::cout << std::endl;
#endif

    // Find pivot based ordering of coef_matrix rows 
    std::vector<unsigned int> sortedRows(n);
    std::iota(sortedRows.begin(), sortedRows.end(), 0);
    std::sort(sortedRows.begin(), sortedRows.end(), [&](unsigned int i, unsigned int j){
        return pivot[i] < pivot[j];
    });

#ifdef DEBUG
    // Debug: print sorted rows
    std::cout << std::endl;
    std::cout << "Sorted rows: ";
    for(unsigned int i = 0; i < n; ++i) {
        std::cout << sortedRows[i] << " ";
    }
    std::cout << std::endl;
#endif

    // Reorder matrix, vector and pivot (TODO: debug only, not mandatory)
    std::vector<double> reord_coef_matrix(n*n);
    std::vector<double> reord_vector(n);
    std::vector<double> pivotReord(n);
    std::vector<unsigned int> sortedRowsReord(n);
    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            reord_coef_matrix[i*n+j] = coef_matrix[sortedRows[i]*n+j];
        }
    }
    for(unsigned int i = 0; i < n; ++i) {
        reord_vector[i] = vector[sortedRows[i]];
    }
    for(unsigned int i = 0; i < n; ++i) {
        reord_vector[i] = vector[sortedRows[i]];
    }
    for(unsigned int i = 0; i < n; ++i) {
        pivotReord[i] = pivot[sortedRows[i]];
    }
    for(unsigned int i = 0; i < n; ++i) {
        sortedRowsReord[i] = sortedRows[i];
    }

#ifdef DEBUG
    std::cout << std::endl;
    std::cout << "Reordered matrix:" << std::endl;
    for(unsigned int i = 0; i < n; ++i) {
        for(unsigned int j = 0; j < n; ++j) {
            std::cout << reord_coef_matrix[i*n + j] << " ";
        }
        std::cout << "| " << reord_vector[i] << std::endl;
    }
    std::cout << std::endl;
#endif

    // Perform Gaussian elimination and transform the matrix in echelon form
    bool rowScaled = true;
    while(rowScaled) {    
        std::cout<<std::endl;
        rowScaled = false;
        for(unsigned int i = n-1; i>0; i--){
            if(pivotReord[i] == pivotReord[i-1]) {
                const double scalingFactor = - reord_coef_matrix[i*n + pivotReord[i]]/reord_coef_matrix[(i-1)*n + pivotReord[i-1]];
                reord_coef_matrix[i*n+pivotReord[i]] = 0.0;
                for(unsigned int j=pivotReord[i]+1; j<n; j++) {
                    reord_coef_matrix[i*n + j] += scalingFactor * reord_coef_matrix[(i-1)*n + j];
                    if(std::abs(reord_coef_matrix[i*n + j]) < eps_zero) {
                        reord_coef_matrix[i*n + j] = 0.0;
                        // pivotReord[i] += 1;// TODO: can be done here?
                    }
                }
                reord_vector[i] += scalingFactor * reord_vector[i-1];
                if(std::abs(reord_vector[i]) < eps_zero) {
                    reord_vector[i] = 0.0;
                }
                rowScaled = true;
            }
        }


        if(rowScaled) {


#ifdef DEBUG
            // Debug: print the matrix after scaling
            std::cout << std::endl;
            std::cout << "Matrix after scaling:" << std::endl;
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = 0; j < n; ++j) {
                    std::cout << reord_coef_matrix[i*n + j] << " ";
                }
                std::cout << "| " << reord_vector[i] << std::endl;
            }
            std::cout<<std::endl;
#endif

            // Update pivot positions after scaling
            for(unsigned int i = 0; i < n; ++i) {
                pivot[i] = n-1; // reset pivot positions
            }
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = 0; j < n; ++j) {
                    if(std::abs(reord_coef_matrix[i*n+j]) > eps_zero) {
                        pivot[i] = j;
                        break;
                    }
                    if(pivot[i] == 0) {
                        std::cerr << "Matrix is singular or not square." << std::endl;
                        return -1;
                    }
                }
            }

#ifdef DEBUG
            std::cout << std::endl;
            std::cout << "Updated pivot positions: ";
            for(unsigned int i = 0; i < n; ++i) {
                std::cout << pivot[i] << " ";
            }
            std::cout << std::endl;
#endif 
            std::iota(sortedRows.begin(), sortedRows.end(), 0);
            std::sort(sortedRows.begin(), sortedRows.end(), [&](unsigned int i, unsigned int j){
                return pivot[i] < pivot[j];
            });

#ifdef DEBUG
            std::cout << std::endl;
            std::cout << "Updated sorted rows: ";
            for(unsigned int i = 0; i < n; ++i) {
                std::cout << sortedRows[i] << " ";
            }
            std::cout << std::endl; 
#endif 

            std::vector<double> reord_coef_matrix_tmp(reord_coef_matrix);
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = 0; j < n; ++j) {
                    reord_coef_matrix[i*n+j] = reord_coef_matrix_tmp[sortedRows[i]*n+j];
                }
            }
            std::vector<double> reord_vector_tmp(reord_vector);
            for(unsigned int i = 0; i < n; ++i) {
                reord_vector[i] = reord_vector_tmp[sortedRows[i]];
            }
            for(unsigned int i = 0; i < n; ++i) {
                pivotReord[i] = pivot[sortedRows[i]];
            }
            // TODO: avoid memory duplication
            std::vector<unsigned int> sortedRowsTmp(sortedRowsReord);
            for(unsigned int i = 0; i < n; ++i) {
                sortedRowsReord[i] = sortedRowsTmp[sortedRows[i]];
            }
        }
    }    

    // Perform back substitution to solve the system
    for(int i = n-1; i >= 0; --i) {
        std::cout << "Solving for solution[" << i << "]" << std::endl;
        solution[i] = reord_vector[i];
        std::cout << "solution[" << i << "] = " << solution[i] << std::endl;
        for(unsigned int j = i+1; j<n; ++j) {
            solution[i] -= reord_coef_matrix[i*n + j]*solution[j];
            std::cout << "solution[" << i << "] = " << solution[i] << std::endl;
        }
        solution[i] /= reord_coef_matrix[i*n + i]; 
        std::cout << "solution[" << i << "] = " << solution[i] << std::endl;

    }

#ifdef DEBUG
    // Print sorted rows reordering
    std::cout << "Sorted Rows Reordering: ";
    for(unsigned int i = 0; i < n; ++i) {
        std::cout << sortedRowsReord[i] << " ";
    }
    std::cout << std::endl; 

    // Print solution
    std::cout << "Reord Solution: ";
    for(unsigned int i = 0; i < n; ++i) {
        std::cout << solution[i] << " ";
    }
    std::cout << std::endl;
#endif



    return 0;
}
