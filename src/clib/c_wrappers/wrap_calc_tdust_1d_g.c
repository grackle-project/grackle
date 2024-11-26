#include "grackle_macros.h"

// TODO: FORTRAN_NAME not needed for unit tests
extern void FORTRAN_NAME(calc_tdust_1d_g)(double* tdust, double* tgas, double* nh, double* gasgr,
                                          double* gamma_isrfa, double* isrf, int* itmask, double* trad,
                                          int* in, int* is, int* ie, int* j, int* k);

extern void FORTRAN_NAME(set_iteration_mask_and_initial_guess)(double *gamma_isrf, double *isrf, double *gamma_isrfa,
                                                               int *itmask, int *nm_itmask, int *bi_itmask,
                                                               double *tdustnow, double *tgas, double *trad,
                                                               double *t_subl, double *radf, double *kgr1,
                                                               double *pert_i, double *pert, int *c_done, int *nm_done,
                                                               int *is, int *ie);

extern void FORTRAN_NAME(calc_tdplus)(int* in, int *is, int *ie, int *nm_itmask,
                                      double *pert, double *tdustnow, double *tdplus);

extern void FORTRAN_NAME(run_newton_method_teration)(int* in, int* is, int* ie, int* nm_itmask, int* bi_itmask,
                                                     double* solplus, double* sol, double* slope, double* tdustnow,
                                                     double* tdustold, double* pert, double* minpert, double* trad,
                                                     double* tol, int* nm_done, int* c_done);

extern void FORTRAN_NAME(run_bisection_method_iteration)(int* in, int* is, int* ie, int* bi_itmask, double* tdustnow,
                                                         double* bi_t_high, double* bi_t_mid, int* iter, double* t_subl,
                                                         double* kgr, double* tgas, double* trad4, double* gasgr,
                                                         double* gamma_isrf, double* nh, double* sol, double* bi_tol,
                                                         int* c_done, int* c_total);

extern void FORTRAN_NAME(calc_gr_balance_g)(double *tdust, double *tgas, double *kgr, double *trad4,
                                            double *gasgr, double *gamma_isrf, double *nh,
                                            int *itmask, double *sol, int *in, int *is, int *ie);

extern void FORTRAN_NAME(calc_kappa_gr_g)(double *tdust, double *kgr, int *itmask,
                                          int *in, int *is, int *ie, double *t_subl);
