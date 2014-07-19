#include <R.h>
#include <Rmath.h>
#include <time.h>

#ifndef max
#  define max(a, b) ((a) > (b) ? (a) : (b))
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

/* trace of a matrix */
double mtrace(const double *A, const int n);

/* Bunch-Kaufman factorization for a symmetric matrix */
extern
void F77_SUB(dsytrf)(const char *uplo, const int *n, double *a, const int *lda,
                     int *ipiv, double *work, const int *lwork, int *info);

/* inverse based on Bunch-Kaufman factorization */
extern
void F77_SUB(dsytri)(const char *uplo, const int *n, double *a, const int *lda,
                     int *ipiv, double *work, int *info);

/* Cholesky factorization for a symmetric positive definite matrix */
extern
void F77_SUB(dpotrf)(const char *uplo, const int *n, double *a, const int *lda,
                     int *info);

/* inverse based on Cholesky factorization */
extern
void F77_SUB(dpotri)(const char *uplo, const int *n, double *a, const int *lda,
                     int *info);
