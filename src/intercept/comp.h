#ifndef CRITTER__INTERCEPT__COMP_H_
#define CRITTER__INTERCEPT__COMP_H_

#ifdef MKL
// Note: this MKL inclusion should be conditional on config.mk
#include "mkl.h"
#else
#define CBLAS_LAYOUT int
#define CBLAS_SIDE int
#define CBLAS_DIAG int
#define CBLAS_TRANSPOSE int
#define CBLAS_UPLO int
#endif /* MKL */

namespace critter{
namespace internal{

// Note: I am only adding what is currently necessary. Interfaces for other BLAS routines can be added later as necessary
// Note: Fortran interface has not been added. It can be added later as necessary.
void _saxpy_(const int n , const float a , const float *x , const int incx , float *y , const int incy);
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy);
void _sscal_(const int n , const float a , float *x , const int incx);
void _dscal_(const int n , const double a , double *x , const int incx);
void _sgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const float alpha , const float *a ,
             const int lda , const float *b , const int ldb , const float beta , float *c , const int ldc);
void _dgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc);
void _strmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb);
void _dtrmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _strsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb);
void _dtrsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _ssyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const float alpha , const float *a , const int lda ,
             const float beta , float *c , const int ldc);
void _dsyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc);
void _sger_(const CBLAS_LAYOUT Layout , const int m , const int n , const float alpha , const float *x , const int incx ,
            const float *y , const int incy , float *a , const int lda);
void _dger_(const CBLAS_LAYOUT Layout , const int m , const int n , const double alpha , const double *x , const int incx , const double *y , const int incy , double *a ,
            const int lda);

// Note: I am only adding what is currently necessary. Interfaces for other LAPACKE routines can be added later as necessary
// Note: Fortran interface has not been added. It can be added later as necessary.
void _sgetrf_(int matrix_layout , int m , int n , float* a , int lda , int* ipiv);
void _dgetrf_(int matrix_layout , int m , int n , double* a , int lda , int* ipiv);
void _spotrf_(int matrix_layout , char uplo , int n , float* a , int lda);
void _dpotrf_(int matrix_layout , char uplo , int n , double* a , int lda);
void _strtri_(int matrix_layout , char uplo , char diag , int n , float* a , int lda);
void _dtrtri_(int matrix_layout , char uplo , char diag , int n , double* a , int lda);
void _sgeqrf_(int matrix_layout , int m , int n , float* a , int lda , float* tau);
void _dgeqrf_(int matrix_layout , int m , int n , double* a , int lda , double* tau);
void _sorgqr_(int matrix_layout , int m , int n , int k , float* a , int lda , const float* tau);
void _dorgqr_(int matrix_layout , int m , int n , int k , double* a , int lda , const double* tau);
void _sormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const float * a , int lda , const float * tau , float * c , int ldc);
void _dormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc);
void _sgetri_(int matrix_layout , int n , float * a , int lda , const int * ipiv);
void _dgetri_(int matrix_layout , int n , double * a , int lda , const int * ipiv);
void _stpqrt_(int matrix_layout , int m , int n , int l , int nb , float * a , int lda , float * b , int ldb , float * t , int ldt);
void _dtpqrt_(int matrix_layout , int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt);
void _stpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const float * v ,
               int ldv , const float * t , int ldt , float * a , int lda , float * b , int ldb);
void _dtpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb);

}
}

#endif /*CRITTER__INTERCEPT__COMP_H_*/
