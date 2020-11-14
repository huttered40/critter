#ifndef CRITTER__INTERCEPT__COMP_H_
#define CRITTER__INTERCEPT__COMP_H_

namespace critter{
namespace internal{

// Note: I am only adding what is currently necessary. Interfaces for other BLAS routines can be added later as necessary
// Note: Fortran interface has not been added. It can be added later as necessary.
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy);
void _dscal_(const int n , const double a , double *x , const int incx);

void _dger_(const int m , const int n , const double alpha , const double *x , const int incx , const double *y , const int incy , double *a ,
            const int lda);
void _dgemv_(const int trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy );
void _dtrmv_(const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx );
void _dgemm_(const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc);
void _dtrmm_(const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _dtrsm_(const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _dsyrk_(const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc);

void __dgemv__(const char trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy );
void __dtrmv__(const char uplo , const char trans , const char diag , const int n , const double *a , const int lda , double *x, const int incx );
void __dgemm__(const char transa , const char transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc);
void __dtrmm__(const char side , const char uplo , const char transa ,
             const char diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void __dtrsm__(const char side , const char uplo , const char transa ,
             const char diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void __dsyrk__(const char uplo, const char trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc);

// Note: I am only adding what is currently necessary. Interfaces for other LAPACKE routines can be added later as necessary
// Note: Fortran interface has not been added. It can be added later as necessary.
void _dgetrf_(int m , int n , double* a , int lda , int* ipiv);
void _dpotrf_(char uplo , int n , double* a , int lda);
void _dtrtri_(char uplo , char diag , int n , double* a , int lda);
void _dgeqrf_(int m , int n , double* a , int lda , double* tau);
void _dorgqr_(int m , int n , int k , double* a , int lda , const double* tau);
void _dormqr_(char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc);
void _dgetri_(int n , double * a , int lda , const int * ipiv);
void _dtpqrt_(int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt);
void _dtpmqrt_(char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb);

void _blk_to_cyc_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim);
void _cyc_to_blk_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim);
}
}

#endif /*CRITTER__INTERCEPT__COMP_H_*/
