#ifndef CRITTER__INTERCEPT__COMP_H_
#define CRITTER__INTERCEPT__COMP_H_

namespace critter{
namespace internal{

// C interface
// BLAS 1
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy);
void _dscal_(const int n , const double a , double *x , const int incx);

// BLAS 2
void _dgbmv_(const int order, const int trans, const int m, const int n, const int kl, const int ku, const double alpha,
             const double *a, const int lda, const double *x, const int incx, const double beta, double *y, const int incy);
void _dgemv_(const int order, const int trans , const int m , const int n, const double alpha , const double *a ,
             const int lda , const double *x, const int incx , const double beta, double *y , const int incy );
void _dger_(const int order, const int m , const int n , const double alpha , const double *x , const int incx ,
            const double *y , const int incy , double *a , const int lda);
void _dsbmv_(const int Layout, const int uplo, const int n, const int k, const double alpha, const double *a,
             const int lda, const double *x, const int incx, const double beta, double *y, const int incy);
void _dspmv_(const int Layout, const int uplo, const int n, const double alpha, const double *ap, const double *x,
             const int incx, const double beta, double *y, const int incy);
void _dspr_(const int Layout, const int uplo, const int n, const double alpha, const double *x,
            const int incx, double *ap);
void _dspr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *ap);
void _dsymv_(const int Layout, const int uplo, const int n, const double alpha, const double *a, const int lda,
            const double *x, const int incx, const double beta, double *y, const int incy);
void _dsyr_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
            double *a, const int lda);
void _dsyr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *a, const int lda);
void _dtrsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *a,
             const int lda, double *x, const int incx);
void _dtrmv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a ,
             const int lda , double *x, const int incx );
void _dtpsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx);
void _dtpmv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx);
void _dtbsv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx);
void _dtbmv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx);

// BLAS 3
void _dgemm_(const int order, const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc);
void _dtrmm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _dtrsm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void _dsyrk_(const int order, const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc);
void _dsyr2k_(const int order, const int uplo, const int trans, const int n, const int k, const double alpha,
              const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
              const int ldc);
void _dsymm_(const int order, const int side, const int uplo, const int m, const int n, const double alpha,
             const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
             const int ldc);

// FORTRAN interface
// BLAS 1
void __daxpy__(const int* n , const double* a , const double *x , const int* incx , double *y , const int* incy);
void __dscal__(const int* n , const double* a , double *x , const int* incx);

// BLAS 2
void __dgbmv__(const char* trans , const int* m , const int* n, const int* kl, const int* ku, const double* alpha ,
               const double *a , const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy );
void __dgemv__(const char* trans , const int* m , const int* n, const double* alpha , const double *a ,
               const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy );
void __dger__(const int* m , const int* n , const double* alpha , const double *x , const int* incx ,
              const double *y , const int* incy , double *a , const int* lda);
void __dsbmv__(const char* uplo, const int* n, const int* k, const double* alpha, const double *a,
               const int* lda, const double *x, const int* incx, const double* beta, double *y, const int* incy);
void __dspmv__(const char* uplo, const int* n, const double* alpha, const double *ap, const double *x,
               const int* incx, const double* beta, double *y, const int* incy);
void __dspr__(const char* uplo, const int* n, const double* alpha, const double *x,
              const int* incx, double *ap);
void __dspr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *ap);
void __dsymv__(const char* uplo, const int* n, const double* alpha, const double *a, const int* lda,
               const double *x, const int* incx, const double* beta, double *y, const int* incy);
void __dsyr__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
              double *a, const int* lda);
void __dsyr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *a, const int* lda);
void __dtrsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *a,
               const int* lda, double *x, const int* incx);
void __dtrmv__(const char* uplo , const char* trans , const char* diag , const int* n , const double *a ,
               const int* lda , double *x, const int* incx );
void __dtpsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx);
void __dtpmv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx);
void __dtbsv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx);
void __dtbmv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx);

// BLAS 3
void __dgemm__(const char* transa , const char* transb ,
               const int* m , const int* n , const int* k , const double* alpha , const double *a ,
               const int* lda , const double *b , const int* ldb , const double* beta , double *c , const int* ldc);
void __dtrmm__(const char* side , const char* uplo , const char* transa ,
               const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
               const int* lda , double *b , const int* ldb);
void __dtrsm__(const char* side , const char* uplo , const char* transa ,
               const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
               const int* lda , double *b , const int* ldb);
void __dsyrk__(const char* uplo, const char* trans ,
               const int* n , const int* k , const double* alpha , const double *a , const int* lda ,
               const double* beta , double *c , const int* ldc);
void __dsyr2k__(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
                const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
                const int* ldc);
void __dsymm__(const char* side, const char* uplo, const int* m, const int* n, const double* alpha,
               const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
               const int* ldc);

// C interface
int _dgetrf_(int matrix_layout, int m , int n , double* a , int lda , int* ipiv);
int _dpotrf_(int matrix_layout, char uplo , int n , double* a , int lda);
int _dtrtri_(int matrix_layout, char uplo , char diag , int n , double* a , int lda);
int _dgeqrf_(int matrix_layout, int m , int n , double* a , int lda , double* tau);
int _dorgqr_(int matrix_layout, int m , int n , int k , double* a , int lda , const double* tau);
int _dormqr_(int matrix_layout, char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc);
int _dgetri_(int matrix_layout, int n , double * a , int lda , const int * ipiv);
int _dtpqrt_(int matrix_layout, int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt);
int _dtpmqrt_(int matrix_layout, char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb);

// FORTRAN interface
void __dgetrf__(const int* m , const int* n , double* a , const int* lda , int* ipiv, int* info);
void __dpotrf__(const char* uplo , const int* n , double* a , const int* lda, int* info);
void __dtrtri__(const char* uplo , const char* diag , const int* n , double* a , const int* lda, int* info);
void __dgeqrf__(const int* m , const int* n , double* a , const int* lda , double* tau, double* work, const int* lwork, int* info);
void __dorgqr__(const int* m , const int* n , const int* k , double* a , const int* lda , const double* tau, double* work, const int* lwork, int* info);
void __dormqr__(const char* side , const char* trans , const int* m , const int* n , const int* k , const double * a , const int* lda , const double * tau ,
                double * c , const int* ldc, double* work, const int* lwork, int* info);
void __dgetri__(const int* n , double * a , const int* lda , const int * ipiv, double* work, const int* lwork, int* info);
void __dtpqrt__(const int* m , const int* n , const int* l , const int* nb , double * a , const int* lda , double* b , const int* ldb , double * t , const int* ldt,
                double* work, int* info);
void __dtpmqrt__(const char* side , const char* trans , const int* m , const int* n , const int* k , const int* l , const int* nb , const double * v ,
                 const int* ldv , const double * t , const int* ldt , double * a , const int* lda , double * b , const int* ldb, double* work, int* info);
}
}

#endif /*CRITTER__INTERCEPT__COMP_H_*/
