#ifndef CRITTER__INTERCEPT__COMP_H_
#define CRITTER__INTERCEPT__COMP_H_

// C interface
// BLAS 1
void critter_daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy);
void critter_dscal_(const int n , const double a , double *x , const int incx);

// BLAS 2
void critter_dgbmv_(const int order, const int trans, const int m, const int n, const int kl, const int ku, const double alpha,
             const double *a, const int lda, const double *x, const int incx, const double beta, double *y, const int incy);
void critter_dgemv_(const int order, const int trans , const int m , const int n, const double alpha , const double *a ,
             const int lda , const double *x, const int incx , const double beta, double *y , const int incy );
void critter_dger_(const int order, const int m , const int n , const double alpha , const double *x , const int incx ,
            const double *y , const int incy , double *a , const int lda);
void critter_dsbmv_(const int Layout, const int uplo, const int n, const int k, const double alpha, const double *a,
             const int lda, const double *x, const int incx, const double beta, double *y, const int incy);
void critter_dspmv_(const int Layout, const int uplo, const int n, const double alpha, const double *ap, const double *x,
             const int incx, const double beta, double *y, const int incy);
void critter_dspr_(const int Layout, const int uplo, const int n, const double alpha, const double *x,
            const int incx, double *ap);
void critter_dspr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *ap);
void critter_dsymv_(const int Layout, const int uplo, const int n, const double alpha, const double *a, const int lda,
            const double *x, const int incx, const double beta, double *y, const int incy);
void critter_dsyr_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
            double *a, const int lda);
void critter_dsyr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *a, const int lda);
void critter_dtrsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *a,
             const int lda, double *x, const int incx);
void critter_dtrmv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a ,
             const int lda , double *x, const int incx );
void critter_dtpsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx);
void critter_dtpmv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx);
void critter_dtbsv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx);
void critter_dtbmv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx);

// BLAS 3
void critter_dgemm_(const int order, const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc);
void critter_dtrmm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void critter_dtrsm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb);
void critter_dsyrk_(const int order, const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc);
void critter_dsyr2k_(const int order, const int uplo, const int trans, const int n, const int k, const double alpha,
              const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
              const int ldc);
void critter_dsymm_(const int order, const int side, const int uplo, const int m, const int n, const double alpha,
             const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
             const int ldc);

// FORTRAN interface
// BLAS 1
void critter__daxpy__(const int* n , const double* a , const double *x , const int* incx , double *y , const int* incy);
void critter__dscal__(const int* n , const double* a , double *x , const int* incx);

// BLAS 2
void critter__dgbmv__(const char* trans , const int* m , const int* n, const int* kl, const int* ku, const double* alpha ,
               const double *a , const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy );
void critter__dgemv__(const char* trans , const int* m , const int* n, const double* alpha , const double *a ,
               const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy );
void critter__dger__(const int* m , const int* n , const double* alpha , const double *x , const int* incx ,
              const double *y , const int* incy , double *a , const int* lda);
void critter__dsbmv__(const char* uplo, const int* n, const int* k, const double* alpha, const double *a,
               const int* lda, const double *x, const int* incx, const double* beta, double *y, const int* incy);
void critter__dspmv__(const char* uplo, const int* n, const double* alpha, const double *ap, const double *x,
               const int* incx, const double* beta, double *y, const int* incy);
void critter__dspr__(const char* uplo, const int* n, const double* alpha, const double *x,
              const int* incx, double *ap);
void critter__dspr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *ap);
void critter__dsymv__(const char* uplo, const int* n, const double* alpha, const double *a, const int* lda,
               const double *x, const int* incx, const double* beta, double *y, const int* incy);
void critter__dsyr__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
              double *a, const int* lda);
void critter__dsyr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *a, const int* lda);
void critter__dtrsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *a,
               const int* lda, double *x, const int* incx);
void critter__dtrmv__(const char* uplo , const char* trans , const char* diag , const int* n , const double *a ,
               const int* lda , double *x, const int* incx );
void critter__dtpsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx);
void critter__dtpmv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx);
void critter__dtbsv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx);
void critter__dtbmv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx);

// BLAS 3
void critter__dgemm__(const char* transa , const char* transb ,
               const int* m , const int* n , const int* k , const double* alpha , const double *a ,
               const int* lda , const double *b , const int* ldb , const double* beta , double *c , const int* ldc);
void critter__dtrmm__(const char* side , const char* uplo , const char* transa ,
               const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
               const int* lda , double *b , const int* ldb);
void critter__dtrsm__(const char* side , const char* uplo , const char* transa ,
               const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
               const int* lda , double *b , const int* ldb);
void critter__dsyrk__(const char* uplo, const char* trans ,
               const int* n , const int* k , const double* alpha , const double *a , const int* lda ,
               const double* beta , double *c , const int* ldc);
void critter__dsyr2k__(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
                const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
                const int* ldc);
void critter__dsymm__(const char* side, const char* uplo, const int* m, const int* n, const double* alpha,
               const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
               const int* ldc);

// C interface
int critter_dgetrf_(int matrix_layout, int m , int n , double* a , int lda , int* ipiv);
int critter_dpotrf_(int matrix_layout, char uplo , int n , double* a , int lda);
int critter_dtrtri_(int matrix_layout, char uplo , char diag , int n , double* a , int lda);
int critter_dgeqrf_(int matrix_layout, int m , int n , double* a , int lda , double* tau);
int critter_dorgqr_(int matrix_layout, int m , int n , int k , double* a , int lda , const double* tau);
int critter_dormqr_(int matrix_layout, char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc);
int critter_dgetri_(int matrix_layout, int n , double * a , int lda , const int * ipiv);
int critter_dtpqrt_(int matrix_layout, int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt);
int critter_dtpmqrt_(int matrix_layout, char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb);

// FORTRAN interface
void critter__dgetrf__(const int* m , const int* n , double* a , const int* lda , int* ipiv, int* info);
void critter__dpotrf__(const char* uplo , const int* n , double* a , const int* lda, int* info);
void critter__dtrtri__(const char* uplo , const char* diag , const int* n , double* a , const int* lda, int* info);
void critter__dgeqrf__(const int* m , const int* n , double* a , const int* lda , double* tau, double* work, const int* lwork, int* info);
void critter__dorgqr__(const int* m , const int* n , const int* k , double* a , const int* lda , const double* tau, double* work, const int* lwork, int* info);
void critter__dormqr__(const char* side , const char* trans , const int* m , const int* n , const int* k , const double * a , const int* lda , const double * tau ,
                double * c , const int* ldc, double* work, const int* lwork, int* info);
void critter__dgetri__(const int* n , double * a , const int* lda , const int * ipiv, double* work, const int* lwork, int* info);
void critter__dtpqrt__(const int* m , const int* n , const int* l , const int* nb , double * a , const int* lda , double* b , const int* ldb , double * t , const int* ldt,
                double* work, int* info);
void critter__dtpmqrt__(const char* side , const char* trans , const int* m , const int* n , const int* k , const int* l , const int* nb , const double * v ,
                 const int* ldv , const double * t , const int* ldt , double * a , const int* lda , double * b , const int* ldb, double* work, int* info);

#endif /*CRITTER__INTERCEPT__COMP_H_*/
