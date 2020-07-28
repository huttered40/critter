#include "comp.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{
namespace internal{

void _saxpy_(const int n , const float a , const float *x , const int incx , float *y , const int incy){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 2.*_n;
    initiate_comp(_BLAS_axpy__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_saxpy(n,a,x,incx,y,incy);
#endif
#endif
    complete_comp(_BLAS_axpy__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_saxpy(n,a,x,incx,y,incy);
#endif
#endif
  }
}
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 2.*_n;
    initiate_comp(_BLAS_axpy__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_daxpy(n,a,x,incx,y,incy);
#endif
#endif
    complete_comp(_BLAS_axpy__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_daxpy(n,a,x,incx,y,incy);
#endif
#endif
  }
}
void _sscal_(const int n , const float a , float *x , const int incx){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1.*_n;
    initiate_comp(_BLAS_scal__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_sscal(n,a,x,incx);
#endif
#endif
    complete_comp(_BLAS_scal__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_sscal(n,a,x,incx);
#endif
#endif
  }
}
void _dscal_(const int n , const double a , double *x , const int incx){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1.*_n;
    initiate_comp(_BLAS_scal__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dscal(n,a,x,incx);
#endif
#endif
    complete_comp(_BLAS_scal__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dscal(n,a,x,incx);
#endif
#endif
  }
}
void _sger_(const CBLAS_LAYOUT Layout , const int m , const int n , const float alpha , const float *x , const int incx ,
            const float *y , const int incy , float *a , const int lda){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = 2.*_m*_n;
    initiate_comp(_BLAS_ger__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_sger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
    complete_comp(_BLAS_ger__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_sger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
  }
}
void _dger_(const CBLAS_LAYOUT Layout , const int m , const int n , const double alpha , const double *x , const int incx , const double *y , const int incy , double *a ,
            const int lda){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = 2.*_m*_n;
    initiate_comp(_BLAS_ger__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
    complete_comp(_BLAS_ger__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
  }
}
void _sgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const float alpha , const float *a ,
             const int lda , const float *b , const int ldb , const float beta , float *c , const int ldc){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m; double _k = k;
    double flops = 2.*_m*_n*_k;
    initiate_comp(_BLAS_gemm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_sgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
    complete_comp(_BLAS_gemm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_sgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
  }
}
void _dgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m; double _k = k;
    double flops = 2.*_m*_n*_k;
    initiate_comp(_BLAS_gemm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
    complete_comp(_BLAS_gemm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
  }
}
void _strmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    initiate_comp(_BLAS_trmm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_strmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(_BLAS_trmm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_strmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
  }
}
void _dtrmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    initiate_comp(_BLAS_trmm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dtrmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(_BLAS_trmm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dtrmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
  }
}
void _strsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    initiate_comp(_BLAS_trsm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_strsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(_BLAS_trsm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_strsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
  }
}
void _dtrsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    initiate_comp(_BLAS_trsm__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dtrsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(_BLAS_trsm__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dtrsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
  }
}
void _ssyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const float alpha , const float *a , const int lda ,
             const float beta , float *c , const int ldc){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _k = k;
    double flops = 2*_k*_n*_n;//Note: the factor of 2 might be wrong
    initiate_comp(_BLAS_syrk__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_ssyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
    complete_comp(_BLAS_syrk__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_ssyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
  }
}
void _dsyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _k = k;
    double flops = 2*_k*_n*_n;//Note: the factor of 2 might be wrong
    initiate_comp(_BLAS_syrk__id,curtime);
#ifdef MKL
#ifdef CBLAS
    cblas_dsyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
    complete_comp(_BLAS_syrk__id,flops);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dsyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
  }
}

void _sgetrf_(int matrix_layout , int m , int n , float* a , int lda , int* ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 2./3. * _n*_n*_n; if (m > n) flops = (1./3.)*_n*_n*(3*_m-_n); if (m < n) flops = (1./3.)*_m*_m*(3*_n-_m);
    initiate_comp(_LAPACK_getrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv);
#endif
#endif
    complete_comp(_LAPACK_getrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv);
#endif
#endif
  }
}
void _dgetrf_(int matrix_layout , int m , int n , double* a , int lda , int* ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 2./3. * _n*_n*_n; if (m > n) flops = (1./3.)*_n*_n*(3*_m-_n); if (m < n) flops = (1./3.)*_m*_m*(3*_n-_m);
    initiate_comp(_LAPACK_getrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv);
#endif
#endif
    complete_comp(_LAPACK_getrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv);
#endif
#endif
  }
}
void _spotrf_(int matrix_layout , char uplo , int n , float* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n;
    initiate_comp(_LAPACK_potrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_spotrf(matrix_layout,uplo,n,a,lda);
#endif
#endif
    complete_comp(_LAPACK_potrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_spotrf(matrix_layout,uplo,n,a,lda);
#endif
#endif
  }
}
void _dpotrf_(int matrix_layout , char uplo , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n;
    initiate_comp(_LAPACK_potrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda);
#endif
#endif
    complete_comp(_LAPACK_potrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda);
#endif
#endif
  }
}
void _strtri_(int matrix_layout , char uplo , char diag , int n , float* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n;
    initiate_comp(_LAPACK_trtri__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda);
#endif
#endif
    complete_comp(_LAPACK_trtri__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda);
#endif
#endif
  }
}
void _dtrtri_(int matrix_layout , char uplo , char diag , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n;
    initiate_comp(_LAPACK_trtri__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda);
#endif
#endif
    complete_comp(_LAPACK_trtri__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda);
#endif
#endif
  }
}
void _sgeqrf_(int matrix_layout , int m , int n , float* a , int lda , float* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 4./3. * _n*_n*_n; if (m > n) flops = (2./3.)*_n*_n*(3*_m-_n); if (m < n) flops = (2./3.)*_m*_m*(3*_n-_m);
    initiate_comp(_LAPACK_geqrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
#endif
#endif
    complete_comp(_LAPACK_geqrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
#endif
#endif
  }
}
void _dgeqrf_(int matrix_layout , int m , int n , double* a , int lda , double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 4./3. * _n*_n*_n; if (m > n) flops = (2./3.)*_n*_n*(3*_m-_n); if (m < n) flops = (2./3.)*_m*_m*(3*_n-_m);
    initiate_comp(_LAPACK_geqrf__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau);
#endif
#endif
    complete_comp(_LAPACK_geqrf__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau);
#endif
#endif
  }
}
void _sorgqr_(int matrix_layout , int m , int n , int k , float* a , int lda , const float* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 4.*_m*_n*_k - 2.*(_m+_n) * _k*_k + (4./3.)*_k*_k*_k;//Note: this routine, which forms an explicit Q, may not even count as flops, rather as overhead
    initiate_comp(_LAPACK_orgqr__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau);
#endif
#endif
    complete_comp(_LAPACK_orgqr__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau);
#endif
#endif
  }
}
void _dorgqr_(int matrix_layout , int m , int n , int k , double* a , int lda , const double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 4.*_m*_n*_k - 2.*(_m+_n) * _k*_k + (4./3.)*_k*_k*_k;//Note: this routine, which forms an explicit Q, may not even count as flops, rather as overhead
    initiate_comp(_LAPACK_orgqr__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau);
#endif
#endif
    complete_comp(_LAPACK_orgqr__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau);
#endif
#endif
  }
}
void _sormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const float * a , int lda , const float * tau , float * c , int ldc){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_ormqr__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
#endif
#endif
    complete_comp(_LAPACK_ormqr__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
#endif
#endif
  }
}
void _dormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_ormqr__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
#endif
#endif
    complete_comp(_LAPACK_ormqr__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
#endif
#endif
  }
}
void _sgetri_(int matrix_layout , int n , float * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 4./3.*_n*_n*_n;
    initiate_comp(_LAPACK_getri__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv);
#endif
#endif
    complete_comp(_LAPACK_getri__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv);
#endif
#endif
  }
}
void _dgetri_(int matrix_layout , int n , double * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 4./3.*_n*_n*_n;
    initiate_comp(_LAPACK_getri__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv);
#endif
#endif
    complete_comp(_LAPACK_getri__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv);
#endif
#endif
  }
}
void _stpqrt_(int matrix_layout , int m , int n , int l , int nb , float * a , int lda , float * b , int ldb , float * t , int ldt){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _l = l;
    double flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_tpqrt__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
#endif
#endif
    complete_comp(_LAPACK_tpqrt__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
#endif
#endif
  }
}
void _dtpqrt_(int matrix_layout , int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _l = l;
    double flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_tpqrt__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
#endif
#endif
    complete_comp(_LAPACK_tpqrt__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
#endif
#endif
  }
}
void _stpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const float * v ,
               int ldv , const float * t , int ldt , float * a , int lda , float * b , int ldb){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_tpmqrt__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
#endif
#endif
    complete_comp(_LAPACK_tpmqrt__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
#endif
#endif
  }
}
void _dtpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    initiate_comp(_LAPACK_tpmqrt__id,curtime);
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
#endif
#endif
    complete_comp(_LAPACK_tpmqrt__id,flops);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
#endif
#endif
  }
}

}
}
