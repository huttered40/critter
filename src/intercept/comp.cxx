#include "comp.h"
#include "../util/util.h"

namespace critter{
namespace internal{

//#ifdef CBLAS
void _saxpy_(const int n , const float a , const float *x , const int incx , float *y , const int incy){
  if (mode && track_blas){
    cblas_saxpy(n,a,x,incx,y,incy);
  } else{
    cblas_saxpy(n,a,x,incx,y,incy);
  }
}
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy){
  if (mode && track_blas){
    cblas_daxpy(n,a,x,incx,y,incy);
  } else{
    cblas_daxpy(n,a,x,incx,y,incy);
  }
}
void _sscal_(const int n , const float a , float *x , const int incx){
  if (mode && track_blas){
    cblas_sscal(n,a,x,incx);
  } else{
    cblas_sscal(n,a,x,incx);
  }
}
void _dscal_(const int n , const double a , double *x , const int incx){
  if (mode && track_blas){
    cblas_dscal(n,a,x,incx);
  } else{
    cblas_dscal(n,a,x,incx);
  }
}
void _sger_(const CBLAS_LAYOUT Layout , const int m , const int n , const float alpha , const float *x , const int incx ,
            const float *y , const int incy , float *a , const int lda){
  if (mode && track_blas){
    cblas_sger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
  } else{
    cblas_sger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
  }
}
void _dger_(const CBLAS_LAYOUT Layout , const int m , const int n , const double alpha , const double *x , const int incx , const double *y , const int incy , double *a ,
            const int lda){
  if (mode && track_blas){
    cblas_dger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
  } else{
    cblas_dger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
  }
}
void _sgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const float alpha , const float *a ,
             const int lda , const float *b , const int ldb , const float beta , float *c , const int ldc){
  if (mode && track_blas){
    cblas_sgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  } else{
    cblas_sgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }
}
void _dgemm_(const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  if (mode && track_blas){
    cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  } else{
    cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }
}
void _strmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb){
  if (mode && track_blas){
    cblas_strmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  } else{
    cblas_strmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _dtrmm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas){
    cblas_dtrmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  } else{
    cblas_dtrmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _strsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const float alpha , const float *a ,
             const int lda , float *b , const int ldb){
  if (mode && track_blas){
    cblas_strsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  } else{
    cblas_strsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _dtrsm_(const CBLAS_LAYOUT Layout , const CBLAS_SIDE side , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE transa ,
             const CBLAS_DIAG diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas){
    cblas_dtrsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  } else{
    cblas_dtrsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _ssyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const float alpha , const float *a , const int lda ,
             const float beta , float *c , const int ldc){
  if (mode && track_blas){
    cblas_ssyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
  } else{
    cblas_ssyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
  }
}
void _dsyrk_(const CBLAS_LAYOUT Layout , const CBLAS_UPLO uplo , const CBLAS_TRANSPOSE trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  if (mode && track_blas){
    cblas_dsyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
  } else{
    cblas_dsyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
  }
}
//#endif /* CBLAS */

//#ifdef LAPACKE
void _sgetrf_(int matrix_layout , int m , int n , float* a , int lda , int* ipiv){
  if (mode && track_lapack){
    LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv);
  } else{
    LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv);
  }
}
void _dgetrf_(int matrix_layout , int m , int n , double* a , int lda , int* ipiv){
  if (mode && track_lapack){
    LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv);
  } else{
    LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv);
  }
}
void _spotrf_(int matrix_layout , char uplo , int n , float* a , int lda){
  if (mode && track_lapack){
    LAPACKE_spotrf(matrix_layout,uplo,n,a,lda);
  } else{
    LAPACKE_spotrf(matrix_layout,uplo,n,a,lda);
  }
}
void _dpotrf_(int matrix_layout , char uplo , int n , double* a , int lda){
  if (mode && track_lapack){
    LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda);
  } else{
    LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda);
  }
}
void _strtri_(int matrix_layout , char uplo , char diag , int n , float* a , int lda){
  if (mode && track_lapack){
    LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda);
  } else{
    LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda);
  }
}
void _dtrtri_(int matrix_layout , char uplo , char diag , int n , double* a , int lda){
  if (mode && track_lapack){
    LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda);
  } else{
    LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda);
  }
}
void _sgeqrf_(int matrix_layout , int m , int n , float* a , int lda , float* tau){
  if (mode && track_lapack){
    LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
  } else{
    LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
  }
}
void _dgeqrf_(int matrix_layout , int m , int n , double* a , int lda , double* tau){
  if (mode && track_lapack){
    LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau);
  } else{
    LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau);
  }
}
void _sorgqr_(int matrix_layout , int m , int n , int k , float* a , int lda , const float* tau){
  if (mode && track_lapack){
    LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau);
  } else{
    LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau);
  }
}
void _dorgqr_(int matrix_layout , int m , int n , int k , double* a , int lda , const double* tau){
  if (mode && track_lapack){
    LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau);
  } else{
    LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau);
  }
}
void _sormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const float * a , int lda , const float * tau , float * c , int ldc){
  if (mode && track_lapack){
    LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
  } else{
    LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
  }
}
void _dormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  if (mode && track_lapack){
    LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
  } else{
    LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
  }
}
void _sgetri_(int matrix_layout , int n , float * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv);
  } else{
    LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv);
  }
}
void _dgetri_(int matrix_layout , int n , double * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv);
  } else{
    LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv);
  }
}
void _stpqrt_(int matrix_layout , int m , int n , int l , int nb , float * a , int lda , float * b , int ldb , float * t , int ldt){
  if (mode && track_lapack){
    LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
  } else{
    LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
  }
}
void _dtpqrt_(int matrix_layout , int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  if (mode && track_lapack){
    LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
  } else{
    LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
  }
}
void _stpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const float * v ,
               int ldv , const float * t , int ldt , float * a , int lda , float * b , int ldb){
  if (mode && track_lapack){
    LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
  } else{
    LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
  }
}
void _dtpmqrt_(int matrix_layout , char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb){
  if (mode && track_lapack){
    LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
  } else{
    LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
  }
}
//#endif /* LAPACKE */

}
}
