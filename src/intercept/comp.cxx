#include "comp.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{
namespace internal{

/*
void _saxpy_(const int n , const float a , const float *x , const int incx , float *y , const int incy){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 2.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x),reinterpret_cast<intptr_t>(y)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_axpy__id,curtime,flops,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_saxpy(n,a,x,incx,y,incy);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_axpy__id,flops,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x),reinterpret_cast<intptr_t>(y)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_axpy__id,curtime,flops,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_daxpy(n,a,x,incx,y,incy);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_axpy__id,flops,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_scal__id,curtime,flops,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_sscal(n,a,x,incx);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_scal__id,flops,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_scal__id,curtime,flops,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dscal(n,a,x,incx);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_scal__id,flops,n);
  } else{
#ifdef MKL
#ifdef CBLAS
    cblas_dscal(n,a,x,incx);
#endif
#endif
  }
}
*/

void _sger_(const CBLAS_LAYOUT Layout , const int m , const int n , const float alpha , const float *x , const int incx ,
            const float *y , const int incy , float *a , const int lda){
  if (mode && track_blas){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = 2.*_m*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x),reinterpret_cast<intptr_t>(y),reinterpret_cast<intptr_t>(a)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_ger__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_sger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_ger__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(x),reinterpret_cast<intptr_t>(y),reinterpret_cast<intptr_t>(a)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_ger__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dger(Layout,m,n,alpha,x,incx,y,incy,a,lda);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_ger__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(c)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_gemm__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_sgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_gemm__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(c)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_BLAS_gemm__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision){
      cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    }
#endif
#endif
    complete_comp(special_time,ptrs,_BLAS_gemm__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_trmm__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_strmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_trmm__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_trmm__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dtrmm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_trmm__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_trsm__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_strsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_trsm__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_trsm__id,curtime,flops,m,n);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dtrsm(Layout,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_trsm__id,flops,m,n);
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
    double flops = _k*_n*(_n+1);
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(c)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_syrk__id,curtime,flops,n,k);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_ssyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_syrk__id,flops,n,k);
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
    double flops = _k*_n*(_n+1);
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(c)};
    bool schedule_decision = initiate_comp(ptrs,_BLAS_syrk__id,curtime,flops,n,k);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dsyrk(Layout,uplo,trans,n,k,alpha,a,lda,beta,c,ldc);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_syrk__id,flops,n,k);
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
    double flops = 0;
    if (m>=n){
      flops = _m*_n*_n - 1./3.*_n*_n*_n - 1./2.*_n*_n + 5./6.*_n;
    } else{
      flops = _n*_m*_m - 1./3.*_m*_m*_m - 1./2.*_m*_m + 5./6.*_m;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_getrf__id,flops,m,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
#endif
#endif
  }
}
void _dgetrf_(int matrix_layout , int m , int n , double* a , int lda , int* ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 0;
    if (m>=n){
      flops = _m*_n*_n - 1./3.*_n*_n*_n - 1./2.*_n*_n + 5./6.*_n;
    } else{
      flops = _n*_m*_m - 1./3.*_m*_m*_m - 1./2.*_m*_m + 5./6.*_m;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_getrf__id,flops,m,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv)==0);
#endif
#endif
  }
}
void _spotrf_(int matrix_layout , char uplo , int n , float* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 1./2.*_n*_n + 1./6.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_potrf__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_spotrf(matrix_layout,uplo,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_spotrf(matrix_layout,uplo,n,a,lda)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_potrf__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_spotrf(matrix_layout,uplo,n,a,lda)==0);
#endif
#endif
  }
}
void _dpotrf_(int matrix_layout , char uplo , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 1./2.*_n*_n + 1./6.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_potrf__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_potrf__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda)==0);
#endif
#endif
  }
}
void _strtri_(int matrix_layout , char uplo , char diag , int n , float* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 2./3.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_trtri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_trtri__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda)==0);
#endif
#endif
  }
}
void _dtrtri_(int matrix_layout , char uplo , char diag , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 2./3.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_trtri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_trtri__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda)==0);
#endif
#endif
  }
}
void _sgeqrf_(int matrix_layout , int m , int n , float* a , int lda , float* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 0;
    if (m>=n){
      flops = 2.*_m*_n*_n - 2./3.*_n*_n*_n + _m*_n + _n*_n + 14./3.*_n;
    } else{
      flops = 2.*_n*_m*_m - 2./3.*_m*_m*_m + 3.*_m*_n - _m*_m + 14./3.*_m;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_geqrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_geqrf__id,flops,m,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau)==0);
#endif
#endif
  }
}
void _dgeqrf_(int matrix_layout , int m , int n , double* a , int lda , double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 0;
    if (m>=n){
      flops = 2.*_m*_n*_n - 2./3.*_n*_n*_n + _m*_n + _n*_n + 14./3.*_n;
    } else{
      flops = 2.*_n*_m*_m - 2./3.*_m*_m*_m + 3.*_m*_n - _m*_m + 14./3.*_m;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_geqrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_geqrf__id,flops,m,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau)==0);
#endif
#endif
  }
}
void _sorgqr_(int matrix_layout , int m , int n , int k , float* a , int lda , const float* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 4.*_m*_n*_k - 2.*(_m+_n)*_k*_k + (4./3.)*_k*_k*_k + 3.*_n*_k - _m*_k - _k*_k - 4./3.*_k;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_orgqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau)==0);
      else{
        special_time = MPI_Wtime();
        float* tau_temp = (float*)tau;
        for (int i=0; i<n; i++){
          memset(a+i*lda,0,m*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
        }
        memset(tau_temp,1,k*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau_temp)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_orgqr__id,flops,m,n,k);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau)==0);
#endif
#endif
  }
}
void _dorgqr_(int matrix_layout , int m , int n , int k , double* a , int lda , const double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 4.*_m*_n*_k - 2.*(_m+_n)*_k*_k + (4./3.)*_k*_k*_k + 3.*_n*_k - _m*_k - _k*_k - 4./3.*_k;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_orgqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau)==0);
      else{
        special_time = MPI_Wtime();
        double* tau_temp = (double*)tau;
        for (int i=0; i<n; i++){
          memset(a+i*lda,0,m*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
        }
        memset(tau_temp,1,k*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau_temp)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_orgqr__id,flops,m,n,k);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau)==0);
#endif
#endif
  }
}
void _sormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const float * a , int lda , const float * tau , float * c , int ldc){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 0;
    if (side == 'L'){
      flops = 4.*_m*_n*_k - 2.*_n*_k*_k + 3.*_n*_k;
    } else{
      flops = 4.*_n*_m*_k - 2.*_m*_k*_k + 2.*m*_k + _n*_k - 1./2.*_k*_k + 1./2.*_k;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_ormqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(c+i*ldc,1,m*sizeof(float));// Assumes column-major
        }
        float* a_temp = (float*)a;
        float* tau_temp = (float*)tau;
        if (side == 'L'){
          for (int i=0; i<k; i++){
            memset(a_temp+i*lda,0,m*sizeof(float));// Assumes column-major
            a_temp[i*lda+i] = 1;
          }
        } else{

          for (int i=0; i<k; i++){
            memset(a_temp+i*lda,0,n*sizeof(float));// Assumes column-major
            a_temp[i*lda+i] = 1;
          }
        }
        memset(tau_temp,1,k*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_ormqr__id,flops,m,n,k);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
#endif
#endif
  }
}
void _dormqr_(int matrix_layout , char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 0;
    if (side == 'L'){
      flops = 4.*_m*_n*_k - 2.*_n*_k*_k + 3.*_n*_k;
    } else{
      flops = 4.*_n*_m*_k - 2.*_m*_k*_k + 2.*m*_k + _n*_k - 1./2.*_k*_k + 1./2.*_k;
    }
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_ormqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(c+i*ldc,1,m*sizeof(double));// Assumes column-major
        }
        double* a_temp = (double*)a;
        double* tau_temp = (double*)tau;
        if (side == 'L'){
          for (int i=0; i<k; i++){
            memset(a_temp+i*lda,0,m*sizeof(double));// Assumes column-major
            a_temp[i*lda+i] = 1;
          }
        } else{
          for (int i=0; i<k; i++){
            memset(a_temp+i*lda,0,n*sizeof(double));// Assumes column-major
            a_temp[i*lda+i] = 1;
          }
        }
        memset(tau_temp,1,k*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_ormqr__id,flops,m,n,k);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
#endif
#endif
  }
}
void _sgetri_(int matrix_layout , int n , float * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 4./3.*_n*_n*_n - _n*_n + 5./3.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv)==0);
      else{
        special_time = MPI_Wtime();
        float* ipiv_temp = (float*)ipiv;
        for (int i=0; i<n; i++){
          memset(a+i*lda,0,n*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
          ipiv_temp[i] = i;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_getri__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv)==0);
#endif
#endif
  }
}
void _dgetri_(int matrix_layout , int n , double * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 4./3.*_n*_n*_n - _n*_n + 5./3.*_n;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv)==0);
      else{
        special_time = MPI_Wtime();
        double* ipiv_temp = (double*)ipiv;
        for (int i=0; i<n; i++){
          memset(a+i*lda,0,n*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
          ipiv_temp[i] = i;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_getri__id,flops,n);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv)==0);
#endif
#endif
  }
}
void _stpqrt_(int matrix_layout , int m , int n , int l , int nb , float * a , int lda , float * b , int ldb , float * t , int ldt){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _l = l;
    double flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(t)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,(i+1)*sizeof(float));// Assumes column-major
          memset(a+i*lda+i+1,0,(n-i-1)*sizeof(float));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,(i+1)*sizeof(float));// Assumes column-major
          memset(b+i*ldb+i+1,0,(std::max(0,(m-l)-i-1))*sizeof(float));// Assumes column-major
          memset(b+i*ldb+(m-l),1,std::min(i+1,l)*sizeof(float));// Assumes column-major
          memset(b+i*ldb+(m-l)+i+1,0,std::max(0,(l-i-1))*sizeof(float));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpqrt__id,flops,m,n,l,nb);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
#endif
#endif
  }
}
void _dtpqrt_(int matrix_layout , int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _l = l;
    double flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(t)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,(i+1)*sizeof(double));// Assumes column-major
          memset(a+i*lda+i+1,0,(n-i-1)*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,(i+1)*sizeof(double));// Assumes column-major
          memset(b+i*ldb+i+1,0,(std::max(0,(m-l)-i-1))*sizeof(double));// Assumes column-major
          memset(b+i*ldb+(m-l),1,std::min(i+1,l)*sizeof(double));// Assumes column-major
          memset(b+i*ldb+(m-l)+i+1,0,std::max(0,(l-i-1))*sizeof(double));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpqrt__id,flops,m,n,l,nb);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(t),reinterpret_cast<intptr_t>(v)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
      else{
        special_time = MPI_Wtime();
        float* v_temp = (float*)v;
        float* t_temp = (float*)t;
        if (side == 'L'){
          for (int i=0; i<k; i++){
            memset(v_temp+i*ldv,1,m*sizeof(float));// Assumes column-major
          }
        } else{
          for (int i=0; i<n; i++){
            memset(v_temp+i*ldv,1,n*sizeof(float));// Assumes column-major
          }
        }
        for (int i=0; i<k; i++){
          memset(t_temp+i*ldt,1,nb*sizeof(float));// Assumes column-major
        }
        if (side=='L'){
          for (int i=0; i<n; i++){
            memset(a+i*lda,1,k*sizeof(float));// Assumes column-major
          }
        } else{
          for (int i=0; i<k; i++){
            memset(a+i*lda,1,m*sizeof(float));// Assumes column-major
          }
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,m*sizeof(float));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,ldv,t_temp,ldt,a,lda,b,ldb)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(t),reinterpret_cast<intptr_t>(v)};
    double special_time=0;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
      else{
        special_time = MPI_Wtime();
        double* v_temp = (double*)v;
        double* t_temp = (double*)t;
        if (side == 'L'){
          for (int i=0; i<k; i++){
            memset(v_temp+i*ldv,1,m*sizeof(double));// Assumes column-major
          }
        } else{
          for (int i=0; i<n; i++){
            memset(v_temp+i*ldv,1,n*sizeof(double));// Assumes column-major
          }
        }
        for (int i=0; i<k; i++){
          memset(t_temp+i*ldt,1,nb*sizeof(double));// Assumes column-major
        }
        if (side=='L'){
          for (int i=0; i<n; i++){
            memset(a+i*lda,1,k*sizeof(double));// Assumes column-major
          }
        } else{
          for (int i=0; i<k; i++){
            memset(a+i*lda,1,m*sizeof(double));// Assumes column-major
          }
        }
        for (int i=0; i<n; i++){
          memset(b+i*ldb,1,m*sizeof(double));// Assumes column-major
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,ldv,t_temp,ldt,a,lda,b,ldb)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
  } else{
#ifdef MKL
#ifdef LAPACKE
    assert(LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
#endif
#endif
  }
}

void _blk_to_cyc_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim){
  if (mode){
    volatile double curtime = MPI_Wtime();
    double flops = 0;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(blocked),reinterpret_cast<intptr_t>(cyclic)};
    bool schedule_decision = initiate_comp(ptrs,_CAPITAL_blktocyc__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
    if (schedule_decision){
      int write_idx = 0; int read_idx = 0;
      int offset = num_rows_local*num_columns_local;
      int num_rows_global = num_rows_local*sliceDim;
      int num_columns_global = num_columns_local*sliceDim;
      for (int i=0; i<num_columns_local; i++){
        for (int j=0; j<sliceDim; j++){
          for (int k=0; k<num_rows_local; k++){
            for (int z=0; z<sliceDim; z++){
              read_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
              cyclic[write_idx++] = blocked[read_idx];
            }
          }
        }
      }
      // Remove scalars that should be zeros
      for (int i=0; i<num_columns_global; i++){
        for (int j=i+1; j<num_rows_global; j++){
          cyclic[i*num_rows_global+j]=0.;
        }
      }
    }
    complete_comp(0,ptrs,_CAPITAL_blktocyc__id,flops,num_rows_local,num_columns_local,sliceDim);
  } else{
    int write_idx = 0; int read_idx = 0;
    int offset = num_rows_local*num_columns_local;
    int num_rows_global = num_rows_local*sliceDim;
    int num_columns_global = num_columns_local*sliceDim;
    for (int i=0; i<num_columns_local; i++){
      for (int j=0; j<sliceDim; j++){
        for (int k=0; k<num_rows_local; k++){
          for (int z=0; z<sliceDim; z++){
            read_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
            cyclic[write_idx++] = blocked[read_idx];
          }
        }
      }
    }
    // Remove scalars that should be zeros
    for (int i=0; i<num_columns_global; i++){
      for (int j=i+1; j<num_rows_global; j++){
        cyclic[i*num_rows_global+j]=0.;
      }
    }
  }
}
void _cyc_to_blk_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim){
  if (mode){
    volatile double curtime = MPI_Wtime();
    double flops = 0;
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(blocked),reinterpret_cast<intptr_t>(cyclic)};
    bool schedule_decision = initiate_comp(ptrs,_CAPITAL_cyctoblk__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
    if (schedule_decision){
      int write_idx = 0; int read_idx = 0; int offset = num_rows_local*num_columns_local;
      for (int i=0; i<num_columns_local; i++){
        for (int j=0; j<sliceDim; j++){
          for (int k=0; k<num_rows_local; k++){
            for (int z=0; z<sliceDim; z++){
              write_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
              blocked[write_idx] = cyclic[read_idx++];
            }
          }
        }
      }
    }
    complete_comp(0,ptrs,_CAPITAL_cyctoblk__id,flops,num_rows_local,num_columns_local,sliceDim);
  } else{
    int write_idx = 0; int read_idx = 0; int offset = num_rows_local*num_columns_local;
    for (int i=0; i<num_columns_local; i++){
      for (int j=0; j<sliceDim; j++){
        for (int k=0; k<num_rows_local; k++){
          for (int z=0; z<sliceDim; z++){
            write_idx = z*offset*sliceDim + k + j*offset + i*num_rows_local;
            blocked[write_idx] = cyclic[read_idx++];
          }
        }
      }
    }
  }
}

}
}
