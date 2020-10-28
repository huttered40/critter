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
    bool schedule_decision = initiate_comp(ptrs,_BLAS_gemm__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef CBLAS
    if (schedule_decision) cblas_dgemm(Layout,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
#endif
#endif
    complete_comp(0,ptrs,_BLAS_gemm__id,flops,m,n,k);
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
    double flops = 2*_k*_n*_n;//Note: the factor of 2 might be wrong
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
    double flops = 2*_k*_n*_n;//Note: the factor of 2 might be wrong
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
    double flops = 2./3. * _n*_n*_n; if (m > n) flops = (1./3.)*_n*_n*(3*_m-_n); if (m < n) flops = (1./3.)*_m*_m*(3*_n-_m);
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_sgetrf(matrix_layout,m,n,a,lda,ipiv);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GETRF -- must force nonsingular
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(float));// assumes column-major layout
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dgetrf(matrix_layout,m,n,a,lda,ipiv);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GETRF -- must force nonsingular
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(double));// assumes column-major layout
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_potrf__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_spotrf(matrix_layout,uplo,n,a,lda);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // POTRF -- must force positive-definiteness
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(float));
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_potrf__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dpotrf(matrix_layout,uplo,n,a,lda);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // POTRF -- must force positive-definiteness
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(double));
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_trtri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_strtri(matrix_layout,uplo,diag,n,a,lda);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TRTRI -- must force positive-definiteness
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(float));
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_trtri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dtrtri(matrix_layout,uplo,diag,n,a,lda);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TRTRI -- must force positive-definiteness
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(double));
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_geqrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
      int ret = LAPACKE_sgeqrf(matrix_layout,m,n,a,lda,tau);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GEQRF -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,1,m*n*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sgeqrf(matrix_layout,m,n,a,/*lda*/m,tau)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_geqrf__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_geqrf__id,curtime,flops,m,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dgeqrf(matrix_layout,m,n,a,lda,tau);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GEQRF -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,1,m*n*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgeqrf(matrix_layout,m,n,a,/*lda*/m,tau)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_geqrf__id,flops,m,n);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_orgqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_sorgqr(matrix_layout,m,n,k,a,lda,tau);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // ORGQR -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        float* tau_temp = (float*)tau;//(float*)malloc(k*sizeof(float));
        memset(a,1,/*lda*/m*n*sizeof(float));// Assumes column-major
        memset(tau_temp,1,k*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_sorgqr(matrix_layout,m,n,k,a,/*lda*/m,tau_temp)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_orgqr__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_orgqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dorgqr(matrix_layout,m,n,k,a,lda,tau);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // ORGQR -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        double* tau_temp = (double*)tau;//(double*)malloc(k*sizeof(double));
        memset(a,1,/*lda*/m*n*sizeof(double));// Assumes column-major
        memset(tau_temp,1,k*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dorgqr(matrix_layout,m,n,k,a,/*lda*/m,tau_temp)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_orgqr__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_ormqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // ORMQR -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        float* a_temp = (float*)a;//malloc(lda*k*sizeof(float));
        float* tau_temp = (float*)tau;//malloc(k*sizeof(float));
        memset(c,1,/*ldc*/m*n*sizeof(float));// Assumes column-major
        if (side == 'L'){
          memset(a_temp,1,/*lda*/m*k*sizeof(float));// Assumes column-major
        } else{
          memset(a_temp,1,/*lda*/n*k*sizeof(float));// Assumes column-major
        }
        memset(tau_temp,1,k*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        if (side == 'L'){
          assert(LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a_temp,/*lda*/m,tau_temp,c,/*ldc*/m)==0);
        } else{
          assert(LAPACKE_sormqr(matrix_layout,side,trans,m,n,k,a_temp,/*lda*/n,tau_temp,c,/*ldc*/m)==0);
        }
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_ormqr__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_ormqr__id,curtime,flops,m,n,k);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // ORMQR -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        double* a_temp = (double*)a;//malloc(lda*k*sizeof(double));
        double* tau_temp = (double*)tau;//malloc(k*sizeof(double));
        memset(c,1,/*ldc*/m*n*sizeof(double));// Assumes column-major
        if (side == 'L'){
          memset(a_temp,1,/*lda*/m*k*sizeof(double));// Assumes column-major
        } else{
          memset(a_temp,1,/*lda*/n*k*sizeof(double));// Assumes column-major
        }
        memset(tau_temp,1,k*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        if (side == 'L'){
          assert(LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a_temp,/*lda*/m,tau_temp,c,/*ldc*/m)==0);
        } else{
          assert(LAPACKE_dormqr(matrix_layout,side,trans,m,n,k,a_temp,/*lda*/n,tau_temp,c,/*ldc*/m)==0);
        }
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_ormqr__id,flops,m,n,k);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_sgetri(matrix_layout,n,a,lda,ipiv);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GETRI -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(float));// Assumes column-major
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
          //ipiv[i] = i;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_getri__id,curtime,flops,n);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dgetri(matrix_layout,n,a,lda,ipiv);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // GETRI -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,0,lda*n*sizeof(double));// Assumes column-major
        for (int i=0; i<n; i++){
          a[i*lda+i] = 1;
          //ipiv[i] = i;
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(t)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TPQRT -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,1,/*lda*/n*n*sizeof(float));// Assumes column-major
        memset(b,1,/*ldb*/m*n*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_stpqrt(matrix_layout,m,n,l,nb,a,/*lda*/n,b,/*ldb*/m,t,ldt)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpqrt__id,flops,m,n,l,nb);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(t)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TPQRT -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        memset(a,1,/*lda*/n*n*sizeof(double));// Assumes column-major
        memset(b,1,/*ldb*/m*n*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dtpqrt(matrix_layout,m,n,l,nb,a,/*lda*/n,b,/*ldb*/m,t,ldt)==0);
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpqrt__id,flops,m,n,l,nb);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(t),reinterpret_cast<intptr_t>(v)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TPMQRT -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        float* v_temp = (float*)v;
        float* t_temp = (float*)t;
        if (side == 'L'){
          memset(v_temp,1,/*ldv*/m*k*sizeof(float));// Assumes column-major
        } else{
          memset(v_temp,1,/*ldv*/n*n*sizeof(float));// Assumes column-major
        }
        memset(t_temp,1,/*ldt*/nb*k*sizeof(float));// Assumes column-major
        if (side=='L') memset(a,1,/*lda*/k*n*sizeof(float));// Assumes column-major
        if (side=='R') memset(a,1,/*lda*/m*k*sizeof(float));// Assumes column-major
        memset(b,1,/*ldb*/m*n*sizeof(float));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        if (side == 'L'){
          assert(LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,/*ldv*/m,t_temp,/*ldt*/nb,a,/*lda*/k,b,/*ldb*/m)==0);
        } else{
          assert(LAPACKE_stpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,/*ldv*/n,t_temp,/*ldt*/nb,a,/*lda*/m,b,/*ldb*/m)==0);
        }
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
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
    std::vector<intptr_t> ptrs = {reinterpret_cast<intptr_t>(a),reinterpret_cast<intptr_t>(b),reinterpret_cast<intptr_t>(t),reinterpret_cast<intptr_t>(v)};
    double special_time;
    bool schedule_decision = initiate_comp(ptrs,_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
#ifdef MKL
#ifdef LAPACKE
    if (schedule_decision){
      int ret = LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
      if ((ret!=0) && (mechanism != 1)) assert(0);
      if (ret != 0){
        // Specialized checks:
        // TPMQRT -- must force values to be legal
        //       -- Safest approach: Overwrite matrix to form identity matrix
        special_time = MPI_Wtime();
        double* v_temp = (double*)v;
        double* t_temp = (double*)t;
        if (side == 'L'){
          memset(v_temp,1,/*ldv*/m*k*sizeof(double));// Assumes column-major
        } else{
          memset(v_temp,1,/*ldv*/n*n*sizeof(double));// Assumes column-major
        }
        memset(t_temp,1,/*ldt*/nb*k*sizeof(double));// Assumes column-major
        if (side=='L') memset(a,1,/*lda*/k*n*sizeof(double));// Assumes column-major
        if (side=='R') memset(a,1,/*lda*/m*k*sizeof(double));// Assumes column-major
        memset(b,1,/*ldb*/m*n*sizeof(double));// Assumes column-major
        special_time = MPI_Wtime() - special_time;
        if (side == 'L'){
          assert(LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,/*ldv*/m,t_temp,/*ldt*/nb,a,/*lda*/k,b,/*ldb*/m)==0);
        } else{
          assert(LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v_temp,/*ldv*/n,t_temp,/*ldt*/nb,a,/*lda*/m,b,/*ldb*/m)==0);
        }
      }
    }
#endif
#endif
    complete_comp(special_time,ptrs,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
  } else{
#ifdef MKL
#ifdef LAPACKE
    LAPACKE_dtpmqrt(matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
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

}
}
