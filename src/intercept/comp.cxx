#include "comp.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

#ifdef MKL
// Note: this MKL inclusion should be conditional on config.mk
#include "mkl.h"
//#include <mkl_cblas.h>
#else
#define CBLAS_LAYOUT int
#define CBLAS_SIDE int
#define CBLAS_DIAG int
#define CBLAS_TRANSPOSE int
#define CBLAS_UPLO int
#endif /* MKL */

namespace critter{
namespace internal{

void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy){
  if (mode && track_blas1){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 2.*_n;
    bool schedule_decision = initiate_comp(_BLAS_axpy__id,curtime,flops,n);
    if (schedule_decision) cblas_daxpy(n,a,x,incx,y,incy);
    complete_comp(0,_BLAS_axpy__id,flops,n);
  } else{
    cblas_daxpy(n,a,x,incx,y,incy);
  }
}
void _dscal_(const int n , const double a , double *x , const int incx){
  if (mode && track_blas1){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1.*_n;
    bool schedule_decision = initiate_comp(_BLAS_scal__id,curtime,flops,n);
    if (schedule_decision) cblas_dscal(n,a,x,incx);
    complete_comp(0,_BLAS_scal__id,flops,n);
  } else{
    cblas_dscal(n,a,x,incx);
  }
}

void _dger_(const int m , const int n , const double alpha , const double *x , const int incx , const double *y , const int incy , double *a ,
            const int lda){
  if (mode && track_blas2){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = 2.*_m*_n+_m+_n;
    bool schedule_decision = initiate_comp(_BLAS_ger__id,curtime,flops,m,n);
    if (schedule_decision) cblas_dger(CblasColMajor,m,n,alpha,x,incx,y,incy,a,lda);
    complete_comp(0,_BLAS_ger__id,flops,m,n);
  } else{
    cblas_dger(CblasColMajor,m,n,alpha,x,incx,y,incy,a,lda);
  }
}
void _dgemv_(const int trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy ){
  if (mode && track_blas2){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = _m*_n+_m+_n;
    bool schedule_decision = initiate_comp(_BLAS_gemv__id,curtime,flops,m,n);
    if (schedule_decision) cblas_dgemv(CblasColMajor,(CBLAS_TRANSPOSE)trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
    complete_comp(0,_BLAS_gemv__id,flops,m,n);
  } else{
    cblas_dgemv(CblasColMajor,(CBLAS_TRANSPOSE)trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
  }
}

void _dtrmv_(const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx ){
  if (mode && track_blas2){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = (0.5*(_n+1)*_n + 2*_n);
    bool schedule_decision = initiate_comp(_BLAS_trmv__id,curtime,flops,n);
    if (schedule_decision) cblas_dtrmv(CblasColMajor,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
    complete_comp(0,_BLAS_trmv__id,flops,n);
  } else{
    cblas_dtrmv(CblasColMajor,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
  }
}

void _dgemm_(const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  if (mode && track_blas3){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m; double _k = k;
    double flops = 2.*_m*_n*_k;
    double special_time=0;
    bool schedule_decision = initiate_comp(_BLAS_gemm__id,curtime,flops,m,n,k);
    if (schedule_decision){
      cblas_dgemm(CblasColMajor,(CBLAS_TRANSPOSE)transa,(CBLAS_TRANSPOSE)transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    }
    complete_comp(special_time,_BLAS_gemm__id,flops,m,n,k);
  } else{
    cblas_dgemm(CblasColMajor,(CBLAS_TRANSPOSE)transa,(CBLAS_TRANSPOSE)transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }
}
void _dtrmm_(const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas3){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = (CBLAS_SIDE)side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    bool schedule_decision = initiate_comp(_BLAS_trmm__id,curtime,flops,m,n);
    if (schedule_decision) cblas_dtrmm(CblasColMajor,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
    complete_comp(0,_BLAS_trmm__id,flops,m,n);
  } else{
    cblas_dtrmm(CblasColMajor,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _dtrsm_(const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  if (mode && track_blas3){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _m = m;
    double flops = (CBLAS_SIDE)side==CblasLeft ? _m*_m*_n : _m*_n*_n;// Note: might want an extra factor of 2.
    bool schedule_decision = initiate_comp(_BLAS_trsm__id,curtime,flops,m,n);
    if (schedule_decision) cblas_dtrsm(CblasColMajor,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
    complete_comp(0,_BLAS_trsm__id,flops,m,n);
  } else{
    cblas_dtrsm(CblasColMajor,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
  }
}
void _dsyrk_(const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  if (mode && track_blas3){
    volatile double curtime = MPI_Wtime();
    double _n = n; double _k = k;
    double flops = _k*_n*(_n+1);
    bool schedule_decision = initiate_comp(_BLAS_syrk__id,curtime,flops,n,k);
    if (schedule_decision) cblas_dsyrk(CblasColMajor,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,beta,c,ldc);
    complete_comp(0,_BLAS_syrk__id,flops,n,k);
  } else{
    cblas_dsyrk(CblasColMajor,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,beta,c,ldc);
  }
}

void __dgemv__(const char trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy ){
  CBLAS_TRANSPOSE _trans;
  if (trans=='T') _trans = CblasTrans; else if (trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dgemv_(_trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
}
void __dtrmv__(const char uplo , const char trans , const char diag , const int n , const double *a , const int lda , double *x, const int incx ){
  CBLAS_TRANSPOSE _trans;
  if (trans=='T') _trans = CblasTrans; else if (trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtrmv_((uplo=='U' ? CblasUpper : CblasLower), _trans, (diag=='U' ? CblasUnit : CblasNonUnit), n,a,lda,x,incx);
}
void __dgemm__(const char transa , const char transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  CBLAS_TRANSPOSE _transa;
  if (transa=='T') _transa = CblasTrans; else if (transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  CBLAS_TRANSPOSE _transb;
  if (transb=='T') _transb = CblasTrans; else if (transb=='N') _transb = CblasNoTrans; else _transb = CblasConjTrans;
  _dgemm_(_transa, _transb, m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void __dtrmm__(const char side , const char uplo , const char transa ,
             const char diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  CBLAS_TRANSPOSE _transa;
  if (transa=='T') _transa = CblasTrans; else if (transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  _dtrmm_((side=='L' ? CblasLeft : CblasRight), (uplo=='U' ? CblasUpper : CblasLower), _transa,
          (diag=='U' ? CblasUnit : CblasNonUnit), m,n,alpha,a,lda,b,ldb);
}
void __dtrsm__(const char side , const char uplo , const char transa ,
             const char diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  CBLAS_TRANSPOSE _transa;
  if (transa=='T') _transa = CblasTrans; else if (transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  _dtrsm_((side=='L' ? CblasLeft : CblasRight), (uplo=='U' ? CblasUpper : CblasLower), _transa,
          (diag=='U' ? CblasUnit : CblasNonUnit), m,n,alpha,a,lda,b,ldb);
}
void __dsyrk__(const char uplo , const char trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  CBLAS_TRANSPOSE _trans;
  if (trans=='T') _trans = CblasTrans; else if (trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dsyrk_((uplo=='U' ? CblasUpper : CblasLower), _trans, n,k,alpha,a,lda,beta,c,ldc);
}

void _dgetrf_(int m , int n , double* a , int lda , int* ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 0;
    if (m>=n){
      flops = _m*_n*_n - 1./3.*_n*_n*_n - 1./2.*_n*_n + 5./6.*_n;
    } else{
      flops = _n*_m*_m - 1./3.*_m*_m*_m - 1./2.*_m*_m + 5./6.*_m;
    }
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_getrf__id,curtime,flops,m,n);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgetrf(LAPACK_COL_MAJOR,m,n,a,lda,ipiv)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgetrf(LAPACK_COL_MAJOR,m,n,a,lda,ipiv)==0);
      }
    }
    complete_comp(special_time,_LAPACK_getrf__id,flops,m,n);
  } else{
    assert(LAPACKE_dgetrf(LAPACK_COL_MAJOR,m,n,a,lda,ipiv)==0);
  }
}
void _dpotrf_(char uplo , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 1./2.*_n*_n + 1./6.*_n;
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_potrf__id,curtime,flops,n);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dpotrf(LAPACK_COL_MAJOR,uplo,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dpotrf(LAPACK_COL_MAJOR,uplo,n,a,lda)==0);
      }
    }
    complete_comp(special_time,_LAPACK_potrf__id,flops,n);
  } else{
    assert(LAPACKE_dpotrf(LAPACK_COL_MAJOR,uplo,n,a,lda)==0);
  }
}
void _dtrtri_(char uplo , char diag , int n , double* a , int lda){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 1./3.*_n*_n*_n + 2./3.*_n;
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_trtri__id,curtime,flops,n);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtrtri(LAPACK_COL_MAJOR,uplo,diag,n,a,lda)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,n*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dtrtri(LAPACK_COL_MAJOR,uplo,diag,n,a,lda)==0);
      }
    }
    complete_comp(special_time,_LAPACK_trtri__id,flops,n);
  } else{
    assert(LAPACKE_dtrtri(LAPACK_COL_MAJOR,uplo,diag,n,a,lda)==0);
  }
}
void _dgeqrf_(int m , int n , double* a , int lda , double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n;
    double flops = 0;
    if (m>=n){
      flops = 2.*_m*_n*_n - 2./3.*_n*_n*_n + _m*_n + _n*_n + 14./3.*_n;
    } else{
      flops = 2.*_n*_m*_m - 2./3.*_m*_m*_m + 3.*_m*_n - _m*_m + 14./3.*_m;
    }
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_geqrf__id,curtime,flops,m,n);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgeqrf(LAPACK_COL_MAJOR,m,n,a,lda,tau)==0);
      else{
        special_time = MPI_Wtime();
        for (int i=0; i<n; i++){
          memset(a+i*lda,1,m*sizeof(double));// Assumes column-major
        }
        for (int i=0; i<n; i++){
          a[i*lda+i] = 4.*n;
        }
        special_time = MPI_Wtime() - special_time;
        assert(LAPACKE_dgeqrf(LAPACK_COL_MAJOR,m,n,a,lda,tau)==0);
      }
    }
    complete_comp(special_time,_LAPACK_geqrf__id,flops,m,n);
  } else{
    assert(LAPACKE_dgeqrf(LAPACK_COL_MAJOR,m,n,a,lda,tau)==0);
  }
}
void _dorgqr_(int m , int n , int k , double* a , int lda , const double* tau){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 4.*_m*_n*_k - 2.*(_m+_n)*_k*_k + (4./3.)*_k*_k*_k + 3.*_n*_k - _m*_k - _k*_k - 4./3.*_k;
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_orgqr__id,curtime,flops,m,n,k);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dorgqr(LAPACK_COL_MAJOR,m,n,k,a,lda,tau)==0);
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
        assert(LAPACKE_dorgqr(LAPACK_COL_MAJOR,m,n,k,a,lda,tau_temp)==0);
      }
    }
    complete_comp(special_time,_LAPACK_orgqr__id,flops,m,n,k);
  } else{
    assert(LAPACKE_dorgqr(LAPACK_COL_MAJOR,m,n,k,a,lda,tau)==0);
  }
}
void _dormqr_(char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 0;
    if (side == 'L'){
      flops = 4.*_m*_n*_k - 2.*_n*_k*_k + 3.*_n*_k;
    } else{
      flops = 4.*_n*_m*_k - 2.*_m*_k*_k + 2.*m*_k + _n*_k - 1./2.*_k*_k + 1./2.*_k;
    }
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_ormqr__id,curtime,flops,m,n,k);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dormqr(LAPACK_COL_MAJOR,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
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
        assert(LAPACKE_dormqr(LAPACK_COL_MAJOR,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
      }
    }
    complete_comp(special_time,_LAPACK_ormqr__id,flops,m,n,k);
  } else{
    assert(LAPACKE_dormqr(LAPACK_COL_MAJOR,side,trans,m,n,k,a,lda,tau,c,ldc)==0);
  }
}
void _dgetri_(int n , double * a , int lda , const int * ipiv){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _n = n;
    double flops = 4./3.*_n*_n*_n - _n*_n + 5./3.*_n;
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_getri__id,curtime,flops,n);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dgetri(LAPACK_COL_MAJOR,n,a,lda,ipiv)==0);
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
        assert(LAPACKE_dgetri(LAPACK_COL_MAJOR,n,a,lda,ipiv)==0);
      }
    }
    complete_comp(special_time,_LAPACK_getri__id,flops,n);
  } else{
    assert(LAPACKE_dgetri(LAPACK_COL_MAJOR,n,a,lda,ipiv)==0);
  }
}
void _dtpqrt_(int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _l = l;
    double flops = 2.*_m*_n*_l;//Note: this is an educated guess. There is no information on this flop count
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_tpqrt__id,curtime,flops,m,n,l,nb);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtpqrt(LAPACK_COL_MAJOR,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
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
        assert(LAPACKE_dtpqrt(LAPACK_COL_MAJOR,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
      }
    }
    complete_comp(special_time,_LAPACK_tpqrt__id,flops,m,n,l,nb);
  } else{
    assert(LAPACKE_dtpqrt(LAPACK_COL_MAJOR,m,n,l,nb,a,lda,b,ldb,t,ldt)==0);
  }
}
void _dtpmqrt_(char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb){
  if (mode && track_lapack){
    volatile double curtime = MPI_Wtime();
    double _m = m; double _n = n; double _k = k;
    double flops = 2.*_m*_n*_k;//Note: this is an educated guess. There is no information on this flop count
    double special_time=0;
    bool schedule_decision = initiate_comp(_LAPACK_tpmqrt__id,curtime,flops,m,n,k,l,nb);
    if (schedule_decision){
      if (mechanism == 0 && autotuning_debug==0) assert(LAPACKE_dtpmqrt(LAPACK_COL_MAJOR,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
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
        assert(LAPACKE_dtpmqrt(LAPACK_COL_MAJOR,side,trans,m,n,k,l,nb,v_temp,ldv,t_temp,ldt,a,lda,b,ldb)==0);
      }
    }
    complete_comp(special_time,_LAPACK_tpmqrt__id,flops,m,n,k,l,nb);
  } else{
    assert(LAPACKE_dtpmqrt(LAPACK_COL_MAJOR,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb)==0);
  }
}

void _blk_to_cyc_rect_(double* blocked, double* cyclic, int num_rows_local, int num_columns_local, int sliceDim){
  if (mode){
    volatile double curtime = MPI_Wtime();
    double flops = 0;
    bool schedule_decision = initiate_comp(_CAPITAL_blktocyc__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
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
    complete_comp(0,_CAPITAL_blktocyc__id,flops,num_rows_local,num_columns_local,sliceDim);
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
    bool schedule_decision = initiate_comp(_CAPITAL_cyctoblk__id,curtime,flops,num_rows_local,num_columns_local,sliceDim);
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
    complete_comp(0,_CAPITAL_cyctoblk__id,flops,num_rows_local,num_columns_local,sliceDim);
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
