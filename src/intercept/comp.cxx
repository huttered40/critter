#include "comp.h"
#include "../util/util.h"
#include "../dispatch/dispatch.h"

#ifdef MKL
#include "mkl.h"
//#include <mkl_cblas.h>
#else
#define CBLAS_LAYOUT int
#define CBLAS_SIDE int
#define CBLAS_DIAG int
#define CBLAS_TRANSPOSE int
#define CBLAS_UPLO int
#endif /* MKL */

#include <tuple>

#include "func_generators.h"

namespace critter{
namespace internal{

/*
	Interfaces for BLAS/LAPACK APIs for both C/Fortran
*/

// BLAS 1
void _daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy){
  conditional_blas_engine(_BLAS_axpy__id,track_blas1,std::make_tuple(n),std::make_tuple(n),IndexPack<0>{},
                          &cblas_daxpy,n,a,x,incx,y,incy);
}
void _dscal_(const int n , const double a , double *x , const int incx){
  conditional_blas_engine(_BLAS_scal__id,track_blas1,std::make_tuple(n),std::make_tuple(n),IndexPack<0>{},
                          &cblas_dscal,n,a,x,incx);
}

// BLAS 2
void _dgbmv_(const int order, const int trans, const int m, const int n, const int kl, const int ku, const double alpha,
             const double *a, const int lda, const double *x, const int incx, const double beta, double *y, const int incy){
  conditional_blas_engine(_BLAS_gbmv__id,track_blas2,std::make_tuple(m,n,kl,ku,trans*1000+(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,kl,ku,trans*1000+(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2,3,4>{},
                          &cblas_dgbmv,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy);
}
void _dgemv_(const int order, const int trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy ){
  conditional_blas_engine(_BLAS_gemv__id,track_blas2,std::make_tuple(m,n,trans*1000+(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,trans*1000+(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2>{},
                          &cblas_dgemv,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
}
void _dger_(const int order, const int m , const int n , const double alpha , const double *x , const int incx ,
            const double *y , const int incy , double *a , const int lda){
  conditional_blas_engine(_BLAS_ger__id,track_blas2,std::make_tuple(m,n,(alpha!=0)),std::make_tuple(m,n,(alpha!=0)),IndexPack<0,1,2>{},
                          &cblas_dger,(CBLAS_ORDER)order,m,n,alpha,x,incx,y,incy,a,lda);
}
void _dsbmv_(const int Layout, const int uplo, const int n, const int k, const double alpha, const double *a,
             const int lda, const double *x, const int incx, const double beta, double *y, const int incy){
  conditional_blas_engine(_BLAS_sbmv__id,track_blas2,std::make_tuple(n,k,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2,3>{},
                          &cblas_dsbmv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,k,alpha,a,lda,x,incx,beta,y,incy);
}
void _dspmv_(const int Layout, const int uplo, const int n, const double alpha, const double *ap, const double *x,
             const int incx, const double beta, double *y, const int incy){
  conditional_blas_engine(_BLAS_spmv__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2>{},
                          &cblas_dspmv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,ap,x,incx,beta,y,incy);
}
void _dspr_(const int Layout, const int uplo, const int n, const double alpha, const double *x,
            const int incx, double *ap){
  conditional_blas_engine(_BLAS_spr__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),IndexPack<0,1,2>{},
                          &cblas_dspr,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,ap);
}
void _dspr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *ap){
  conditional_blas_engine(_BLAS_spr2__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),IndexPack<0,1,2>{},
                          &cblas_dspr2,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,y,incy,ap);
}
void _dsymv_(const int Layout, const int uplo, const int n, const double alpha, const double *a, const int lda,
            const double *x, const int incx, const double beta, double *y, const int incy){
  conditional_blas_engine(_BLAS_symv__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2>{},
                          &cblas_dsymv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,a,lda,x,incx,beta,y,incy);
}
void _dsyr_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
            double *a, const int lda){
  conditional_blas_engine(_BLAS_syr__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),IndexPack<0,1,2>{},
                          &cblas_dsyr,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,a,lda);
}
void _dsyr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *a, const int lda){
  conditional_blas_engine(_BLAS_syr2__id,track_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),IndexPack<0,1,2>{},
                          &cblas_dsyr2,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,y,incy,a,lda);
}
void _dtrsv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx ){
  conditional_blas_engine(_BLAS_trsv__id,track_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),IndexPack<0,1,2,3>{},
                          &cblas_dtrsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
}
void _dtrmv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx ){
  conditional_blas_engine(_BLAS_trmv__id,track_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),IndexPack<0,1,2,3>{},
                          &cblas_dtrmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
}
void _dtpsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx){
  conditional_blas_engine(_BLAS_tpsv__id,track_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),IndexPack<0,1,2,3>{},
                          &cblas_dtpsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,ap,x,incx);
}
void _dtpmv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx){
  conditional_blas_engine(_BLAS_tpmv__id,track_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),IndexPack<0,1,2,3>{},
                          &cblas_dtpmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,ap,x,incx);
}
void _dtbsv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx){
  conditional_blas_engine(_BLAS_tbsv__id,track_blas2,std::make_tuple(n,k,uplo,trans,diag),std::make_tuple(n,k,uplo,trans,diag),IndexPack<0,1,2,3,4>{},
                          &cblas_dtbsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,k,a,lda,x,incx);
}
void _dtbmv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx){
  conditional_blas_engine(_BLAS_tbmv__id,track_blas2,std::make_tuple(n,k,uplo,trans,diag),std::make_tuple(n,k,uplo,trans,diag),IndexPack<0,1,2,3,4>{},
                          &cblas_dtbmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,k,a,lda,x,incx);
}

// BLAS 3
void _dgemm_(const int order, const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  conditional_blas_engine(_BLAS_gemm__id,track_blas3,std::make_tuple(m,n,k,transa*1000+transb,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,k,transa*1000+transb,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2,3,4>{},
                          &cblas_dgemm,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)transa,(CBLAS_TRANSPOSE)transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void _dtrmm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  conditional_blas_engine(_BLAS_trmm__id,track_blas3,std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)),
                          std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)),IndexPack<0,1,2,3,4>{},
                          &cblas_dtrmm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
}
void _dtrsm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  conditional_blas_engine(_BLAS_trsm__id,track_blas3,
                          (CBLAS_SIDE)side==CblasLeft ? std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)) : std::make_tuple(n,m,side*1000+uplo,transa*1000+diag,(alpha!=0)),
                          std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)), IndexPack<0,1,2,3,4>{},
                          &cblas_dtrsm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
}
void _dsyrk_(const int order, const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  conditional_blas_engine(_BLAS_syrk__id,track_blas3,std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2,3>{},
                          &cblas_dsyrk,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,beta,c,ldc);
}
void _dsyr2k_(const int order, const int uplo, const int trans, const int n, const int k, const double alpha,
              const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
              const int ldc){
  conditional_blas_engine(_BLAS_syr2k__id,track_blas3,std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),IndexPack<0,1,2,3>{},
                          &cblas_dsyr2k,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void _dsymm_(const int order, const int side, const int uplo, const int m, const int n, const double alpha,
             const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
             const int ldc){
  conditional_blas_engine(_BLAS_symm__id,track_blas3,
                          (CBLAS_SIDE)side==CblasLeft ? std::make_tuple(m,n,side*1000+uplo,(alpha!=0)*2+(beta!=0)) : std::make_tuple(n,m,side*1000+uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,side*1000+uplo,(alpha!=0)*2+(beta!=0)), IndexPack<0,1,2,3>{},
                          &cblas_dsymm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc);
}

// **********************************************************************************************************************************
// BLAS 1
void __daxpy__(const int* n , const double* a , const double *x , const int* incx , double *y , const int* incy){
  _daxpy_(*n,*a,x,*incx,y,*incy);
}
void __dscal__(const int* n , const double* a , double *x , const int* incx){
  _dscal_(*n,*a,x,*incx);
}

// BLAS 2
void __dgbmv__(const char* trans , const int* m , const int* n, const int* kl, const int* ku, const double* alpha ,
               const double *a , const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dgbmv_(CblasColMajor,_trans,*m,*n,*kl,*ku,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void __dgemv__(const char* trans , const int* m , const int* n, const double* alpha , const double *a , const int* lda ,
               const double *x, const int* incx , const double* beta, double *y , const int* incy ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dgemv_(CblasColMajor,_trans,*m,*n,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void __dger__(const int* m , const int* n , const double* alpha , const double *x , const int* incx , const double *y ,
              const int* incy , double *a ,
            const int* lda){
  _dger_(CblasColMajor,*m,*n,*alpha,x,*incx,y,*incy,a,*lda);
}
void __dsbmv__(const char* uplo, const int* n, const int* k, const double* alpha, const double *a,
               const int* lda, const double *x, const int* incx, const double* beta, double *y, const int* incy){
  _dsbmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*k,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void __dspmv__(const char* uplo, const int* n, const double* alpha, const double *ap, const double *x,
               const int* incx, const double* beta, double *y, const int* incy){
  _dspmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,ap,x,*incx,*beta,y,*incy);
}
void __dspr__(const char* uplo, const int* n, const double* alpha, const double *x,
              const int* incx, double *ap){
  _dspr_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,ap);
}
void __dspr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *ap){
  _dspr2_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,y,*incy,ap);
}
void __dsymv__(const char* uplo, const int* n, const double* alpha, const double *a, const int* lda,
               const double *x, const int* incx, const double* beta, double *y, const int* incy){
  _dsymv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void __dsyr__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
              double *a, const int* lda){
  _dsyr_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,a,*lda);
}
void __dsyr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *a, const int* lda){
  _dsyr2_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,y,*incy,a,*lda);
}
void __dtrsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *a,
               const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtrsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,a,*lda,x,*incx);
}
void __dtrmv__(const char* uplo , const char* trans , const char* diag , const int* n , const double *a ,
               const int* lda , double *x, const int* incx ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtrmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,a,*lda,x,*incx);
}
void __dtpsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtpsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,ap,x,*incx);
}
void __dtpmv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtpmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,ap,x,*incx);
}
void __dtbsv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtbsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,*k,a,*lda,x,*incx);
}
void __dtbmv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dtbmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,*k,a,*lda,x,*incx);
}

// BLAS 3
void __dgemm__(const char* transa , const char* transb ,
             const int* m , const int* n , const int* k , const double* alpha , const double *a ,
             const int* lda , const double *b , const int* ldb , const double* beta , double *c , const int* ldc){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  CBLAS_TRANSPOSE _transb;
  if (*transb=='T') _transb = CblasTrans; else if (*transb=='N') _transb = CblasNoTrans; else _transb = CblasConjTrans;
  _dgemm_(CblasColMajor,_transa, _transb, *m,*n,*k,*alpha,a,*lda,b,*ldb,*beta,c,*ldc);
}
void __dtrmm__(const char* side , const char* uplo , const char* transa ,
             const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
             const int* lda , double *b , const int* ldb){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  _dtrmm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower), _transa,
          (*diag=='U' ? CblasUnit : CblasNonUnit), *m,*n,*alpha,a,*lda,b,*ldb);
}
void __dtrsm__(const char* side , const char* uplo , const char* transa ,
             const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
             const int* lda , double *b , const int* ldb){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  _dtrsm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower), _transa,
          (*diag=='U' ? CblasUnit : CblasNonUnit), *m,*n,*alpha,a,*lda,b,*ldb);
}
void __dsyrk__(const char* uplo, const char* trans ,
             const int* n , const int* k , const double* alpha , const double *a , const int* lda ,
             const double* beta , double *c , const int* ldc){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dsyrk_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, *n,*k,*alpha,a,*lda,*beta,c,*ldc);
}
void __dsyr2k__(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
                const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
                const int* ldc){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  _dsyr2k_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans,*n,*k,*alpha,a,*lda,b,*ldb,*beta,c,*ldc);
}
void __dsymm__(const char* side, const char* uplo, const int* m, const int* n, const double* alpha,
               const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
               const int* ldc){
  _dsymm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower),*m,*n,*alpha,
          a,*lda,b,*ldb,*beta,c,*ldc);
}

// **********************************************************************************************************************************
// C LAPACK interface
int _dgetrf_(int matrix_layout, int m , int n , double* a , int lda , int* ipiv){
  return conditional_lapack_engine(_LAPACK_getrf__id,track_lapack,
                                 (m>=n ? std::make_tuple(m,n) : std::make_pair(n,m)), std::make_tuple(m,n),IndexPack<0,1>{},
                                 &LAPACKE_dgetrf,std::make_tuple(matrix_layout,m,n,a,lda,ipiv),IndexPack<0,1,2,3,4,5>{},
                                std::make_tuple(a,m,n,lda,1,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 4.*n;}}));
}
int _dpotrf_(int matrix_layout, char uplo , int n , double* a , int lda){
  return conditional_lapack_engine(_LAPACK_potrf__id,track_lapack,
                                   std::make_tuple(n,uplo), std::make_tuple(n,uplo),IndexPack<0,1>{},
                                   &LAPACKE_dpotrf,std::make_tuple(matrix_layout,uplo,n,a,lda),IndexPack<0,1,2,3,4>{},
                                   std::make_tuple(a,n,n,lda,1,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 4.*n;}}));
}
int _dtrtri_(int matrix_layout, char uplo , char diag , int n , double* a , int lda){
  return conditional_lapack_engine(_LAPACK_trtri__id,track_lapack,
                                   std::make_tuple(n,uplo,diag), std::make_tuple(n,uplo,diag),IndexPack<0,1,2>{},
                                   &LAPACKE_dtrtri,std::make_tuple(matrix_layout,uplo,diag,n,a,lda),IndexPack<0,1,2,3,4,5>{},
                                   std::make_tuple(a,n,n,lda,1,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 4.*n;}}));
}

int _dgeqrf_(int matrix_layout, int m , int n , double* a , int lda , double* tau){
  return conditional_lapack_engine(_LAPACK_geqrf__id,track_lapack,
                                   (m>=n ? std::make_tuple(m,n) : std::make_tuple(n,m)),std::make_tuple(m,n),IndexPack<0,1>{},
                                   &LAPACKE_dgeqrf,std::make_tuple(matrix_layout,m,n,a,lda,tau),IndexPack<0,1,2,3,4,5>{},
                                   std::make_tuple(a,m,n,lda,1,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 4.*n;}}));
}
int _dorgqr_(int matrix_layout, int m , int n , int k , double* a , int lda , const double* tau){
  return conditional_lapack_engine(_LAPACK_orgqr__id,track_lapack,
                                   std::make_tuple(m,n,k),std::make_tuple(m,n,k),IndexPack<0,1,2>{},
                                   &LAPACKE_dorgqr,std::make_tuple(matrix_layout,m,n,k,a,lda,tau),IndexPack<0,1,2,3,4,5,6>{},
                                   std::make_tuple(a,m,n,lda,0,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 1.;}}),
                                   std::make_tuple((double*)tau,k,1,k,1,[](double* a, int m, int n, int lda){}));
}
int _dormqr_(int matrix_layout, char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  return conditional_lapack_engine(_LAPACK_ormqr__id,track_lapack,
                                   (side == 'L' ? std::make_tuple(m,n,k,side,trans) : std::make_tuple(n,m,k,side,trans)),
                                   std::make_tuple(m,n,k,side,trans),IndexPack<0,1,2,3,4>{},
                                   &LAPACKE_dormqr,std::make_tuple(matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc),IndexPack<0,1,2,3,4,5,6,7,8,9,10>{},
                                   std::make_tuple(c,m,n,ldc,1,[](double* a, int m, int n, int lda){}),
                                   std::make_tuple((double*)a,(side=='L'?m:n),k,lda,0,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 1.;}}),
                                   std::make_tuple((double*)tau,k,1,k,1,[](double* a, int m, int n, int lda){}));
}
int _dgetri_(int matrix_layout, int n , double * a , int lda , const int * ipiv){
  return conditional_lapack_engine(_LAPACK_getri__id,track_lapack,
                                   std::make_tuple(n),std::make_tuple(n),IndexPack<0>{},
                                   &LAPACKE_dgetri,std::make_tuple(matrix_layout,n,a,lda,ipiv),IndexPack<0,1,2,3,4>{},
                                   std::make_tuple(a,n,n,lda,0,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i*lda+i] = 1.;}}),
                                   std::make_tuple((double*)ipiv,n,1,n,1,[](double* a, int m, int n, int lda){for (int i=0; i<n; i++){a[i] = i;}}));
}
int _dtpqrt_(int matrix_layout, int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  return conditional_lapack_engine_tpqrt_(&LAPACKE_dtpqrt,matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
}
int _dtpmqrt_(int matrix_layout, char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb){
  return conditional_lapack_engine_tpmqrt_(&LAPACKE_dtpmqrt,matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
}

// FORTRAN LAPACK interface -> just call the C interface.
void __dgetrf__(const int* m , const int* n , double* a , const int* lda , int* ipiv, int* info){
  *info = _dgetrf_(LAPACK_COL_MAJOR,*m,*n,a,*lda,ipiv);
}
void __dpotrf__(const char* uplo , const int* n , double* a , const int* lda, int* info){
  *info = _dpotrf_(LAPACK_COL_MAJOR,*uplo,*n,a,*lda);
}
void __dtrtri__(const char* uplo , const char* diag , const int* n , double* a , const int* lda, int* info){
  *info = _dtrtri_(LAPACK_COL_MAJOR,*uplo,*diag,*n,a,*lda);
}
void __dgeqrf__(const int* m , const int* n , double* a , const int* lda , double* tau, double* work, const int* lwork, int* info){
  *info = _dgeqrf_(LAPACK_COL_MAJOR,*m,*n,a,*lda,tau);
}
void __dorgqr__(const int* m , const int* n , const int* k , double* a , const int* lda , const double* tau, double* work, const int* lwork, int* info){
  *info = _dorgqr_(LAPACK_COL_MAJOR,*m,*n,*k,a,*lda,tau);
}
void __dormqr__(const char* side , const char* trans , const int* m , const int* n , const int* k , const double * a , const int* lda , const double * tau ,
                double * c , const int* ldc, double* work, const int* lwork, int* info){
  *info = _dormqr_(LAPACK_COL_MAJOR,*side,*trans,*m,*n,*k,a,*lda,tau,c,*ldc);
}
void __dgetri__(const int* n , double * a , const int* lda , const int * ipiv, double* work, const int* lwork, int* info){
  *info = _dgetri_(LAPACK_COL_MAJOR,*n,a,*lda,ipiv);
}
void __dtpqrt__(const int* m , const int* n , const int* l , const int* nb , double * a , const int* lda , double* b , const int* ldb , double * t , const int* ldt,
                double* work, int* info){
  *info = _dtpqrt_(LAPACK_COL_MAJOR,*m,*n,*l,*nb,a,*lda,b,*ldb,t,*ldt);
}
void __dtpmqrt__(const char* side , const char* trans , const int* m , const int* n , const int* k , const int* l , const int* nb , const double * v ,
                 const int* ldv , const double * t , const int* ldt , double * a , const int* lda , double * b , const int* ldb, double* work, int* info){
  *info = _dtpmqrt_(LAPACK_COL_MAJOR,*side,*trans,*m,*n,*k,*l,*nb,v,*ldv,t,*ldt,a,*lda,b,*ldb);
}

}
}
