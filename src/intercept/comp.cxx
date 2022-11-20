#include "comp.h"
#include "../util.h"

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

#include <tuple>

#include "func_generators.h"

// BLAS 1
void critter_daxpy_(const int n , const double a , const double *x , const int incx , double *y , const int incy){
  internal::engine(0,internal::profile_blas1,std::make_tuple(n),std::make_tuple(n),internal::IndexPack<0>{},
                          &cblas_daxpy,n,a,x,incx,y,incy);
}
void critter_dscal_(const int n , const double a , double *x , const int incx){
  internal::engine(1,internal::profile_blas1,std::make_tuple(n),std::make_tuple(n),internal::IndexPack<0>{},
                          &cblas_dscal,n,a,x,incx);
}

// BLAS 2
void critter_dgbmv_(const int order, const int trans, const int m, const int n, const int kl, const int ku, const double alpha,
             const double *a, const int lda, const double *x, const int incx, const double beta, double *y, const int incy){
  internal::engine(20,internal::profile_blas2,std::make_tuple(m,n,kl,ku,trans*1000+(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,kl,ku,trans*1000+(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dgbmv,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy);
}
void critter_dgemv_(const int order, const int trans , const int m , const int n, const double alpha , const double *a , const int lda , const double *x, const int incx ,
             const double beta, double *y , const int incy ){
  internal::engine(21,internal::profile_blas2,std::make_tuple(m,n,trans*1000+(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,trans*1000+(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dgemv,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
}
void critter_dger_(const int order, const int m , const int n , const double alpha , const double *x , const int incx ,
            const double *y , const int incy , double *a , const int lda){
  internal::engine(22,internal::profile_blas2,std::make_tuple(m,n,(alpha!=0)),std::make_tuple(m,n,(alpha!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dger,(CBLAS_ORDER)order,m,n,alpha,x,incx,y,incy,a,lda);
}
void critter_dsbmv_(const int Layout, const int uplo, const int n, const int k, const double alpha, const double *a,
             const int lda, const double *x, const int incx, const double beta, double *y, const int incy){
  internal::engine(23,internal::profile_blas2,std::make_tuple(n,k,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2,3>{},
                          &cblas_dsbmv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,k,alpha,a,lda,x,incx,beta,y,incy);
}
void critter_dspmv_(const int Layout, const int uplo, const int n, const double alpha, const double *ap, const double *x,
             const int incx, const double beta, double *y, const int incy){
  internal::engine(24,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dspmv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,ap,x,incx,beta,y,incy);
}
void critter_dspr_(const int Layout, const int uplo, const int n, const double alpha, const double *x,
            const int incx, double *ap){
  internal::engine(25,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dspr,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,ap);
}
void critter_dspr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *ap){
  internal::engine(26,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dspr2,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,y,incy,ap);
}
void critter_dsymv_(const int Layout, const int uplo, const int n, const double alpha, const double *a, const int lda,
            const double *x, const int incx, const double beta, double *y, const int incy){
  internal::engine(27,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,uplo,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dsymv,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,a,lda,x,incx,beta,y,incy);
}
void critter_dsyr_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
            double *a, const int lda){
  internal::engine(28,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dsyr,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,a,lda);
}
void critter_dsyr2_(const int Layout, const int uplo, const int n, const double alpha, const double *x, const int incx,
             const double *y, const int incy, double *a, const int lda){
  internal::engine(29,internal::profile_blas2,std::make_tuple(n,uplo,(alpha!=0)),std::make_tuple(n,uplo,(alpha!=0)),internal::IndexPack<0,1,2>{},
                          &cblas_dsyr2,(CBLAS_ORDER)Layout,(CBLAS_UPLO)uplo,n,alpha,x,incx,y,incy,a,lda);
}
void critter_dtrsv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx ){
  internal::engine(30,internal::profile_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),internal::IndexPack<0,1,2,3>{},
                          &cblas_dtrsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
}
void critter_dtrmv_(const int order, const int uplo , const int trans , const int diag , const int n , const double *a , const int lda , double *x, const int incx ){
  internal::engine(31,internal::profile_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),internal::IndexPack<0,1,2,3>{},
                          &cblas_dtrmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,a,lda,x,incx);
}
void critter_dtpsv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx){
  internal::engine(32,internal::profile_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),internal::IndexPack<0,1,2,3>{},
                          &cblas_dtpsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,ap,x,incx);
}
void critter_dtpmv_(const int order, const int uplo, const int trans, const int diag, const int n, const double *ap,
             double *x, const int incx){
  internal::engine(33,internal::profile_blas2,std::make_tuple(n,uplo,trans,diag),std::make_tuple(n,uplo,trans,diag),internal::IndexPack<0,1,2,3>{},
                          &cblas_dtpmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,ap,x,incx);
}
void critter_dtbsv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx){
  internal::engine(34,internal::profile_blas2,std::make_tuple(n,k,uplo,trans,diag),std::make_tuple(n,k,uplo,trans,diag),internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dtbsv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,k,a,lda,x,incx);
}
void critter_dtbmv_(const int order, const int uplo, const int trans, const int diag, const int n, const int k,
             const double *a, const int lda, double *x, const int incx){
  internal::engine(35,internal::profile_blas2,std::make_tuple(n,k,uplo,trans,diag),std::make_tuple(n,k,uplo,trans,diag),internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dtbmv,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,(CBLAS_DIAG)diag,n,k,a,lda,x,incx);
}

// BLAS 3
void critter_dgemm_(const int order, const int transa , const int transb ,
             const int m , const int n , const int k , const double alpha , const double *a ,
             const int lda , const double *b , const int ldb , const double beta , double *c , const int ldc){
  internal::engine(50,internal::profile_blas3,std::make_tuple(m,n,k,transa*1000+transb,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,k,transa*1000+transb,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dgemm,(CBLAS_ORDER)order,(CBLAS_TRANSPOSE)transa,(CBLAS_TRANSPOSE)transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void critter_dtrmm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  internal::engine(51,internal::profile_blas3,std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)),
                          std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)),internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dtrmm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
}
void critter_dtrsm_(const int order, const int side , const int uplo , const int transa ,
             const int diag , const int m , const int n , const double alpha , const double *a ,
             const int lda , double *b , const int ldb){
  internal::engine(52,internal::profile_blas3,
                          (CBLAS_SIDE)side==CblasLeft ? std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)) : std::make_tuple(n,m,side*1000+uplo,transa*1000+diag,(alpha!=0)),
                          std::make_tuple(m,n,side*1000+uplo,transa*1000+diag,(alpha!=0)), internal::IndexPack<0,1,2,3,4>{},
                          &cblas_dtrsm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)transa,(CBLAS_DIAG)diag,m,n,alpha,a,lda,b,ldb);
}
void critter_dsyrk_(const int order, const int uplo , const int trans ,
             const int n , const int k , const double alpha , const double *a , const int lda ,
             const double beta , double *c , const int ldc){
  internal::engine(53,internal::profile_blas3,std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2,3>{},
                          &cblas_dsyrk,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,beta,c,ldc);
}
void critter_dsyr2k_(const int order, const int uplo, const int trans, const int n, const int k, const double alpha,
              const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
              const int ldc){
  internal::engine(54,internal::profile_blas3,std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(n,k,uplo*1000+trans,(alpha!=0)*2+(beta!=0)),internal::IndexPack<0,1,2,3>{},
                          &cblas_dsyr2k,(CBLAS_ORDER)order,(CBLAS_UPLO)uplo,(CBLAS_TRANSPOSE)trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void critter_dsymm_(const int order, const int side, const int uplo, const int m, const int n, const double alpha,
             const double *a, const int lda, const double *b, const int ldb, const double beta, double *c,
             const int ldc){
  internal::engine(55,internal::profile_blas3,
                          (CBLAS_SIDE)side==CblasLeft ? std::make_tuple(m,n,side*1000+uplo,(alpha!=0)*2+(beta!=0)) : std::make_tuple(n,m,side*1000+uplo,(alpha!=0)*2+(beta!=0)),
                          std::make_tuple(m,n,side*1000+uplo,(alpha!=0)*2+(beta!=0)), internal::IndexPack<0,1,2,3>{},
                          &cblas_dsymm,(CBLAS_ORDER)order,(CBLAS_SIDE)side,(CBLAS_UPLO)uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc);
}

// **********************************************************************************************************************************
// BLAS 1
void critter__daxpy__(const int* n , const double* a , const double *x , const int* incx , double *y , const int* incy){
  critter_daxpy_(*n,*a,x,*incx,y,*incy);
}
void critter__dscal__(const int* n , const double* a , double *x , const int* incx){
  critter_dscal_(*n,*a,x,*incx);
}

// BLAS 2
void critter__dgbmv__(const char* trans , const int* m , const int* n, const int* kl, const int* ku, const double* alpha ,
               const double *a , const int* lda , const double *x, const int* incx ,
               const double* beta, double *y , const int* incy ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dgbmv_(CblasColMajor,_trans,*m,*n,*kl,*ku,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void critter__dgemv__(const char* trans , const int* m , const int* n, const double* alpha , const double *a , const int* lda ,
               const double *x, const int* incx , const double* beta, double *y , const int* incy ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dgemv_(CblasColMajor,_trans,*m,*n,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void critter__dger__(const int* m , const int* n , const double* alpha , const double *x , const int* incx , const double *y ,
              const int* incy , double *a ,
            const int* lda){
  critter_dger_(CblasColMajor,*m,*n,*alpha,x,*incx,y,*incy,a,*lda);
}
void critter__dsbmv__(const char* uplo, const int* n, const int* k, const double* alpha, const double *a,
               const int* lda, const double *x, const int* incx, const double* beta, double *y, const int* incy){
  critter_dsbmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*k,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void critter__dspmv__(const char* uplo, const int* n, const double* alpha, const double *ap, const double *x,
               const int* incx, const double* beta, double *y, const int* incy){
  critter_dspmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,ap,x,*incx,*beta,y,*incy);
}
void critter__dspr__(const char* uplo, const int* n, const double* alpha, const double *x,
              const int* incx, double *ap){
  critter_dspr_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,ap);
}
void critter__dspr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *ap){
  critter_dspr2_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,y,*incy,ap);
}
void critter__dsymv__(const char* uplo, const int* n, const double* alpha, const double *a, const int* lda,
               const double *x, const int* incx, const double* beta, double *y, const int* incy){
  critter_dsymv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,a,*lda,x,*incx,*beta,y,*incy);
}
void critter__dsyr__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
              double *a, const int* lda){
  critter_dsyr_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,a,*lda);
}
void critter__dsyr2__(const char* uplo, const int* n, const double* alpha, const double *x, const int* incx,
               const double *y, const int* incy, double *a, const int* lda){
  critter_dsyr2_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower),*n,*alpha,x,*incx,y,*incy,a,*lda);
}
void critter__dtrsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *a,
               const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtrsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,a,*lda,x,*incx);
}
void critter__dtrmv__(const char* uplo , const char* trans , const char* diag , const int* n , const double *a ,
               const int* lda , double *x, const int* incx ){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtrmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,a,*lda,x,*incx);
}
void critter__dtpsv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtpsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,ap,x,*incx);
}
void critter__dtpmv__(const char* uplo, const char* trans, const char* diag, const int* n, const double *ap,
               double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtpmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,ap,x,*incx);
}
void critter__dtbsv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtbsv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,*k,a,*lda,x,*incx);
}
void critter__dtbmv__(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
               const double *a, const int* lda, double *x, const int* incx){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dtbmv_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, (*diag=='U' ? CblasUnit : CblasNonUnit), *n,*k,a,*lda,x,*incx);
}

// BLAS 3
void critter__dgemm__(const char* transa , const char* transb ,
             const int* m , const int* n , const int* k , const double* alpha , const double *a ,
             const int* lda , const double *b , const int* ldb , const double* beta , double *c , const int* ldc){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  CBLAS_TRANSPOSE _transb;
  if (*transb=='T') _transb = CblasTrans; else if (*transb=='N') _transb = CblasNoTrans; else _transb = CblasConjTrans;
  critter_dgemm_(CblasColMajor,_transa, _transb, *m,*n,*k,*alpha,a,*lda,b,*ldb,*beta,c,*ldc);
}
void critter__dtrmm__(const char* side , const char* uplo , const char* transa ,
             const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
             const int* lda , double *b , const int* ldb){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  critter_dtrmm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower), _transa,
          (*diag=='U' ? CblasUnit : CblasNonUnit), *m,*n,*alpha,a,*lda,b,*ldb);
}
void critter__dtrsm__(const char* side , const char* uplo , const char* transa ,
             const char* diag , const int* m , const int* n , const double* alpha , const double *a ,
             const int* lda , double *b , const int* ldb){
  CBLAS_TRANSPOSE _transa;
  if (*transa=='T') _transa = CblasTrans; else if (*transa=='N') _transa = CblasNoTrans; else _transa = CblasConjTrans;
  critter_dtrsm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower), _transa,
          (*diag=='U' ? CblasUnit : CblasNonUnit), *m,*n,*alpha,a,*lda,b,*ldb);
}
void critter__dsyrk__(const char* uplo, const char* trans ,
             const int* n , const int* k , const double* alpha , const double *a , const int* lda ,
             const double* beta , double *c , const int* ldc){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dsyrk_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans, *n,*k,*alpha,a,*lda,*beta,c,*ldc);
}
void critter__dsyr2k__(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
                const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
                const int* ldc){
  CBLAS_TRANSPOSE _trans;
  if (*trans=='T') _trans = CblasTrans; else if (*trans=='N') _trans = CblasNoTrans; else _trans = CblasConjTrans;
  critter_dsyr2k_(CblasColMajor,(*uplo=='U' ? CblasUpper : CblasLower), _trans,*n,*k,*alpha,a,*lda,b,*ldb,*beta,c,*ldc);
}
void critter__dsymm__(const char* side, const char* uplo, const int* m, const int* n, const double* alpha,
               const double *a, const int* lda, const double *b, const int* ldb, const double* beta, double *c,
               const int* ldc){
  critter_dsymm_(CblasColMajor,(*side=='L' ? CblasLeft : CblasRight), (*uplo=='U' ? CblasUpper : CblasLower),*m,*n,*alpha,
          a,*lda,b,*ldb,*beta,c,*ldc);
}

// **********************************************************************************************************************************
// C interface
int critter_dgetrf_(int matrix_layout, int m , int n , double* a , int lda , int* ipiv){
  return internal::engine(100,internal::profile_lapack,
                                 (m>=n ? std::make_tuple(m,n) : std::make_pair(n,m)), std::make_tuple(m,n),internal::IndexPack<0,1>{},
                                 &LAPACKE_dgetrf,matrix_layout,m,n,a,lda,ipiv);
}
int critter_dpotrf_(int matrix_layout, char uplo , int n , double* a , int lda){
  return internal::engine(101,internal::profile_lapack,
                                   std::make_tuple(n,uplo), std::make_tuple(n,uplo),internal::IndexPack<0,1>{},
                                   &LAPACKE_dpotrf,matrix_layout,uplo,n,a,lda);
}
int critter_dtrtri_(int matrix_layout, char uplo , char diag , int n , double* a , int lda){
  return internal::engine(102,internal::profile_lapack,
                                   std::make_tuple(n,uplo,diag), std::make_tuple(n,uplo,diag),internal::IndexPack<0,1,2>{},
                                   &LAPACKE_dtrtri,matrix_layout,uplo,diag,n,a,lda);
}

int critter_dgeqrf_(int matrix_layout, int m , int n , double* a , int lda , double* tau){
  return internal::engine(103,internal::profile_lapack,
                                   (m>=n ? std::make_tuple(m,n) : std::make_tuple(n,m)),std::make_tuple(m,n),internal::IndexPack<0,1>{},
                                   &LAPACKE_dgeqrf,matrix_layout,m,n,a,lda,tau);
}
int critter_dorgqr_(int matrix_layout, int m , int n , int k , double* a , int lda , const double* tau){
  return internal::engine(104,internal::profile_lapack,
                                   std::make_tuple(m,n,k),std::make_tuple(m,n,k),internal::IndexPack<0,1,2>{},
                                   &LAPACKE_dorgqr,matrix_layout,m,n,k,a,lda,tau);
}
int critter_dormqr_(int matrix_layout, char side , char trans , int m , int n , int k , const double * a , int lda , const double * tau , double * c , int ldc){
  return internal::engine(105,internal::profile_lapack,
                                   (side == 'L' ? std::make_tuple(m,n,k,side,trans) : std::make_tuple(n,m,k,side,trans)),
                                   std::make_tuple(m,n,k,side,trans),internal::IndexPack<0,1,2,3,4>{},
                                   &LAPACKE_dormqr,matrix_layout,side,trans,m,n,k,a,lda,tau,c,ldc);
}
int critter_dgetri_(int matrix_layout, int n , double * a , int lda , const int * ipiv){
  return internal::engine(106,internal::profile_lapack,
                                   std::make_tuple(n),std::make_tuple(n),internal::IndexPack<0>{},
                                   &LAPACKE_dgetri,matrix_layout,n,a,lda,ipiv);
}
int critter_dtpqrt_(int matrix_layout, int m , int n , int l , int nb , double * a , int lda , double * b , int ldb , double * t , int ldt){
  return internal::engine(107,internal::profile_lapack,
                                   std::make_tuple(m,n,l,nb),std::make_tuple(m,n,l,nb),internal::IndexPack<0,1,2,3>{},
                                   &LAPACKE_dtpqrt,matrix_layout,m,n,l,nb,a,lda,b,ldb,t,ldt);
}
int critter_dtpmqrt_(int matrix_layout, char side , char trans , int m , int n , int k , int l , int nb , const double * v ,
               int ldv , const double * t , int ldt , double * a , int lda , double * b , int ldb){
  return internal::engine(108,internal::profile_lapack,
                                   std::make_tuple(m,n,k,l,nb),std::make_tuple(m,n,k,l,nb),internal::IndexPack<0,1,2,3,4>{},
                                   &LAPACKE_dtpmqrt,matrix_layout,side,trans,m,n,k,l,nb,v,ldv,t,ldt,a,lda,b,ldb);
}

// FORTRAN interface
void critter__dgetrf__(const int* m , const int* n , double* a , const int* lda , int* ipiv, int* info){
  *info = critter_dgetrf_(LAPACK_COL_MAJOR,*m,*n,a,*lda,ipiv);
}
void critter__dpotrf__(const char* uplo , const int* n , double* a , const int* lda, int* info){
  *info = critter_dpotrf_(LAPACK_COL_MAJOR,*uplo,*n,a,*lda);
}
void critter__dtrtri__(const char* uplo , const char* diag , const int* n , double* a , const int* lda, int* info){
  *info = critter_dtrtri_(LAPACK_COL_MAJOR,*uplo,*diag,*n,a,*lda);
}
void critter__dgeqrf__(const int* m , const int* n , double* a , const int* lda , double* tau, double* work, const int* lwork, int* info){
  *info = critter_dgeqrf_(LAPACK_COL_MAJOR,*m,*n,a,*lda,tau);
}
void critter__dorgqr__(const int* m , const int* n , const int* k , double* a , const int* lda , const double* tau, double* work, const int* lwork, int* info){
  *info = critter_dorgqr_(LAPACK_COL_MAJOR,*m,*n,*k,a,*lda,tau);
}
void critter__dormqr__(const char* side , const char* trans , const int* m , const int* n , const int* k , const double * a , const int* lda , const double * tau ,
                double * c , const int* ldc, double* work, const int* lwork, int* info){
  *info = critter_dormqr_(LAPACK_COL_MAJOR,*side,*trans,*m,*n,*k,a,*lda,tau,c,*ldc);
}
void critter__dgetri__(const int* n , double * a , const int* lda , const int * ipiv, double* work, const int* lwork, int* info){
  *info = critter_dgetri_(LAPACK_COL_MAJOR,*n,a,*lda,ipiv);
}
void critter__dtpqrt__(const int* m , const int* n , const int* l , const int* nb , double * a , const int* lda , double* b , const int* ldb , double * t , const int* ldt,
                double* work, int* info){
  *info = critter_dtpqrt_(LAPACK_COL_MAJOR,*m,*n,*l,*nb,a,*lda,b,*ldb,t,*ldt);
}
void critter__dtpmqrt__(const char* side , const char* trans , const int* m , const int* n , const int* k , const int* l , const int* nb , const double * v ,
                 const int* ldv , const double * t , const int* ldt , double * a , const int* lda , double * b , const int* ldb, double* work, int* info){
  *info = critter_dtpmqrt_(LAPACK_COL_MAJOR,*side,*trans,*m,*n,*k,*l,*nb,v,*ldv,t,*ldt,a,*lda,b,*ldb);
}
