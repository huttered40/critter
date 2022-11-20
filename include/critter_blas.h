#ifndef __CRITTER_BLAS_H__
#define __CRITTER_BLAS_H__

#include "../src/intercept/comp.h"

// *****************************************************************************************************************************************************************
// Note: these are defined specially for 'capital', which abstracts the call to blas routines. Double arguents are always used.

// C interface

// BLAS 1
#define cblas_daxpy(a,b,c,d,e,f)\
    critter_daxpy_(a,b,c,d,e,f)

#define cblas_dscal(a,b,c,d)\
    critter_dscal_(a,b,c,d)

// BLAS 2
#define cblas_dgbmv(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
    critter_dgbmv_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)

#define cblas_dgemv(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter_dgemv_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dger(a,b,c,d,e,f,g,h,i,j)\
    critter_dger_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dsbmv(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter_dsbmv_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dspmv(a,b,c,d,e,f,g,h,i,j)\
    critter_dspmv_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dspr(a,b,c,d,e,f,g)\
    critter_dspr_(a,b,c,d,e,f,g)

#define cblas_dspr2(a,b,c,d,e,f,g,h,i)\
    critter_dspr2_(a,b,c,d,e,f,g,h,i)

#define cblas_dsymv(a,b,c,d,e,f,g,h,i,j,k)\
    critter_dsymv_(a,b,c,d,e,f,g,h,i,j,k)

#define cblas_dsyr(a,b,c,d,e,f,g,h)\
    critter_dsyr_(a,b,c,d,e,f,g,h)

#define cblas_dsyr2(a,b,c,d,e,f,g,h,i,j)\
    critter_dsyr2_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dtrsv(a,b,c,d,e,f,g,h,i)\
    critter_dtrsv_(a,b,c,d,e,f,g,h,i)

#define cblas_dtrmv(a,b,c,d,e,f,g,h,i)\
    critter_dtrmv_(a,b,c,d,e,f,g,h,i)

#define cblas_dtpsv(a,b,c,d,e,f,g,h)\
    critter_dtpsv_(a,b,c,d,e,f,g,h)

#define cblas_dtpmv(a,b,c,d,e,f,g,h)\
    critter_dtpmv_(a,b,c,d,e,f,g,h)

#define cblas_dtbsv(a,b,c,d,e,f,g,h,i,j)\
    critter_dtbsv_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dtbmv(a,b,c,d,e,f,g,h,i,j)\
    critter_dtbmv_(a,b,c,d,e,f,g,h,i,j)

// BLAS 3
#define cblas_dtrmm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter_dtrmm_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dtrsm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter_dtrsm_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
    critter_dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)

#define cblas_dsyrk(a,b,c,d,e,f,g,h,i,j,k)\
    critter_dsyrk_(a,b,c,d,e,f,g,h,i,j,k)

#define cblas_dsyr2k(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter_dsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define cblas_dsymm(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter_dsymm_(a,b,c,d,e,f,g,h,i,j,k,l,m)

// FORTRAN interface
#define daxpy_(a,b,c,d,e,f)\
    critter__daxpy__(a,b,c,d,e,f)
#define daxpy(a,b,c,d,e,f)\
    critter__daxpy__(a,b,c,d,e,f)
#define DAXPY_(a,b,c,d,e,f)\
    critter__daxpy__(a,b,c,d,e,f)
#define DAXPY(a,b,c,d,e,f)\
    critter__daxpy__(a,b,c,d,e,f)

#define dscal_(a,b,c,d)\
    critter__dscal__(a,b,c,d)
#define dscal(a,b,c,d)\
    critter__dscal__(a,b,c,d)
#define DSCAL_(a,b,c,d)\
    critter__dscal__(a,b,c,d)
#define DSCAL(a,b,c,d)\
    critter__dscal__(a,b,c,d)

// BLAS 2
#define dgbmv_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define dgbmv(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGBMV_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGBMV(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define dgemv_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define dgemv(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define DGEMV_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define DGEMV(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dgemv__(a,b,c,d,e,f,g,h,i,j,k)

#define dger_(a,b,c,d,e,f,g,h,i)\
    critter__dger__(a,b,c,d,e,f,g,h,i)
#define dger(a,b,c,d,e,f,g,h,i)\
    critter__dger__(a,b,c,d,e,f,g,h,i)
#define DGER_(a,b,c,d,e,f,g,h,i)\
    critter__dger__(a,b,c,d,e,f,g,h,i)
#define DGER(a,b,c,d,e,f,g,h,i)\
    critter__dger__(a,b,c,d,e,f,g,h,i)

#define dsbmv_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dsbmv__(a,b,c,d,e,f,g,h,i,j,k)
#define dsbmv(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dsbmv__(a,b,c,d,e,f,g,h,i,j,k)
#define DSBMV_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dsbmv__(a,b,c,d,e,f,g,h,i,j,k)
#define DSBMV(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dsbmv__(a,b,c,d,e,f,g,h,i,j,k)

#define dspmv_(a,b,c,d,e,f,g,h,i)\
    critter__dspmv__(a,b,c,d,e,f,g,h,i)
#define dspmv(a,b,c,d,e,f,g,h,i)\
    critter__dspmv__(a,b,c,d,e,f,g,h,i)
#define DSPMV_(a,b,c,d,e,f,g,h,i)\
    critter__dspmv__(a,b,c,d,e,f,g,h,i)
#define DSPMV(a,b,c,d,e,f,g,h,i)\
    critter__dspmv__(a,b,c,d,e,f,g,h,i)

#define dspr_(a,b,c,d,e,f)\
    critter__dspr__(a,b,c,d,e,f)
#define dspr(a,b,c,d,e,f)\
    critter__dspr__(a,b,c,d,e,f)
#define DSPR_(a,b,c,d,e,f)\
    critter__dspr__(a,b,c,d,e,f)
#define DSPR(a,b,c,d,e,f)\
    critter__dspr__(a,b,c,d,e,f)

#define dspr2_(a,b,c,d,e,f,g,h)\
    critter__dspr2__(a,b,c,d,e,f,g,h)
#define dspr2(a,b,c,d,e,f,g,h)\
    critter__dspr2__(a,b,c,d,e,f,g,h)
#define DSPR2_(a,b,c,d,e,f,g,h)\
    critter__dspr2__(a,b,c,d,e,f,g,h)
#define DSPR2(a,b,c,d,e,f,g,h)\
    critter__dspr2__(a,b,c,d,e,f,g,h)

#define dsymv_(a,b,c,d,e,f,g,h,i,j)\
    critter__dsymv__(a,b,c,d,e,f,g,h,i,j)
#define dsymv(a,b,c,d,e,f,g,h,i,j)\
    critter__dsymv__(a,b,c,d,e,f,g,h,i,j)
#define DSYMV_(a,b,c,d,e,f,g,h,i,j)\
    critter__dsymv__(a,b,c,d,e,f,g,h,i,j)
#define DSYMV(a,b,c,d,e,f,g,h,i,j)\
    critter__dsymv__(a,b,c,d,e,f,g,h,i,j)

#define dsyr_(a,b,c,d,e,f,g)\
    critter__dsyr__(a,b,c,d,e,f,g)
#define dsyr(a,b,c,d,e,f,g)\
    critter__dsyr__(a,b,c,d,e,f,g)
#define DSYR_(a,b,c,d,e,f,g)\
    critter__dsyr__(a,b,c,d,e,f,g)
#define DSYR(a,b,c,d,e,f,g)\
    critter__dsyr__(a,b,c,d,e,f,g)

#define dsyr2_(a,b,c,d,e,f,g,h,i)\
    critter__dsyr2__(a,b,c,d,e,f,g,h,i)
#define dsyr2(a,b,c,d,e,f,g,h,i)\
    critter__dsyr2__(a,b,c,d,e,f,g,h,i)
#define DSYR2_(a,b,c,d,e,f,g,h,i)\
    critter__dsyr2__(a,b,c,d,e,f,g,h,i)
#define DSYR2(a,b,c,d,e,f,g,h,i)\
    critter__dsyr2__(a,b,c,d,e,f,g,h,i)

#define dtrsv_(a,b,c,d,e,f,g,h)\
    critter__dtrsv__(a,b,c,d,e,f,g,h)
#define dtrsv(a,b,c,d,e,f,g,h)\
    critter__dtrsv__(a,b,c,d,e,f,g,h)
#define DTRSV_(a,b,c,d,e,f,g,h)\
    critter__dtrsv__(a,b,c,d,e,f,g,h)
#define DTRSV(a,b,c,d,e,f,g,h)\
    critter__dtrsv__(a,b,c,d,e,f,g,h)

#define dtrmv_(a,b,c,d,e,f,g,h)\
    critter__dtrmv__(a,b,c,d,e,f,g,h)
#define dtrmv(a,b,c,d,e,f,g,h)\
    critter__dtrmv__(a,b,c,d,e,f,g,h)
#define DTRMV_(a,b,c,d,e,f,g,h)\
    critter__dtrmv__(a,b,c,d,e,f,g,h)
#define DTRMV(a,b,c,d,e,f,g,h)\
    critter__dtrmv__(a,b,c,d,e,f,g,h)

#define dtpsv_(a,b,c,d,e,f,g)\
    critter__dtpsv__(a,b,c,d,e,f,g)
#define dtpsv(a,b,c,d,e,f,g)\
    critter__dtpsv__(a,b,c,d,e,f,g)
#define DTPSV_(a,b,c,d,e,f,g)\
    critter__dtpsv__(a,b,c,d,e,f,g)
#define DTPSV(a,b,c,d,e,f,g)\
    critter__dtpsv__(a,b,c,d,e,f,g)

#define dtpmv_(a,b,c,d,e,f,g)\
    critter__dtpmv__(a,b,c,d,e,f,g)
#define dtpmv(a,b,c,d,e,f,g)\
    critter__dtpmv__(a,b,c,d,e,f,g)
#define DTPMV_(a,b,c,d,e,f,g)\
    critter__dtpmv__(a,b,c,d,e,f,g)
#define DTPMV(a,b,c,d,e,f,g)\
    critter__dtpmv__(a,b,c,d,e,f,g)

#define dtbsv_(a,b,c,d,e,f,g,h,i)\
    critter__dtbsv__(a,b,c,d,e,f,g,h,i)
#define dtbsv(a,b,c,d,e,f,g,h,i)\
    critter__dtbsv__(a,b,c,d,e,f,g,h,i)
#define DTBSV_(a,b,c,d,e,f,g,h,i)\
    critter__dtbsv__(a,b,c,d,e,f,g,h,i)
#define DTBSV(a,b,c,d,e,f,g,h,i)\
    critter__dtbsv__(a,b,c,d,e,f,g,h,i)

#define dtbmv_(a,b,c,d,e,f,g,h,i)\
    critter__dtbmv__(a,b,c,d,e,f,g,h,i)
#define dtbmv(a,b,c,d,e,f,g,h,i)\
    critter__dtbmv__(a,b,c,d,e,f,g,h,i)
#define DTBMV_(a,b,c,d,e,f,g,h,i)\
    critter__dtbmv__(a,b,c,d,e,f,g,h,i)
#define DTBMV(a,b,c,d,e,f,g,h,i)\
    critter__dtbmv__(a,b,c,d,e,f,g,h,i)

// BLAS 3
#define dtrmm_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define dtrmm(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRMM_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRMM(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)

#define dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGEMM_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGEMM(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define dtrsm_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define dtrsm(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRSM_(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRSM(a,b,c,d,e,f,g,h,i,j,k)\
    critter__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)

#define dsyrk_(a,b,c,d,e,f,g,h,i,j)\
    critter__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define dsyrk(a,b,c,d,e,f,g,h,i,j)\
    critter__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define DSYRK_(a,b,c,d,e,f,g,h,i,j)\
    critter__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define DSYRK(a,b,c,d,e,f,g,h,i,j)\
    critter__dsyrk__(a,b,c,d,e,f,g,h,i,j)

#define dsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define dsyr2k(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYR2K_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYR2K(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)

#define dsymm_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define dsymm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYMM_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYMM(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)

#endif /*CRITTER_BLAS_H_*/
