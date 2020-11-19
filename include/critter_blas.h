#ifndef __CRITTER_BLAS_H__
#define __CRITTER_BLAS_H__

#include "../src/intercept/comp.h"

// *****************************************************************************************************************************************************************
// Note: these are defined specially for 'capital', which abstracts the call to blas routines. Double arguents are always used.

// C interface

// BLAS 1
#define cblas_daxpy(a,b,c,d,e,f)\
    critter::internal::_daxpy_(a,b,c,d,e,f)

#define cblas_dscal(a,b,c,d)\
    critter::internal::_dscal_(a,b,c,d)

// BLAS 2
#define cblas_dgbmv(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
    critter::internal::_dgbmv_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)

#define cblas_dgemv(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::_dgemv_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dger(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::_dger_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dtrsv(a,b,c,d,e,f,g,h,i)\
    critter::internal::_dtrsv_(a,b,c,d,e,f,g,h,i)

#define cblas_dtrmv(a,b,c,d,e,f,g,h,i)\
    critter::internal::_dtrmv_(a,b,c,d,e,f,g,h,i)

#define cblas_dtpsv(a,b,c,d,e,f,g,h)\
    critter::internal::_dtpsv_(a,b,c,d,e,f,g,h)

#define cblas_dtpmv(a,b,c,d,e,f,g,h)\
    critter::internal::_dtpmv_(a,b,c,d,e,f,g,h)

#define cblas_dtbsv(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::_dtbsv_(a,b,c,d,e,f,g,h,i,j)

#define cblas_dtbmv(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::_dtbmv_(a,b,c,d,e,f,g,h,i,j)

// BLAS 3
#define cblas_dtrmm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::_dtrmm_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dtrsm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::_dtrsm_(a,b,c,d,e,f,g,h,i,j,k,l)

#define cblas_dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
    critter::internal::_dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)

#define cblas_dsyrk(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::_dsyrk_(a,b,c,d,e,f,g,h,i,j,k)

#define cblas_dsyr2k(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::_dsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define cblas_dsymm(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::_dsymm_(a,b,c,d,e,f,g,h,i,j,k,l,m)

// FORTRAN interface
#define daxpy_(a,b,c,d,e,f)\
    critter::internal::__daxpy__(a,b,c,d,e,f)
#define daxpy(a,b,c,d,e,f)\
    critter::internal::__daxpy__(a,b,c,d,e,f)
#define DAXPY_(a,b,c,d,e,f)\
    critter::internal::__daxpy__(a,b,c,d,e,f)
#define DAXPY(a,b,c,d,e,f)\
    critter::internal::__daxpy__(a,b,c,d,e,f)

#define dscal_(a,b,c,d)\
    critter::internal::__dscal__(a,b,c,d)
#define dscal(a,b,c,d)\
    critter::internal::__dscal__(a,b,c,d)
#define DSCAL_(a,b,c,d)\
    critter::internal::__dscal__(a,b,c,d)
#define DSCAL(a,b,c,d)\
    critter::internal::__dscal__(a,b,c,d)

// BLAS 2
#define dgbmv_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define dgbmv(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGBMV_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGBMV(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgbmv__(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define dgemv_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define dgemv(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define DGEMV_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dgemv__(a,b,c,d,e,f,g,h,i,j,k)
#define DGEMV(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dgemv__(a,b,c,d,e,f,g,h,i,j,k)

#define dger_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dger__(a,b,c,d,e,f,g,h,i)
#define dger(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dger__(a,b,c,d,e,f,g,h,i)
#define DGER_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dger__(a,b,c,d,e,f,g,h,i)
#define DGER(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dger__(a,b,c,d,e,f,g,h,i)

#define dtrsv_(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrsv__(a,b,c,d,e,f,g,h)
#define dtrsv(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrsv__(a,b,c,d,e,f,g,h)
#define DTRSV_(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrsv__(a,b,c,d,e,f,g,h)
#define DTRSV(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrsv__(a,b,c,d,e,f,g,h)

#define dtrmv_(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrmv__(a,b,c,d,e,f,g,h)
#define dtrmv(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrmv__(a,b,c,d,e,f,g,h)
#define DTRMV_(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrmv__(a,b,c,d,e,f,g,h)
#define DTRMV(a,b,c,d,e,f,g,h)\
    critter::internal::__dtrmv__(a,b,c,d,e,f,g,h)

#define dtpsv_(a,b,c,d,e,f,g)\
    critter::internal::__dtpsv__(a,b,c,d,e,f,g)
#define dtpsv(a,b,c,d,e,f,g)\
    critter::internal::__dtpsv__(a,b,c,d,e,f,g)
#define DTPSV_(a,b,c,d,e,f,g)\
    critter::internal::__dtpsv__(a,b,c,d,e,f,g)
#define DTPSV(a,b,c,d,e,f,g)\
    critter::internal::__dtpsv__(a,b,c,d,e,f,g)

#define dtpmv_(a,b,c,d,e,f,g)\
    critter::internal::__dtpmv__(a,b,c,d,e,f,g)
#define dtpmv(a,b,c,d,e,f,g)\
    critter::internal::__dtpmv__(a,b,c,d,e,f,g)
#define DTPMV_(a,b,c,d,e,f,g)\
    critter::internal::__dtpmv__(a,b,c,d,e,f,g)
#define DTPMV(a,b,c,d,e,f,g)\
    critter::internal::__dtpmv__(a,b,c,d,e,f,g)

#define dtbsv_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbsv__(a,b,c,d,e,f,g,h,i)
#define dtbsv(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbsv__(a,b,c,d,e,f,g,h,i)
#define DTBSV_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbsv__(a,b,c,d,e,f,g,h,i)
#define DTBSV(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbsv__(a,b,c,d,e,f,g,h,i)

#define dtbmv_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbmv__(a,b,c,d,e,f,g,h,i)
#define dtbmv(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbmv__(a,b,c,d,e,f,g,h,i)
#define DTBMV_(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbmv__(a,b,c,d,e,f,g,h,i)
#define DTBMV(a,b,c,d,e,f,g,h,i)\
    critter::internal::__dtbmv__(a,b,c,d,e,f,g,h,i)

// BLAS 3
#define dtrmm_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define dtrmm(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRMM_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRMM(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrmm__(a,b,c,d,e,f,g,h,i,j,k)

#define dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGEMM_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DGEMM(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter::internal::__dgemm__(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define dtrsm_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define dtrsm(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRSM_(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)
#define DTRSM(a,b,c,d,e,f,g,h,i,j,k)\
    critter::internal::__dtrsm__(a,b,c,d,e,f,g,h,i,j,k)

#define dsyrk_(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define dsyrk(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define DSYRK_(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::__dsyrk__(a,b,c,d,e,f,g,h,i,j)
#define DSYRK(a,b,c,d,e,f,g,h,i,j)\
    critter::internal::__dsyrk__(a,b,c,d,e,f,g,h,i,j)

#define dsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define dsyr2k(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYR2K_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYR2K(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsyr2k__(a,b,c,d,e,f,g,h,i,j,k,l)

#define dsymm_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define dsymm(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYMM_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DSYMM(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter::internal::__dsymm__(a,b,c,d,e,f,g,h,i,j,k,l)

#endif /*CRITTER_BLAS_H_*/
