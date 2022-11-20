#ifndef __CRITTER_LAPACK_H__
#define __CRITTER_LAPACK_H__

#include "../src/intercept/comp.h"

// C interface
#define LAPACKE_dgetrf(a,b,c,d,e,f)\
    critter_dgetrf_(a,b,c,d,e,f)

#define LAPACKE_dpotrf(a,b,c,d,e)\
    critter_dpotrf_(a,b,c,d,e)

#define LAPACKE_dtrtri(a,b,c,d,e,f)\
    critter_dtrtri_(a,b,c,d,e,f)

#define LAPACKE_dgeqrf(a,b,c,d,e,f)\
    critter_dgeqrf_(a,b,c,d,e,f)

#define LAPACKE_dorgqr(a,b,c,d,e,f,g)\
    critter_dorgqr_(a,b,c,d,e,f,g)

#define LAPACKE_dormqr(a,b,c,d,e,f,g,h,i,j,k)\
    critter_dormqr_(a,b,c,d,e,f,g,h,i,j,k)

#define LAPACKE_dgetri(a,b,c,d,e)\
    critter_dgetri_(a,b,c,d,e)

#define LAPACKE_dtpmqrt(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\
    critter_dtpmqrt_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)

#define LAPACKE_dtpqrt(a,b,c,d,e,f,g,h,i,j,k)\
    critter_dtpqrt_(a,b,c,d,e,f,g,h,i,j,k)


// FORTRAN interface
#define dgetrf_(a,b,c,d,e,f)\
    critter__dgetrf__(a,b,c,d,e,f)
#define dgetrf(a,b,c,d,e,f)\
    critter__dgetrf__(a,b,c,d,e,f)
#define DGETRF_(a,b,c,d,e,f)\
    critter__dgetrf__(a,b,c,d,e,f)
#define DGETRF(a,b,c,d,e,f)\
    critter__dgetrf__(a,b,c,d,e,f)

#define dpotrf_(a,b,c,d,e)\
    critter__dpotrf__(a,b,c,d,e)
#define dpotrf(a,b,c,d,e)\
    critter__dpotrf__(a,b,c,d,e)
#define DPOTRF_(a,b,c,d,e)\
    critter__dpotrf__(a,b,c,d,e)
#define DPOTRF(a,b,c,d,e)\
    critter__dpotrf__(a,b,c,d,e)

#define dtrtri_(a,b,c,d,e,f)\
    critter__dtrtri__(a,b,c,d,e,f)
#define dtrtri(a,b,c,d,e,f)\
    critter__dtrtri__(a,b,c,d,e,f)
#define DTRTRI_(a,b,c,d,e,f)\
    critter__dtrtri__(a,b,c,d,e,f)
#define DTRTRI(a,b,c,d,e,f)\
    critter__dtrtri__(a,b,c,d,e,f)

#define dgeqrf_(a,b,c,d,e,f,g,h)\
    critter__dgeqrf__(a,b,c,d,e,f,g,h)
#define dgeqrf(a,b,c,d,e,f,g,h)\
    critter__dgeqrf__(a,b,c,d,e,f,g,h)
#define DGEQRF_(a,b,c,d,e,f,g,h)\
    critter__dgeqrf__(a,b,c,d,e,f,g,h)
#define DGEQRF(a,b,c,d,e,f,g,h)\
    critter__dgeqrf__(a,b,c,d,e,f,g,h)

#define dorgqr_(a,b,c,d,e,f,g,h,i)\
    critter__dorgqr__(a,b,c,d,e,f,g,h,i)
#define dorgqr(a,b,c,d,e,f,g,h,i)\
    critter__dorgqr__(a,b,c,d,e,f,g,h,i)
#define DORGQR_(a,b,c,d,e,f,g,h,i)\
    critter__dorgqr__(a,b,c,d,e,f,g,h,i)
#define DORGQR(a,b,c,d,e,f,g,h,i)\
    critter__dorgqr__(a,b,c,d,e,f,g,h,i)

#define dormqr_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dormqr__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define dormqr(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dormqr__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DORMQR_(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dormqr__(a,b,c,d,e,f,g,h,i,j,k,l,m)
#define DORMQR(a,b,c,d,e,f,g,h,i,j,k,l,m)\
    critter__dormqr__(a,b,c,d,e,f,g,h,i,j,k,l,m)

#define dgetri_(a,b,c,d,e,f,g)\
    critter__dgetri__(a,b,c,d,e,f,g)
#define dgetri(a,b,c,d,e,f,g)\
    critter__dgetri__(a,b,c,d,e,f,g)
#define DGETRI_(a,b,c,d,e,f,g)\
    critter__dgetri__(a,b,c,d,e,f,g)
#define DGETRI(a,b,c,d,e,f,g)\
    critter__dgetri__(a,b,c,d,e,f,g)

#define dtpmqrt_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)\
    critter__dtpmqrt__(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
#define dtpmqrt(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)\
    critter__dtpmqrt__(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
#define DTPMQRT_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)\
    critter__dtpmqrt__(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
#define DTPMQRT(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)\
    critter__dtpmqrt__(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)

#define dtpqrt_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dtpqrt__(a,b,c,d,e,f,g,h,i,j,k,l)
#define dtpqrt(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dtpqrt__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DTPQRT_(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dtpqrt__(a,b,c,d,e,f,g,h,i,j,k,l)
#define DTPQRT(a,b,c,d,e,f,g,h,i,j,k,l)\
    critter__dtpqrt__(a,b,c,d,e,f,g,h,i,j,k,l)

#endif /*CRITTER_LAPACK_H_*/
