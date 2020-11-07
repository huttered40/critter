#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "../src/intercept/comp.h"
#include "../src/intercept/comm.h"
#include "../src/intercept/symbol.h"

// *****************************************************************************************************************************************************************
#define MPI_Init(argc, argv)\
  do {\
     critter::internal::init(argc,argv);\
  } while (0)

#define MPI_Init_thread(argc, argv, required, provided)\
  do{\
     critter::internal::init_thread(argc,argv,required,provided);\
   } while (0)

#define MPI_Finalize()\
  do {\
    critter::internal::finalize();\
    } while (0)

#define MPI_Barrier(cm)\
  do {\
    critter::internal::barrier(cm);\
  } while (0)

#define MPI_Comm_split(comm,color,key,newcomm)\
  do {\
    critter::internal::comm_split(comm,color,key,newcomm);\
  } while (0)

#define MPI_Comm_free(cm)\
  do {\
    critter::internal::comm_free(cm);\
  } while (0)

#define MPI_Bcast(buf, nelem, t, root, cm)\
  do {\
    critter::internal::bcast(buf,nelem,t,root,cm);\
  } while (0)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
  do {\
    critter::internal::reduce(sbuf,rbuf,nelem,t,op,root,cm);\
  } while (0)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
  do {\
    critter::internal::allreduce(sbuf,rbuf,nelem,t,op,cm);\
  } while (0)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    critter::internal::gather(sbuf,scount,st,rbuf,rcount,rt,root,cm);\
  } while (0)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    critter::internal::allgather(sbuf,scount,st,rbuf,rcount,rt,cm);\
  } while (0)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
  do {\
    critter::internal::scatter(sbuf,scount,st,rbuf,rcount,rt,root,cm);\
  } while (0)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
  do {\
    critter::internal::reduce_scatter(sbuf,rbuf,rcounts,t,op,cm);\
  } while (0)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
  do {\
    critter::internal::alltoall(sbuf,scount,st,rbuf,rcount,rt,cm);\
  } while (0)

#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
  do {\
    critter::internal::gatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm);\
  } while (0)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    critter::internal::allgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm);\
  } while (0)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)\
  do {\
    critter::internal::scatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm);\
  } while (0)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
  do {\
    critter::internal::alltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm);\
  } while (0)

#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
  do {\
    critter::internal::sendrecv(sbuf,scnt,st,dest,stag,rbuf,rcnt,rt,src,rtag,cm,status);\
  } while (0)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
  do {\
    critter::internal::sendrecv_replace(sbuf,scnt,st,dest,stag,src,rtag,cm,status);\
  } while (0)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
  do {\
    critter::internal::ssend(buf,nelem,t,dest,tag,cm);\
  } while (0)

#define MPI_Bsend(buf, nelem, t, dest, tag, cm)\
  do {\
    critter::internal::bsend(buf,nelem,t,dest,tag,cm);\
  } while (0)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
  do {\
    critter::internal::send(buf,nelem,t,dest,tag,cm);\
  } while (0)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
  do {\
    critter::internal::recv(buf,nelem,t,src,tag,cm,status);\
  } while (0)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
  do {\
    critter::internal::isend(buf,nelem,t,dest,tag,cm,req);\
  } while (0)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
  do {\
    critter::internal::irecv(buf, nelem, t, src, tag, cm, req);\
  } while (0)

#define MPI_Ibcast(buf, nelem, t, root, cm, req)\
  do {\
    critter::internal::ibcast(buf,nelem,t,root,cm,req);\
  } while (0)

#define MPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req)\
  do {\
    critter::internal::iallreduce(sbuf,rbuf,nelem,t,op,cm,req);\
  } while (0)

#define MPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req)\
  do {\
    critter::internal::ireduce(sbuf,rbuf,nelem,t,op,root,cm,req);\
  } while (0)

#define MPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    critter::internal::igather(sbuf,scount,st,rbuf,rcount,rt,root,cm,req);\
  } while (0)

#define MPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req)\
  do {\
    critter::internal::igatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm,req);\
  } while (0)

#define MPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    critter::internal::iallgather(sbuf,scount,st,rbuf,rcount,rt,cm,req);\
  } while (0)

#define MPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    critter::internal::iallgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm,req);\
  } while (0)

#define MPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    critter::internal::iscatter(sbuf,scount,st,rbuf,rcount,rt,root,cm,req);\
  } while (0)

#define MPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req)\
  do {\
    critter::internal::iscatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm,req);\
  } while (0)

#define MPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req)\
  do {\
    critter::internal::ireduce_scatter(sbuf,rbuf,rcounts,t,op,cm,req);\
  } while (0)

#define MPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
  do {\
    critter::internal::ialltoall(sbuf,scount,st,rbuf,rcount,rt,cm,req);\
  } while (0)

#define MPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req)\
  do {\
    critter::internal::ialltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm,req);\
  } while (0)

#define MPI_Wait(req, stat)\
  do {\
    critter::internal::wait(req,stat);\
  } while (0)

#define MPI_Waitany(cnt, reqs, indx, stat)\
  do {\
    critter::internal::waitany(cnt, reqs, indx, stat);\
  } while (0)

#define MPI_Waitsome(incnt, reqs, outcnt, indices, stats)\
  do {\
    critter::internal::waitsome(incnt,reqs,outcnt,indices,stats);\
  } while (0)

#define MPI_Waitall(cnt, reqs, stats)\
  do {\
    critter::internal::waitall(cnt,reqs,stats);\
  } while (0)

// *****************************************************************************************************************************************************************

#define CRITTER_START(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define CRITTER_STOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

#define TAU_START(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define TAU_STOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

#define TAU_FSTART(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define TAU_FSTOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

// *****************************************************************************************************************************************************************
// Note: these are defined specially for 'capital', which abstracts the call to blas routines. Double arguents are always used.

/*
#define _axpy_(a,b,c,d,e,f)\
  do{\
    critter::internal::_daxpy_(a,b,c,d,e,f);\
  } while (0);

#define _scal_(a,b,c,d)\
  do{\
    critter::internal::_dscal_(a,b,c,d);\
  } while (0);

#define _ger_(a,b,c,d,e,f,g,h,i,j)\
  do{\
    critter::internal::_dger_(a,b,c,d,e,f,g,h,i,j);\
  } while (0);
*/

#define _trmm_(a,b,c,d,e,f,g,h,i,j,k,l)\
  do{\
    critter::internal::_dtrmm_(a,b,c,d,e,f,g,h,i,j,k,l);\
  } while (0);

#define _trsm_(a,b,c,d,e,f,g,h,i,j,k,l)\
  do{\
    critter::internal::_dtrsm_(a,b,c,d,e,f,g,h,i,j,k,l);\
  } while (0);

#define _gemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
  do{\
    critter::internal::_dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n);\
  } while (0);

#define _syrk_(a,b,c,d,e,f,g,h,i,j,k)\
  do{\
    critter::internal::_dsyrk_(a,b,c,d,e,f,g,h,i,j,k);\
  } while (0);

// *****************************************************************************************************************************************************************
// Note: as we are testing on Stampede2, we are using cblas interface.
/*
#define cblas_strmm(a,b,c,d,e,f,g,h,i,j,k,l)\
  do{\
    critter::internal::_strmm_(a,b,c,d,e,f,g,h,i,j,k,l);\
  } while (0);

#define cblas_dtrmm(a,b,c,d,e,f,g,h,i,j,k,l)\
  do{\
    critter::internal::_dtrmm_(a,b,c,d,e,f,g,h,i,j,k,l);\
  } while (0);

#define cblas_sgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
  do{\
    critter::internal::_sgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n);\
  } while (0);

#define cblas_dgemm(a,b,c,d,e,f,g,h,i,j,k,l,m,n)\
  do{\
    critter::internal::_dgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m,n);\
  } while (0);

#define cblas_ssyrk(a,b,c,d,e,f,g,h,i,j,k)\
  do{\
    critter::internal::_ssyrk_(a,b,c,d,e,f,g,h,i,j,k);\
  } while (0);

#define cblas_dsyrk(a,b,c,d,e,f,g,h,i,j,k)\
  do{\
    critter::internal::_dsyrk_(a,b,c,d,e,f,g,h,i,j,k);\
  } while (0);
*/
// *****************************************************************************************************************************************************************
// Note: these are defined specially for 'capital', which abstracts the call to lapack routines. Double arguents are always used.

#define _getrf_(a,b,c,d,e,f)\
  do{\
    critter::internal::_dgetrf_(a,b,c,d,e,f);\
  } while (0);

#define _potrf_(a,b,c,d,e)\
  do{\
    critter::internal::_dpotrf_(a,b,c,d,e);\
  } while (0);

#define _trtri_(a,b,c,d,e,f)\
  do{\
    critter::internal::_dtrtri_(a,b,c,d,e,f);\
  } while (0);

#define _geqrf_(a,b,c,d,e,f)\
  do{\
    critter::internal::_dgeqrf_(a,b,c,d,e,f);\
  } while (0);

#define _orgqr_(a,b,c,d,e,f,g)\
  do{\
    critter::internal::_dorgqr_(a,b,c,d,e,f,g);\
  } while (0);

#define _ormqr_(a,b,c,d,e,f,g,h,i,j,k)\
  do{\
    critter::internal::_dormqr_(a,b,c,d,e,f,g,h,i,j,k);\
  } while (0);

#define _getri_(a,b,c,d,e)\
  do{\
    critter::internal::_dgetri_(a,b,c,d,e);\
  } while (0);

#define _tpmqrt_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)\
  do{\
    critter::internal::_dtpmqrt_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);\
  } while (0);

#define _tpqrt_(a,b,c,d,e,f,g,h,i,j,k)\
  do{\
    critter::internal::_dtpqrt_(a,b,c,d,e,f,g,h,i,j,k);\
  } while (0);
// *****************************************************************************************************************************************************************
// Note: as we are testing on Stampede2, we are using LAPACKE interface.
/*
#define LAPACKE_spotrf(a,b,c,d,e)\
  do{\
    critter::internal::_spotrf_(a,b,c,d,e);\
  } while (0);

#define LAPACKE_dpotrf(a,b,c,d,e)\
  do{\
    critter::internal::_dpotrf_(a,b,c,d,e);\
  } while (0);

#define LAPACKE_strtri(a,b,c,d,e,f)\
  do{\
    critter::internal::_strtri_(a,b,c,d,e,f);\
  } while (0);

#define LAPACKE_dtrtri(a,b,c,d,e,f)\
  do{\
    critter::internal::_dtrtri_(a,b,c,d,e,f);\
  } while (0);

#define LAPACKE_sgeqrf_(a,b,c,d,e,f)\
  do{\
    critter::internal::_sgeqrf_(a,b,c,d,e,f);\
  } while (0);

#define LAPACKE_dgeqrf_(a,b,c,d,e,f)\
  do{\
    critter::internal::_dgeqrf_(a,b,c,d,e,f);\
  } while (0);

#define LAPACKE_sorgqr(a,b,c,d,e,f,g)\
  do{\
    critter::internal::_sorgqr_(a,b,c,d,e,f,g);\
  } while (0);

#define LAPACKE_dorgqr(a,b,c,d,e,f,g)\
  do{\
    critter::internal::_dorgqr_(a,b,c,d,e,f,g);\
  } while (0);
*/

#define blk_to_cyc_rect(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)\
  do{\
    critter::internal::_blk_to_cyc_rect_(blocked,cyclic,num_rows_local,num_columns_local,sliceDim);\
  } while (0);

#define cyc_to_blk_rect(blocked,cyclic,num_rows_local,num_columns_local,sliceDim)\
  do{\
    critter::internal::_cyc_to_blk_rect_(blocked,cyclic,num_rows_local,num_columns_local,sliceDim);\
  } while (0);

#endif /*CRITTER_H_*/
