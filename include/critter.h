#ifndef __CRITTER_H__
#define __CRITTER_H__

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

#endif /*CRITTER_H_*/
