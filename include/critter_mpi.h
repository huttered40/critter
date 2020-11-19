#ifndef __CRITTER_MPI_H__
#define __CRITTER_MPI_H__

#include "../src/intercept/comm.h"

#define MPI_Init(argc, argv)\
     critter::internal::init(argc,argv)

#define MPI_Init_thread(argc, argv, required, provided)\
     critter::internal::init_thread(argc,argv,required,provided)

#define MPI_Finalize()\
    critter::internal::finalize()

#define MPI_Barrier(cm)\
    critter::internal::barrier(cm)

#define MPI_Comm_split(comm,color,key,newcomm)\
    critter::internal::comm_split(comm,color,key,newcomm)

#define MPI_Comm_dup(comm,newcomm)\
    critter::internal::comm_dup(comm,newcomm)

#define MPI_Comm_free(cm)\
    critter::internal::comm_free(cm)

#define MPI_Get_count(status,dt,count)\
    critter::internal::get_count(status,dt,count)

#define MPI_Bcast(buf, nelem, t, root, cm)\
    critter::internal::bcast(buf,nelem,t,root,cm)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
    critter::internal::reduce(sbuf,rbuf,nelem,t,op,root,cm)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
    critter::internal::allreduce(sbuf,rbuf,nelem,t,op,cm)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
    critter::internal::gather(sbuf,scount,st,rbuf,rcount,rt,root,cm)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
    critter::internal::allgather(sbuf,scount,st,rbuf,rcount,rt,cm)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
    critter::internal::scatter(sbuf,scount,st,rbuf,rcount,rt,root,cm)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
    critter::internal::reduce_scatter(sbuf,rbuf,rcounts,t,op,cm)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
    critter::internal::alltoall(sbuf,scount,st,rbuf,rcount,rt,cm)

#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
    critter::internal::gatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
    critter::internal::allgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)\
    critter::internal::scatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
    critter::internal::alltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm)

#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
    critter::internal::sendrecv(sbuf,scnt,st,dest,stag,rbuf,rcnt,rt,src,rtag,cm,status)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
    critter::internal::sendrecv_replace(sbuf,scnt,st,dest,stag,src,rtag,cm,status)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
    critter::internal::ssend(buf,nelem,t,dest,tag,cm)

#define MPI_Bsend(buf, nelem, t, dest, tag, cm)\
    critter::internal::bsend(buf,nelem,t,dest,tag,cm)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
    critter::internal::send(buf,nelem,t,dest,tag,cm)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
    critter::internal::recv(buf,nelem,t,src,tag,cm,status)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
    critter::internal::isend(buf,nelem,t,dest,tag,cm,req)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
    critter::internal::irecv(buf, nelem, t, src, tag, cm, req)

#define MPI_Ibcast(buf, nelem, t, root, cm, req)\
    critter::internal::ibcast(buf,nelem,t,root,cm,req)

#define MPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req)\
    critter::internal::iallreduce(sbuf,rbuf,nelem,t,op,cm,req)

#define MPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req)\
    critter::internal::ireduce(sbuf,rbuf,nelem,t,op,root,cm,req)

#define MPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
    critter::internal::igather(sbuf,scount,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req)\
    critter::internal::igatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm,req)

#define MPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
    critter::internal::iallgather(sbuf,scount,st,rbuf,rcount,rt,cm,req)

#define MPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req)\
    critter::internal::iallgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm,req)

#define MPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
    critter::internal::iscatter(sbuf,scount,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req)\
    critter::internal::iscatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req)\
    critter::internal::ireduce_scatter(sbuf,rbuf,rcounts,t,op,cm,req)

#define MPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
    critter::internal::ialltoall(sbuf,scount,st,rbuf,rcount,rt,cm,req)

#define MPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req)\
    critter::internal::ialltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm,req)

#define MPI_Wait(req, stat)\
    critter::internal::wait(req,stat)

#define MPI_Waitany(cnt, reqs, indx, stat)\
    critter::internal::waitany(cnt, reqs, indx, stat)

#define MPI_Waitsome(incnt, reqs, outcnt, indices, stats)\
    critter::internal::waitsome(incnt,reqs,outcnt,indices,stats)

#define MPI_Waitall(cnt, reqs, stats)\
    critter::internal::waitall(cnt,reqs,stats)

#define MPI_Test(req,flag,st)\
    critter::internal::test(req,flag,st)

#define MPI_Probe(src,tag,cm,st)\
    critter::internal::probe(src,tag,cm,st)

#define MPI_Iprobe(src,tag,cm,fl,st)\
    critter::internal::iprobe(src,tag,cm,fl,st)

#endif /*CRITTER_MPI_H_*/
