#ifndef __CRITTER_H__
#define __CRITTER_H__

#include "../src/interface.h"
#include "../src/intercept/comm.h"

#define CRITTER_START(ARG)\
    critter_start_timer(#ARG);

#define CRITTER_STOP(ARG)\
    critter_stop_timer(#ARG);

#define MPI_Init(argc, argv)\
     critter_init(argc,argv)

#define MPI_Init_thread(argc, argv, required, provided)\
     critter_init_thread(argc,argv,required,provided)

#define MPI_Finalize()\
    critter_finalize()

#define MPI_Barrier(cm)\
    critter_barrier(cm)

#define MPI_Comm_split(comm,color,key,newcomm)\
    critter_comm_split(comm,color,key,newcomm)

#define MPI_Comm_dup(comm,newcomm)\
    critter_comm_dup(comm,newcomm)

#define MPI_Comm_free(cm)\
    critter_comm_free(cm)

#define MPI_Bcast(buf, nelem, t, root, cm)\
    critter_bcast(buf,nelem,t,root,cm)

#define MPI_Reduce(sbuf, rbuf, nelem, t, op, root, cm)\
    critter_reduce(sbuf,rbuf,nelem,t,op,root,cm)

#define MPI_Allreduce(sbuf, rbuf, nelem, t, op, cm)\
    critter_allreduce(sbuf,rbuf,nelem,t,op,cm)

#define MPI_Gather(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
    critter_gather(sbuf,scount,st,rbuf,rcount,rt,root,cm)

#define MPI_Allgather(sbuf, scount, st, rbuf, rcount, rt, cm)\
    critter_allgather(sbuf,scount,st,rbuf,rcount,rt,cm)

#define MPI_Scatter(sbuf, scount, st, rbuf, rcount, rt, root, cm)\
    critter_scatter(sbuf,scount,st,rbuf,rcount,rt,root,cm)

#define MPI_Reduce_scatter(sbuf, rbuf, rcounts, t, op, cm)\
    critter_reduce_scatter(sbuf,rbuf,rcounts,t,op,cm)

#define MPI_Alltoall(sbuf, scount, st, rbuf, rcount, rt, cm)\
    critter_alltoall(sbuf,scount,st,rbuf,rcount,rt,cm)

#define MPI_Gatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm)\
    critter_gatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm)

#define MPI_Allgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm)\
    critter_allgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm)

#define MPI_Scatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm)\
    critter_scatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm)

#define MPI_Alltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm)\
    critter_alltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm)

#define MPI_Sendrecv(sbuf, scnt, st, dest, stag, rbuf, rcnt, rt, src, rtag, cm, status)\
    critter_sendrecv(sbuf,scnt,st,dest,stag,rbuf,rcnt,rt,src,rtag,cm,status)

#define MPI_Sendrecv_replace(sbuf, scnt, st, dest, stag, src, rtag, cm, status)\
    critter_sendrecv_replace(sbuf,scnt,st,dest,stag,src,rtag,cm,status)

#define MPI_Ssend(buf, nelem, t, dest, tag, cm)\
    critter_ssend(buf,nelem,t,dest,tag,cm)

#define MPI_Bsend(buf, nelem, t, dest, tag, cm)\
    critter_bsend(buf,nelem,t,dest,tag,cm)

#define MPI_Send(buf, nelem, t, dest, tag, cm)\
    critter_send(buf,nelem,t,dest,tag,cm)

#define MPI_Recv(buf, nelem, t, src, tag, cm, status)\
    critter_recv(buf,nelem,t,src,tag,cm,status)

#define MPI_Isend(buf, nelem, t, dest, tag, cm, req)\
    critter_isend(buf,nelem,t,dest,tag,cm,req)

#define MPI_Irecv(buf, nelem, t, src, tag, cm, req)\
    critter_irecv(buf, nelem, t, src, tag, cm, req)

#define MPI_Ibcast(buf, nelem, t, root, cm, req)\
    critter_ibcast(buf,nelem,t,root,cm,req)

#define MPI_Iallreduce(sbuf, rbuf, nelem, t, op, cm, req)\
    critter_iallreduce(sbuf,rbuf,nelem,t,op,cm,req)

#define MPI_Ireduce(sbuf, rbuf, nelem, t, op, root, cm, req)\
    critter_ireduce(sbuf,rbuf,nelem,t,op,root,cm,req)

#define MPI_Igather(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
    critter_igather(sbuf,scount,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Igatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, root, cm, req)\
    critter_igatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,root,cm,req)

#define MPI_Iallgather(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
    critter_iallgather(sbuf,scount,st,rbuf,rcount,rt,cm,req)

#define MPI_Iallgatherv(sbuf, scount, st, rbuf, rcounts, rdispsls, rt, cm, req)\
    critter_iallgatherv(sbuf,scount,st,rbuf,rcounts,rdispsls,rt,cm,req)

#define MPI_Iscatter(sbuf, scount, st, rbuf, rcount, rt, root, cm, req)\
    critter_iscatter(sbuf,scount,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Iscatterv(sbuf, scounts, sdispls, st, rbuf, rcount, rt, root, cm, req)\
    critter_iscatterv(sbuf,scounts,sdispls,st,rbuf,rcount,rt,root,cm,req)

#define MPI_Ireduce_scatter(sbuf, rbuf, rcounts, t, op, cm, req)\
    critter_ireduce_scatter(sbuf,rbuf,rcounts,t,op,cm,req)

#define MPI_Ialltoall(sbuf, scount, st, rbuf, rcount, rt, cm, req)\
    critter_ialltoall(sbuf,scount,st,rbuf,rcount,rt,cm,req)

#define MPI_Ialltoallv(sbuf, scounts, sdispls, st, rbuf, rcounts, rdispsls, rt, cm, req)\
    critter_ialltoallv(sbuf,scounts,sdispls,st,rbuf,rcounts,rdispsls,rt,cm,req)

#define MPI_Wait(req, stat)\
    critter_wait(req,stat)

#define MPI_Waitany(cnt, reqs, indx, stat)\
    critter_waitany(cnt, reqs, indx, stat)

#define MPI_Waitsome(incnt, reqs, outcnt, indices, stats)\
    critter_waitsome(incnt,reqs,outcnt,indices,stats)

#define MPI_Waitall(cnt, reqs, stats)\
    critter_waitall(cnt,reqs,stats)

#define MPI_Test(req, flag, st)\
    critter_test(req,flag,st)

#define MPI_Testany(cnt, reqs, indx, flag, st)\
    critter_testany(cnt,reqs,indx,flag,st)

#define MPI_Testsome(incnt, reqs, outcnt, indices, stats)\
    critter_testsome(incnt,reqs,outcnt,indices,stats)

#define MPI_Testall(cnt, reqs, flag, stats)\
    critter_testall(cnt,reqs,flag,stats)

#endif /*CRITTER_H_*/
