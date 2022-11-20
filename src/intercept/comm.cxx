#include "comm.h"
#include "../util.h"
#include "../interface.h"
#include "../profiler.h"
#include "../container/comm_tracker.h"

int critter_init(int* argc, char*** argv){
  int ret = PMPI_Init(argc,argv);
  _init(argc, argv);
  return ret;
}

int critter_init_thread(int* argc, char*** argv, int required, int* provided){
  assert(required == MPI_THREAD_SINGLE);
  int ret = PMPI_Init_thread(argc,argv,required,provided);
  _init(argc, argv);
  return ret;
}

int critter_finalize(){
  if (internal::auto_capture) { critter_stop(); critter_record(); }
  return PMPI_Finalize();
}

int critter_comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[33],curtime, 0, MPI_CHAR, comm);
    ret = PMPI_Comm_split(comm,color,key,newcomm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[33]);
  }
  else{
    ret = PMPI_Comm_split(comm,color,key,newcomm);
  }
  return ret;
}

int critter_comm_dup(MPI_Comm comm, MPI_Comm* newcomm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[34],curtime, 0, MPI_CHAR, comm);
    ret = PMPI_Comm_dup(comm,newcomm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[34]);
  }
  else{
    ret = PMPI_Comm_dup(comm,newcomm);
  }
  return ret;
}

int critter_comm_free(MPI_Comm* comm){
  return PMPI_Comm_free(comm);
}

int critter_barrier(MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[0],curtime, 0, MPI_DOUBLE, comm);
    /*if (schedule_decision)*/ ret = PMPI_Barrier(comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[0]);
  }
  else{
    ret = PMPI_Barrier(comm);
  }
  return ret;
}

int critter_bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[1], curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Bcast(buffer, count, datatype, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[1]);
  }
  else{
    ret = PMPI_Bcast(buffer, count, datatype, root, comm);
  }
  return ret;
}

int critter_reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[2],curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[2]);
  }
  else{
    ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
  return ret;
}

int critter_allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[3],curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[3]);
  }
  else{
    ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }
  return ret;
}

int critter_gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
           MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[4],curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[4]);
  }
  else{
    ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int critter_allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[5],curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[5]);
  }
  else{
    ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  return ret;
}

int critter_scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
            MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[6],curtime, sendbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[6]);
  }
  else{
    ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int critter_reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[7],curtime, tot_recv, datatype, comm);
    if (schedule_decision) ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[7]);
  }
  else{
    ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }
  return ret;
}

int critter_alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[8],curtime,recvbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[8]);
  }
  else{
    ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  return ret;
}

int critter_gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[9],curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[9]);
   }
   else{
    ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  }
  return ret;
}

int critter_allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[10],curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[10]);
  }
  else{
    ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
  return ret;
}

int critter_scatterv(const void* sendbuf, const int* sendcounts, const int* displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0; int comm_size;MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[11],curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    if (schedule_decision) ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[11]);
  }
  else{
    ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int critter_alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int* recvcounts, const int* rdispls, MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[12],curtime, std::max(tot_send,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[12]);
  }
  else{
    ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  }
  return ret;
}

int critter_sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[13],curtime, std::max(sendcount,recvcount), sendtype, comm, true, dest, sendtag, source, recvtag);
    if (schedule_decision) ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount,
                                               recvtype, (source==MPI_ANY_SOURCE ? internal::save_wildcard_id : source), recvtag, comm, status);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[13]);
  }
  else{
    ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
  return ret;
}

int critter_sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[14],curtime, count, datatype, comm, true, dest, sendtag, source, recvtag);
    if (schedule_decision) ret = PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, (source==MPI_ANY_SOURCE ? internal::save_wildcard_id : source), recvtag, comm, status);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[14]);
   }
  else{
    ret = PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
  return ret;
}

int critter_ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    // Check for conflicting communication tag(s)
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[15],curtime, count, datatype, comm, true, dest, tag);
    if (schedule_decision) ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[15]);
  }
  else{
    ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int critter_bsend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    // Check for conflicting communication tag(s)
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[16],curtime, count, datatype, comm, true, dest, tag);
    if (schedule_decision) ret = PMPI_Bsend(buf, count, datatype, dest, tag, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[16]);
  }
  else{
    ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int critter_send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[17],curtime, count, datatype, comm, true, dest, tag);
    if (schedule_decision) ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[17]);
  }
  else{
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int critter_recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::initiate_comm(*(internal::blocking*)internal::list[18],curtime, count, datatype, comm, false, source, tag);
    if (schedule_decision) ret = PMPI_Recv(buf, count, datatype, (source==MPI_ANY_SOURCE ? internal::save_wildcard_id : source), tag, comm, status);
    internal::profiler::complete_comm(*(internal::blocking*)internal::list[18]);
  }
  else{
    ret = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  }
  return ret;
}

int critter_isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[19], curtime, count, datatype, comm, true, dest, tag);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[19], itime, count, datatype, comm, request, true, dest, tag);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }
  return ret;
}

int critter_irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[20], curtime, count, datatype, comm, false, source, tag);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[20], itime, count, datatype, comm, request, false, source, tag);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }
  return ret;
}

int critter_ibcast(void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[21], curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ibcast(buf, count, datatype, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[21], itime, count, datatype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Ibcast(buf, count, datatype, root, comm, request);
  }
  return ret;
}

int critter_iallreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[22], curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[22], itime, count, datatype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
  }
  return ret;
}

int critter_ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[23], curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[23], itime, count, datatype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
  }
  return ret;
}

int critter_igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[24], curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[24], itime, recvbuf_size, sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int critter_igatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_rank,comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[25], curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[25], itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
     ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
  }
  return ret;
}

int critter_iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size); int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[26], curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[26], itime, recvbuf_size, sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
  return ret;
}

int critter_iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[27], curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[27], itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
  }
  return ret;
}

int critter_iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[28], curtime, sendbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[28], itime, sendbuf_size, sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int critter_iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0;
    int comm_rank, comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[29], curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[29], itime, std::max(tot_send,(int64_t)recvcount), sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int critter_ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[30], curtime, tot_recv, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[30], itime, tot_recv, datatype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
  }
  return ret;
}

int critter_ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[31], curtime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[31], itime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
  return ret;
}

int critter_ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    bool schedule_decision = internal::profiler::inspect_comm(*(internal::nonblocking*)internal::list[32], curtime, std::max(tot_send,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      internal::profiler::initiate_comm(*(internal::nonblocking*)internal::list[32], itime, std::max(tot_send,tot_recv), sendtype, comm, request);
    } else{
      *request = internal::request_id++;
    }
  }
  else{
    ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
  }
  return ret;
}

int critter_wait(MPI_Request* request, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime,request, status);
  }
  else{
    ret = PMPI_Wait(request, status);
  }
  return ret;
}

int critter_waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime, count, array_of_requests, indx, status);
  }
  else{
    ret = PMPI_Waitany(count, array_of_requests, indx, status);
  }
  return ret;
}

int critter_waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime, incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  else{
    ret = PMPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  return ret;
}

int critter_waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime,count,array_of_requests,array_of_statuses);
  }
  else{
    ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
  }
  return ret;
}

int critter_test(MPI_Request* request, int* flag, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime,request, status,1,flag);
  }
  else{
    ret = PMPI_Test(request, flag, status);
  }
  return ret;
}

int critter_testany(int count, MPI_Request array_of_requests[], int* indx, int* flag, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime, count, array_of_requests, indx, status,1,flag);
  }
  else{
    ret = PMPI_Testany(count, array_of_requests, indx, flag, status);
  }
  return ret;
}

int critter_testsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime, incount, array_of_requests, outcount, array_of_indices, array_of_statuses,1);
  }
  else{
    ret = PMPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  return ret;
}

int critter_testall(int count, MPI_Request array_of_requests[], int* flag, MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (internal::mode && internal::profile_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = internal::profiler::complete_comm(curtime,count,array_of_requests,array_of_statuses,1,flag);
  }
  else{
    ret = PMPI_Testall(count, array_of_requests, flag, array_of_statuses);
  }
  return ret;
}
