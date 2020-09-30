#include "comm.h"
#include "../util/util.h"
#include "../discretization/util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{

void start(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 1;
  }
  internal::stack_id++;
  //if (internal::stack_id>1) { return; }
  assert(internal::internal_comm_info.size() == 0);
  internal::wait_id=true;
  internal::reset(schedule_kernels_override,force_steady_statistical_data_overide);

  // Barrier used to make as certain as possible that 'computation_timer' starts in synch.
  PMPI_Barrier(MPI_COMM_WORLD);
  internal::computation_timer = MPI_Wtime();
  internal::wall_timer.push_back((double)internal::computation_timer);
}

void stop(){
  volatile double last_time = MPI_Wtime();
  assert(internal::wall_timer.size()>0);
  internal::wall_timer[internal::wall_timer.size()-1] = last_time - internal::wall_timer[internal::wall_timer.size()-1];
  internal::stack_id--; 
  //if (internal::stack_id>0) { return; }
  PMPI_Barrier(MPI_COMM_WORLD);

  assert(internal::internal_comm_info.size() == 0);
  internal::final_accumulate(MPI_COMM_WORLD,last_time); 
  internal::propagate(MPI_COMM_WORLD);
  internal::collect(MPI_COMM_WORLD);
  internal::mode = 0; internal::wait_id=false;
}

void record(int variantID, int print_mode, double overhead_time){
  internal::record(std::cout,variantID,print_mode,overhead_time);
  if (internal::flag) {internal::record(internal::stream,variantID,print_mode,overhead_time);}
  internal::is_first_iter = false;
  assert(internal::wall_timer.size()>0);
  internal::wall_timer.pop_back();
}

void clear(){
  internal::mode = 0; internal::wait_id=false;
  internal::clear();
}

void set_mode(int input_mode){
  if (input_mode != -1) { internal::mode = input_mode; }
  else{
    if (std::getenv("CRITTER_MODE") != NULL){
      internal::mode = atoi(std::getenv("CRITTER_MODE"));
    } else{
      internal::mode = 1;
    }
  }
}

namespace internal{

// These routines aim to achieve agnosticity to mechanism.

void _init(int* argc, char*** argv){
  mode=0;
  stack_id=0;
  internal_tag = 31133;
  internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2;
  internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4;
  internal_tag5 = internal_tag+5;
  delete_comm = 1;
  flag = 0;
  file_name="";
  stream_name="";
  if (std::getenv("CRITTER_MECHANISM") != NULL){
    mechanism = atoi(std::getenv("CRITTER_MECHANISM"));
  } else{
    mechanism = 0;
  }
  if (std::getenv("CRITTER_AUTO") != NULL){
    auto_capture = atoi(std::getenv("CRITTER_AUTO"));
  } else{
    auto_capture = 0;
  }
  if (std::getenv("CRITTER_TRACK_BLAS") != NULL){
    track_blas = atoi(std::getenv("CRITTER_TRACK_BLAS"));
  } else{
    track_blas = 1;
  }
  if (std::getenv("CRITTER_TRACK_LAPACK") != NULL){
    track_lapack = atoi(std::getenv("CRITTER_TRACK_LAPACK"));
  } else{
    track_lapack = 1;
  }
  if (std::getenv("CRITTER_TRACK_COLLECTIVE") != NULL){
    track_collective = atoi(std::getenv("CRITTER_TRACK_COLLECTIVE"));
  } else{
    track_collective = 1;
  }
  if (std::getenv("CRITTER_TRACK_P2P") != NULL){
    track_p2p = atoi(std::getenv("CRITTER_TRACK_P2P"));
  } else{
    track_p2p = 1;
  }
  if (std::getenv("CRITTER_TRACK_P2P_IDLE") != NULL){
    track_p2p_idle = atoi(std::getenv("CRITTER_TRACK_P2P_IDLE"));
  } else{
    track_p2p_idle = 1;
  }
  if (std::getenv("CRITTER_EAGER_P2P") != NULL){
    eager_p2p = atoi(std::getenv("CRITTER_EAGER_P2P"));
  } else{
    eager_p2p = 0;
  }
  if (std::getenv("CRITTER_DELETE_COMM") != NULL){
    delete_comm = atoi(std::getenv("CRITTER_DELETE_COMM"));
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    file_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name = file_name + ".txt";
    flag=1;
  }
  is_first_iter = true;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (world_rank == 0){ is_world_root = true; }
  else                 { is_world_root = false; }
  if (flag == 1){
    if (is_world_root){
      stream.open(stream_name.c_str(),std::ofstream::app);
    }
  }

  _MPI_Barrier__id = 0;
  _MPI_Bcast__id = 1;
  _MPI_Reduce__id = 2;
  _MPI_Allreduce__id = 3;
  _MPI_Gather__id = 4;
  _MPI_Allgather__id = 5;
  _MPI_Scatter__id = 6;
  _MPI_Reduce_scatter__id = 7;
  _MPI_Alltoall__id = 8;
  _MPI_Gatherv__id = 9;
  _MPI_Allgatherv__id = 10;
  _MPI_Scatterv__id = 11;
  _MPI_Alltoallv__id = 12;
  _MPI_Sendrecv__id = 13;
  _MPI_Sendrecv_replace__id = 14;
  _MPI_Ssend__id = 15;
  _MPI_Send__id = 16;
  _MPI_Recv__id = 17;
  _MPI_Isend__id = 18;
  _MPI_Irecv__id = 19;
  _MPI_Ibcast__id = 20;
  _MPI_Iallreduce__id = 21;
  _MPI_Ireduce__id = 22;
  _MPI_Igather__id = 23;
  _MPI_Igatherv__id = 24;
  _MPI_Iallgather__id = 25;
  _MPI_Iallgatherv__id = 26;
  _MPI_Iscatter__id = 27;
  _MPI_Iscatterv__id = 28;
  _MPI_Ireduce_scatter__id = 29;
  _MPI_Ialltoall__id = 30;
  _MPI_Ialltoallv__id = 31;
  _MPI_Bsend__id = 32;

  _BLAS_axpy__id = 0;
  _BLAS_scal__id = 1;
  _BLAS_ger__id = 2;
  _BLAS_gemm__id = 3;
  _BLAS_trmm__id = 4;
  _BLAS_trsm__id = 5;
  _BLAS_syrk__id = 6;

  _LAPACK_getrf__id = 100;
  _LAPACK_potrf__id = 101;
  _LAPACK_trtri__id = 102;
  _LAPACK_geqrf__id = 103;
  _LAPACK_orgqr__id = 104;
  _LAPACK_ormqr__id = 105;
  _LAPACK_getri__id = 106;
  _LAPACK_tpqrt__id = 107;
  _LAPACK_tpmqrt__id = 108;

  allocate(MPI_COMM_WORLD);
  if (auto_capture) start();
}


void init(int* argc, char*** argv){
  PMPI_Init(argc,argv);
  _init(argc, argv);
}

void init_thread(int* argc, char*** argv, int required, int* provided){
  assert(required == MPI_THREAD_SINGLE);
  PMPI_Init_thread(argc,argv,required,provided);
  _init(argc, argv);
}

void barrier(MPI_Comm comm){
  if (mode){
    volatile double curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Barrier__id,curtime, 0, MPI_CHAR, comm);
    if (schedule_decision) PMPI_Barrier(comm);// For now, barriers will not be conditionally scheduled
    complete_comm(_MPI_Barrier__id);
  }
  else{
    PMPI_Barrier(comm);
  }
}

void comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm){
  if (mode){
    volatile double curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Barrier__id,curtime, 0, MPI_CHAR, comm);
    PMPI_Comm_split(comm,color,key,newcomm);
    complete_comm(_MPI_Barrier__id);
  }
  else{
    PMPI_Comm_split(comm,color,key,newcomm);
  }
  exchange_communicators(comm,*newcomm);
}

void comm_free(MPI_Comm* comm){
  if (mode){
    if (delete_comm){
      PMPI_Comm_free(comm);
    }
  }
  else{
    PMPI_Comm_free(comm);
  }
}

void bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Bcast__id,curtime, count, datatype, comm);
    if (schedule_decision) PMPI_Bcast(buffer, count, datatype, root, comm);
    complete_comm(_MPI_Bcast__id);
  }
  else{
    PMPI_Bcast(buffer, count, datatype, root, comm);
  }
}

void reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Reduce__id,curtime, count, datatype, comm);
    if (schedule_decision) PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    complete_comm(_MPI_Reduce__id);
  }
  else{
    PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
}

void allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Allreduce__id,curtime, count, datatype, comm);
    if (schedule_decision) PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    complete_comm(_MPI_Allreduce__id);
  }
  else{
    PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }
}

void gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Gather__id,curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision) PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Gather__id);
  }
  else{
    PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Allgather__id,curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision) PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    complete_comm(_MPI_Allgather__id);
  }
  else{
    PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
}

void scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Scatter__id,curtime, sendbuf_size, sendtype, comm);
    if (schedule_decision) PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Scatter__id);
  }
  else{
    PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Reduce_scatter__id,curtime, tot_recv, datatype, comm);
    if (schedule_decision) PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    complete_comm(_MPI_Reduce_scatter__id);
  }
  else{
    PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }
}

void alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Alltoall__id,curtime,recvbuf_size, sendtype, comm);
    if (schedule_decision) PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    complete_comm(_MPI_Alltoall__id);
  }
  else{
    PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
}

void gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    bool schedule_decision = initiate_comm(_MPI_Gatherv__id,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
    complete_comm(_MPI_Gatherv__id);
   }
   else{
    PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  }
}

void allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Allgatherv__id,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
    complete_comm(_MPI_Allgatherv__id);
  }
  else{
    PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
}

void scatterv(const void* sendbuf, const int* sendcounts, const int* displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0; int comm_size;MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    bool schedule_decision = initiate_comm(_MPI_Scatterv__id,curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    if (schedule_decision) PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Scatterv__id);
  }
  else{
    PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int* recvcounts, const int* rdispls, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Alltoallv__id,curtime, std::max(tot_send,tot_recv), sendtype, comm);
    if (schedule_decision) PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
    complete_comm(_MPI_Alltoallv__id);
  }
  else{
    PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  }
}

void sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Sendrecv__id,curtime, std::max(sendcount,recvcount), sendtype, comm, true, dest, source);
    if (schedule_decision) PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    complete_comm(_MPI_Sendrecv__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
  }
  else{
    PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
}

void sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Sendrecv_replace__id,curtime, count, datatype, comm, true, dest, source);
    if (schedule_decision) PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
    complete_comm(_MPI_Sendrecv_replace__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
   }
  else{
    PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
}

void ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Ssend__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) PMPI_Ssend(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Ssend__id);
  }
  else{
    PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
}

void bsend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Bsend__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) PMPI_Bsend(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Bsend__id);
  }
  else{
    PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
}

void send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Send__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) PMPI_Send(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Send__id);
  }
  else{
    PMPI_Send(buf, count, datatype, dest, tag, comm);
  }
}

void recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Recv__id,curtime, count, datatype, comm, false, source);
    if (schedule_decision) PMPI_Recv(buf, count, datatype, source, tag, comm, status);
    complete_comm(_MPI_Recv__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
  }
  else{
    PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  }
}

void isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    volatile double itime = MPI_Wtime();
    PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Isend__id,curtime, itime, count, datatype, comm, request, true, dest);
  }
  else{
    PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }
}

void irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    volatile double itime = MPI_Wtime();
    PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Irecv__id,curtime, itime, count, datatype, comm, request, false, source);
  }
  else{
    PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }
}

void ibcast(void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Ibcast(buf, count, datatype, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Ibcast__id,curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Ibcast(buf, count, datatype, root, comm, request);
  }
}

void iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Iallreduce__id, curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
  }
}

void ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Iallreduce__id,curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
  }
}

void igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Igather__id, curtime, itime, recvbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_rank,comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Igatherv__id, curtime, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
  }
  else{
     PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
  }
}

void iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size); int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
     bool schedule_decision = initiate_comm(_MPI_Iallgather__id, curtime, itime, recvbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
}

void iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Iallgatherv__id, curtime, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
  }
  else{
    PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
  }
}

void iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Iscatter__id,curtime, itime, sendbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0;
    int comm_rank, comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    volatile double itime = MPI_Wtime();
    PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Iscatterv__id, curtime, itime, std::max(tot_send,(int64_t)recvcount), sendtype, comm, request);
  }
  else{
    PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Ireduce_scatter__id, curtime, itime, tot_recv, datatype, comm, request);
  }
  else{
    PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
  }
}

void ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    volatile double itime = MPI_Wtime();
    PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Ialltoall__id, curtime, itime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm, request);
  }
  else{
    PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
}

void ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  if (mode && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    bool schedule_decision = initiate_comm(_MPI_Ialltoallv__id, curtime, itime, std::max(tot_send,tot_recv), sendtype, comm, request);
  }
  else{
    PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
  }
}

void wait(MPI_Request* request, MPI_Status* status){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete_comm(curtime,request, status);
  }
  else{
    PMPI_Wait(request, status);
  }
}

void waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete_comm(curtime, count, array_of_requests, indx, status);
  }
  else{
    PMPI_Waitany(count, array_of_requests, indx, status);
  }
}

void waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete_comm(curtime, incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  else{
    PMPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
}

void waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  if (mode && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete_comm(curtime,count,array_of_requests,array_of_statuses);
  }
  else{
    PMPI_Waitall(count, array_of_requests, array_of_statuses);
  }
}

void finalize(){
  if (auto_capture) stop();
  if (is_world_root){
    if (flag == 1){
      stream.close();
    }
  }
  internal::_finalize();
  PMPI_Finalize();
}

}
}
