#include "comm.h"
#include "../util.h"
#include "../record/record.h"
#include "../mechanism/path/dispatch.h"
#include "../mechanism/volumetric/volumetric.h"
#include "../mechanism/per-process/per-process.h"
#include "../mechanism/path/dispatch.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{

void start(){
  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 0;
  }
  assert(internal::mode>=0 && internal::mode<=3);
  internal::stack_id++;
  if (internal::stack_id>1) { return; }
  assert(internal::internal_comm_info.size() == 0);
  internal::wait_id=true; internal::waitall_id=false;
  // TODO: How to allow different number of cost models. Perhaps just put an assert that both cost models must be on? Or don't use these altogether?
  if (internal::is_world_root){
    if (!internal::is_first_iter){
      if (internal::flag) {internal::stream << "\n";} else {std::cout << "\n";}
    }
  }

  for (auto i=0; i<internal::list_size; i++){ internal::list[i]->init(); }
  for (auto i=0; i<internal::critical_path_costs.size(); i++){ internal::critical_path_costs[i]=0.; }
  for (auto i=0; i<internal::max_per_process_costs.size(); i++){ internal::max_per_process_costs[i]=0.; }
  for (auto i=0; i<internal::volume_costs.size(); i++){ internal::volume_costs[i]=0.; }

  PMPI_Barrier(MPI_COMM_WORLD);
  internal::computation_timer=MPI_Wtime();
}

void stop(){
  volatile double last_time = MPI_Wtime();
  internal::stack_id--; 
  if (internal::stack_id>0) { return; }
  assert(internal::internal_comm_info.size() == 0);
  internal::critical_path_costs[internal::num_critical_path_measures-2]+=(last_time-internal::computation_timer);	// update critical path computation time
  internal::critical_path_costs[internal::num_critical_path_measures-1]+=(last_time-internal::computation_timer);	// update critical path runtime
  internal::volume_costs[internal::num_volume_measures-2]+=(last_time-internal::computation_timer);			// update computation time volume
  internal::volume_costs[internal::num_volume_measures-1]+=(last_time-internal::computation_timer);			// update runtime volume
  for (size_t i=0; i<internal::breakdown_size; i++){ internal::critical_path_costs[internal::critical_path_costs_size-1-i] += (last_time-internal::computation_timer); }
  PMPI_Barrier(MPI_COMM_WORLD);
  internal::propagate(MPI_COMM_WORLD);
  internal::per_process::collect(MPI_COMM_WORLD);
  internal::volumetric::collect(MPI_COMM_WORLD);

  internal::record(std::cout);
  if (internal::flag) {internal::record(internal::stream);}

  internal::mode = 0; internal::wait_id=false; internal::waitall_id=false; internal::is_first_iter = false;
  internal::save_info.clear(); internal::comm_pattern_table1.clear(); internal::p2p_table.clear(); internal::comm_pattern_seq.clear();
  internal::need_new_line=false; internal::symbol_timers.clear();
}


namespace internal{

// The goal with these routines is to have them be agnostic to mechanism. This might take some tweaking on a few of them.

void _init(int* argc, char*** argv){
  PMPI_Init(argc,argv);\
  mode=0;
  stack_id=0;
  num_ftimer_measures = 2;
  mode_1_width = 25;
  mode_2_width = 15;
  internal_tag = 31133;
  internal_tag1 = internal_tag+1;
  internal_tag2 = internal_tag+2;
  internal_tag3 = internal_tag+3;
  internal_tag4 = internal_tag+4;
  internal_tag5 = internal_tag+5;
  flag = 0;
  file_name="";
  stream_name="";
  pattern_stream_name="";
  if (std::getenv("CRITTER_MECHANISM") != NULL){
    mechanism = atoi(std::getenv("CRITTER_MECHANISM"));
  } else{
    mechanism = 0;
  }
  if (std::getenv("CRITTER_MODEL_SELECT") != NULL){
    _cost_models_ = std::getenv("CRITTER_MODEL_SELECT");
  } else{
    _cost_models_ = "11";
  }
  if (std::getenv("CRITTER_PATH_SELECT") != NULL){
    _breakdown_ = std::getenv("CRITTER_PATH_SELECT");
  } else{
    _breakdown_ = "000000000";
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    flag = 1;
    file_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name = file_name + ".txt";
    pattern_stream_name = file_name + "_pattern.txt";
  }
  mode = 0;
  if (std::getenv("CRITTER_MAX_NUM_SYMBOLS") != NULL){
    max_num_symbols = atoi(std::getenv("CRITTER_MAX_NUM_SYMBOLS"));
  } else{
    max_num_symbols = 40;
  }
  if (std::getenv("CRITTER_MAX_SYMBOL_LENGTH") != NULL){
    max_timer_name_length = atoi(std::getenv("CRITTER_MAX_SYMBOL_LENGTH"));
  } else{
    max_timer_name_length = 40;
  }
  if (std::getenv("CRITTER_AUTO") != NULL){
    auto_capture = atoi(std::getenv("CRITTER_AUTO"));
  } else{
    auto_capture = 0;
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
  assert(mode>=0 && mode<=3);
  cost_model_size=0; breakdown_size=0;
  is_first_iter = true;
  need_new_line = false;
  int _critter_rank,_critter_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&_critter_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&_critter_size);
  synch_pad_send.resize(_critter_size);
  synch_pad_recv.resize(_critter_size);
  barrier_pad_send.resize(_critter_size);
  barrier_pad_recv.resize(_critter_size);
  if (_critter_rank == 0){
    is_world_root = true;
  } else {is_world_root=false;}
  if (flag == 1){
    if (_critter_rank==0){
      stream.open(stream_name.c_str());
      pattern_stream.open(pattern_stream_name.c_str());
    }
  }
  for (auto i=0; i<2; i++){
    if (_cost_models_[i] == '1'){ cost_model_size++; }
    cost_models.push_back(_cost_models_[i]);
  } 
  for (auto i=0; i<9; i++){
    if (_breakdown_[i] == '1'){ breakdown_size++; }
    breakdown.push_back(_breakdown_[i]);
  } 

  num_critical_path_measures 		= 5+2*cost_model_size;
  num_per_process_measures 		= 6+2*cost_model_size;
  num_volume_measures 		= 6+2*cost_model_size;
  num_tracker_critical_path_measures 	= 3+2*cost_model_size;
  num_tracker_per_process_measures 	= 3+2*cost_model_size;
  num_tracker_volume_measures 	= 3+2*cost_model_size;
  critical_path_costs_size            = num_critical_path_measures+num_tracker_critical_path_measures*breakdown_size*list_size+2*breakdown_size;
  per_process_costs_size              = num_per_process_measures+num_tracker_per_process_measures*breakdown_size*list_size+2*breakdown_size;
  volume_costs_size                   = num_volume_measures+num_tracker_volume_measures*list_size;

  decisions.resize(breakdown_size);
  critical_path_costs.resize(critical_path_costs_size);
  max_per_process_costs.resize(per_process_costs_size);
  volume_costs.resize(volume_costs_size);
  new_cs.resize(critical_path_costs_size);  
  symbol_pad_cp.resize(max_timer_name_length*max_num_symbols);
  symbol_pad_ncp.resize(max_timer_name_length*max_num_symbols);
  symbol_len_pad_cp.resize(max_num_symbols);
  symbol_len_pad_ncp.resize(max_num_symbols);
  symbol_timer_pad_local_cp.resize((num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols);
  symbol_timer_pad_global_cp.resize((num_ftimer_measures*num_critical_path_measures+1)*max_num_symbols);
  symbol_timer_pad_local_pp.resize((num_ftimer_measures*num_per_process_measures+1)*max_num_symbols);
  symbol_timer_pad_global_pp.resize((num_ftimer_measures*num_per_process_measures+1)*max_num_symbols);
  symbol_timer_pad_local_vol.resize((num_ftimer_measures*num_volume_measures+1)*max_num_symbols);
  symbol_timer_pad_global_vol.resize((num_ftimer_measures*num_volume_measures+1)*max_num_symbols);
  symbol_order.resize(max_num_symbols);
  info_sender.resize(num_critical_path_measures);
  info_receiver.resize(num_critical_path_measures);

  if (auto_capture) start();
}


void init(int* argc, char*** argv){
  _init(argc, argv);
}

void init_thread(int* argc, char*** argv, int required, int* provided){
  _init(argc, argv);
}

void barrier(MPI_Comm comm){
  if (mode>=1){
    volatile double curtime = MPI_Wtime();
    initiate(_MPI_Barrier,curtime, 0, MPI_CHAR, comm);
    PMPI_Barrier(comm);
    complete(_MPI_Barrier);
  }
  else{
    PMPI_Barrier(comm);
  }
}

void bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    initiate(_MPI_Bcast,curtime, count, datatype, comm);
    PMPI_Bcast(buffer, count, datatype, root, comm);
    complete(_MPI_Bcast);
  }
  else{
    PMPI_Bcast(buffer, count, datatype, root, comm);
  }
}

void reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    initiate(_MPI_Reduce,curtime, count, datatype, comm);
    PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    complete(_MPI_Reduce);
  }
  else{
    PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
}

void allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    initiate(_MPI_Allreduce,curtime, count, datatype, comm);
    PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    complete(_MPI_Allreduce);
  }
  else{
    PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }
}

void gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    initiate(_MPI_Gather,curtime, recvbuf_size, sendtype, comm);
    PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete(_MPI_Gather);
  }
  else{
    PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    initiate(_MPI_Allgather,curtime, recvbuf_size, sendtype, comm);
    PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    complete(_MPI_Allgather);
  }
  else{
    PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
}

void scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){\
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    initiate(_MPI_Scatter,curtime, sendbuf_size, sendtype, comm);
    PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete(_MPI_Scatter);
  }
  else{
    PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    initiate(_MPI_Reduce_scatter,curtime, tot_recv, datatype, comm);
    PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    complete(_MPI_Reduce_scatter);
  }
  else{
    PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }
}

void alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    initiate(_MPI_Alltoall,curtime,recvbuf_size, sendtype, comm);
    PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    complete(_MPI_Alltoall);
  }
  else{
    PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
}

void gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    initiate(_MPI_Gatherv,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
    complete(_MPI_Gatherv);
   }
   else{
    PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  }
}

void allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    initiate(_MPI_Allgatherv,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
    complete(_MPI_Allgatherv);
  }
  else{
    PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
}

void scatterv(const void* sendbuf, const int* sendcounts, const int* displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0; int comm_size;MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    initiate(_MPI_Scatterv,curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete(_MPI_Scatterv);
  }
  else{
    PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
}

void alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int* recvcounts, const int* rdispls, MPI_Datatype recvtype, MPI_Comm comm){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    initiate(_MPI_Alltoallv,curtime, std::max(tot_send,tot_recv), sendtype, comm);
    PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
    complete(_MPI_Alltoallv);
  }
  else{
    PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  }
}

void sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    initiate(_MPI_Sendrecv,curtime, std::max(sendcount,recvcount), sendtype, comm, true, dest, source);
    PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    complete(_MPI_Sendrecv);
  }
  else{
    PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
}

void sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    initiate(_MPI_Sendrecv_replace,curtime, count, datatype, comm, true, dest, source);
    PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
    complete(_MPI_Sendrecv_replace);
   }
  else{
    PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
}

void ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    initiate(_MPI_Ssend,curtime, count, datatype, comm, true, dest);
    PMPI_Ssend(buf, count, datatype, dest, tag, comm);
    complete(_MPI_Ssend);
  }
  else{
    PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
}

void send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    initiate(_MPI_Send,curtime, count, datatype, comm, true, dest);
    PMPI_Send(buf, count, datatype, dest, tag, comm);
    complete(_MPI_Send);
  }
  else{
    PMPI_Send(buf, count, datatype, dest, tag, comm);
  }
}

void recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    initiate(_MPI_Recv,curtime, count, datatype, comm, false, source);
    PMPI_Recv(buf, count, datatype, source, tag, comm, status);
    complete(_MPI_Recv);
  }
  else{
    PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  }
}

void isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    volatile double itime = MPI_Wtime();
    PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Isend,curtime, itime, count, datatype, comm, request, true, dest);
  }
  else{
    PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }
}

void irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    assert(tag != internal_tag);
    volatile double itime = MPI_Wtime();
    PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Irecv,curtime, itime, count, datatype, comm, request, false, source);
  }
  else{
    PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }
}

void ibcast(void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Ibcast(buf, count, datatype, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Ibcast,curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Ibcast(buf, count, datatype, root, comm, request);
  }
}

void iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Iallreduce, curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
  }
}

void ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    volatile double itime = MPI_Wtime();
    PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Iallreduce,curtime, itime, count, datatype, comm, request);
  }
  else{
    PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
  }
}

void igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Igather, curtime, itime, recvbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_rank,comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Igatherv, curtime, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
  }
  else{
     PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
  }
}

void iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size); int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
     initiate(_MPI_Iallgather, curtime, itime, recvbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
}

void iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Iallgatherv, curtime, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
  }
  else{
    PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
  }
}

void iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    volatile double itime = MPI_Wtime();
    PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Iscatter,curtime, itime, sendbuf_size, sendtype, comm, request);
  }
  else{
    PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0;
    int comm_rank, comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    volatile double itime = MPI_Wtime();
    PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Iscatterv, curtime, itime, std::max(tot_send,(int64_t)recvcount), sendtype, comm, request);
  }
  else{
    PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
}

void ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Ireduce_scatter, curtime, itime, tot_recv, datatype, comm, request);
  }
  else{
    PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
  }
}

void ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    volatile double itime = MPI_Wtime();
    PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Ialltoall, curtime, itime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm, request);
  }
  else{
    PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
}

void ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  if (mode>=1 && track_collective){
    volatile double curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    volatile double itime = MPI_Wtime();
    PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
    itime = MPI_Wtime()-itime;
    initiate(_MPI_Ialltoallv, curtime, itime, std::max(tot_send,tot_recv), sendtype, comm, request);
  }
  else{
    PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
  }
}

void wait(MPI_Request* request, MPI_Status* status){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete(curtime,request, status);
  }
  else{
    PMPI_Wait(request, status);
  }
}

void waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete(curtime, count, array_of_requests, indx, status);
  }
  else{
    PMPI_Waitany(count, array_of_requests, indx, status);
  }
}

void waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete(curtime, incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  else{
    PMPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
}

void waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  if (mode>=1 && track_p2p){
    volatile double curtime = MPI_Wtime();
    complete(curtime,count,array_of_requests,array_of_statuses);
  }
  else{
    PMPI_Waitall(count, array_of_requests, array_of_statuses);
  }
}

void finalize(){
  if (auto_capture) stop();
  cost_models.clear();
  breakdown.clear();
  decisions.clear();
  critical_path_costs.clear();
  max_per_process_costs.clear();
  volume_costs.clear();
  new_cs.clear();
  symbol_pad_cp.clear();
  symbol_pad_ncp.clear();
  symbol_len_pad_cp.clear();
  symbol_len_pad_ncp.clear();
  symbol_timer_pad_local_cp.clear();
  symbol_timer_pad_global_cp.clear();
  symbol_timer_pad_local_pp.clear();
  symbol_timer_pad_global_pp.clear();
  symbol_timer_pad_local_vol.clear();
  symbol_timer_pad_global_vol.clear();
  symbol_order.clear();
  info_sender.clear();
  info_receiver.clear();
  if (is_world_root){
    if (flag == 1){
      stream.close();
      pattern_stream.close();
    }
  }
  PMPI_Finalize();
}

}
}
