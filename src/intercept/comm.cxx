#include "comm.h"
#include "../util/util.h"
#include "../discretization/util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{

void init(std::vector<std::string>& symbols){
  internal::init_symbol(symbols);
}

void start(bool schedule_kernels_override, bool force_steady_statistical_data_overide){
  internal::stack_id++;
  //if (internal::stack_id>1) { return; }
  internal::reset(schedule_kernels_override,force_steady_statistical_data_overide);
  internal::reset_counter++;

  // Barrier used to make as certain as possible that 'computation_timer' starts in synch.
  PMPI_Barrier(MPI_COMM_WORLD);
  internal::computation_timer = MPI_Wtime();
  internal::wall_timer.push_back((double)internal::computation_timer);
}

void stop(){
  volatile auto last_time = MPI_Wtime();
  assert(internal::wall_timer.size()>0);
  internal::wall_timer[internal::wall_timer.size()-1] = last_time - internal::wall_timer[internal::wall_timer.size()-1];
  internal::stack_id--; 
  //if (internal::stack_id>0) { return; }
  PMPI_Barrier(MPI_COMM_WORLD);

  internal::final_accumulate(MPI_COMM_WORLD,last_time); 
  internal::propagate(MPI_COMM_WORLD);
  internal::collect(MPI_COMM_WORLD);
  assert(internal::wall_timer.size()>0);
  internal::wall_timer.pop_back();
  if (internal::stack_id==0) internal::mode = 0;
}

void record(int variantID, int print_mode, double overhead_time){
  internal::print(variantID,print_mode,overhead_time);
  internal::write_file(variantID,print_mode,overhead_time);
}

void clear(int mode, int tag_count, int* distribution_tags){
  if (mode==1){
    internal::clear_aggregates();
    internal::generate_initial_aggregate();
  } else{
    internal::clear_counter++;
    internal::clear(tag_count, distribution_tags);
  }
}

void set_mode(int input_mode){
  if (input_mode != -1) { internal::mode = input_mode; }
  else{
    if (std::getenv("CRITTER_MODE") != NULL){
      internal::mode = atoi(std::getenv("CRITTER_MODE"));
      assert(internal::mode >=0 && internal::mode <=1);
    } else{
      internal::mode = 1;
    }
  }
}

void set_debug(int debug_mode){
  internal::autotuning_debug = debug_mode;
  if (debug_mode != 0){
    if (internal::mechanism == 1){
      internal::set_reference_values();
    }
  }
  else{
     internal::save_reference_values();
  }
}

void set_mechanism(int input_mechanism){
  if (input_mechanism != -1) { internal::mechanism = input_mechanism; }
  else{
    if (std::getenv("CRITTER_MECHANISM") != NULL){
      internal::mechanism = atoi(std::getenv("CRITTER_MECHANISM"));
      assert(internal::mechanism >=0 && internal::mechanism <=2);
    } else{
      internal::mechanism = 1;
    }
  }
}

int get_critical_path_costs(){
  return internal::get_critical_path_costs();
}
void get_critical_path_costs(float* costs){
  return internal::get_critical_path_costs(costs);
}
int get_max_per_process_costs(){
  return internal::get_max_per_process_costs();
}
void get_max_per_process_costs(float* costs){
  return internal::get_max_per_process_costs(costs);
}
int get_volumetric_costs(){
  return internal::get_volumetric_costs();
}
void get_volumetric_costs(float* costs){
  return internal::get_volumetric_costs(costs);
}

namespace internal{

// These routines aim to achieve agnosticity to mechanism.

void _init(int* argc, char*** argv){
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  is_world_root = false;
  if (std::getenv("CRITTER_DEBUG_RANK") != NULL){
    debug_rank = atof(std::getenv("CRITTER_DEBUG_RANK"));
  } else{
    debug_rank = 0;
  }
  if (world_rank == debug_rank){ is_world_root = true; }

  mode=0;
  stack_id=0;
  delete_comm = 1;
  request_id = 100;
  track_p2p = 1;
  track_collective = 1;
  if (std::getenv("CRITTER_AUTO") != NULL){
    auto_capture = atoi(std::getenv("CRITTER_AUTO"));
    assert(auto_capture >= 0 && auto_capture <= 1);
  } else{
    auto_capture = 0;
  }
  if (std::getenv("CRITTER_TRACK_BLAS1") != NULL){
    track_blas1 = atoi(std::getenv("CRITTER_TRACK_BLAS1"));
    assert(track_blas1 >=0 && track_blas1 <= 1);
  } else{
    track_blas1 = 0;
  }
  if (std::getenv("CRITTER_TRACK_BLAS2") != NULL){
    track_blas2 = atoi(std::getenv("CRITTER_TRACK_BLAS2"));
    assert(track_blas2 >=0 && track_blas2 <= 1);
  } else{
    track_blas2 = 0;
  }
  if (std::getenv("CRITTER_TRACK_BLAS3") != NULL){
    track_blas3 = atoi(std::getenv("CRITTER_TRACK_BLAS3"));
    assert(track_blas3 >=0 && track_blas3 <= 1);
  } else{
    track_blas3 = 1;
  }
  if (std::getenv("CRITTER_TRACK_LAPACK") != NULL){
    track_lapack = atoi(std::getenv("CRITTER_TRACK_LAPACK"));
    assert(track_lapack >=0 && track_lapack <= 1);
  } else{
    track_lapack = 1;
  }
  if (std::getenv("CRITTER_EAGER_LIMIT") != NULL){
    eager_limit = atoi(std::getenv("CRITTER_EAGER_LIMIT"));
    assert(eager_limit >= 0);
  } else{
    eager_limit = 32768;
  }
  if (std::getenv("CRITTER_DELETE_COMM") != NULL){
    delete_comm = atoi(std::getenv("CRITTER_DELETE_COMM"));
    assert(delete_comm >= 0 && delete_comm <= 1);
  }
  if (std::getenv("CRITTER_RESET_MATRIX") != NULL){
    reset_matrix = atof(std::getenv("CRITTER_RESET_MATRIX"));
    assert(reset_matrix >= 0 && reset_matrix <= 1); 
  } else{
    reset_matrix = 1;
  }
  if (std::getenv("CRITTER_TRACK_NUM_COMM_KERNELS") != NULL){
    comm_kernel_select_count = atof(std::getenv("CRITTER_TRACK_NUM_COMM_KERNELS"));
    assert(comm_kernel_select_count>=0);
  } else{
    comm_kernel_select_count = 0;
  }
  if (std::getenv("CRITTER_TRACK_NUM_COMP_KERNELS") != NULL){
    comp_kernel_select_count = atof(std::getenv("CRITTER_TRACK_NUM_COMP_KERNELS"));
    assert(comp_kernel_select_count>=0);
  } else{
    comp_kernel_select_count = 0;
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
  _MPI_Bsend__id = 16;
  _MPI_Send__id = 17;
  _MPI_Recv__id = 18;
  _MPI_Isend__id = 19;
  _MPI_Irecv__id = 20;
  _MPI_Ibcast__id = 21;
  _MPI_Iallreduce__id = 22;
  _MPI_Ireduce__id = 23;
  _MPI_Igather__id = 24;
  _MPI_Igatherv__id = 25;
  _MPI_Iallgather__id = 26;
  _MPI_Iallgatherv__id = 27;
  _MPI_Iscatter__id = 28;
  _MPI_Iscatterv__id = 29;
  _MPI_Ireduce_scatter__id = 30;
  _MPI_Ialltoall__id = 31;
  _MPI_Ialltoallv__id = 32;

  _BLAS_axpy__id = 100;
  _BLAS_scal__id = 101;
  _BLAS_gbmv__id = 120;
  _BLAS_gemv__id = 121;
  _BLAS_ger__id = 123;
  _BLAS_sbmv__id = 124;
  _BLAS_spmv__id = 125;
  _BLAS_spr__id = 126;
  _BLAS_spr2__id = 127;
  _BLAS_symv__id = 128;
  _BLAS_syr__id = 129;
  _BLAS_syr2__id = 130;
  _BLAS_trsv__id = 131;
  _BLAS_trmv__id = 132;
  _BLAS_tpsv__id = 133;
  _BLAS_tpmv__id = 134;
  _BLAS_tbsv__id = 135;
  _BLAS_tbmv__id = 136;
  _BLAS_gemm__id = 150;
  _BLAS_trmm__id = 151;
  _BLAS_trsm__id = 152;
  _BLAS_syrk__id = 153;
  _BLAS_syr2k__id = 154;
  _BLAS_symm__id = 155;

  _LAPACK_getrf__id = 200;
  _LAPACK_potrf__id = 201;
  _LAPACK_trtri__id = 202;
  _LAPACK_geqrf__id = 203;
  _LAPACK_orgqr__id = 204;
  _LAPACK_ormqr__id = 205;
  _LAPACK_getri__id = 206;
  _LAPACK_tpqrt__id = 207;
  _LAPACK_tpmqrt__id = 208;

  autotuning_debug = 0;
  reset_counter = 0;
  comp_kernel_counter = 0;
  comm_kernel_counter = 0;
  clear_counter = 0;
  symbol_id_count = 300;

  comp_kernel_key ex_1;
  MPI_Datatype comp_kernel_key_internal_type[2] = { MPI_INT, MPI_FLOAT };
  int comp_kernel_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comp_kernel_key_internal_type_disp[2] = { (char*)&ex_1.tag-(char*)&ex_1, (char*)&ex_1.flops-(char*)&ex_1 };
  PMPI_Type_create_struct(2,comp_kernel_key_internal_type_block_len,comp_kernel_key_internal_type_disp,comp_kernel_key_internal_type,&comp_kernel_key_type);
  PMPI_Type_commit(&comp_kernel_key_type);

  comm_kernel_key ex_2;
  MPI_Datatype comm_kernel_key_internal_type[2] = { MPI_INT, MPI_FLOAT };
  int comm_kernel_key_internal_type_block_len[2] = { 7,1 };
  MPI_Aint comm_kernel_key_internal_type_disp[2] = { (char*)&ex_2.tag-(char*)&ex_2, (char*)&ex_2.msg_size-(char*)&ex_2 };
  PMPI_Type_create_struct(2,comm_kernel_key_internal_type_block_len,comm_kernel_key_internal_type_disp,comm_kernel_key_internal_type,&comm_kernel_key_type);
  PMPI_Type_commit(&comm_kernel_key_type);


  mechanism=0;
  allocate(MPI_COMM_WORLD);
  mechanism=1;
  allocate(MPI_COMM_WORLD);
  mechanism=2;
  allocate(MPI_COMM_WORLD);

  if (std::getenv("CRITTER_MECHANISM") != NULL){
    mechanism = atoi(std::getenv("CRITTER_MECHANISM"));
    assert(mechanism >=0 && mechanism <=2);
  } else{
    mechanism = 0;
  }
  if (auto_capture) start();
}


int init(int* argc, char*** argv){
  int ret = PMPI_Init(argc,argv);
  _init(argc, argv);
  return ret;
}

int init_thread(int* argc, char*** argv, int required, int* provided){
  assert(required == MPI_THREAD_SINGLE);
  int ret = PMPI_Init_thread(argc,argv,required,provided);
  _init(argc, argv);
  return ret;
}

int finalize(){
  if (auto_capture) stop();
  internal::_finalize();
  return PMPI_Finalize();
}

int comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm){
  int ret = MPI_SUCCESS;
  if (mode){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Barrier__id,curtime, 0, MPI_CHAR, comm);
    ret = PMPI_Comm_split(comm,color,key,newcomm);
    complete_comm(_MPI_Barrier__id);
  }
  else{
    ret = PMPI_Comm_split(comm,color,key,newcomm);
  }
  exchange_communicators(comm,*newcomm);
  return ret;
}

int comm_dup(MPI_Comm comm, MPI_Comm* newcomm){
  int ret = MPI_SUCCESS;
  if (mode){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Barrier__id,curtime, 0, MPI_CHAR, comm);
    ret = PMPI_Comm_dup(comm,newcomm);
    complete_comm(_MPI_Barrier__id);
  }
  else{
    ret = PMPI_Comm_dup(comm,newcomm);
  }
  exchange_communicators(comm,*newcomm);
  return ret;
}

int comm_free(MPI_Comm* comm){
  int ret = MPI_SUCCESS;
  if (mode){
    if (delete_comm){
      ret = PMPI_Comm_free(comm);
    }
  }
  else{
    ret = PMPI_Comm_free(comm);
  }
  return ret;
}

int get_count(MPI_Status* status, MPI_Datatype, int* count){
  //TODO: Not implemented yet.
  assert(0);
  return 0;
}

int barrier(MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Barrier__id,curtime, 0, MPI_CHAR, comm);
    if (schedule_decision) ret = PMPI_Barrier(comm);
    complete_comm(_MPI_Barrier__id);
  }
  else{
    ret = PMPI_Barrier(comm);
  }
  return ret;
}

int bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Bcast__id,curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Bcast(buffer, count, datatype, root, comm);
    complete_comm(_MPI_Bcast__id);
  }
  else{
    ret = PMPI_Bcast(buffer, count, datatype, root, comm);
  }
  return ret;
}

int reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Reduce__id,curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    complete_comm(_MPI_Reduce__id);
  }
  else{
    ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
  return ret;
}

int allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = initiate_comm(_MPI_Allreduce__id,curtime, count, datatype, comm);
    if (schedule_decision) ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    complete_comm(_MPI_Allreduce__id);
  }
  else{
    ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }
  return ret;
}

int gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
           MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Gather__id,curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Gather__id);
  }
  else{
    ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int allgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Allgather__id,curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision){
      ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    }
    complete_comm(_MPI_Allgather__id);
  }
  else{
    ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  return ret;
}

int scatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
            MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Scatter__id,curtime, sendbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Scatter__id);
  }
  else{
    ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int reduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Reduce_scatter__id,curtime, tot_recv, datatype, comm);
    if (schedule_decision) ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    complete_comm(_MPI_Reduce_scatter__id);
  }
  else{
    ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
  }
  return ret;
}

int alltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount,
             MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = initiate_comm(_MPI_Alltoall__id,curtime,recvbuf_size, sendtype, comm);
    if (schedule_decision) ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
    complete_comm(_MPI_Alltoall__id);
  }
  else{
    ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
  }
  return ret;
}

int gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    bool schedule_decision = initiate_comm(_MPI_Gatherv__id,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
    complete_comm(_MPI_Gatherv__id);
   }
   else{
    ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
  }
  return ret;
}

int allgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int* recvcounts, const int* displs,
             MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Allgatherv__id,curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
    complete_comm(_MPI_Allgatherv__id);
  }
  else{
    ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
  return ret;
}

int scatterv(const void* sendbuf, const int* sendcounts, const int* displs, MPI_Datatype sendtype,
              void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0; int comm_size;MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    bool schedule_decision = initiate_comm(_MPI_Scatterv__id,curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    if (schedule_decision) ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
    complete_comm(_MPI_Scatterv__id);
  }
  else{
    ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm);
  }
  return ret;
}

int alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, MPI_Datatype sendtype, void* recvbuf,
               const int* recvcounts, const int* rdispls, MPI_Datatype recvtype, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    bool schedule_decision = initiate_comm(_MPI_Alltoallv__id,curtime, std::max(tot_send,tot_recv), sendtype, comm);
    if (schedule_decision) ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
    complete_comm(_MPI_Alltoallv__id);
  }
  else{
    ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
  }
  return ret;
}

int sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void* recvbuf, int recvcount,
              MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Sendrecv__id,curtime, std::max(sendcount,recvcount), sendtype, comm, true, dest, source);
    if (schedule_decision) ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount,
                                               recvtype, source, recvtag, comm, status);
    complete_comm(_MPI_Sendrecv__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
  }
  else{
    ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
  return ret;
}

int sendrecv_replace(void* buf, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag,
                      MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(sendtag != internal_tag); assert(recvtag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Sendrecv_replace__id,curtime, count, datatype, comm, true, dest, source);
    if (schedule_decision) ret = PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
    complete_comm(_MPI_Sendrecv_replace__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
   }
  else{
    ret = PMPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
  }
  return ret;
}

int ssend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Ssend__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Ssend__id);
  }
  else{
    ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int bsend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Bsend__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) ret = PMPI_Bsend(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Bsend__id);
  }
  else{
    ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int send(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Send__id,curtime, count, datatype, comm, true, dest);
    if (schedule_decision) ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    complete_comm(_MPI_Send__id);
  }
  else{
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
  }
  return ret;
}

int recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = initiate_comm(_MPI_Recv__id,curtime, count, datatype, comm, false, source);
    if (schedule_decision) ret = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
    complete_comm(_MPI_Recv__id,(source==MPI_ANY_SOURCE ? status->MPI_SOURCE : -1));
  }
  else{
    ret = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
  }
  return ret;
}

int isend(const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = inspect_comm(_MPI_Isend__id, curtime, count, datatype, comm, tag, true, dest);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Isend__id, itime, count, datatype, comm, request, tag, true, dest);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }
  return ret;
}

int irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    //assert(tag != internal_tag);
    bool schedule_decision = inspect_comm(_MPI_Irecv__id, curtime, count, datatype, comm, tag, false, source);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Irecv__id, itime, count, datatype, comm, request, tag, false, source);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }
  return ret;
}

int ibcast(void* buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = inspect_comm(_MPI_Ibcast__id, curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ibcast(buf, count, datatype, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Ibcast__id, itime, count, datatype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Ibcast(buf, count, datatype, root, comm, request);
  }
  return ret;
}

int iallreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                MPI_Request *request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = inspect_comm(_MPI_Iallreduce__id, curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iallreduce__id, itime, count, datatype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
  }
  return ret;
}

int ireduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    bool schedule_decision = inspect_comm(_MPI_Iallreduce__id, curtime, count, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iallreduce__id, itime, count, datatype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
  }
  return ret;
}

int igather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
             int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = inspect_comm(_MPI_Igather__id, curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Igather__id, itime, recvbuf_size, sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int igatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
              MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_rank,comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_recv += ((int*)recvcounts)[i]; }
    bool schedule_decision = inspect_comm(_MPI_Igatherv__id, curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Igatherv__id, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
     ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, request);
  }
  return ret;
}

int iallgather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
                MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size); int64_t recvbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = inspect_comm(_MPI_Iallgather__id, curtime, recvbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iallgather__id, itime, recvbuf_size, sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
  return ret;
}

int iallgatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[],
                 MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0; int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = inspect_comm(_MPI_Iallgatherv__id, curtime, std::max((int64_t)sendcount,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iallgatherv__id, itime, std::max((int64_t)sendcount,tot_recv), sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, request);
  }
  return ret;
}

int iscatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
              MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    int64_t sendbuf_size = std::max((int64_t)sendcount,(int64_t)recvcount) * comm_size;
    bool schedule_decision = inspect_comm(_MPI_Iscatter__id, curtime, sendbuf_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iscatter__id, itime, sendbuf_size, sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int iscatterv(const void* sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void* recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0;
    int comm_rank, comm_size; MPI_Comm_rank(comm, &comm_rank); MPI_Comm_size(comm, &comm_size);
    if (comm_rank == root) for (int i=0; i<comm_size; i++){ tot_send += ((int*)sendcounts)[i]; } 
    bool schedule_decision = inspect_comm(_MPI_Iscatterv__id, curtime, std::max(tot_send,(int64_t)recvcount), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Iscatterv__id, itime, std::max(tot_send,(int64_t)recvcount), sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
  }
  return ret;
}

int ireduce_scatter(const void* sendbuf, void* recvbuf, const int recvcounts[], MPI_Datatype datatype, MPI_Op op,
                     MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_recv += recvcounts[i]; }
    bool schedule_decision = inspect_comm(_MPI_Ireduce_scatter__id, curtime, tot_recv, datatype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Ireduce_scatter__id, itime, tot_recv, datatype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
  }
  return ret;
}

int ialltoall(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype,
               MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int comm_size; MPI_Comm_size(comm, &comm_size);
    bool schedule_decision = inspect_comm(_MPI_Ialltoall__id, curtime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Ialltoall__id, itime, std::max((int64_t)sendcount,(int64_t)recvcount)*comm_size, sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request);
  }
  return ret;
}

int ialltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf,
                const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm, MPI_Request* request){
  int ret = MPI_SUCCESS;
  if (mode && track_collective){
    volatile auto curtime = MPI_Wtime();
    int64_t tot_send=0, tot_recv=0;
    int comm_size; MPI_Comm_size(comm, &comm_size);
    for (int i=0; i<comm_size; i++){ tot_send += sendcounts[i]; tot_recv += recvcounts[i]; }
    bool schedule_decision = inspect_comm(_MPI_Ialltoallv__id, curtime, std::max(tot_send,tot_recv), sendtype, comm);
    if (schedule_decision){
      volatile auto itime = MPI_Wtime();
      ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
      itime = MPI_Wtime()-itime;
      initiate_comm(_MPI_Ialltoallv__id, itime, std::max(tot_send,tot_recv), sendtype, comm, request);
    } else{
      *request = request_id++;
    }
  }
  else{
    ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, request);
  }
  return ret;
}

int wait(MPI_Request* request, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = complete_comm(curtime,request, status);
  }
  else{
    ret = PMPI_Wait(request, status);
  }
  return ret;
}

int waitany(int count, MPI_Request array_of_requests[], int* indx, MPI_Status* status){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = complete_comm(curtime, count, array_of_requests, indx, status);
  }
  else{
    ret = PMPI_Waitany(count, array_of_requests, indx, status);
  }
  return ret;
}

int waitsome(int incount, MPI_Request array_of_requests[], int* outcount, int array_of_indices[], MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = complete_comm(curtime, incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  else{
    ret = PMPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
  }
  return ret;
}

int waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
  int ret = MPI_SUCCESS;
  if (mode && track_p2p){
    volatile auto curtime = MPI_Wtime();
    ret = complete_comm(curtime,count,array_of_requests,array_of_statuses);
  }
  else{
    ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
  }
  return ret;
}

int test(MPI_Request* request, int* flag, MPI_Status* status){
  assert(0);//TODO
  return 0;
}

int probe(int source, int tag, MPI_Comm comm, MPI_Status* status){
  assert(0);//TODO
  return 0;
}

int iprobe(int source, int tag, MPI_Comm comm, int* flag, MPI_Status* status){
  assert(0);//TODO
  return 0;
}

}
}
