#include "comm.h"
#include "../util/util.h"
#include "../discretization/util/util.h"
#include "../dispatch/dispatch.h"

namespace critter{

void start(bool track_statistical_data_override, bool clear_statistical_data, bool schedule_kernels_override){
  if (std::getenv("CRITTER_MODE") != NULL){
    internal::mode = atoi(std::getenv("CRITTER_MODE"));
  } else{
    internal::mode = 1;
  }
  internal::stack_id++;
  if (internal::stack_id>1) { return; }
  assert(internal::internal_comm_info.size() == 0);
  internal::wait_id=true;
  internal::reset();

  if (clear_statistical_data){
    // I don't see any reason to clear the communicator map. In fact, doing so would be harmful
    // Below could be moved to reset, but its basically harmless here
    internal::comm_pattern_param1_map.clear();
    internal::comp_pattern_param1_map.clear();
    internal::steady_state_comm_pattern_keys.clear();
    internal::active_comm_pattern_keys.clear();
    internal::steady_state_comp_pattern_keys.clear();
    internal::active_comp_pattern_keys.clear();
    internal::steady_state_patterns.clear();
    internal::active_patterns.clear();
  }
  // Reset this global variable, as we are updating it using function arguments for convenience
  if (std::getenv("CRITTER_PATTERN_PARAM") != NULL){
    internal::pattern_param = atoi(std::getenv("CRITTER_PATTERN_PARAM"));
  } else{
    internal::pattern_param = 0;
  }
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    internal::schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
  } else{
    internal::schedule_kernels = 1;
  }
  if (internal::pattern_param>0){ internal::pattern_param = track_statistical_data_override ? internal::pattern_param : 0; }
  if (internal::schedule_kernels==1){ internal::schedule_kernels = schedule_kernels_override ? internal::schedule_kernels : 0; }

  // Barrier used to make as certain as possible that 'computation_timer' starts in synch.
  PMPI_Barrier(MPI_COMM_WORLD);
  internal::computation_timer=MPI_Wtime();
}

void stop(double* data, bool track_statistical_data_override, bool clear_statistical_data, bool print_statistical_data){
  volatile double last_time = MPI_Wtime();
  internal::stack_id--; 
  if (internal::stack_id>0) { return; }
  PMPI_Barrier(MPI_COMM_WORLD);

  assert(internal::internal_comm_info.size() == 0);
  internal::final_accumulate(last_time); 
  internal::propagate(MPI_COMM_WORLD);
  internal::collect(MPI_COMM_WORLD);

  if (data != nullptr) *data = internal::critical_path_costs[internal::num_critical_path_measures-1];

  // Lets iterate over the map to create two counters, then reduce them to get a global idea:
  //   Another idea is to cache this list over the critical path, but that might be too much.
  if (internal::pattern_param>0){
    int patterns[4] = {0,0,0,0};
    double communications[4] = {0,0,0,0};
    if (internal::pattern_param==1){
      for (auto& it : internal::comm_pattern_param1_map){
        auto& pattern_list = it.second.is_active == true ? internal::active_patterns : internal::steady_state_patterns;
        patterns[0] += pattern_list[it.second.val_index].num_schedules;
        patterns[1] += pattern_list[it.second.val_index].num_non_schedules;
        communications[0] += pattern_list[it.second.val_index].num_scheduled_units;
        communications[1] += pattern_list[it.second.val_index].num_non_scheduled_units;
      }
      for (auto& it : internal::comp_pattern_param1_map){
        auto& pattern_list = it.second.is_active == true ? internal::active_patterns : internal::steady_state_patterns;
        patterns[2] += pattern_list[it.second.val_index].num_schedules;
        patterns[3] += pattern_list[it.second.val_index].num_non_schedules;
        communications[2] += pattern_list[it.second.val_index].num_scheduled_units;
        communications[3] += pattern_list[it.second.val_index].num_non_scheduled_units;
      }
    }
    PMPI_Allreduce(MPI_IN_PLACE,&patterns[0],4,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    PMPI_Allreduce(MPI_IN_PLACE,&communications[0],4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (print_statistical_data){
    if (rank==0){
      std::cout << internal::pattern_count_limit << " " << internal::pattern_count_limit << " " << internal::pattern_error_limit << std::endl;
      if (internal::pattern_param==1){
        for (auto& it : internal::comm_pattern_param1_map){
          auto& pattern_list = it.second.is_active == true ? internal::active_patterns : internal::steady_state_patterns;
          auto& key_list = it.second.is_active == true ? internal::active_comm_pattern_keys : internal::steady_state_comm_pattern_keys;
          std::cout << "Rank 0 Communication pattern (" << key_list[it.second.key_index].tag
                    << "," << key_list[it.second.key_index].comm_size
                    << "," << key_list[it.second.key_index].comm_color
                    << "," << key_list[it.second.key_index].msg_size
                    << "," << key_list[it.second.key_index].partner_offset
                    << ") - with byte-count " << key_list[it.second.key_index].msg_size  << "\n";
          std::cout << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
                    << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
                    << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
                    << ", NumScheduledBytes - " << pattern_list[it.second.val_index].num_scheduled_units
                    << ", NumSkippedBytes - " << pattern_list[it.second.val_index].num_non_scheduled_units << std::endl;
          std::cout << "\t\tArithmeticMean - " << internal::discretization::get_arithmetic_mean(it.second)
                    << ", StdDev - " << internal::discretization::get_std_dev(it.second)
                    << ", StdError - " << internal::discretization::get_std_error(it.second)
                    << ", 95% confidence interval len - " << internal::discretization::get_confidence_interval(it.second)
                    << ", Stopping criterion - " << internal::discretization::get_confidence_interval(it.second)/(2*internal::discretization::get_arithmetic_mean(it.second)) << std::endl;
/*
          for (auto k=0; k<it.second.save_comm_times.size(); k++){
            std::cout << "\t\t\tCommTime - " << it.second.save_comm_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        }
        std::cout << std::endl;
        std::cout << std::endl;
        for (auto& it : internal::comp_pattern_param1_map){
          auto& pattern_list = it.second.is_active == true ? internal::active_patterns : internal::steady_state_patterns;
          auto& key_list = it.second.is_active == true ? internal::active_comp_pattern_keys : internal::steady_state_comp_pattern_keys;
          std::cout << "Rank 0 Computation pattern (" << it.first.tag
                    << "," << key_list[it.second.key_index].param1
                    << "," << key_list[it.second.key_index].param2
                    << "," << key_list[it.second.key_index].param3
                    << "," << key_list[it.second.key_index].param4
                    << "," << key_list[it.second.key_index].param5
                    << ") - with flop-count "
                    << it.first.flops << "\n";
          std::cout << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
                    << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
                    << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
                    << ", NumScheduledBytes - " << pattern_list[it.second.val_index].num_scheduled_units
                    << ", NumSkippedBytes - " << pattern_list[it.second.val_index].num_non_scheduled_units << std::endl;
          std::cout << "\t\tArithmeticMean - " << internal::discretization::get_arithmetic_mean(it.second)
                    << ", StdDev - " << internal::discretization::get_std_dev(it.second)
                    << ", StdError - " << internal::discretization::get_std_error(it.second)
                    << ", 95% confidence interval len - " << internal::discretization::get_confidence_interval(it.second)
                    << ", Stopping criterion - " << internal::discretization::get_confidence_interval(it.second)/(2*internal::discretization::get_arithmetic_mean(it.second)) << std::endl;
/*
          for (auto k=0; k<it.second.save_comp_times.size(); k++){
            std::cout << "\t\t\tCompTime - " << it.second.save_comp_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        }
      }
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "Execution path parameterization #" << internal::pattern_param << ": volumetric communication:\n";
      std::cout << "\tNum scheduled patterns - " << patterns[0] << std::endl;
      std::cout << "\tNum total patterns - " << patterns[0]+patterns[1] << std::endl;
      std::cout << "\tPattern hit ratio - " << 1.-(patterns[0] * 1. / (patterns[0]+patterns[1])) << std::endl;
      std::cout << "\tNum scheduled bytes - " << communications[0] << std::endl;
      std::cout << "\tNum total bytes - " << communications[0]+communications[1] << std::endl;
      std::cout << "\tCommunication byte hit ratio - " << 1. - (communications[0] * 1. / (communications[0]+communications[1])) << std::endl;
      std::cout << "Execution path parameterization #" << internal::pattern_param << ": volumetric computation:\n";
      std::cout << "\tNum scheduled patterns - " << patterns[2] << std::endl;
      std::cout << "\tNum total patterns - " << patterns[2]+patterns[3] << std::endl;
      std::cout << "\tPattern hit ratio - " << 1.-(patterns[2] * 1. / (patterns[2]+patterns[3])) << std::endl;
      std::cout << "\tNum scheduled flops - " << communications[2] << std::endl;
      std::cout << "\tNum total flops - " << communications[2]+communications[3] << std::endl;
      std::cout << "\tComputation flop hit ratio - " << 1. - (communications[2] * 1. / (communications[2]+communications[3])) << std::endl;
    }
  }
  }


  internal::record(std::cout);
  if (internal::flag) {internal::record(internal::stream);}
  internal::mode = 0; internal::wait_id=false; internal::is_first_iter = false;
  internal::clear();
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
  if (std::getenv("CRITTER_OPT") != NULL){
    opt = atoi(std::getenv("CRITTER_OPT"));
    delete_comm = 0;
  } else{
    opt = 0;
  }
  if (std::getenv("CRITTER_OPT_MAX_ITER") != NULL){
    opt_max_iter = atoi(std::getenv("CRITTER_OPT_MAX_ITER"));
  } else{
    opt_max_iter = 5;
  }
  if (std::getenv("CRITTER_OPT_NUM_GRADIENT_POINTS") != NULL){
    num_gradient_points = atoi(std::getenv("CRITTER_OPT_NUM_GRADIENT_POINTS"));
  } else{
    num_gradient_points = 10;
  }
  if (std::getenv("CRITTER_OPT_GRADIENT_JUMP_SIZE") != NULL){
    gradient_jump_size = atoi(std::getenv("CRITTER_OPT_GRADIENT_JUMP_SIZE"));
  } else{
    gradient_jump_size = 5;// signifies 5%
  }
  if (std::getenv("CRITTER_MODEL_SELECT") != NULL){
    _cost_models_ = std::getenv("CRITTER_MODEL_SELECT");
  } else{
    _cost_models_ = "11";
  }
  if (std::getenv("CRITTER_SYMBOL_PATH_SELECT") != NULL){
    _symbol_path_select_ = std::getenv("CRITTER_SYMBOL_PATH_SELECT");
  } else{
    _symbol_path_select_ = "000000000";
  }
  if (std::getenv("CRITTER_COMM_PATH_SELECT") != NULL){
    _comm_path_select_ = std::getenv("CRITTER_COMM_PATH_SELECT");
  } else{
    _comm_path_select_ = "000000000";
  }
  if (std::getenv("CRITTER_VIZ_FILE") != NULL){
    flag = 1;
    file_name = std::getenv("CRITTER_VIZ_FILE");
    stream_name = file_name + ".txt";
  }
  if (std::getenv("CRITTER_MAX_NUM_SYMBOLS") != NULL){
    max_num_symbols = atoi(std::getenv("CRITTER_MAX_NUM_SYMBOLS"));
  } else{
    max_num_symbols = 15;
  }
  if (std::getenv("CRITTER_MAX_SYMBOL_LENGTH") != NULL){
    max_timer_name_length = atoi(std::getenv("CRITTER_MAX_SYMBOL_LENGTH"));
  } else{
    max_timer_name_length = 25;
  }
  if (std::getenv("CRITTER_AUTO") != NULL){
    auto_capture = atoi(std::getenv("CRITTER_AUTO"));
  } else{
    auto_capture = 0;
  }
  if (std::getenv("CRITTER_PATTERN_PARAM") != NULL){
    pattern_param = atoi(std::getenv("CRITTER_PATTERN_PARAM"));
  } else{
    pattern_param = 0;
  }
  if (std::getenv("CRITTER_PATTERN_COUNT_LIMIT") != NULL){
    pattern_count_limit = atoi(std::getenv("CRITTER_PATTERN_COUNT_LIMIT"));
  } else{
    pattern_count_limit = 1;
  }
  if (std::getenv("CRITTER_PATTERN_TIME_LIMIT") != NULL){
    pattern_time_limit = atof(std::getenv("CRITTER_PATTERN_TIME_LIMIT"));
  } else{
    pattern_time_limit = .001;
  }
  if (std::getenv("CRITTER_PATTERN_ERROR_LIMIT") != NULL){
    pattern_error_limit = atof(std::getenv("CRITTER_PATTERN_ERROR_LIMIT"));
  } else{
    pattern_error_limit = .5;
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
  if (std::getenv("CRITTER_SCHEDULE_KERNELS") != NULL){
    schedule_kernels = atoi(std::getenv("CRITTER_SCHEDULE_KERNELS"));
  } else{
    schedule_kernels = 1;
  }
  if (std::getenv("CRITTER_AUTOTUNING") != NULL){
    is_autotuning = atoi(std::getenv("CRITTER_AUTOTUNING"));
  } else{
    is_autotuning = 0;
  }
  assert(_cost_models_.size()==2);
  assert(_comm_path_select_.size()==9);
  assert(_symbol_path_select_.size()==9);
  is_first_iter = true;
  int _world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&_world_rank);
  if (_world_rank == 0){ is_world_root = true; }
  else                 { is_world_root = false; }
  if (flag == 1){
    if (_world_rank==0){
      stream.open(stream_name.c_str());
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
  PMPI_Finalize();
}

}
}
