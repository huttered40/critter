#include "interface.h"
#include "util.h"
#include "record.h"
#include "container/comm_tracker.h"
#include "profiler.h"

#include <string>

void _init(int* argc, char*** argv){
  MPI_Comm_rank(MPI_COMM_WORLD,&internal::world_rank);
  internal::is_world_root = false;
  internal::debug_rank = 0;
  if (internal::world_rank == internal::debug_rank){ internal::is_world_root = true; }

  internal::mode=0;
  internal::stack_id=0;
  internal::request_id = 100;
  if (std::getenv("CRITTER_PROPAGATE_P2P") != NULL){
    internal::propagate_p2p = atoi(std::getenv("CRITTER_PROPAGATE_P2P"));
    assert(internal::propagate_p2p >= 0 && internal::propagate_p2p <= 1);
  } else{
    internal::propagate_p2p = 1;
  }
  if (std::getenv("CRITTER_PROPAGATE_COLLECTIVE") != NULL){
    internal::propagate_collective = atoi(std::getenv("CRITTER_PROPAGATE_COLLECTIVE"));
    assert(internal::propagate_collective >= 0 && internal::propagate_collective <= 1);
  } else{
    internal::propagate_collective = 1;
  }
  if (std::getenv("CRITTER_PROFILE_P2P") != NULL){
    internal::profile_p2p = atoi(std::getenv("CRITTER_PROFILE_P2P"));
    assert(internal::profile_p2p >= 0 && internal::profile_p2p <= 1);
  } else{
    internal::profile_p2p = 1;
  }
  if (std::getenv("CRITTER_PROFILE_COLLECTIVE") != NULL){
    internal::profile_collective = atoi(std::getenv("CRITTER_PROFILE_COLLECTIVE"));
    assert(internal::profile_collective >= 0 && internal::profile_collective <= 1);
  } else{
    internal::profile_collective = 1;
  }
  if (std::getenv("CRITTER_AUTO_PROFILE") != NULL){
    internal::auto_capture = atoi(std::getenv("CRITTER_AUTO_PROFILE"));
    assert(internal::auto_capture >= 0 && internal::auto_capture <= 1);
  } else{
    internal::auto_capture = 0;
  }
  if (std::getenv("CRITTER_PROFILE_BLAS1") != NULL){
    internal::profile_blas1 = atoi(std::getenv("CRITTER_PROFILE_BLAS1"));
    assert(internal::profile_blas1 >=0 && internal::profile_blas1 <= 1);
  } else{
    internal::profile_blas1 = 0;
  }
  if (std::getenv("CRITTER_PROFILE_BLAS2") != NULL){
    internal::profile_blas2 = atoi(std::getenv("CRITTER_PROFILE_BLAS2"));
    assert(internal::profile_blas2 >=0 && internal::profile_blas2 <= 1);
  } else{
    internal::profile_blas2 = 0;
  }
  if (std::getenv("CRITTER_PROFILE_BLAS3") != NULL){
    internal::profile_blas3 = atoi(std::getenv("CRITTER_PROFILE_BLAS3"));
    assert(internal::profile_blas3 >=0 && internal::profile_blas3 <= 1);
  } else{
    internal::profile_blas3 = 1;
  }
  if (std::getenv("CRITTER_PROFILE_LAPACK") != NULL){
    internal::profile_lapack = atoi(std::getenv("CRITTER_PROFILE_LAPACK"));
    assert(internal::profile_lapack >=0 && internal::profile_lapack <= 1);
  } else{
    internal::profile_lapack = 1;
  }
  if (std::getenv("CRITTER_EAGER_LIMIT") != NULL){
    internal::eager_limit = atoi(std::getenv("CRITTER_EAGER_LIMIT"));
    assert(internal::eager_limit >= 0);
  } else{
    internal::eager_limit = 32768;
  }
  if (std::getenv("CRITTER_DEBUG") != NULL){
    internal::debug = atoi(std::getenv("CRITTER_DEBUG"));
    assert(internal::debug >= 0 && internal::debug <= 1);
  } else{
    internal::debug = 0;
  }
  if (std::getenv("CRITTER_EXECUTE_KERNELS") != NULL){
    internal::execute_kernels = atoi(std::getenv("CRITTER_EXECUTE_KERNELS"));
    assert(internal::execute_kernels >= 0 && internal::execute_kernels <= 1);
  } else{
    internal::execute_kernels = 1;
  }
  if (std::getenv("CRITTER_EXECUTE_KERNELS_MAX_MESSAGE_SIZE") != NULL){
    internal::execute_kernels_max_message_size = atoi(std::getenv("CRITTER_EXECUTE_KERNELS_MAX_MESSAGE_SIZE"));
    assert(internal::execute_kernels_max_message_size >= 0);
  } else{
    internal::execute_kernels_max_message_size = 100;
  }
  if (internal::execute_kernels == 1) internal::execute_kernels_max_message_size = 0;

  internal::propagate_within_timer = true;
  internal::initialize(MPI_COMM_WORLD);
  if (internal::is_world_root && internal::debug==1){
    std::cout << "internal::propagate_p2p -- " << internal::propagate_p2p << std::endl;
    std::cout << "internal::propagate_collective -- " << internal::propagate_collective << std::endl;
    std::cout << "internal::auto_capture -- " << internal::auto_capture << std::endl;
    std::cout << "internal::profile_p2p -- " << internal::profile_p2p << std::endl;
    std::cout << "internal::profile_collective -- " << internal::profile_collective << std::endl;
    std::cout << "internal::profile_blas1 -- " << internal::profile_blas1 << std::endl;
    std::cout << "internal::profile_blas2 -- " << internal::profile_blas2 << std::endl;
    std::cout << "internal::profile_blas3 -- " << internal::profile_blas3 << std::endl;
    std::cout << "internal::profile_lapack -- " << internal::profile_lapack << std::endl;
    std::cout << "internal::eager_limit -- " << internal::eager_limit << std::endl;
    std::cout << "internal::path_decomposition -- " << internal::path_decomposition << std::endl;
    std::cout << "internal::path_count -- " << internal::path_count << std::endl;
    std::cout << "internal::cost_model -- " << internal::cost_model << std::endl;
    std::cout << "internal::path_select -- " << internal::path_select << std::endl;
    std::cout << "internal::path_measure_select -- " << internal::path_measure_select << std::endl;
    std::cout << "internal::cp_costs size -- " << internal::cp_costs.size() << std::endl;
    std::cout << "internal::max_num_tracked_kernels -- " << internal::max_num_tracked_kernels << std::endl;
  }
  if (internal::auto_capture) critter_start(MPI_COMM_WORLD);
}


void critter_start(MPI_Comm cm){
  internal::stack_id++;
  internal::reset();

  // Barrier used to make as certain as possible that 'computation_timer' starts in synch.
  PMPI_Barrier(cm);
  internal::computation_timer = MPI_Wtime();
  internal::wall_timer.push_back((double)internal::computation_timer);
}

void critter_stop(MPI_Comm cm){
  volatile auto last_time = MPI_Wtime();
  assert(internal::wall_timer.size()>0);
  internal::wall_timer[internal::wall_timer.size()-1] = last_time - internal::wall_timer[internal::wall_timer.size()-1];
  internal::stack_id--; 
  PMPI_Barrier(cm);

  internal::update_time(last_time); 
  internal::_MPI_Barrier.comm = cm;
  internal::_MPI_Barrier.partner1 = -1;
  internal::_MPI_Barrier.partner2 = -1;
  internal::profiler::propagate(internal::_MPI_Barrier);
  internal::collect_volumetric_statistics(cm);

  internal::wall_timer.pop_back();
  if (internal::stack_id==0 && internal::auto_capture==0) internal::mode = 0;// Shut off critter if outside of a window
}

void critter_record(int variantID){
  internal::print(variantID);
  internal::write_file(variantID);
}

// Profiling interface
void critter_register_timer(const char* timer_name){
  //NOTE: this routine used to specify a fixed ordering of timers
  //      across all processors (so as to enable propagation of
  //      kernel information with less overhead and correctness)
  //NOTE: mode can be 0 here.
  if (internal::path_decomposition==2) internal::register_timer(timer_name);
}

void critter_start_timer(const char* timer_name, bool propagate_within, MPI_Comm cm){
  if (internal::mode==1 && internal::path_decomposition==2){
    volatile auto save_time = MPI_Wtime();
    if (cm != MPI_COMM_NULL){
      internal::update_time(save_time);
      internal::_MPI_Barrier.comm = cm;
      internal::_MPI_Barrier.partner1 = -1;
      internal::_MPI_Barrier.partner2 = -1;
      internal::profiler::propagate(internal::_MPI_Barrier);
    }
    internal::__start_timer__(timer_name,save_time,propagate_within);
  }
}

void critter_stop_timer(const char* timer_name, MPI_Comm cm){
  if (internal::mode==1 && internal::path_decomposition==2){
    volatile auto save_time = MPI_Wtime();
    internal::__stop_timer__(timer_name,save_time);
    if (cm != MPI_COMM_NULL){
      //NOTE: Don't need to call update_time(save_time) here because stop_timer(...) has already been called.
      internal::_MPI_Barrier.comm = cm;
      internal::_MPI_Barrier.partner1 = -1;
      internal::_MPI_Barrier.partner2 = -1;
      internal::profiler::propagate(internal::_MPI_Barrier);
    }
  }
}

int critter_get_critical_path_costs(){
  return internal::num_cp_measures;
}
void critter_get_critical_path_costs(float* costs){
  std::memcpy(costs,&internal::cp_costs[0],sizeof(float)*internal::num_cp_measures);
}
int critter_get_max_per_process_costs(){
  return internal::num_pp_measures;
}
void critter_get_max_per_process_costs(float* costs){
  std::memcpy(costs,&internal::max_pp_costs[0],sizeof(float)*internal::num_pp_measures);
}
int critter_get_volumetric_costs(){
  return internal::num_vol_measures;
}
void critter_get_volumetric_costs(float* costs){
  std::memcpy(costs,&internal::vol_costs[0],sizeof(float)*internal::num_vol_measures);
}

// Auxiliary routines
void critter_set_mode(int input_mode){
  assert(input_mode >=0 && input_mode <= 1);
  internal::mode = input_mode;
}
