#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

void record::invoke(std::ofstream& Stream){}

void record::invoke(std::ostream& Stream, double* data, bool track_statistical_data_override, bool clear_statistical_data, bool print_statistical_data, bool save_statistical_data){

  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (data != nullptr) *data = critical_path_costs[num_critical_path_measures-1];

  // Lets iterate over the map to create two counters, then reduce them to get a global idea:
  //   Another idea is to cache this list over the critical path, but that might be too much.
  if (autotuning_mode>0){
    int patterns[4] = {0,0,0,0};
    double communications[4] = {0,0,0,0};
    int propagations[2] = {0,0};
    for (auto& it : comm_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      patterns[0] += pattern_list[it.second.val_index].num_schedules;
      patterns[1] += pattern_list[it.second.val_index].num_non_schedules;
      communications[0] += pattern_list[it.second.val_index].num_scheduled_units;
      communications[1] += pattern_list[it.second.val_index].num_non_scheduled_units;
      propagations[0] += pattern_list[it.second.val_index].num_propagations;
      propagations[1] += pattern_list[it.second.val_index].num_non_propagations;
    }
    for (auto& it : comp_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      patterns[2] += pattern_list[it.second.val_index].num_schedules;
      patterns[3] += pattern_list[it.second.val_index].num_non_schedules;
      communications[2] += pattern_list[it.second.val_index].num_scheduled_units;
      communications[3] += pattern_list[it.second.val_index].num_non_scheduled_units;
    }
    PMPI_Allreduce(MPI_IN_PLACE,&patterns[0],4,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    PMPI_Allreduce(MPI_IN_PLACE,&communications[0],4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    PMPI_Allreduce(MPI_IN_PLACE,&propagations[0],2,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (save_statistical_data && data != nullptr){
      data[1] = patterns[0];
      data[2] = patterns[1];
      data[3] = communications[0];
      data[4] = communications[1];
      data[5] = patterns[2];
      data[6] = patterns[3];
      data[7] = communications[2];
      data[8] = communications[3];
      //data[9] = propagations[0];
      //data[10] = propagations[1];
    }
    if (print_statistical_data){
      if (rank==0) { Stream << pattern_count_limit << " " << pattern_count_limit << " " << pattern_error_limit << std::endl; }

      double total_scheduled_comm_time=0;
      double total_scheduled_comp_time=0;
      for (auto& it : comm_pattern_map){
        auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
        auto& key_list = it.second.is_active == true ? active_comm_pattern_keys : steady_state_comm_pattern_keys;
        if (rank==0) {
          Stream << "Rank 0 Communication pattern (" << key_list[it.second.key_index].tag
                 << "," << key_list[it.second.key_index].comm_size
                 << "," << key_list[it.second.key_index].comm_color
                 << "," << key_list[it.second.key_index].msg_size
                 << "," << key_list[it.second.key_index].partner_offset
                 << ") - with byte-count " << key_list[it.second.key_index].msg_size
                 << std::endl;
          Stream << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
                 << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
                 << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
                 << ", NumScheduledBytes - " << pattern_list[it.second.val_index].num_scheduled_units
                 << ", NumSkippedBytes - " << pattern_list[it.second.val_index].num_non_scheduled_units
                 << ", NumPropagations - " << pattern_list[it.second.val_index].num_propagations
                 << ", NumSkippedPropagations - " << pattern_list[it.second.val_index].num_non_propagations
                 << ", M1 - " << pattern_list[it.second.val_index].M1
                 << ", M2 - " << pattern_list[it.second.val_index].M2
                 << std::endl;
          Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comm_pattern_param)
                 << ", StdDev - " << discretization::get_std_dev(it.second,comm_pattern_param)
                 << ", StdError - " << discretization::get_std_error(it.second,comm_pattern_param)
                 << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comm_pattern_param)
                 << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comm_pattern_param)/(2*discretization::get_estimate(it.second,comm_pattern_param))
                 << std::endl;
        }
/*
          for (auto k=0; k<it.second.save_comm_times.size(); k++){
            Stream << "\t\t\tCommTime - " << it.second.save_comm_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        total_scheduled_comm_time += pattern_list[it.second.val_index].total_exec_time;
      }
      if (rank==0) { Stream << std::endl << std::endl; }
      for (auto& it : comp_pattern_map){
        auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
        auto& key_list = it.second.is_active == true ? active_comp_pattern_keys : steady_state_comp_pattern_keys;
        if (rank==0) {
           Stream << "Rank 0 Computation pattern (" << it.first.tag
                  << "," << key_list[it.second.key_index].param1
                  << "," << key_list[it.second.key_index].param2
                  << "," << key_list[it.second.key_index].param3
                  << "," << key_list[it.second.key_index].param4
                  << "," << key_list[it.second.key_index].param5
                  << ") - with flop-count "
                  << it.first.flops
                  << std::endl;
           Stream << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
                  << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
                  << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
                  << ", NumScheduledFlops - " << pattern_list[it.second.val_index].num_scheduled_units
                  << ", NumSkippedFlops - " << pattern_list[it.second.val_index].num_non_scheduled_units
                  << ", M1 - " << pattern_list[it.second.val_index].M1
                  << ", M2 - " << pattern_list[it.second.val_index].M2
                  << std::endl;
           Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comp_pattern_param,comp_pattern_param)
                  << ", StdDev - " << discretization::get_std_dev(it.second,comp_pattern_param)
                  << ", StdError - " << discretization::get_std_error(it.second,comp_pattern_param)
                  << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comp_pattern_param)
                  << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comp_pattern_param)/(2*discretization::get_estimate(it.second,comp_pattern_param))
                  << std::endl;
         }
/*
          for (auto k=0; k<it.second.save_comp_times.size(); k++){
            Stream << "\t\t\tCompTime - " << it.second.save_comp_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        total_scheduled_comp_time += pattern_list[it.second.val_index].total_exec_time;
      }
      std::vector<double> total_schedled_comm_time_gather(world_size);
      std::vector<double> total_schedled_comp_time_gather(world_size);
      PMPI_Gather(&total_scheduled_comm_time,1,MPI_DOUBLE,&total_schedled_comm_time_gather[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      PMPI_Gather(&total_scheduled_comp_time,1,MPI_DOUBLE,&total_schedled_comp_time_gather[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      if (rank==0){
        Stream << std::endl;
        Stream << std::endl;
        Stream << "Execution path parameterization #" << comm_pattern_param << " " << comp_pattern_param << ": volumetric communication:\n";
        Stream << "\tNum scheduled patterns - " << patterns[0] << std::endl;
        Stream << "\tNum total patterns - " << patterns[0]+patterns[1] << std::endl;
        Stream << "\tPattern hit ratio - " << 1.-(patterns[0] * 1. / (patterns[0]+patterns[1])) << std::endl;
        Stream << "\tNum scheduled propagations - " << propagations[0] << std::endl;
        Stream << "\tNum total propagations - " << propagations[0]+propagations[1] << std::endl;
        Stream << "\tPropagation hit ratio - " << 1.-(propagations[0] * 1. / (propagations[0]+propagations[1])) << std::endl;
        Stream << "\tNum scheduled bytes - " << communications[0] << std::endl;
        Stream << "\tNum total bytes - " << communications[0]+communications[1] << std::endl;
        Stream << "\tCommunication byte hit ratio - " << 1. - (communications[0] * 1. / (communications[0]+communications[1])) << std::endl;
        Stream << "Execution path parameterization #" << comm_pattern_param << " " << comp_pattern_param << ": volumetric computation:\n";
        Stream << "\tNum scheduled patterns - " << patterns[2] << std::endl;
        Stream << "\tNum total patterns - " << patterns[2]+patterns[3] << std::endl;
        Stream << "\tPattern hit ratio - " << 1.-(patterns[2] * 1. / (patterns[2]+patterns[3])) << std::endl;
        Stream << "\tNum scheduled flops - " << communications[2] << std::endl;
        Stream << "\tNum total flops - " << communications[2]+communications[3] << std::endl;
        Stream << "\tComputation flop hit ratio - " << 1. - (communications[2] * 1. / (communications[2]+communications[3])) << std::endl;
        for (auto i=0; i< total_schedled_comm_time_gather.size(); i++){
          Stream << "\tScheduled communication time on rank " << i << ": " << total_schedled_comm_time_gather[i] << std::endl;
        }
        for (auto i=0; i< total_schedled_comp_time_gather.size(); i++){
          Stream << "\tScheduled computation time on rank " << i << ": " << total_schedled_comp_time_gather[i] << std::endl;
        }
      }
    }
  }
  if (rank==0){
    std::cout << "Max per-process computational overhead produced by critter - " << comp_intercept_overhead << std::endl;
    std::cout << "Max per-process stage 1 communication overhead produced by critter - " << comm_intercept_overhead_stage1 << std::endl;
    std::cout << "Max per-process stage 2 communication overhead produced by critter - " << comm_intercept_overhead_stage2 << std::endl;
    std::cout << "Max per-process stage 3 communication overhead produced by critter - " << comm_intercept_overhead_stage3 << std::endl;
    std::cout << "Max per-process stage 4 communication overhead produced by critter - " << comm_intercept_overhead_stage4 << std::endl;
  }
}

}
}
}
