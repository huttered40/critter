#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

void record::invoke(std::ofstream& Stream){}

void record::invoke(std::ostream& Stream, double* data, bool track_statistical_data_override, bool clear_statistical_data, bool print_statistical_data, bool save_statistical_data){
  if (data != nullptr) *data = critical_path_costs[num_critical_path_measures-1];

  // Lets iterate over the map to create two counters, then reduce them to get a global idea:
  //   Another idea is to cache this list over the critical path, but that might be too much.
  if (autotuning_mode>0){
    int patterns[4] = {0,0,0,0};
    double communications[4] = {0,0,0,0};
    for (auto& it : comm_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      patterns[0] += pattern_list[it.second.val_index].num_schedules;
      patterns[1] += pattern_list[it.second.val_index].num_non_schedules;
      communications[0] += pattern_list[it.second.val_index].num_scheduled_units;
      communications[1] += pattern_list[it.second.val_index].num_non_scheduled_units;
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
    }
    if (print_statistical_data){
      if (rank==0){
        Stream << pattern_count_limit << " " << pattern_count_limit << " " << pattern_error_limit << std::endl;
        for (auto& it : comm_pattern_map){
          auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
          auto& key_list = it.second.is_active == true ? active_comm_pattern_keys : steady_state_comm_pattern_keys;
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
                    << std::endl;
          Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comm_pattern_param)
                    << ", StdDev - " << discretization::get_std_dev(it.second,comm_pattern_param)
                    << ", StdError - " << discretization::get_std_error(it.second,comm_pattern_param)
                    << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comm_pattern_param)
                    << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comm_pattern_param)/(2*discretization::get_estimate(it.second,comm_pattern_param))
                    << std::endl;
/*
          for (auto k=0; k<it.second.save_comm_times.size(); k++){
            Stream << "\t\t\tCommTime - " << it.second.save_comm_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        }
        Stream << std::endl;
        Stream << std::endl;
        for (auto& it : comp_pattern_map){
          auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
          auto& key_list = it.second.is_active == true ? active_comp_pattern_keys : steady_state_comp_pattern_keys;
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
                    << ", NumScheduledBytes - " << pattern_list[it.second.val_index].num_scheduled_units
                    << ", NumSkippedBytes - " << pattern_list[it.second.val_index].num_non_scheduled_units
                    << std::endl;
          Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comp_pattern_param,comp_pattern_param)
                    << ", StdDev - " << discretization::get_std_dev(it.second,comp_pattern_param)
                    << ", StdError - " << discretization::get_std_error(it.second,comp_pattern_param)
                    << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comp_pattern_param)
                    << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comp_pattern_param)/(2*discretization::get_estimate(it.second,comp_pattern_param))
                    << std::endl;
/*
          for (auto k=0; k<it.second.save_comp_times.size(); k++){
            Stream << "\t\t\tCompTime - " << it.second.save_comp_times[k] << ", Arithmetic mean - " << it.second.save_arithmetic_means[k] << ", StdDev - " << it.second.save_std_dev[k] << ", StdError - " << it.second.save_std_error[k]
                      << ", 95% confidence interval len - " << it.second.save_confidence_interval[k] << ", Stopping criterion len-  " << it.second.save_confidence_interval[k]/(2.*it.second.save_arithmetic_means[k]) << std::endl;
          }
*/
        }
        Stream << std::endl;
        Stream << std::endl;
        Stream << "Execution path parameterization #" << comm_pattern_param << " " << comp_pattern_param << ": volumetric communication:\n";
        Stream << "\tNum scheduled patterns - " << patterns[0] << std::endl;
        Stream << "\tNum total patterns - " << patterns[0]+patterns[1] << std::endl;
        Stream << "\tPattern hit ratio - " << 1.-(patterns[0] * 1. / (patterns[0]+patterns[1])) << std::endl;
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
      }
    }
  }
}

}
}
}
