#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

void record::invoke(std::ofstream& Stream, double* data, bool print_statistical_data, bool save_statistical_data){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (data != nullptr){
    for (auto i=0; i<critical_path_costs.size(); i++){
      data[critical_path_costs.size()-1-i] = critical_path_costs[i];
    }
  }

  // Lets iterate over the map to create two counters, then reduce them to get a global idea:
  //   Another idea is to cache this list over the critical path, but that might be too much.
  int max_patterns[6] = {0,0,0,0,0,0};
  double max_units[6] = {0,0,0,0,0,0};
  double max_exec_times[4] = {0,0,0,0};
  int vol_patterns[6] = {0,0,0,0,0,0};
  double vol_units[6] = {0,0,0,0,0,0};
  double vol_exec_times[4] = {0,0,0,0};
  for (auto& it : comm_pattern_map){
    auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
    max_patterns[0] += pattern_list[it.second.val_index].num_schedules;
    max_patterns[1] += pattern_list[it.second.val_index].num_local_schedules;
    max_patterns[2] += pattern_list[it.second.val_index].num_non_schedules;
    max_units[0] += pattern_list[it.second.val_index].num_scheduled_units;
    max_units[1] += pattern_list[it.second.val_index].num_local_scheduled_units;
    max_units[2] += pattern_list[it.second.val_index].num_non_scheduled_units;
    max_exec_times[0] += pattern_list[it.second.val_index].total_exec_time;
    max_exec_times[1] += pattern_list[it.second.val_index].total_local_exec_time;
    vol_patterns[0] += pattern_list[it.second.val_index].num_schedules;
    vol_patterns[1] += pattern_list[it.second.val_index].num_local_schedules;
    vol_patterns[2] += pattern_list[it.second.val_index].num_non_schedules;
    vol_units[0] += pattern_list[it.second.val_index].num_scheduled_units;
    vol_units[1] += pattern_list[it.second.val_index].num_local_scheduled_units;
    vol_units[2] += pattern_list[it.second.val_index].num_non_scheduled_units;
    vol_exec_times[0] += pattern_list[it.second.val_index].total_exec_time;
    vol_exec_times[1] += pattern_list[it.second.val_index].total_local_exec_time;
  }
  for (auto& it : comp_pattern_map){
    auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
    max_patterns[3] += pattern_list[it.second.val_index].num_schedules;
    max_patterns[4] += pattern_list[it.second.val_index].num_local_schedules;
    max_patterns[5] += pattern_list[it.second.val_index].num_non_schedules;
    max_units[3] += pattern_list[it.second.val_index].num_scheduled_units;
    max_units[4] += pattern_list[it.second.val_index].num_local_scheduled_units;
    max_units[5] += pattern_list[it.second.val_index].num_non_scheduled_units;
    max_exec_times[2] += pattern_list[it.second.val_index].total_exec_time;
    max_exec_times[3] += pattern_list[it.second.val_index].total_local_exec_time;
    vol_patterns[3] += pattern_list[it.second.val_index].num_schedules;
    vol_patterns[4] += pattern_list[it.second.val_index].num_local_schedules;
    vol_patterns[5] += pattern_list[it.second.val_index].num_non_schedules;
    vol_units[3] += pattern_list[it.second.val_index].num_scheduled_units;
    vol_units[4] += pattern_list[it.second.val_index].num_local_scheduled_units;
    vol_units[5] += pattern_list[it.second.val_index].num_non_scheduled_units;
    vol_exec_times[2] += pattern_list[it.second.val_index].total_exec_time;
    vol_exec_times[3] += pattern_list[it.second.val_index].total_local_exec_time;
  }
  PMPI_Allreduce(MPI_IN_PLACE,&max_patterns[0],6,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  PMPI_Allreduce(MPI_IN_PLACE,&max_units[0],6,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  PMPI_Allreduce(MPI_IN_PLACE,&max_exec_times[0],4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  PMPI_Allreduce(MPI_IN_PLACE,&vol_patterns[0],6,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  PMPI_Allreduce(MPI_IN_PLACE,&vol_units[0],6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  PMPI_Allreduce(MPI_IN_PLACE,&vol_exec_times[0],4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (save_statistical_data && data != nullptr){
    data[0] = critical_path_costs[num_critical_path_measures-1];
    data[1] = critical_path_costs[num_critical_path_measures-2];
    data[2] = critical_path_costs[num_critical_path_measures-3];
    data[3] = critical_path_costs[num_critical_path_measures-4];
    data[4] = max_patterns[0];
    data[5] = max_patterns[1];
    data[6] = max_patterns[2];
    data[7] = max_units[0];
    data[8] = max_units[1];
    data[9] = max_units[2];
    data[10] = max_exec_times[0];
    data[11] = max_exec_times[1];
    data[12] = max_patterns[3];
    data[13] = max_patterns[4];
    data[14] = max_patterns[5];
    data[15] = max_units[3];
    data[16] = max_units[4];
    data[17] = max_units[5];
    data[18] = max_exec_times[2];
    data[19] = max_exec_times[3];
    data[20] = vol_patterns[0];
    data[21] = vol_patterns[1];
    data[22] = vol_patterns[2];
    data[23] = vol_units[0];
    data[24] = vol_units[1];
    data[25] = vol_units[2];
    data[26] = vol_exec_times[0];
    data[27] = vol_exec_times[1];
    data[28] = vol_patterns[3];
    data[29] = vol_patterns[4];
    data[30] = vol_patterns[5];
    data[31] = vol_units[3];
    data[32] = vol_units[4];
    data[33] = vol_units[5];
    data[34] = vol_exec_times[2];
    data[35] = vol_exec_times[3];
    data[36] = comp_intercept_overhead;
    data[37] = comm_intercept_overhead_stage1;
    data[38] = comm_intercept_overhead_stage2;
  }
  if (print_statistical_data){
    if (world_rank==0) { Stream << pattern_count_limit << " " << pattern_count_limit << " " << pattern_error_limit << std::endl; }

    double total_scheduled_comm_time=0;
    double total_scheduled_comp_time=0;
    for (auto& it : comm_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      auto& key_list = it.second.is_active == true ? active_comm_pattern_keys : steady_state_comm_pattern_keys;
      if (world_rank==0) {
        Stream << "Rank 0 Communication pattern (" << key_list[it.second.key_index].tag
               << ",(" << key_list[it.second.key_index].dim_sizes[0] << "," << key_list[it.second.key_index].dim_sizes[1] << "," << key_list[it.second.key_index].dim_sizes[2] << ")"
               << ",(" << key_list[it.second.key_index].dim_strides[0] << "," << key_list[it.second.key_index].dim_strides[1] << "," << key_list[it.second.key_index].dim_strides[2] << ")"
               << "," << key_list[it.second.key_index].msg_size
               << "," << key_list[it.second.key_index].partner_offset
               << ") - with byte-count " << key_list[it.second.key_index].msg_size
               << std::endl;
        Stream << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
               << "\tLocalScheduledTime - " << pattern_list[it.second.val_index].total_local_exec_time
               << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
               << ", NumLocalSchedules - " << pattern_list[it.second.val_index].num_local_schedules
               << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
               << ", NumScheduledBytes - " << pattern_list[it.second.val_index].num_scheduled_units
               << ", NumLocalScheduledBytes - " << pattern_list[it.second.val_index].num_local_scheduled_units
               << ", NumSkippedBytes - " << pattern_list[it.second.val_index].num_non_scheduled_units
               << ", M1 - " << pattern_list[it.second.val_index].M1
               << ", M2 - " << pattern_list[it.second.val_index].M2
               << std::endl;
        Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comm_analysis_param)
               << ", StdDev - " << discretization::get_std_dev(it.second,comm_analysis_param)
               << ", StdError - " << discretization::get_std_error(it.second,comm_analysis_param)
               << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comm_analysis_param)
               << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comm_analysis_param)/(2*discretization::get_estimate(it.second,comm_analysis_param))
               << std::endl;
      }
      total_scheduled_comm_time += pattern_list[it.second.val_index].total_local_exec_time;
    }
    if (world_rank==0) { Stream << std::endl << std::endl; }
    for (auto& it : comp_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      auto& key_list = it.second.is_active == true ? active_comp_pattern_keys : steady_state_comp_pattern_keys;
      if (world_rank==0) {
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
                << "\tLocalScheduledTime - " << pattern_list[it.second.val_index].total_local_exec_time
                << ", NumSchedules - " << pattern_list[it.second.val_index].num_schedules
                << ", NumLocalSchedules - " << pattern_list[it.second.val_index].num_local_schedules
                << ", NumScheduleSkips - " << pattern_list[it.second.val_index].num_non_schedules
                << ", NumScheduledFlops - " << pattern_list[it.second.val_index].num_scheduled_units
                << ", NumLocalScheduledFlops - " << pattern_list[it.second.val_index].num_local_scheduled_units
                << ", NumSkippedFlops - " << pattern_list[it.second.val_index].num_non_scheduled_units
                << ", M1 - " << pattern_list[it.second.val_index].M1
                << ", M2 - " << pattern_list[it.second.val_index].M2
                << std::endl;
         Stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comp_analysis_param,comp_analysis_param)
                << ", StdDev - " << discretization::get_std_dev(it.second,comp_analysis_param)
                << ", StdError - " << discretization::get_std_error(it.second,comp_analysis_param)
                << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comp_analysis_param)
                << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comp_analysis_param)/(2*discretization::get_estimate(it.second,comp_analysis_param))
                << std::endl;
       }
      total_scheduled_comp_time += pattern_list[it.second.val_index].total_local_exec_time;
    }
/*
    std::vector<double> total_schedled_comm_time_gather(world_size);
    std::vector<double> total_schedled_comp_time_gather(world_size);
    PMPI_Gather(&total_scheduled_comm_time,1,MPI_DOUBLE,&total_schedled_comm_time_gather[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    PMPI_Gather(&total_scheduled_comp_time,1,MPI_DOUBLE,&total_schedled_comp_time_gather[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
*/
    if (world_rank==0){
      Stream << std::endl;
      Stream << "Max -- communication\n";
      Stream << "\tExecution path parameterization #" << comm_envelope_param << " " << comm_unit_param << " " << comm_analysis_param << ":\n";
      Stream << "\t\tNum scheduled max_patterns - " << max_patterns[1] << std::endl;
      Stream << "\t\tNum total max_patterns - " << max_patterns[1]+max_patterns[2] << std::endl;
      Stream << "\t\tPattern hit ratio - " << 1.-(max_patterns[1] * 1. / (max_patterns[1]+max_patterns[2])) << std::endl;
      Stream << "\t\tNum scheduled bytes - " << max_units[1] << std::endl;
      Stream << "\t\tNum total bytes - " << max_units[1]+max_units[2] << std::endl;
      Stream << "\t\tCommunication byte hit ratio - " << 1. - (max_units[1] * 1. / (max_units[1]+max_units[2])) << std::endl;
      Stream << "Max -- computation\n";
      Stream << "\tExecution path parameterization #" << comp_envelope_param << " " << comp_unit_param << " " << comp_analysis_param << ":\n";
      Stream << "\t\tNum scheduled max_patterns - " << max_patterns[4] << std::endl;
      Stream << "\t\tNum total max_patterns - " << max_patterns[4]+max_patterns[5] << std::endl;
      Stream << "\t\tPattern hit ratio - " << 1.-(max_patterns[4] * 1. / (max_patterns[4]+max_patterns[5])) << std::endl;
      Stream << "\t\tNum scheduled flops - " << max_units[4] << std::endl;
      Stream << "\t\tNum total flops - " << max_units[4]+max_units[5] << std::endl;
      Stream << "\t\tComputation flop hit ratio - " << 1. - (max_units[4] * 1. / (max_units[4]+max_units[5])) << std::endl;
      Stream << "Vol -- communication\n";
      Stream << "\tExecution path parameterization #" << comm_envelope_param << " " << comm_unit_param << " " << comm_analysis_param << ":\n";
      Stream << "\t\tNum scheduled vol_patterns - " << vol_patterns[1] << std::endl;
      Stream << "\t\tNum total vol_patterns - " << vol_patterns[1]+vol_patterns[2] << std::endl;
      Stream << "\t\tPattern hit ratio - " << 1.-(vol_patterns[1] * 1. / (vol_patterns[1]+vol_patterns[2])) << std::endl;
      Stream << "\t\tNum scheduled bytes - " << vol_units[1] << std::endl;
      Stream << "\t\tNum total bytes - " << vol_units[1]+vol_units[2] << std::endl;
      Stream << "\t\tCommunication byte hit ratio - " << 1. - (vol_units[1] * 1. / (vol_units[1]+vol_units[2])) << std::endl;
      Stream << "Vol -- computation\n";
      Stream << "\tExecution path parameterization #" << comp_envelope_param << " " << comp_unit_param << " " << comp_analysis_param << ":\n";
      Stream << "\t\tNum scheduled vol_patterns - " << vol_patterns[4] << std::endl;
      Stream << "\t\tNum total vol_patterns - " << vol_patterns[4]+vol_patterns[5] << std::endl;
      Stream << "\t\tPattern hit ratio - " << 1.-(vol_patterns[4] * 1. / (vol_patterns[4]+vol_patterns[5])) << std::endl;
      Stream << "\t\tNum scheduled flops - " << vol_units[4] << std::endl;
      Stream << "\t\tNum total flops - " << vol_units[4]+vol_units[5] << std::endl;
      Stream << "\t\tComputation flop hit ratio - " << 1. - (vol_units[4] * 1. / (vol_units[4]+vol_units[5])) << std::endl;
      Stream << "Critter computational overheads:\n";
      Stream << "\tMax per-process: " << comp_intercept_overhead << std::endl;
      Stream << "Critter communication overheads:\n";
      Stream << "\tStage1, Max per-process: " << comm_intercept_overhead_stage1 << std::endl;
      Stream << "\tStage2, Max per-process: " << comm_intercept_overhead_stage2 << std::endl;
/*
      for (auto i=0; i< total_schedled_comm_time_gather.size(); i++){
        Stream << "\tScheduled communication time on world_rank " << i << ": " << total_schedled_comm_time_gather[i] << std::endl;
      }
      for (auto i=0; i< total_schedled_comp_time_gather.size(); i++){
        Stream << "\tScheduled computation time on world_rank " << i << ": " << total_schedled_comp_time_gather[i] << std::endl;
      }
*/
    }
  }
}

void record::invoke(std::ostream& Stream, double* data, bool print_statistical_data, bool save_statistical_data){}

}
}
}
