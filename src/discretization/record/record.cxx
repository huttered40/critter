#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

std::vector<double> record::set_tuning_statistics(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
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

    if (world_rank==0) { stream << pattern_count_limit << " " << pattern_count_limit << " " << pattern_error_limit << std::endl; }

    double total_scheduled_comm_time=0;
    double total_scheduled_comp_time=0;
    for (auto& it : comm_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      auto& key_list = it.second.is_active == true ? active_comm_pattern_keys : steady_state_comm_pattern_keys;
      if (world_rank==0) {
        stream << "Rank 0 Communication pattern (" << key_list[it.second.key_index].tag
               << ",(" << key_list[it.second.key_index].dim_sizes[0] << "," << key_list[it.second.key_index].dim_sizes[1] << "," << key_list[it.second.key_index].dim_sizes[2] << ")"
               << ",(" << key_list[it.second.key_index].dim_strides[0] << "," << key_list[it.second.key_index].dim_strides[1] << "," << key_list[it.second.key_index].dim_strides[2] << ")"
               << "," << key_list[it.second.key_index].msg_size
               << "," << key_list[it.second.key_index].partner_offset
               << ") - with byte-count " << key_list[it.second.key_index].msg_size
               << std::endl;
        stream << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
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
        stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comm_analysis_param)
               << ", StdDev - " << discretization::get_std_dev(it.second,comm_analysis_param)
               << ", StdError - " << discretization::get_std_error(it.second,comm_analysis_param)
               << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.second,comm_analysis_param)
               << ", Stopping criterion - " << discretization::get_confidence_interval(it.second,comm_analysis_param)/(2*discretization::get_estimate(it.second,comm_analysis_param))
               << std::endl;
      }
      total_scheduled_comm_time += pattern_list[it.second.val_index].total_local_exec_time;
    }
    if (world_rank==0) { stream << std::endl << std::endl; }
    for (auto& it : comp_pattern_map){
      auto& pattern_list = it.second.is_active == true ? active_patterns : steady_state_patterns;
      auto& key_list = it.second.is_active == true ? active_comp_pattern_keys : steady_state_comp_pattern_keys;
      if (world_rank==0) {
         stream << "Rank 0 Computation pattern (" << it.first.tag
                << "," << key_list[it.second.key_index].param1
                << "," << key_list[it.second.key_index].param2
                << "," << key_list[it.second.key_index].param3
                << "," << key_list[it.second.key_index].param4
                << "," << key_list[it.second.key_index].param5
                << ") - with flop-count "
                << it.first.flops
                << std::endl;
         stream << "\tScheduledTime - " << pattern_list[it.second.val_index].total_exec_time
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
         stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comp_analysis_param,comp_analysis_param)
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
      stream << std::endl;
      stream << "Max -- communication\n";
      stream << "\tExecution path parameterization #" << comm_envelope_param << " " << comm_unit_param << " " << comm_analysis_param << ":\n";
      stream << "\t\tNum scheduled max_patterns - " << max_patterns[1] << std::endl;
      stream << "\t\tNum total max_patterns - " << max_patterns[1]+max_patterns[2] << std::endl;
      stream << "\t\tPattern hit ratio - " << 1.-(max_patterns[1] * 1. / (max_patterns[1]+max_patterns[2])) << std::endl;
      stream << "\t\tNum scheduled bytes - " << max_units[1] << std::endl;
      stream << "\t\tNum total bytes - " << max_units[1]+max_units[2] << std::endl;
      stream << "\t\tCommunication byte hit ratio - " << 1. - (max_units[1] * 1. / (max_units[1]+max_units[2])) << std::endl;
      stream << "Max -- computation\n";
      stream << "\tExecution path parameterization #" << comp_envelope_param << " " << comp_unit_param << " " << comp_analysis_param << ":\n";
      stream << "\t\tNum scheduled max_patterns - " << max_patterns[4] << std::endl;
      stream << "\t\tNum total max_patterns - " << max_patterns[4]+max_patterns[5] << std::endl;
      stream << "\t\tPattern hit ratio - " << 1.-(max_patterns[4] * 1. / (max_patterns[4]+max_patterns[5])) << std::endl;
      stream << "\t\tNum scheduled flops - " << max_units[4] << std::endl;
      stream << "\t\tNum total flops - " << max_units[4]+max_units[5] << std::endl;
      stream << "\t\tComputation flop hit ratio - " << 1. - (max_units[4] * 1. / (max_units[4]+max_units[5])) << std::endl;
      stream << "Vol -- communication\n";
      stream << "\tExecution path parameterization #" << comm_envelope_param << " " << comm_unit_param << " " << comm_analysis_param << ":\n";
      stream << "\t\tNum scheduled vol_patterns - " << vol_patterns[1] << std::endl;
      stream << "\t\tNum total vol_patterns - " << vol_patterns[1]+vol_patterns[2] << std::endl;
      stream << "\t\tPattern hit ratio - " << 1.-(vol_patterns[1] * 1. / (vol_patterns[1]+vol_patterns[2])) << std::endl;
      stream << "\t\tNum scheduled bytes - " << vol_units[1] << std::endl;
      stream << "\t\tNum total bytes - " << vol_units[1]+vol_units[2] << std::endl;
      stream << "\t\tCommunication byte hit ratio - " << 1. - (vol_units[1] * 1. / (vol_units[1]+vol_units[2])) << std::endl;
      stream << "Vol -- computation\n";
      stream << "\tExecution path parameterization #" << comp_envelope_param << " " << comp_unit_param << " " << comp_analysis_param << ":\n";
      stream << "\t\tNum scheduled vol_patterns - " << vol_patterns[4] << std::endl;
      stream << "\t\tNum total vol_patterns - " << vol_patterns[4]+vol_patterns[5] << std::endl;
      stream << "\t\tPattern hit ratio - " << 1.-(vol_patterns[4] * 1. / (vol_patterns[4]+vol_patterns[5])) << std::endl;
      stream << "\t\tNum scheduled flops - " << vol_units[4] << std::endl;
      stream << "\t\tNum total flops - " << vol_units[4]+vol_units[5] << std::endl;
      stream << "\t\tComputation flop hit ratio - " << 1. - (vol_units[4] * 1. / (vol_units[4]+vol_units[5])) << std::endl;
      stream << "Critter computational overheads:\n";
      stream << "\tMax per-process: " << comp_intercept_overhead << std::endl;
      stream << "Critter communication overheads:\n";
      stream << "\tStage1, Max per-process: " << comm_intercept_overhead_stage1 << std::endl;
      stream << "\tStage2, Max per-process: " << comm_intercept_overhead_stage2 << std::endl;
/*
      for (auto i=0; i< total_schedled_comm_time_gather.size(); i++){
        stream << "\tScheduled communication time on world_rank " << i << ": " << total_schedled_comm_time_gather[i] << std::endl;
      }
      for (auto i=0; i< total_schedled_comp_time_gather.size(); i++){
        stream << "\tScheduled computation time on world_rank " << i << ": " << total_schedled_comp_time_gather[i] << std::endl;
      }
*/
    }

  std::vector<double> tuning_stats(35);
  tuning_stats[0] = max_patterns[0];
  tuning_stats[1] = max_patterns[1];
  tuning_stats[2] = max_patterns[2];
  tuning_stats[3] = max_units[0];
  tuning_stats[4] = max_units[1];
  tuning_stats[5] = max_units[2];
  tuning_stats[6] = max_exec_times[0];
  tuning_stats[7] = max_exec_times[1];
  tuning_stats[8] = max_patterns[3];
  tuning_stats[9] = max_patterns[4];
  tuning_stats[10] = max_patterns[5];
  tuning_stats[11] = max_units[3];
  tuning_stats[12] = max_units[4];
  tuning_stats[13] = max_units[5];
  tuning_stats[14] = max_exec_times[2];
  tuning_stats[15] = max_exec_times[3];
  tuning_stats[16] = vol_patterns[0];
  tuning_stats[17] = vol_patterns[1];
  tuning_stats[18] = vol_patterns[2];
  tuning_stats[19] = vol_units[0];
  tuning_stats[20] = vol_units[1];
  tuning_stats[21] = vol_units[2];
  tuning_stats[22] = vol_exec_times[0];
  tuning_stats[23] = vol_exec_times[1];
  tuning_stats[24] = vol_patterns[3];
  tuning_stats[25] = vol_patterns[4];
  tuning_stats[26] = vol_patterns[5];
  tuning_stats[27] = vol_units[3];
  tuning_stats[28] = vol_units[4];
  tuning_stats[29] = vol_units[5];
  tuning_stats[30] = vol_exec_times[2];
  tuning_stats[31] = vol_exec_times[3];
  tuning_stats[32] = comp_intercept_overhead;
  tuning_stats[33] = comm_intercept_overhead_stage1;
  tuning_stats[34] = comm_intercept_overhead_stage2;
  return tuning_stats;
}

void record::write_file(int variantID, int print_mode, double overhead_time){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

  double _wall_time = wall_timer[wall_timer.size()-1];
  //if (print_mode==1){ _wall_time = MPI_Wtime() - _wall_time; }
  _wall_time -= overhead_time;
  PMPI_Allreduce(MPI_IN_PLACE,&_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (world_rank==0){
    if (is_first_iter){
      if (autotuning_test_id==2) stream_tune << std::left << std::setw(mode_1_width) << "ID";
      stream_tune << std::left << std::setw(mode_1_width) << "TestID";
      stream_tune << std::left << std::setw(mode_1_width) << "IsOptimized";
      stream_tune << std::left << std::setw(mode_1_width) << "AggregationMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CommEnvelopeParam";
      stream_tune << std::left << std::setw(mode_1_width) << "CommUnitParam";
      stream_tune << std::left << std::setw(mode_1_width) << "CommAnalysisParam";
      stream_tune << std::left << std::setw(mode_1_width) << "CompEnvelopeParam";
      stream_tune << std::left << std::setw(mode_1_width) << "CompUnitParam";
      stream_tune << std::left << std::setw(mode_1_width) << "CompAnalysisParam";
      stream_tune << std::left << std::setw(mode_1_width) << "PatternCountLimit";
      stream_tune << std::left << std::setw(mode_1_width) << "PatternTimeLimit";
      stream_tune << std::left << std::setw(mode_1_width) << "PatternErrorLimit";
      stream_tune << std::left << std::setw(mode_1_width) << "TuningTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedComm";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalComm";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSkipComm";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSkipBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedCommTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalCommTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedComp";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalComp";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSkipComp";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSkipFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedCompTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxSchedLocalCompTime";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedComm";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalComm";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSkipComm";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSkipBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedCommTime";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalCommTime";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedComp";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalComp";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSkipComp";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSkipFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedCompTime";
      stream_tune << std::left << std::setw(mode_1_width) << "VolSchedLocalCompTime";
      stream_tune << std::left << std::setw(mode_1_width) << "CompOverhead";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead1";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead2";
      stream_tune << std::endl;

      stream_reconstruct << std::left << std::setw(mode_1_width) << "ID";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "TestID";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "Delta";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "AggregationMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommEnvelopeParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommUnitParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommAnalysisParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompEnvelopeParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompUnitParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompAnalysisParam";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "PatternCountLimit";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "PatternTimeLimit";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "PatternErrorLimit";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ReconstructTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpEstET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpEstCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpEstCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpEstCommKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppEstET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppEstCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppEstCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppEstCommKTime";
      stream_reconstruct << std::endl;
    }
  }

  // To simplify interface, 'update_analysis' is used to differentiate between printing tuning data and reconstruct data
  if (print_mode==0){
    auto tuning_data = set_tuning_statistics();
    if (world_rank == 0){ 
      if (autotuning_test_id==2) stream_tune << std::left << std::setw(mode_1_width) << variantID;
      stream_tune << std::left << std::setw(mode_1_width) << autotuning_test_id;
      stream_tune << std::left << std::setw(mode_1_width) << tuning_delta;
      stream_tune << std::left << std::setw(mode_1_width) << aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comm_envelope_param;
      stream_tune << std::left << std::setw(mode_1_width) << comm_unit_param;
      stream_tune << std::left << std::setw(mode_1_width) << comm_analysis_param;
      stream_tune << std::left << std::setw(mode_1_width) << comp_envelope_param;
      stream_tune << std::left << std::setw(mode_1_width) << comp_unit_param;
      stream_tune << std::left << std::setw(mode_1_width) << comp_analysis_param;
      stream_tune << std::left << std::setw(mode_1_width) << pattern_count_limit;
      stream_tune << std::left << std::setw(mode_1_width) << pattern_time_limit;
      stream_tune << std::left << std::setw(mode_1_width) << pattern_error_limit;
      stream_tune << std::left << std::setw(mode_1_width) << _wall_time;
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[0];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[1];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[2];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[3];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[4];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[5];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[6];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[7];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[8];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[9];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[10];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[11];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[12];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[13];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[14];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[15];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[16];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[17];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[18];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[19];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[20];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[21];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[22];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[23];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[24];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[25];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[26];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[27];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[28];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[29];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[30];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[31];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[32];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[33];
      stream_tune << std::left << std::setw(mode_1_width) << tuning_data[34];
      stream_tune << std::endl;
    }
  }
  else if (print_mode==1){
    if (world_rank == 0){
      stream_reconstruct << std::left << std::setw(mode_1_width) << variantID;
      stream_reconstruct << std::left << std::setw(mode_1_width) << autotuning_test_id;
      stream_reconstruct << std::left << std::setw(mode_1_width) << tuning_delta;
      stream_reconstruct << std::left << std::setw(mode_1_width) << aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_envelope_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_unit_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_analysis_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_envelope_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_unit_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_analysis_param;
      stream_reconstruct << std::left << std::setw(mode_1_width) << pattern_count_limit;
      stream_reconstruct << std::left << std::setw(mode_1_width) << pattern_time_limit;
      stream_reconstruct << std::left << std::setw(mode_1_width) << pattern_error_limit;
      stream_reconstruct << std::left << std::setw(mode_1_width) << _wall_time;
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-4];
      stream_reconstruct << std::endl;
    }
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, int print_mode, double overhead_time){}

}
}
}
