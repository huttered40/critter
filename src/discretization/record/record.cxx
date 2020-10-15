#include "record.h"
#include "../util/util.h"
#include "../../decomposition/util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

std::vector<double> record::set_tuning_statistics(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  if (world_rank==0) { stream << pattern_count_limit << " " << pattern_count_limit << " " << pattern_error_limit << std::endl; }
  {
    double total_scheduled_comm_time=0;
    double total_scheduled_comp_time=0;
    for (auto& it : comm_pattern_map){
      auto& pattern_list = active_patterns;
      auto& key_list = active_comm_pattern_keys;
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
      auto& pattern_list = active_patterns;
      auto& key_list = active_comp_pattern_keys;
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
    if (world_rank==0){
      stream << std::endl;
/*
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
*/
      stream << "Critter computational overheads:\n";
      stream << "\tMax per-process: " << global_intercept_overhead[0] << std::endl;
      stream << "Critter communication overheads:\n";
      stream << "\tStage1, Max per-process: " << global_intercept_overhead[1] << std::endl;
      stream << "\tStage2, Max per-process: " << global_intercept_overhead[2] << std::endl;
    }
  }

  std::vector<double> tuning_stats(19);
  tuning_stats[0] = global_comm_kernel_stats[0];
  tuning_stats[1] = global_comm_kernel_stats[1];
  tuning_stats[2] = global_comm_kernel_stats[2];
  tuning_stats[3] = global_comm_kernel_stats[3];
  tuning_stats[4] = global_comm_kernel_stats[4];
  tuning_stats[5] = global_comm_kernel_stats[5];
  tuning_stats[6] = global_comm_kernel_stats[6];
  tuning_stats[7] = global_comm_kernel_stats[7];
  tuning_stats[8] = global_comp_kernel_stats[0];
  tuning_stats[9] = global_comp_kernel_stats[1];
  tuning_stats[10] = global_comp_kernel_stats[2];
  tuning_stats[11] = global_comp_kernel_stats[3];
  tuning_stats[12] = global_comp_kernel_stats[4];
  tuning_stats[13] = global_comp_kernel_stats[5];
  tuning_stats[14] = global_comp_kernel_stats[6];
  tuning_stats[15] = global_comp_kernel_stats[7];
  tuning_stats[16] = global_intercept_overhead[0];
  tuning_stats[17] = global_intercept_overhead[1];
  tuning_stats[18] = global_intercept_overhead[2];
  return tuning_stats;
}

void record::write_file(int variantID, int print_mode, double overhead_time){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

  double _wall_time = wall_timer[wall_timer.size()-1];
  _wall_time -= overhead_time;
  PMPI_Allreduce(MPI_IN_PLACE,&_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (world_rank==0){
    if (is_first_iter){
      if (autotuning_test_id==2) stream_tune << std::left << std::setw(mode_1_width) << "ID";
      stream_tune << std::left << std::setw(mode_1_width) << "TestID";
      stream_tune << std::left << std::setw(mode_1_width) << "IsOptimized";
      stream_tune << std::left << std::setw(mode_1_width) << "CompSampleAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CompStateAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CommSampleAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CommStateAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "SampleConstraintMode";
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
      stream_tune << std::left << std::setw(mode_1_width) << "CompOverhead";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead1";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead2";
      stream_tune << std::endl;

      stream_reconstruct << std::left << std::setw(mode_1_width) << "ID";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "TestID";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "Delta";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompSampleAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompStateAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommSampleAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommStateAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "SampleConstraintMode";
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
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volEstET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volEstCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volEstCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volEstCommKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "cpCommKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ppCommKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volET";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volCompTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "volCommKTime";
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
      stream_tune << std::left << std::setw(mode_1_width) << comp_sample_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comp_state_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comm_sample_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comm_state_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << sample_constraint_mode;
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
      stream_tune << std::endl;
    }
  }
  else if (print_mode==1){
    if (world_rank == 0){
      stream_reconstruct << std::left << std::setw(mode_1_width) << variantID;
      stream_reconstruct << std::left << std::setw(mode_1_width) << autotuning_test_id;
      stream_reconstruct << std::left << std::setw(mode_1_width) << tuning_delta;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_sample_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_state_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_sample_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_state_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << sample_constraint_mode;
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
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures-4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs[num_volume_measures-4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::critical_path_costs[critter::internal::decomposition::num_critical_path_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::critical_path_costs[critter::internal::decomposition::num_critical_path_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::critical_path_costs[critter::internal::decomposition::num_critical_path_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::critical_path_costs[critter::internal::decomposition::num_critical_path_measures-5];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::max_per_process_costs[critter::internal::decomposition::num_per_process_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::max_per_process_costs[critter::internal::decomposition::num_per_process_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::max_per_process_costs[critter::internal::decomposition::num_per_process_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::max_per_process_costs[critter::internal::decomposition::num_per_process_measures-5];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::volume_costs[critter::internal::decomposition::num_volume_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::volume_costs[critter::internal::decomposition::num_volume_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::volume_costs[critter::internal::decomposition::num_volume_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critter::internal::decomposition::volume_costs[critter::internal::decomposition::num_volume_measures-5];
      stream_reconstruct << std::endl;
    }
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, int print_mode, double overhead_time){}

}
}
}
