#include "record.h"
#include "../util/util.h"
#include "../../decomposition/util/util.h"
#include "../../skeletonization/util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

std::vector<double> record::set_tuning_statistics(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  {
    double total_scheduled_comm_time=0;
    double total_scheduled_comp_time=0;
    for (auto& it : comm_kernel_map){
      auto& kernel_list = active_kernels;
      auto& key_list = active_comm_kernel_keys;
      auto& skel_kernel_list = skeletonization::active_kernels;
      auto& skel_key_list = skeletonization::active_comm_kernel_keys;
      // Don't bother printing if num_local_schedules == 0
      if (kernel_list[it.second.val_index].num_local_schedules == 0) continue;
      int skel_count = -1;
      if (skeletonization::comm_kernel_map.find(it.first) != skeletonization::comm_kernel_map.end()){
        skel_count = skel_kernel_list[skeletonization::comm_kernel_map[it.first].val_index];
      }
      if (world_rank==0) {
        stream << "Rank 0 Communication kernel (" << key_list[it.second.key_index].tag
               << ",(" << key_list[it.second.key_index].dim_sizes[0] << "," << key_list[it.second.key_index].dim_sizes[1] << ")"
               << ",(" << key_list[it.second.key_index].dim_strides[0] << "," << key_list[it.second.key_index].dim_strides[1] << ")"
               << "," << key_list[it.second.key_index].msg_size
               << "," << key_list[it.second.key_index].partner_offset
               << ") - with byte-count " << key_list[it.second.key_index].msg_size
               << std::endl;
        stream << "\tScheduledTime - " << kernel_list[it.second.val_index].total_exec_time
               << "\tLocalScheduledTime - " << kernel_list[it.second.val_index].total_local_exec_time
               << ", NumSchedulesAlongCP - " << skel_count
               << ", NumSchedules - " << kernel_list[it.second.val_index].num_schedules
               << ", NumLocalSchedules - " << kernel_list[it.second.val_index].num_local_schedules
               << ", NumScheduleSkips - " << kernel_list[it.second.val_index].num_non_schedules
               << ", NumScheduledBytes - " << kernel_list[it.second.val_index].num_scheduled_units
               << ", NumLocalScheduledBytes - " << kernel_list[it.second.val_index].num_local_scheduled_units
               << ", NumSkippedBytes - " << kernel_list[it.second.val_index].num_non_scheduled_units
               << ", M1 - " << kernel_list[it.second.val_index].M1
               << ", M2 - " << kernel_list[it.second.val_index].M2
               << std::endl;
        stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comm_analysis_param)
               << ", StdDev - " << discretization::get_std_dev(it.second,comm_analysis_param)
               << ", StdError - " << discretization::get_std_error(it.first,it.second,comm_analysis_param)
               << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param)
               << ", Stopping criterion - " << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param)/(2*discretization::get_estimate(it.second,comm_analysis_param))
               << std::endl;
      }
      total_scheduled_comm_time += kernel_list[it.second.val_index].total_local_exec_time;
    }
    if (world_rank==0) { stream << std::endl << std::endl; }
    for (auto& it : comp_kernel_map){
      auto& kernel_list = active_kernels;
      auto& key_list = active_comp_kernel_keys;
      auto& skel_kernel_list = skeletonization::active_kernels;
      auto& skel_key_list = skeletonization::active_comp_kernel_keys;
      // Don't bother printing if num_local_schedules == 0
      if (kernel_list[it.second.val_index].num_local_schedules == 0) continue;
      int skel_count = -1;
      if (skeletonization::comp_kernel_map.find(it.first) != skeletonization::comp_kernel_map.end()){
        skel_count = skel_kernel_list[skeletonization::comp_kernel_map[it.first].val_index];
      }
      if (world_rank==0) {
         stream << "Rank 0 Computation kernel (" << it.first.tag
                << "," << key_list[it.second.key_index].param1
                << "," << key_list[it.second.key_index].param2
                << "," << key_list[it.second.key_index].param3
                << "," << key_list[it.second.key_index].param4
                << "," << key_list[it.second.key_index].param5
                << ") - with flop-count "
                << it.first.flops
                << std::endl;
         stream << "\tScheduledTime - " << kernel_list[it.second.val_index].total_exec_time
                << "\tLocalScheduledTime - " << kernel_list[it.second.val_index].total_local_exec_time
                << ", NumSchedulesAlongCP - " << skel_count
                << ", NumSchedules - " << kernel_list[it.second.val_index].num_schedules
                << ", NumLocalSchedules - " << kernel_list[it.second.val_index].num_local_schedules
                << ", NumScheduleSkips - " << kernel_list[it.second.val_index].num_non_schedules
                << ", NumScheduledFlops - " << kernel_list[it.second.val_index].num_scheduled_units
                << ", NumLocalScheduledFlops - " << kernel_list[it.second.val_index].num_local_scheduled_units
                << ", NumSkippedFlops - " << kernel_list[it.second.val_index].num_non_scheduled_units
                << ", M1 - " << kernel_list[it.second.val_index].M1
                << ", M2 - " << kernel_list[it.second.val_index].M2
                << std::endl;
         stream << "\t\tEstimate - " << discretization::get_estimate(it.second,comp_analysis_param,comp_analysis_param)
                << ", StdDev - " << discretization::get_std_dev(it.second,comp_analysis_param)
                << ", StdError - " << discretization::get_std_error(it.first,it.second,comp_analysis_param)
                << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param)
                << ", Stopping criterion - " << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param)/(2*discretization::get_estimate(it.second,comp_analysis_param))
                << std::endl;
       }
      total_scheduled_comp_time += kernel_list[it.second.val_index].total_local_exec_time;
    }
    if (world_rank==0){
      stream << std::endl;
      stream << "Critter computational overheads:\n";
      stream << "\tMax per-process: " << global_intercept_overhead[0] << std::endl;
      stream << "Critter communication overheads:\n";
      stream << "\tStage1, Max per-process: " << global_intercept_overhead[1] << std::endl;
      stream << "\tStage2, Max per-process: " << global_intercept_overhead[2] << std::endl;
    }
  }

  std::vector<double> tuning_stats(13);
  tuning_stats[0] = global_comp_kernel_stats[0];
  tuning_stats[1] = global_comp_kernel_stats[1];
  tuning_stats[2] = global_comp_kernel_stats[2];
  tuning_stats[3] = global_comp_kernel_stats[3];
  tuning_stats[4] = global_comp_kernel_stats[4];
  tuning_stats[5] = global_comm_kernel_stats[0];
  tuning_stats[6] = global_comm_kernel_stats[1];
  tuning_stats[7] = global_comm_kernel_stats[2];
  tuning_stats[8] = global_comm_kernel_stats[3];
  tuning_stats[9] = global_comm_kernel_stats[4];
  tuning_stats[10] = global_intercept_overhead[0];
  tuning_stats[11] = global_intercept_overhead[1];
  tuning_stats[12] = global_intercept_overhead[2];
  return tuning_stats;
}

void record::write_file(int variantID, int print_mode, double overhead_time){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

  _wall_time -= overhead_time;
  PMPI_Allreduce(MPI_IN_PLACE,&_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (world_rank==0){
    if (is_first_iter){
      stream_tune << std::left << std::setw(mode_1_width) << "Delta";
      stream_tune << std::left << std::setw(mode_1_width) << "CompSampleAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CompStateAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CommSampleAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CommStateAggMode";
      stream_tune << std::left << std::setw(mode_1_width) << "SampleConstraintMode";
      stream_tune << std::left << std::setw(mode_1_width) << "CompKernelTransfer";
      stream_tune << std::left << std::setw(mode_1_width) << "CommKernelTransfer";
      stream_tune << std::left << std::setw(mode_1_width) << "CompKernelBuffer";
      stream_tune << std::left << std::setw(mode_1_width) << "CompReplace";
      stream_tune << std::left << std::setw(mode_1_width) << "CommReplace";
      stream_tune << std::left << std::setw(mode_1_width) << "ErrorLimit";
      stream_tune << std::left << std::setw(mode_1_width) << "TuningTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCompKs";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCompKSkips";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCompKFlops";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCompKFlopSkips";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCompKTime";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCommKs";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCommKSkips";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCommKBytes";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCommKByteSkips";
      stream_tune << std::left << std::setw(mode_1_width) << "MaxCommKTime";
      stream_tune << std::left << std::setw(mode_1_width) << "CompOverhead";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead1";
      stream_tune << std::left << std::setw(mode_1_width) << "CommOverhead2";
      stream_tune << std::endl;

      stream_reconstruct << std::left << std::setw(mode_1_width) << "ID";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "Delta";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompSampleAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompStateAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommSampleAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommStateAggMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "SampleConstraintMode";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompKernelTransfer";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommKernelTransfer";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompKernelBuffer";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CompReplace";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "CommReplace";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "ErrorLimit";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCompKs";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCompKSkips";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCompKFlops";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCompKFlopSkips";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCompKTime";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCommKs";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCommKSkips";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCommKBytes";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCommKByteSkips";
      stream_reconstruct << std::left << std::setw(mode_1_width) << "MaxCommKTime";
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
      stream_tune << std::left << std::setw(mode_1_width) << tuning_delta;
      stream_tune << std::left << std::setw(mode_1_width) << comp_sample_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comp_state_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comm_sample_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comm_state_aggregation_mode;
      stream_tune << std::left << std::setw(mode_1_width) << sample_constraint_mode;
      stream_tune << std::left << std::setw(mode_1_width) << comp_kernel_transfer_id;
      stream_tune << std::left << std::setw(mode_1_width) << comm_kernel_transfer_id;
      stream_tune << std::left << std::setw(mode_1_width) << comp_kernel_buffer_id;
      stream_tune << std::left << std::setw(mode_1_width) << decomposition::replace_comp;
      stream_tune << std::left << std::setw(mode_1_width) << decomposition::replace_comm;
      stream_tune << std::left << std::setw(mode_1_width) << kernel_error_limit;
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
      stream_tune << std::endl;
    }
  }
  else if (print_mode==1){
    if (world_rank == 0){
      stream_reconstruct << std::left << std::setw(mode_1_width) << variantID;
      stream_reconstruct << std::left << std::setw(mode_1_width) << tuning_delta;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_sample_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_state_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_sample_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_state_aggregation_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << sample_constraint_mode;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_kernel_transfer_id;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comm_kernel_transfer_id;
      stream_reconstruct << std::left << std::setw(mode_1_width) << comp_kernel_buffer_id;
      stream_reconstruct << std::left << std::setw(mode_1_width) << decomposition::replace_comp;
      stream_reconstruct << std::left << std::setw(mode_1_width) << decomposition::replace_comm;
      stream_reconstruct << std::left << std::setw(mode_1_width) << kernel_error_limit;
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comp_kernel_stats[0];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comp_kernel_stats[1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comp_kernel_stats[2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comp_kernel_stats[3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comp_kernel_stats[4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comm_kernel_stats[0];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comm_kernel_stats[1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comm_kernel_stats[2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comm_kernel_stats[3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << local_comm_kernel_stats[4];
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
  else if (print_mode==2){
    set_tuning_statistics();
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, int print_mode, double overhead_time){}

}
}
}
