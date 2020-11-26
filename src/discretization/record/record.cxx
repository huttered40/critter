#include "record.h"
#include "../util/util.h"
#include "../../decomposition/util/util.h"
#include "../../skeletonization/util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace discretization{

void write_comm_kernel_header(std::ofstream& _stream_, bool print_key = true){
  if (print_key){
   _stream_ << std::left << std::setw(mode_1_width) << "reset_count";
   _stream_ << std::left << std::setw(mode_1_width) << "reset_count";
   _stream_ << std::left << std::setw(mode_1_width) << "tag";
   _stream_ << std::left << std::setw(mode_1_width) << "cm_size1";
   _stream_ << std::left << std::setw(mode_1_width) << "cm_size2";
   _stream_ << std::left << std::setw(mode_1_width) << "cm_stride1";
   _stream_ << std::left << std::setw(mode_1_width) << "cm_stride2";
   _stream_ << std::left << std::setw(mode_1_width) << "msg_size";
   _stream_ << std::left << std::setw(mode_1_width) << "partner_offset";
  }
   _stream_ << std::left << std::setw(mode_1_width) << "skel_count";
   _stream_ << std::left << std::setw(mode_1_width) << "total_exec_time";
   _stream_ << std::left << std::setw(mode_1_width) << "local_exec_time";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_loc_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_non_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_sched_bytes";
   _stream_ << std::left << std::setw(mode_1_width) << "num_loc_sched_bytes";
   _stream_ << std::left << std::setw(mode_1_width) << "num_non_sched_bytes";
   _stream_ << std::left << std::setw(mode_1_width) << "M1";
   _stream_ << std::left << std::setw(mode_1_width) << "M2";
   _stream_ << std::left << std::setw(mode_1_width) << "est";
   _stream_ << std::left << std::setw(mode_1_width) << "std_dev";
   _stream_ << std::left << std::setw(mode_1_width) << "std_err";
   _stream_ << std::left << std::setw(mode_1_width) << "ci";
   _stream_ << std::left << std::setw(mode_1_width) << "ci_rel";
  if (print_key){
   _stream_ << std::endl;
  }
}
void write_comp_kernel_header(std::ofstream& _stream_, bool print_key = true){
  if (print_key){
   _stream_ << std::left << std::setw(mode_1_width) << "reset_count";
   _stream_ << std::left << std::setw(mode_1_width) << "reset_count";
   _stream_ << std::left << std::setw(mode_1_width) << "tag";
   _stream_ << std::left << std::setw(mode_1_width) << "param1";
   _stream_ << std::left << std::setw(mode_1_width) << "param2";
   _stream_ << std::left << std::setw(mode_1_width) << "param3";
   _stream_ << std::left << std::setw(mode_1_width) << "param4";
   _stream_ << std::left << std::setw(mode_1_width) << "param5";
   _stream_ << std::left << std::setw(mode_1_width) << "flops";
  }
   _stream_ << std::left << std::setw(mode_1_width) << "skel_count";
   _stream_ << std::left << std::setw(mode_1_width) << "total_exec_time";
   _stream_ << std::left << std::setw(mode_1_width) << "local_exec_time";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_loc_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_k_non_sched";
   _stream_ << std::left << std::setw(mode_1_width) << "num_sched_flops";
   _stream_ << std::left << std::setw(mode_1_width) << "num_loc_sched_flops";
   _stream_ << std::left << std::setw(mode_1_width) << "num_non_sched_flops";
   _stream_ << std::left << std::setw(mode_1_width) << "M1";
   _stream_ << std::left << std::setw(mode_1_width) << "M2";
   _stream_ << std::left << std::setw(mode_1_width) << "est";
   _stream_ << std::left << std::setw(mode_1_width) << "std_dev";
   _stream_ << std::left << std::setw(mode_1_width) << "std_err";
   _stream_ << std::left << std::setw(mode_1_width) << "ci";
   _stream_ << std::left << std::setw(mode_1_width) << "ci_rel";
  if (print_key){
   _stream_ << std::endl;
  }
}

void record::set_kernel_statistics(){
  if (is_world_root){
    int kk_count = 0;
    int kk_total = 0;
    for (auto& it : comm_kernel_map){
      if (active_kernels[it.second.val_index].num_schedules > 0) kk_total++;
    }
    for (auto& it : comm_kernel_map){
      auto& kernel_list = active_kernels;
      auto& key_list = active_comm_kernel_keys;
      auto& skel_kernel_list = skeletonization::active_kernels;
      auto& skel_key_list = skeletonization::active_comm_kernel_keys;
      if (kernel_list[it.second.val_index].num_schedules == 0) continue;
      int skel_count = -1;
      if (skeletonization::comm_kernel_map.find(it.first) != skeletonization::comm_kernel_map.end()){
        skel_count = skel_kernel_list[skeletonization::comm_kernel_map[it.first].val_index];
      }
      stream_comm_kernel << std::left << std::setw(mode_1_width) << comm_kernel_counter++;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << clear_counter + (kk_count*1. / kk_total);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].tag;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].dim_sizes[0];
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].dim_sizes[1];
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].dim_strides[0];
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].dim_strides[1];
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].msg_size;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].partner_offset;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << skel_count;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].total_exec_time;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].total_local_exec_time;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_local_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_non_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_local_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_non_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].M1;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].M2;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(it.second,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(it.second,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,it.second,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param)/discretization::get_estimate(it.second,comm_analysis_param);
      assert(comm_kernel_ref_map.find(it.first) != comm_kernel_ref_map.end());
      auto& ref_p = comm_kernel_ref_map[it.first];
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.total_exec_time;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.total_local_exec_time;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_local_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_non_schedules;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_local_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.num_non_scheduled_units;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.M1;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << ref_p.M2;
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(ref_p,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(ref_p,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,ref_p,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,ref_p,comm_analysis_param);
      stream_comm_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,ref_p,comm_analysis_param)/discretization::get_estimate(ref_p,comm_analysis_param);
      stream_comm_kernel << std::endl;
      kk_count++;
    }
    kk_count = 0;
    kk_total = 0;
    for (auto& it : comp_kernel_map){
      if (active_kernels[it.second.val_index].num_schedules > 0) kk_total++;
    }
    for (auto& it : comp_kernel_map){
      auto& kernel_list = active_kernels;
      auto& key_list = active_comp_kernel_keys;
      auto& skel_kernel_list = skeletonization::active_kernels;
      auto& skel_key_list = skeletonization::active_comp_kernel_keys;
      if (kernel_list[it.second.val_index].num_schedules == 0) continue;
      int skel_count = -1;
      if (skeletonization::comp_kernel_map.find(it.first) != skeletonization::comp_kernel_map.end()){
        skel_count = skel_kernel_list[skeletonization::comp_kernel_map[it.first].val_index];
      }
      stream_comp_kernel << std::left << std::setw(mode_1_width) << comp_kernel_counter++;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << clear_counter + (kk_count*1. / kk_total);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].tag;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].param1;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].param2;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].param3;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].param4;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].param5;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << key_list[it.second.key_index].flops;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << skel_count;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].total_exec_time;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].total_local_exec_time;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_local_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_non_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_local_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].num_non_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].M1;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << kernel_list[it.second.val_index].M2;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(it.second,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(it.second,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,it.second,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param)/discretization::get_estimate(it.second,comp_analysis_param);
      assert(comp_kernel_ref_map.find(it.first) != comp_kernel_ref_map.end());
      auto& ref_p = comp_kernel_ref_map[it.first];
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.total_exec_time;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.total_local_exec_time;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_local_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_non_schedules;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_local_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.num_non_scheduled_units;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.M1;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << ref_p.M2;
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(ref_p,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(ref_p,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,ref_p,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,ref_p,comp_analysis_param);
      stream_comp_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,ref_p,comp_analysis_param)/discretization::get_estimate(ref_p,comp_analysis_param);
      stream_comp_kernel << std::endl;
      kk_count++;
    }

    if (stop_criterion_mode==2){
      for (auto& it : comm_kernel_map){
        std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
        stream_name += "_--" + std::to_string(it.first.tag)
                    + "," + std::to_string(it.first.dim_sizes[0])
                    + "," + std::to_string(it.first.dim_sizes[1])
                    + "," + std::to_string(it.first.dim_strides[0])
                    + "," + std::to_string(it.first.dim_strides[1])
                    + "," + std::to_string(int(it.first.msg_size))
                    + "," + std::to_string(it.first.partner_offset)
                    +"--.txt"; 
        std::ofstream stream_kernel;
        stream_kernel.open(stream_name.c_str(),std::ofstream::app);
        stream_kernel << std::left << std::setw(mode_1_width) << "count";
        write_comm_kernel_header(stream_kernel,false);
        stream_kernel << std::left << std::setw(mode_1_width) << "pred_acc";
        stream_kernel << std::left << std::setw(mode_1_width) << "elasticity";
        stream_kernel << std::endl;

        int skel_count = -1;
        if (skeletonization::comm_kernel_map.find(it.first) != skeletonization::comm_kernel_map.end()){
          skel_count = skeletonization::active_kernels[skeletonization::comm_kernel_map[it.first].val_index];
        }
        double decomp_time = comm_kernel_list[it.first][comm_kernel_list[it.first].size()-1].M1;
        int decomp_num_schedules = comm_kernel_list[it.first][comm_kernel_list[it.first].size()-1].num_schedules;
        // Iterate over all entries in comm_kernel_list
        int s_count=0;
        double elasticity=0.;
        double prev_num_samples = -1;
        double prev_prediction_accuracy = -1;
        for (auto& iitt : comm_kernel_list[it.first]){
          if (s_count>0){
            double curr_num_samples = iitt.num_schedules;
            double _fix_ = iitt.M1 - decomp_time;
            if (_fix_ < 0.) _fix_ *= -1;
            double curr_prediction_accuracy = 1. - (_fix_/decomp_time);
            curr_prediction_accuracy = std::max(curr_prediction_accuracy,0.);
            if (prev_prediction_accuracy == 0){
              elasticity = curr_prediction_accuracy / ((curr_num_samples - prev_num_samples)/prev_num_samples);
            } else{
              elasticity = ((curr_prediction_accuracy - prev_prediction_accuracy)/prev_prediction_accuracy) /
                           ((curr_num_samples - prev_num_samples)/prev_num_samples);
            }
          }
          prev_num_samples = iitt.num_schedules;
          double _fix_ = iitt.M1 - decomp_time;
          if (_fix_ < 0.) _fix_ *= -1;
          prev_prediction_accuracy = 1. - (_fix_/decomp_time);
          prev_prediction_accuracy = std::max(prev_prediction_accuracy,0.);

          stream_kernel << std::left << std::setw(mode_1_width) << s_count++;
          stream_kernel << std::left << std::setw(mode_1_width) << skel_count;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.total_exec_time;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.total_local_exec_time;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_local_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_non_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_local_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_non_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.M1;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.M2;
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(iitt,comm_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(iitt,comm_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,iitt,comm_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,iitt,comm_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,iitt,comm_analysis_param)/discretization::get_estimate(iitt,comm_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << prev_prediction_accuracy;
          stream_kernel << std::left << std::setw(mode_1_width) << elasticity;
          stream_kernel << std::endl;

        }
        stream_kernel.close();
      }
      for (auto& it : comp_kernel_map){
        std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
        stream_name += "_--" + std::to_string(it.first.tag)
                    + "," + std::to_string(it.first.param1)
                    + "," + std::to_string(it.first.param2)
                    + "," + std::to_string(it.first.param3)
                    + "," + std::to_string(it.first.param4)
                    + "," + std::to_string(it.first.param5)
                    +"--.txt"; 
        std::ofstream stream_kernel;
        stream_kernel.open(stream_name.c_str(),std::ofstream::app);
        stream_kernel << std::left << std::setw(mode_1_width) << "count";
        write_comp_kernel_header(stream_kernel,false);
        stream_kernel << std::left << std::setw(mode_1_width) << "pred_acc";
        stream_kernel << std::left << std::setw(mode_1_width) << "elasticity";
        stream_kernel << std::endl;

        int skel_count = -1;
        if (skeletonization::comp_kernel_map.find(it.first) != skeletonization::comp_kernel_map.end()){
          skel_count = skeletonization::active_kernels[skeletonization::comp_kernel_map[it.first].val_index];
        }
        double decomp_time = comp_kernel_list[it.first][comp_kernel_list[it.first].size()-1].M1;
        int decomp_num_schedules = comp_kernel_list[it.first][comp_kernel_list[it.first].size()-1].num_schedules;
        // Iterate over all entries in comp_kernel_list
        int s_count=0;
        double elasticity=0.;
        double prev_num_samples = -1;
        double prev_prediction_accuracy = -1;
        for (auto& iitt : comp_kernel_list[it.first]){
          if (s_count>0){
            double curr_num_samples = iitt.num_schedules;
            double _fix_ = iitt.M1 - decomp_time;
            if (_fix_ < 0.) _fix_ *= -1;
            double curr_prediction_accuracy = 1. - (_fix_/decomp_time);
            if (prev_prediction_accuracy == 0){
              elasticity = curr_prediction_accuracy / ((curr_num_samples - prev_num_samples)/prev_num_samples);
            } else{
              elasticity = ((curr_prediction_accuracy - prev_prediction_accuracy)/prev_prediction_accuracy) /
                         ((curr_num_samples - prev_num_samples)/prev_num_samples);
            }
          }
          prev_num_samples = iitt.num_schedules;
          double _fix_ = iitt.M1 - decomp_time;
          if (_fix_ < 0.) _fix_ *= -1;
          prev_prediction_accuracy = 1. - (_fix_/decomp_time);
          prev_prediction_accuracy = std::max(prev_prediction_accuracy,0.);

          stream_kernel << std::left << std::setw(mode_1_width) << s_count++;
          stream_kernel << std::left << std::setw(mode_1_width) << skel_count;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.total_exec_time;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.total_local_exec_time;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_local_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_non_schedules;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_local_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.num_non_scheduled_units;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.M1;
          stream_kernel << std::left << std::setw(mode_1_width) << iitt.M2;
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_estimate(iitt,comp_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_dev(iitt,comp_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_std_error(it.first,iitt,comp_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,iitt,comp_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << discretization::get_confidence_interval(it.first,iitt,comp_analysis_param)/discretization::get_estimate(iitt,comp_analysis_param);
          stream_kernel << std::left << std::setw(mode_1_width) << prev_prediction_accuracy;
          stream_kernel << std::left << std::setw(mode_1_width) << elasticity;
          stream_kernel << std::endl;

        }
        stream_kernel.close();
      }
    }
  }
}

void record::set_tuning_statistics(){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  int world_rank; MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  {
    double total_scheduled_comm_time=0;
    if (decomposition::replace_comm == 1){
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
        double decomp_time = -1;
        if (decomposition::replace_comm_map_local.find(it.first) != decomposition::replace_comm_map_local.end()){
          decomp_time = decomposition::replace_comm_map_local[it.first].second / decomposition::replace_comm_map_local[it.first].first;
        } else{
          decomp_time = decomposition::replace_comm_map_global[it.first].second; // decomposition::replace_comm_map_global[it.first].first>0 ? decomposition::replace_comm_map_global[it.first].second : -1;
        }
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
	       << ", DecompTime - " << decomp_time
               << ", StdDev - " << discretization::get_std_dev(it.second,comm_analysis_param)
               << ", StdError - " << discretization::get_std_error(it.first,it.second,comm_analysis_param)
               << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param)
               << ", Stopping criterion - " << discretization::get_confidence_interval(it.first,it.second,comm_analysis_param)/discretization::get_estimate(it.second,comm_analysis_param)
               << std::endl;
      }
      total_scheduled_comm_time += kernel_list[it.second.val_index].total_local_exec_time;
    }
    }
    if (world_rank==0) { stream << std::endl << std::endl; }
    double total_scheduled_comp_time=0;
    if (decomposition::replace_comp == 1){
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
        double decomp_time = -1;
        if (decomposition::replace_comp_map_local.find(it.first) != decomposition::replace_comp_map_local.end()){
          decomp_time = decomposition::replace_comp_map_local[it.first].second / decomposition::replace_comp_map_local[it.first].first;
        } else{
          decomp_time = decomposition::replace_comp_map_global[it.first].second;//decomposition::replace_comp_map_global[it.first].first>0 ? decomposition::replace_comp_map_global[it.first].second : -1;
        }
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
	        << ", DecompTime - " << decomp_time
                << ", StdDev - " << discretization::get_std_dev(it.second,comp_analysis_param)
                << ", StdError - " << discretization::get_std_error(it.first,it.second,comp_analysis_param)
                << ", 95% confidence interval len - " << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param)
                << ", Stopping criterion - " << discretization::get_confidence_interval(it.first,it.second,comp_analysis_param)/discretization::get_estimate(it.second,comp_analysis_param)
                << std::endl;
       }
      total_scheduled_comp_time += kernel_list[it.second.val_index].total_local_exec_time;
    }
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

     write_comm_kernel_header(stream_comm_kernel);
     write_comp_kernel_header(stream_comp_kernel);

    }
  }

  // To simplify interface, 'update_analysis' is used to differentiate between printing tuning data and reconstruct data
  if (print_mode==0){
    std::vector<double> tuning_data(13);
    tuning_data[0] = global_comp_kernel_stats[0];
    tuning_data[1] = global_comp_kernel_stats[1];
    tuning_data[2] = global_comp_kernel_stats[2];
    tuning_data[3] = global_comp_kernel_stats[3];
    tuning_data[4] = global_comp_kernel_stats[4];
    tuning_data[5] = global_comm_kernel_stats[0];
    tuning_data[6] = global_comm_kernel_stats[1];
    tuning_data[7] = global_comm_kernel_stats[2];
    tuning_data[8] = global_comm_kernel_stats[3];
    tuning_data[9] = global_comm_kernel_stats[4];
    tuning_data[10] = global_intercept_overhead[0];
    tuning_data[11] = global_intercept_overhead[1];
    tuning_data[12] = global_intercept_overhead[2];
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
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs_ref[num_critical_path_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs_ref[num_critical_path_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs_ref[num_critical_path_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << critical_path_costs_ref[num_critical_path_measures-4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs_ref[num_per_process_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs_ref[num_per_process_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs_ref[num_per_process_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << max_per_process_costs_ref[num_per_process_measures-4];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs_ref[num_volume_measures-1];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs_ref[num_volume_measures-2];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs_ref[num_volume_measures-3];
      stream_reconstruct << std::left << std::setw(mode_1_width) << volume_costs_ref[num_volume_measures-4];
      stream_reconstruct << std::endl;
    }
  }
  else if (print_mode==2){
    set_tuning_statistics();
  }
  else if (print_mode==3){
    set_kernel_statistics();
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, int print_mode, double overhead_time){}

}
}
}
