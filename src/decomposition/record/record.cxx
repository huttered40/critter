#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

std::vector<std::string> parse_file_string(){
  std::string stream_name = std::getenv("CRITTER_VIZ_FILE");
  std::vector<std::string> inputs;
  auto prev=0;
  auto First=false;
  for (auto i=0; i<stream_name.size(); i++){
    if (stream_name[i]=='+'){
      if (First){
        inputs.emplace_back(stream_name.substr(prev,i-prev));
      }
      else{
        First=true;
      }
      prev=i+1;
    }
  }
  return inputs;
}

void print_inputs(int np, std::vector<std::string>& inputs){
  stream << np;
  for (auto input_str : inputs){
    stream << "\t" << input_str;
  }
}

std::string get_measure_title(size_t idx){
  std::vector<std::string> measure_titles(num_per_process_measures);
  measure_titles[0] = "est_comm_cost_bsp";
  measure_titles[1] = "est_comm_cost_ab";
  measure_titles[2] = "est_synch_cost_bsp";
  measure_titles[3] = "est_synch_cost_ab";
  measure_titles[4] = "est_comp_cost";
  measure_titles[5] = "idle time";
  measure_titles[6] = "comm time";
  measure_titles[7] = "synch time";
  measure_titles[8] = "comp time";
  measure_titles[9] = "comp kernel time";
  measure_titles[10] = "exec time";
  if ((cost_models[0]=='0') && (cost_models[1]=='0')){ return measure_titles[idx+2*cost_model_size]; }
  if ((cost_models[0]=='0') && (cost_models[1]=='1')){ if (idx==0) return measure_titles[0]; if (idx==1) return measure_titles[2]; return measure_titles[2+idx];}
  if ((cost_models[0]=='1') && (cost_models[1]=='0')){ if (idx==0) return measure_titles[1]; if (idx==1) return measure_titles[3]; return measure_titles[2+idx];}
  if ((cost_models[0]=='1') && (cost_models[1]=='1')){ return measure_titles[idx]; }
}

void print_cost_model_header(){
  if (cost_models[0]=='1'){
    std::cout << std::left << std::setw(mode_1_width) << "BSPCommCost";
    std::cout << std::left << std::setw(mode_1_width) << "BSPSynchCost";
  }
  if (cost_models[1]=='1'){
    std::cout << std::left << std::setw(mode_1_width) << "ABCommCost";
    std::cout << std::left << std::setw(mode_1_width) << "ABSynchCost";
  }
}

void print_cost_model_header_file(){
  if ((cost_models[0]=='1') && (cost_models[1]=='1')){
    stream << "\tBSPCommCost";
    stream << "\tABCommCost";
    stream << "\tBSPSynchCost";
    stream << "\tABSynchCost";
  } else{
    if (cost_models[0]=='1'){
      stream << "\tBSPCommCost";
      stream << "\tBSPSynchCost";
    }
    if (cost_models[1]=='1'){
      stream << "\tABCommCost";
      stream << "\tABSynchCost";
    }
  }
}

void print_header(size_t num_inputs){
  for (size_t idx = 0; idx < (num_inputs+1); idx++){
    if (idx != 0){
      stream << "\t";
    }
    stream << "Input";
  }
  print_cost_model_header_file();
  stream << "\tCompCost\tCommunicationTime\tSynchronizationTime\tComputationTime\tComputationalKernelTime\tRunTime";// critical path
  print_cost_model_header_file();
  stream << "\tCompCost\tIdleTime\tCommunicationTime\tSynchronizationTime\tComputationTime\tComputationalKernelTime\tRunTime";// per-process
  print_cost_model_header_file();
  stream << "\tCompCost\tIdleTime\tCommunicationTime\tSynchronizationTime\tComputationTime\tComputationalKernelTime\tRunTime";// volume
  for (auto i=0; i<num_tracker_critical_path_measures*comm_path_select_size+num_tracker_per_process_measures*comm_path_select_size+num_tracker_volume_measures;i++){
    for (auto& it : save_info){
     stream << "\t" << it.first;
    }
  }
}

void record::write_file(int variantID, float overhead_time){
  int temp_mode=1;
  if (std::getenv("CRITTER_MODE") != NULL){
    temp_mode = atoi(std::getenv("CRITTER_MODE"));
  }
  if (temp_mode==0){
    float _wall_time = wall_timer[wall_timer.size()-1];
    _wall_time -= overhead_time;
    PMPI_Allreduce(MPI_IN_PLACE,&_wall_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (is_world_root){
      if (is_first_iter){
        stream << std::left << std::setw(mode_1_width) << "ID";
        stream << std::left << std::setw(mode_1_width) << "Execution time";
        stream << std::left << std::setw(mode_1_width) << "Wall time";
        stream << "\n";
      }
      stream << std::left << std::setw(mode_1_width) << variantID;
      stream << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-1];
      stream << std::left << std::setw(mode_1_width) << _wall_time;
      stream << "\n";
    }
  }
  else if (temp_mode){
    auto np=0; MPI_Comm_size(MPI_COMM_WORLD,&np);
    if (is_world_root){
      auto inputs = parse_file_string();
      for (int i=0; i<list_size; i++){
        list[i]->set_header();
      }
      if (is_first_iter){
        if (variantID != -1) stream << "ID\t";
        print_header(inputs.size());
        stream << "\n";
      }
      if (variantID != -1) stream << variantID << "\t";
      print_inputs(np,inputs);
      for (size_t i=0; i<num_critical_path_measures; i++){
        stream << "\t" << critical_path_costs[i];
      }
      for (size_t i=0; i<num_per_process_measures; i++){
        stream << "\t" << max_per_process_costs[i];
      }
      for (size_t i=0; i<num_volume_measures; i++){
        stream << "\t" << volume_costs[i];
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      for (size_t j=0; j<num_tracker_volume_measures; j++){
        for (auto& it : save_info){
          stream << "\t" << it.second[j];
        }
      }
      size_t breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){	// no idle time
        if (comm_path_select[i]=='0') continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_per_process_measures; j++){
          for (auto& it : save_info){
            stream << "\t" << it.second[j];
          }
        }
        stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-1];// comp kernel time
        stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-2];// comp time
        stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-3];// idle time
        stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-4];// comp cost
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        stream << "\t" << critical_path_costs[critical_path_costs_size-comm_path_select_size+breakdown_idx];// comp kernel time
        stream << "\t" << critical_path_costs[critical_path_costs_size-2*comm_path_select_size+breakdown_idx];// comp time
        stream << "\t" << critical_path_costs[critical_path_costs_size-3*comm_path_select_size+breakdown_idx];// idle time
        stream << "\t" << critical_path_costs[critical_path_costs_size-4*comm_path_select_size+breakdown_idx];// comp cost
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_critical_path_measures; j++){
          for (auto& it : save_info){
            stream << "\t" << it.second[j];
          }
        }
        breakdown_idx++;
      }
      stream << "\n";
    }
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, float overhead_time){
  int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int temp_mode=1;
  if (std::getenv("CRITTER_MODE") != NULL){
    temp_mode = atoi(std::getenv("CRITTER_MODE"));
  }
  if (temp_mode==0){
    if (is_world_root){
      std::cout << std::left << std::setw(mode_1_width) << "Execution time:";
      std::cout << "\n";
      std::cout << std::left << std::setw(mode_1_width) << "                  ";
      std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-1];
      std::cout << "\n\n";
    }
  }
  if (temp_mode){
    if (is_world_root){
      std::cout << "\n\n";
      std::cout << std::left << std::setw(mode_1_width) << "Critical path:";
      print_cost_model_header();
      std::cout << std::left << std::setw(mode_1_width) << "CompCost";
      std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
      std::cout << std::left << std::setw(mode_1_width) << "CommTime";
      std::cout << std::left << std::setw(mode_1_width) << "SynchTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompKernelTime";
      std::cout << std::left << std::setw(mode_1_width) << "RunTime";
      std::cout << "\n";
      std::cout << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_critical_path_measures+1; i++){//+1 for idle time (which is not present in 'num_critical_path_measures'
        if (i==(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_1_width) << "0";
        else if ((i<(2*cost_model_size)) && (i%2==0)) std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[(i-1)/2+cost_model_size];
        else if (i==(2*cost_model_size)) std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[i];
        else std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[i-1];
      }
      std::cout << "\n\n";

      std::cout << std::left << std::setw(mode_1_width) << "Per-process max:";
      print_cost_model_header();
      std::cout << std::left << std::setw(mode_1_width) << "CompCost";
      std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
      std::cout << std::left << std::setw(mode_1_width) << "CommTime";
      std::cout << std::left << std::setw(mode_1_width) << "SynchTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompKernelTime";
      std::cout << std::left << std::setw(mode_1_width) << "RunTime";
      std::cout << "\n";
      std::cout << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[(i-1)/2+cost_model_size];
        else std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[i];
      }
      std::cout << "\n\n";

      std::cout << std::left << std::setw(mode_1_width) << "Volumetric avg:";
      print_cost_model_header();
      std::cout << std::left << std::setw(mode_1_width) << "CompCost";
      std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
      std::cout << std::left << std::setw(mode_1_width) << "CommTime";
      std::cout << std::left << std::setw(mode_1_width) << "SynchTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompTime";
      std::cout << std::left << std::setw(mode_1_width) << "CompKernelTime";
      std::cout << std::left << std::setw(mode_1_width) << "RunTime";
      std::cout << "\n";
      std::cout << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) std::cout << std::left << std::setw(mode_1_width) << volume_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) std::cout << std::left << std::setw(mode_1_width) << volume_costs[(i-1)/2+cost_model_size];
        else std::cout << std::left << std::setw(mode_1_width) << volume_costs[i];
      }
      std::cout << "\n\n";

      size_t breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        if (i==0){
          std::cout << std::left << std::setw(mode_1_width) << "BSPCommCost max:";
        } else if (i==1){
          std::cout << std::left << std::setw(mode_1_width) << "ABCommCost max:";
        } else if (i==2){
          std::cout << std::left << std::setw(mode_1_width) << "BSPSynchCost max:";
        } else if (i==3){
          std::cout << std::left << std::setw(mode_1_width) << "ABSynchCost max:";
        } else if (i==4){
          std::cout << std::left << std::setw(mode_1_width) << "CompCost max:";
        } else if (i==5){
          std::cout << std::left << std::setw(mode_1_width) << "CommTime max:";
        } else if (i==6){
          std::cout << std::left << std::setw(mode_1_width) << "SynchTime max:";
        } else if (i==7){
          std::cout << std::left << std::setw(mode_1_width) << "CompTime max:";
        } else if (i==8){
          std::cout << std::left << std::setw(mode_1_width) << "CompKernelTime max:";
        } else if (i==9){
          std::cout << std::left << std::setw(mode_1_width) << "RunTime max:";
        }
        std::cout << std::left << std::setw(mode_1_width) << "MeasureType";
        std::cout << std::left << std::setw(mode_1_width) << "CompKernelTime";
        std::cout << std::left << std::setw(mode_1_width) << "CompTime";
        std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
        std::cout << std::left << std::setw(mode_1_width) << "CompCost";
        print_cost_model_header();
        std::cout << std::left << std::setw(mode_1_width) << "CommTime";
        std::cout << std::left << std::setw(mode_1_width) << "SynchTime";
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationKernelTime";
        std::cout << std::left << std::setw(mode_1_width) << "path";
        std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-comm_path_select_size+breakdown_idx];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationTime";
        std::cout << std::left << std::setw(mode_1_width) << "path";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-2*comm_path_select_size+breakdown_idx];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
        std::cout << std::left << std::setw(mode_1_width) << "path";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-3*comm_path_select_size+breakdown_idx];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationCost";
        std::cout << std::left << std::setw(mode_1_width) << "path";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-4*comm_path_select_size+breakdown_idx];
        for (int j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          std::cout << "\n";
          std::cout << std::left << std::setw(mode_1_width) << it.first;
          std::cout << std::left << std::setw(mode_1_width) << "path";
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_critical_path_measures; j++){
            std::cout << std::left << std::setw(mode_1_width) << it.second[j];
          }
        }
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationKernelTime";
        std::cout << std::left << std::setw(mode_1_width) << "per-process";
        std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-1];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationTime";
        std::cout << std::left << std::setw(mode_1_width) << "per-process";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-2];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "IdleTime";
        std::cout << std::left << std::setw(mode_1_width) << "per-process";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-3];
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << "ComputationCost";
        std::cout << std::left << std::setw(mode_1_width) << "per-process";
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << 0.0;
        std::cout << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+4)-4];
        for (int j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          std::cout << "\n";
          std::cout << std::left << std::setw(mode_1_width) << it.first;
          std::cout << std::left << std::setw(mode_1_width) << "per-process";
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          std::cout << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_per_process_measures; j++){
            std::cout << std::left << std::setw(mode_1_width) << it.second[j];
          }
        }
        breakdown_idx++;
        std::cout << "\n\n";
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      std::cout << std::left << std::setw(mode_1_width) << "Volumetric avg:";
      print_cost_model_header();
      std::cout << std::left << std::setw(mode_1_width) << "CommTime";
      std::cout << std::left << std::setw(mode_1_width) << "SynchTime";
      for (auto& it : save_info){
        std::cout << "\n";
        std::cout << std::left << std::setw(mode_1_width) << it.first;
        for (size_t j=0; j<num_tracker_volume_measures; j++){
          std::cout << std::left << std::setw(mode_1_width) << it.second[j];
        }
      }
      std::cout << "\n";
    }
  }
  if ((temp_mode) && (symbol_path_select_size>0) && (symbol_timers.size()>0)){
    if (is_world_root){
      for (auto z=0; z<symbol_path_select_size; z++){
        std::cout << "***********************************************************************************************************************";
        std::vector<std::pair<std::string,std::array<float,6>>> sort_info(symbol_timers.size());
        for (int i=symbol_path_select.size(); i>=0; i--){// We just iterate over all measures regardless of whether they are set or not.
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          // Reset symbol timers and sort
          size_t j=0;
          float cp_ref,pp_ref,vol_ref;
          for (auto& it : symbol_timers){
            if (it.second.start_timer.size() != 0) { std::cout << "Symbol " << it.first << " is not handled properly\n"; assert(it.second.start_timer.size() == 0); }
              sort_info[j++] = std::make_pair(it.second.name,std::array<float,6>{it.second.cp_numcalls[z][0],it.second.cp_excl_measure[z][i],*it.second.pp_numcalls,it.second.pp_excl_measure[i],*it.second.vol_numcalls/world_size,it.second.vol_excl_measure[i]/world_size});
          }
          std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<float,6>>& vec1, std::pair<std::string,std::array<float,6>>& vec2){return vec1.second[1] > vec2.second[1];});
          if (i==2*cost_model_size){
            cp_ref = 1;	//TODO: A bit of a struggle if comm_path_select != symbol_path_select.
          } else if (i>2*cost_model_size){
            cp_ref = critical_path_costs[i-1];
          } else{
            cp_ref = critical_path_costs[i];
          }
          pp_ref = max_per_process_costs[i];
          vol_ref = volume_costs[i];
          if (cp_ref==0.) cp_ref=1.; if (pp_ref==0.) pp_ref=1.; if (vol_ref==0.) vol_ref=1.;

          // Exclusive
          std::cout << "\n\n\n\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
          std::cout << std::left << std::setw(mode_2_width) << "cp-#calls";
          std::cout << std::left << std::setw(mode_2_width) << "pp-#calls";
          std::cout << std::left << std::setw(mode_2_width) << "vol-#calls";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "cp-excl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "cp-excl";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "pp-excl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "pp-excl";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "vol-excl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "vol-excl";
          std::cout << std::left << std::setw(mode_2_width) << "cp-excl (%)";
          std::cout << std::left << std::setw(mode_2_width) << "pp-excl (%)";
          std::cout << std::left << std::setw(mode_2_width) << "vol-excl (%)";
          float cp_total_exclusive = 0.;
          float pp_total_exclusive = 0.;
          float vol_total_exclusive = 0.;
          for (auto& it : sort_info){
            std::cout << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
            std::cout << std::left << std::setw(mode_2_width) << it.second[0];
            std::cout << std::left << std::setw(mode_2_width) << it.second[2];
            std::cout << std::left << std::setw(mode_2_width) << it.second[4];
            std::cout << std::left << std::setw(mode_2_width) << it.second[1];
            std::cout << std::left << std::setw(mode_2_width) << it.second[3];
            std::cout << std::left << std::setw(mode_2_width) << it.second[5];
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_exclusive += it.second[1];
            pp_total_exclusive += it.second[3];
            vol_total_exclusive += it.second[5];
          }
          std::cout << "\n" << std::left << std::setw(max_timer_name_length) << "total";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << cp_total_exclusive;
          std::cout << std::left << std::setw(mode_2_width) << pp_total_exclusive;
          std::cout << std::left << std::setw(mode_2_width) << vol_total_exclusive;
          std::cout << std::left << std::setw(mode_2_width) << 100.*cp_total_exclusive/cp_ref;
          std::cout << std::left << std::setw(mode_2_width) << 100.*pp_total_exclusive/pp_ref;
          std::cout << std::left << std::setw(mode_2_width) << 100.*vol_total_exclusive/vol_ref;
          std::cout << "\n";

          // Reset symbol timers and sort
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          j=0;
          for (auto& it : symbol_timers){
            assert(it.second.start_timer.size() == 0);
            sort_info[j++] = std::make_pair(it.second.name,std::array<float,6>{it.second.cp_numcalls[z][0],it.second.cp_incl_measure[z][i],*it.second.pp_numcalls,it.second.pp_incl_measure[i],*it.second.vol_numcalls/world_size,it.second.vol_incl_measure[i]/world_size});
          }
          std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<float,6>>& vec1, std::pair<std::string,std::array<float,6>>& vec2){return vec1.second[1] > vec2.second[1];});
          if (i==2*cost_model_size){
            cp_ref = 100.;
          } else if (i>2*cost_model_size){
            cp_ref = critical_path_costs[i-1];
          } else{
            cp_ref = critical_path_costs[i];
          }
          pp_ref = max_per_process_costs[i];
          vol_ref = volume_costs[i];
          if (cp_ref==0.) cp_ref=1.; if (pp_ref==0.) pp_ref=1.; if (vol_ref==0.) vol_ref=1.;

          // Inclusive
          std::cout << "\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
          std::cout << std::left << std::setw(mode_2_width) << "cp-#calls";
          std::cout << std::left << std::setw(mode_2_width) << "pp-#calls";
          std::cout << std::left << std::setw(mode_2_width) << "vol-#calls";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "cp-incl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "cp-incl";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "pp-incl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "pp-incl";
          if (i>=(2*cost_model_size+1)) std::cout << std::left << std::setw(mode_2_width) << "vol-incl (s)";
          else std::cout << std::left << std::setw(mode_2_width) << "vol-incl";
          std::cout << std::left << std::setw(mode_2_width) << "cp-incl (%)";
          std::cout << std::left << std::setw(mode_2_width) << "pp-incl (%)";
          std::cout << std::left << std::setw(mode_2_width) << "vol-incl (%)";
          float cp_total_inclusive = 0.;
          float pp_total_inclusive = 0.;
          float vol_total_inclusive = 0.;
          for (auto& it : sort_info){
            std::cout << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
            std::cout << std::left << std::setw(mode_2_width) << it.second[0];
            std::cout << std::left << std::setw(mode_2_width) << it.second[2];
            std::cout << std::left << std::setw(mode_2_width) << it.second[4];
            std::cout << std::left << std::setw(mode_2_width) << it.second[1];
            std::cout << std::left << std::setw(mode_2_width) << it.second[3];
            std::cout << std::left << std::setw(mode_2_width) << it.second[5];
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            std::cout << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_inclusive = std::max(it.second[1],cp_total_inclusive);
            pp_total_inclusive = std::max(it.second[3],pp_total_inclusive);
            vol_total_inclusive = std::max(it.second[5],vol_total_inclusive);
          }
          std::cout << "\n" << std::left << std::setw(max_timer_name_length) << "total";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << "";
          std::cout << std::left << std::setw(mode_2_width) << cp_total_inclusive;
          std::cout << std::left << std::setw(mode_2_width) << pp_total_inclusive;
          std::cout << std::left << std::setw(mode_2_width) << vol_total_inclusive;
          std::cout << std::left << std::setw(mode_2_width) << 100.*cp_total_inclusive/cp_ref;
          std::cout << std::left << std::setw(mode_2_width) << 100.*pp_total_inclusive/pp_ref;
          std::cout << std::left << std::setw(mode_2_width) << 100.*vol_total_inclusive/vol_ref;
          std::cout << "\n";
        }
      }
    }
  }
}

}
}
}
