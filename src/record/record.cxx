#include "record.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{

std::vector<std::string> parse_file_string(){
  std::vector<std::string> inputs;
  auto prev=0;
  auto First=false;
  for (auto i=0; i<file_name.size(); i++){
    if (file_name[i]=='+'){
      if (First){
        inputs.emplace_back(file_name.substr(prev,i-prev));
      }
      else{
        First=true;
      }
      prev=i+1;
    }
  }
  return inputs;
}

void print_inputs(std::ofstream& Stream, int np, std::vector<std::string>& inputs){
  Stream << np;
  for (auto input_str : inputs){
    Stream << "\t" << input_str;
  }
}

std::string get_measure_title(size_t idx){
  std::vector<std::string> measure_titles(num_per_process_measures);
  measure_titles[0] = "est_comm_bsp";
  measure_titles[1] = "est_comm_ab";
  measure_titles[2] = "est_synch_bsp";
  measure_titles[3] = "est_synch_ab";
  measure_titles[4] = "idle time";
  measure_titles[5] = "comm time";
  measure_titles[6] = "synch time";
  measure_titles[7] = "datamvt time";
  measure_titles[8] = "comp time";
  measure_titles[9] = "runtime";
  if ((cost_models[0]=='0') && (cost_models[1]=='0')){ return measure_titles[idx+2*cost_model_size]; }
  if ((cost_models[0]=='0') && (cost_models[1]=='1')){ if (idx==0) return measure_titles[0]; if (idx==1) return measure_titles[2]; return measure_titles[2+idx];}
  if ((cost_models[0]=='1') && (cost_models[1]=='0')){ if (idx==0) return measure_titles[1]; if (idx==1) return measure_titles[3]; return measure_titles[2+idx];}
  if ((cost_models[0]=='1') && (cost_models[1]=='1')){ return measure_titles[idx]; }
}

void print_cost_model_header(std::ostream& Stream){
  if (cost_models[0]=='1'){
    Stream << std::left << std::setw(mode_1_width) << "BSPCommCost";
    Stream << std::left << std::setw(mode_1_width) << "BSPSynchCost";
  }
  if (cost_models[1]=='1'){
    Stream << std::left << std::setw(mode_1_width) << "ABCommCost";
    Stream << std::left << std::setw(mode_1_width) << "ABSynchCost";
  }
}

void print_cost_model_header_file(std::ofstream& Stream){
  if (cost_models[0]=='1'){
    Stream << "\tBSPCommCost";
    Stream << "\tBSPSynchCost";
  }
  if (cost_models[1]=='1'){
    Stream << "\tABCommCost";
    Stream << "\tABSynchCost";
  }
}

void print_header(std::ofstream& Stream, size_t num_inputs){
  for (size_t idx = 0; idx < (num_inputs+1); idx++){
    if (idx != 0){
      Stream << "\t";
    }
    Stream << "Input";
  }
  print_cost_model_header_file(Stream);
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// critical path
  print_cost_model_header_file(Stream);
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// per-process
  print_cost_model_header_file(Stream);
  Stream << "\tIdleTime\tCommunicationTime\tSynchronizationTime\tDataMvtTime\tComputationTime\tRunTime";// volume
  for (auto i=0; i<num_tracker_critical_path_measures*comm_path_select_size+num_tracker_per_process_measures*comm_path_select_size+num_tracker_volume_measures;i++){
    for (auto& it : save_info){
     Stream << "\t" << it.first;
    }
  }
}

void record(std::ofstream& Stream){
  assert(internal_comm_info.size() == 0);
  if (mode){
    auto np=0; MPI_Comm_size(MPI_COMM_WORLD,&np);
    if (is_world_root){
      auto inputs = parse_file_string();
      for (int i=0; i<list_size; i++){
        list[i]->set_header();
      }
      if (is_first_iter){
        print_header(Stream,inputs.size());
        Stream << "\n";
      }
      print_inputs(Stream,np,inputs);
      for (size_t i=0; i<num_critical_path_measures; i++){
        Stream << "\t" << critical_path_costs[i];
      }
      for (size_t i=0; i<num_per_process_measures; i++){
        Stream << "\t" << max_per_process_costs[i];
      }
      for (size_t i=0; i<num_volume_measures; i++){
        Stream << "\t" << volume_costs[i];
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      for (size_t j=0; j<num_tracker_volume_measures; j++){
        for (auto& it : save_info){
          Stream << "\t" << it.second[j];
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
            Stream << "\t" << it.second[j];
          }
        }
        Stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-2];// comp time
        Stream << "\t" << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-1];// idle time
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        Stream << "\t" << critical_path_costs[critical_path_costs_size-comm_path_select_size+breakdown_idx];// comp time
        Stream << "\t" << critical_path_costs[critical_path_costs_size-2*comm_path_select_size+breakdown_idx];// idle time
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
            Stream << "\t" << it.second[j];
          }
        }
        breakdown_idx++;
      }
    }
  }
}

void record(std::ostream& Stream){
  assert(internal_comm_info.size() == 0);
  int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (mode==0){
    if (is_world_root){
      Stream << std::left << std::setw(mode_1_width) << "Runtime:";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      Stream << std::left << std::setw(mode_1_width) << critical_path_costs[num_critical_path_measures-1];
      Stream << "\n\n";
    }
  }
  if (mode){
    if (is_world_root){
      Stream << "\n\n";
      Stream << std::left << std::setw(mode_1_width) << "Critical path:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "IdleTime";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_critical_path_measures+1; i++){//+1 for idle time (which is not present in 'num_critical_path_measures'
        if (i==(2*cost_model_size)) Stream << std::left << std::setw(mode_1_width) << "0";
        else if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << critical_path_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << critical_path_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << critical_path_costs[i-1];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(mode_1_width) << "Per-process max:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "IdleTime";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << max_per_process_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << max_per_process_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << max_per_process_costs[i];
      }
      Stream << "\n\n";

      Stream << std::left << std::setw(mode_1_width) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "IdleTime";
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      Stream << std::left << std::setw(mode_1_width) << "CompTime";
      Stream << std::left << std::setw(mode_1_width) << "RunTime";
      Stream << "\n";
      Stream << std::left << std::setw(mode_1_width) << "                  ";
      for (size_t i=0; i<num_volume_measures; i++){
        if ((i<(2*cost_model_size)) && (i%2==0)) Stream << std::left << std::setw(mode_1_width) << volume_costs[i/2];
        else if ((i<(2*cost_model_size)) && (i%2==1)) Stream << std::left << std::setw(mode_1_width) << volume_costs[(i-1)/2+cost_model_size];
        else Stream << std::left << std::setw(mode_1_width) << volume_costs[i];
      }
      Stream << "\n\n";

      size_t breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        if (i==0){
          Stream << std::left << std::setw(mode_1_width) << "BSPCommCost max:";
        } else if (i==1){
          Stream << std::left << std::setw(mode_1_width) << "ABCommCost max:";
        } else if (i==2){
          Stream << std::left << std::setw(mode_1_width) << "BSPSynchCost max:";
        } else if (i==3){
          Stream << std::left << std::setw(mode_1_width) << "ABSynchCost max:";
        } else if (i==4){
          Stream << std::left << std::setw(mode_1_width) << "CommTime max:";
        } else if (i==5){
          Stream << std::left << std::setw(mode_1_width) << "SynchTime max:";
        } else if (i==6){
          Stream << std::left << std::setw(mode_1_width) << "DataMvtTime max:";
        } else if (i==7){
          Stream << std::left << std::setw(mode_1_width) << "CompTime max:";
        } else if (i==8){
          Stream << std::left << std::setw(mode_1_width) << "RunTime max:";
        }
        Stream << std::left << std::setw(mode_1_width) << "MeasureType";
        Stream << std::left << std::setw(mode_1_width) << "CompTime";
        Stream << std::left << std::setw(mode_1_width) << "IdleTime";
        print_cost_model_header(Stream);
        Stream << std::left << std::setw(mode_1_width) << "CommTime";
        Stream << std::left << std::setw(mode_1_width) << "SynchTime";
        Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Computation";
        Stream << std::left << std::setw(mode_1_width) << "path";
        Stream << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-comm_path_select_size+breakdown_idx];
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Idle";
        Stream << std::left << std::setw(mode_1_width) << "path";
        Stream << std::left << std::setw(mode_1_width) << 0.0;
        Stream << std::left << std::setw(mode_1_width) << critical_path_costs[critical_path_costs_size-2*comm_path_select_size+breakdown_idx];
        for (int j=0; j<list_size; j++){
          list[j]->set_critical_path_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(mode_1_width) << it.first;
          Stream << std::left << std::setw(mode_1_width) << "path";
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_critical_path_measures; j++){
            Stream << std::left << std::setw(mode_1_width) << it.second[j];
          }
        }
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Computation";
        Stream << std::left << std::setw(mode_1_width) << "per-process";
        Stream << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-2];
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << "Idle";
        Stream << std::left << std::setw(mode_1_width) << "per-process";
        Stream << std::left << std::setw(mode_1_width) << 0.0;
        Stream << std::left << std::setw(mode_1_width) << max_per_process_costs[num_per_process_measures+(breakdown_idx+1)*(num_tracker_per_process_measures*list_size+2)-1];
        for (int j=0; j<list_size; j++){
          list[j]->set_per_process_costs(breakdown_idx);
        }
        for (auto& it : save_info){
          Stream << "\n";
          Stream << std::left << std::setw(mode_1_width) << it.first;
          Stream << std::left << std::setw(mode_1_width) << "per-process";
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          Stream << std::left << std::setw(mode_1_width) << 0.0;
          for (size_t j=0; j<num_tracker_per_process_measures; j++){
            Stream << std::left << std::setw(mode_1_width) << it.second[j];
          }
        }
        breakdown_idx++;
        Stream << "\n\n";
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_volume_costs();
      }
      Stream << std::left << std::setw(mode_1_width) << "Volume:";
      print_cost_model_header(Stream);
      Stream << std::left << std::setw(mode_1_width) << "CommTime";
      Stream << std::left << std::setw(mode_1_width) << "SynchTime";
      Stream << std::left << std::setw(mode_1_width) << "DataMvtTime";
      for (auto& it : save_info){
        Stream << "\n";
        Stream << std::left << std::setw(mode_1_width) << it.first;
        for (size_t j=0; j<num_tracker_volume_measures; j++){
          Stream << std::left << std::setw(mode_1_width) << it.second[j];
        }
      }
      Stream << "\n";
    }
  }
  if ((mode) && (symbol_path_select_size>0) && (symbol_timers.size()>0)){
    if (is_world_root){
      for (auto z=0; z<symbol_path_select_size; z++){
        Stream << "***********************************************************************************************************************";
        std::vector<std::pair<std::string,std::array<double,6>>> sort_info(symbol_timers.size());
        for (int i=symbol_path_select.size(); i>=0; i--){// We just iterate over all measures regardless of whether they are set or not.
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          // Reset symbol timers and sort
          size_t j=0;
          double cp_ref,pp_ref,vol_ref;
          for (auto& it : symbol_timers){
            if (it.second.start_timer.size() != 0) { std::cout << "Symbol " << it.first << " is not handled properly\n"; assert(it.second.start_timer.size() == 0); }
              sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{it.second.cp_numcalls[z][0],it.second.cp_excl_measure[z][i],*it.second.pp_numcalls,it.second.pp_excl_measure[i],*it.second.vol_numcalls/world_size,it.second.vol_excl_measure[i]/world_size});
          }
          std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
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
          Stream << "\n\n\n\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
          Stream << std::left << std::setw(mode_2_width) << "cp-#calls";
          Stream << std::left << std::setw(mode_2_width) << "pp-#calls";
          Stream << std::left << std::setw(mode_2_width) << "vol-#calls";
          if (i>=2*cost_model_size) Stream << std::left << std::setw(mode_2_width) << "cp-excl (s)";
          else Stream << std::left << std::setw(mode_2_width) << "cp-excl";
          if (i>=2*cost_model_size) Stream << std::left << std::setw(mode_2_width) << "pp-excl (s)";
          else Stream << std::left << std::setw(mode_2_width) << "pp-excl";
          Stream << std::left << std::setw(mode_2_width) << "vol-excl (s)";
          Stream << std::left << std::setw(mode_2_width) << "cp-excl (%)";
          Stream << std::left << std::setw(mode_2_width) << "pp-excl (%)";
          Stream << std::left << std::setw(mode_2_width) << "vol-excl (%)";
          double cp_total_exclusive = 0.;
          double pp_total_exclusive = 0.;
          double vol_total_exclusive = 0.;
          for (auto& it : sort_info){
            Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
            Stream << std::left << std::setw(mode_2_width) << it.second[0];
            Stream << std::left << std::setw(mode_2_width) << it.second[2];
            Stream << std::left << std::setw(mode_2_width) << it.second[4];
            Stream << std::left << std::setw(mode_2_width) << it.second[1];
            Stream << std::left << std::setw(mode_2_width) << it.second[3];
            Stream << std::left << std::setw(mode_2_width) << it.second[5];
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_exclusive += it.second[1];
            pp_total_exclusive += it.second[3];
            vol_total_exclusive += it.second[5];
          }
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << cp_total_exclusive;
          Stream << std::left << std::setw(mode_2_width) << pp_total_exclusive;
          Stream << std::left << std::setw(mode_2_width) << vol_total_exclusive;
          Stream << std::left << std::setw(mode_2_width) << 100.*cp_total_exclusive/cp_ref;
          Stream << std::left << std::setw(mode_2_width) << 100.*pp_total_exclusive/pp_ref;
          Stream << std::left << std::setw(mode_2_width) << 100.*vol_total_exclusive/vol_ref;
          Stream << "\n";

          // Reset symbol timers and sort
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          j=0;
          for (auto& it : symbol_timers){
            assert(it.second.start_timer.size() == 0);
            sort_info[j++] = std::make_pair(it.second.name,std::array<double,6>{it.second.cp_numcalls[z][0],it.second.cp_incl_measure[z][i],*it.second.pp_numcalls,it.second.pp_incl_measure[i],*it.second.vol_numcalls/world_size,it.second.vol_incl_measure[i]/world_size});
          }
          std::sort(sort_info.begin(),sort_info.end(),[](std::pair<std::string,std::array<double,6>>& vec1, std::pair<std::string,std::array<double,6>>& vec2){return vec1.second[1] > vec2.second[1];});
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
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << get_measure_title(i);
          Stream << std::left << std::setw(mode_2_width) << "cp-#calls";
          Stream << std::left << std::setw(mode_2_width) << "pp-#calls";
          Stream << std::left << std::setw(mode_2_width) << "vol-#calls";
          if (i>=2*cost_model_size) Stream << std::left << std::setw(mode_2_width) << "cp-incl (s)";
          else Stream << std::left << std::setw(mode_2_width) << "cp-incl";
          if (i>=2*cost_model_size) Stream << std::left << std::setw(mode_2_width) << "pp-incl (s)";
          else Stream << std::left << std::setw(mode_2_width) << "pp-incl";
          Stream << std::left << std::setw(mode_2_width) << "vol-incl (s)";
          Stream << std::left << std::setw(mode_2_width) << "cp-incl (%)";
          Stream << std::left << std::setw(mode_2_width) << "pp-incl (%)";
          Stream << std::left << std::setw(mode_2_width) << "vol-incl (%)";
          double cp_total_inclusive = 0.;
          double pp_total_inclusive = 0.;
          double vol_total_inclusive = 0.;
          for (auto& it : sort_info){
            Stream << "\n" << std::left << std::setw(max_timer_name_length) << it.first;
            Stream << std::left << std::setw(mode_2_width) << it.second[0];
            Stream << std::left << std::setw(mode_2_width) << it.second[2];
            Stream << std::left << std::setw(mode_2_width) << it.second[4];
            Stream << std::left << std::setw(mode_2_width) << it.second[1];
            Stream << std::left << std::setw(mode_2_width) << it.second[3];
            Stream << std::left << std::setw(mode_2_width) << it.second[5];
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            Stream << std::left << std::setw(mode_2_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_inclusive = std::max(it.second[1],cp_total_inclusive);
            pp_total_inclusive = std::max(it.second[3],pp_total_inclusive);
            vol_total_inclusive = std::max(it.second[5],vol_total_inclusive);
          }
          Stream << "\n" << std::left << std::setw(max_timer_name_length) << "total";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << "";
          Stream << std::left << std::setw(mode_2_width) << cp_total_inclusive;
          Stream << std::left << std::setw(mode_2_width) << pp_total_inclusive;
          Stream << std::left << std::setw(mode_2_width) << vol_total_inclusive;
          Stream << std::left << std::setw(mode_2_width) << 100.*cp_total_inclusive/cp_ref;
          Stream << std::left << std::setw(mode_2_width) << 100.*pp_total_inclusive/pp_ref;
          Stream << std::left << std::setw(mode_2_width) << 100.*vol_total_inclusive/vol_ref;
          Stream << "\n";
        }
      }
    }
  }
}

}
}
