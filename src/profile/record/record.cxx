#include "record.h"
#include "../util/util.h"
#include "../container/comm_tracker.h"
#include "../container/symbol_tracker.h"

namespace critter{
namespace internal{
namespace profile{

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
  std::vector<std::string> measure_titles(num_pp_measures);
  measure_titles[0] = "comm_cost";
  measure_titles[1] = "synch_cost";
  measure_titles[2] = "comp_cost";
  measure_titles[3] = "comm time";
  measure_titles[4] = "exec time";
  return measure_titles[idx];
}

void print_cost_model_tracker_header(){
  std::cout << std::left << std::setw(text_width) << "CommCost";
  std::cout << std::left << std::setw(text_width) << "SynchCost";
  std::cout << std::left << std::setw(text_width) << "CommTime";
}

void print_cost_model_header(){
  std::cout << std::left << std::setw(text_width) << "CommCost";
  std::cout << std::left << std::setw(text_width) << "SynchCost";
  std::cout << std::left << std::setw(text_width) << "CompCost";
  std::cout << std::left << std::setw(text_width) << "CommTime";
  std::cout << std::left << std::setw(text_width) << "ExecTime";
  std::cout << "\n";
  std::cout << std::left << std::setw(text_width) << "      ";
}

void print_header(size_t num_inputs, std::map<std::string,std::vector<float>>& save_info){
  for (size_t idx = 0; idx < (num_inputs+1); idx++){
    if (idx != 0){
      stream << "\t";
    }
    stream << "Input";
  }
  stream << "\tCommCost\tSynchCost\tCompCost\tCommunicationTime\tExecTime";// critical path
  stream << "\tCommCost\tSynchCost\tCompCost\tCommunicationTime\tExecTime";// per-process
  stream << "\tCommCost\tSynchCost\tCompCost\tCommunicationTime\tExecTime";// vol
  for (auto i=0; i<num_decomp_cp_measures*path_count+num_decomp_pp_measures*path_count+num_decomp_vol_measures;i++){
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
    PMPI_Allreduce(MPI_IN_PLACE,&_wall_time,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
    if (is_world_root){
      if (is_first_iter){
        stream << std::left << std::setw(text_width) << "ID";
        stream << std::left << std::setw(text_width) << "Execution time";
        stream << std::left << std::setw(text_width) << "Wall time";
        stream << "\n";
      }
      stream << std::left << std::setw(text_width) << variantID;
      stream << std::left << std::setw(text_width) << cp_costs[num_cp_measures-1];
      stream << std::left << std::setw(text_width) << _wall_time;
      stream << "\n";
    }
  }
  else if (temp_mode){
    auto np=0; MPI_Comm_size(MPI_COMM_WORLD,&np);
    if (is_world_root){
      std::map<std::string,std::vector<float>> save_info;
      auto inputs = parse_file_string();
      for (int i=0; i<list_size; i++){
        list[i]->set_header(save_info);
      }
      if (is_first_iter){
        if (variantID != -1) stream << "ID\t";
        print_header(inputs.size(),save_info);
        stream << "\n";
      }
      save_info.clear();
      if (variantID != -1) stream << variantID << "\t";
      print_inputs(np,inputs);
      for (size_t i=0; i<num_cp_measures; i++){
        stream << "\t" << cp_costs[i];
      }
      for (size_t i=0; i<num_pp_measures; i++){
        stream << "\t" << max_pp_costs[i];
      }
      for (size_t i=0; i<num_vol_measures; i++){
        stream << "\t" << vol_costs[i];
      }
      for (int i=0; i<list_size; i++){
        list[i]->set_vol_costs(save_info);
      }
      for (size_t j=0; j<num_decomp_vol_measures; j++){
        for (auto& it : save_info){
          stream << "\t" << it.second[j];
        }
      }
      save_info.clear();
/*
      size_t breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){	// no idle time
        if (comm_path_select[i]=='0') continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_pp_costs(save_info,breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_pp_measures; j++){
          for (auto& it : save_info){
            stream << "\t" << it.second[j];
          }
        }
        stream << "\t" << max_pp_costs[num_pp_measures+(breakdown_idx+1)*(num_tracker_pp_measures*list_size+4)-1];// comp kernel time
        stream << "\t" << max_pp_costs[num_pp_measures+(breakdown_idx+1)*(num_tracker_pp_measures*list_size+4)-2];// comp time
        stream << "\t" << max_pp_costs[num_pp_measures+(breakdown_idx+1)*(num_tracker_pp_measures*list_size+4)-3];// idle time
        stream << "\t" << max_pp_costs[num_pp_measures+(breakdown_idx+1)*(num_tracker_pp_measures*list_size+4)-4];// comp cost
        breakdown_idx++;
      }
      save_info.clear();
      breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        stream << "\t" << cp_costs[cp_costs_size-comm_path_select_size+breakdown_idx];// comp kernel time
        stream << "\t" << cp_costs[cp_costs_size-2*comm_path_select_size+breakdown_idx];// comp time
        stream << "\t" << cp_costs[cp_costs_size-3*comm_path_select_size+breakdown_idx];// idle time
        stream << "\t" << cp_costs[cp_costs_size-4*comm_path_select_size+breakdown_idx];// comp cost
        breakdown_idx++;
      }
      breakdown_idx=0;
      for (auto i=0; i<comm_path_select.size(); i++){
        if (comm_path_select[i]=='0') continue;
        // Save the critter information before printing
        for (size_t j=0; j<list_size; j++){
          list[j]->set_cp_costs(save_info,breakdown_idx);
        }
        for (size_t j=0; j<num_tracker_cp_measures; j++){
          for (auto& it : save_info){
            stream << "\t" << it.second[j];
          }
        }
        breakdown_idx++;
      }
*/
      stream << "\n";
    }
  }
  is_first_iter = false;// set here only beause this routine is called directly after 'invoke' on std::ostream
}

void record::print(int variantID, float overhead_time){
  //int world_size; MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int temp_mode=1;
  if (std::getenv("CRITTER_MODE") != NULL){
    temp_mode = atoi(std::getenv("CRITTER_MODE"));
  }
  if (temp_mode==0){
    if (is_world_root){
      std::cout << std::left << std::setw(text_width) << "Execution time:";
      std::cout << "\n";
      std::cout << std::left << std::setw(text_width) << "              ";
      std::cout << std::left << std::setw(text_width) << cp_costs[num_cp_measures-1];
      std::cout << "\n\n";
    }
  }
  if (temp_mode){
    if (is_world_root){
      std::map<std::string,std::vector<float>> save_info;
      std::cout << "\n\n";
      std::cout << std::left << std::setw(text_width) << "Critical path:";
      print_cost_model_header();
      for (size_t i=0; i<num_cp_measures; i++){
        std::cout << std::left << std::setw(text_width) << cp_costs[i];
      }
      std::cout << "\n\n";

      std::cout << std::left << std::setw(text_width) << "Process max:";
      print_cost_model_header();
      for (size_t i=0; i<num_vol_measures; i++){
        std::cout << std::left << std::setw(text_width) << max_pp_costs[i];
      }
      std::cout << "\n\n";

      std::cout << std::left << std::setw(text_width) << "Volume avg:";
      print_cost_model_header();
      for (size_t i=0; i<num_vol_measures; i++){
        std::cout << std::left << std::setw(text_width) << vol_costs[i];
      }
      std::cout << "\n\n";

      if (path_decomposition == 1){
      for (auto i=0; i<path_index.size(); i++){
        if (path_index[i]==0){
          std::cout << std::left << std::setw(decomp_text_width+1) << "CommCost:";
        } else if (path_index[i]==1){
          std::cout << std::left << std::setw(decomp_text_width+1) << "SynchCost:";
        } else if (path_index[i]==2){
          std::cout << std::left << std::setw(decomp_text_width+1) << "CompCost:";
        } else if (path_index[i]==3){
          std::cout << std::left << std::setw(decomp_text_width+1) << "CommTime:";
        } else if (path_index[i]==4){
          std::cout << std::left << std::setw(decomp_text_width+1) << "ExecTime:";
        }
        std::cout << std::left << std::setw(text_width) << "MeasureType";
        std::cout << std::left << std::setw(text_width) << "CompCost";
        std::cout << std::left << std::setw(text_width) << "CommCost";
        std::cout << std::left << std::setw(text_width) << "SynchCost";
        std::cout << std::left << std::setw(text_width) << "CommTime";
        std::cout << "\n";
        std::cout << std::left << std::setw(decomp_text_width+1) << "CompCost";
        std::cout << std::left << std::setw(text_width) << "path";
        std::cout << std::left << std::setw(text_width) << cp_costs[cp_costs_size-1*path_count+i];//CompCost
        save_info.clear();
        for (int j=0; j<list_size; j++){
          list[j]->set_cp_costs(save_info,i);
        }
        for (auto& it : save_info){
          std::cout << "\n";
          std::cout << std::left << std::setw(decomp_text_width+1) << it.first;
          std::cout << std::left << std::setw(text_width) << "path";
          std::cout << std::left << std::setw(text_width) << 0.0;
          // Print out {CommCost,SynchCost,CommTime}
          for (size_t j=0; j<3; j++){
            std::cout << std::left << std::setw(text_width) << it.second[j];
          }
        }
        std::cout << "\n";
        std::cout << std::left << std::setw(decomp_text_width+1) << "CompCost";
        std::cout << std::left << std::setw(text_width) << "per-process";
        std::cout << std::left << std::setw(text_width) << max_pp_costs[num_pp_measures+(i+1)*(num_decomp_pp_measures*list_size+2)-2];
        save_info.clear();
        for (int j=0; j<list_size; j++){
          list[j]->set_pp_costs(save_info,i);
        }
        for (auto& it : save_info){
          std::cout << "\n";
          std::cout << std::left << std::setw(decomp_text_width+1) << it.first;
          std::cout << std::left << std::setw(text_width) << "per-process";
          std::cout << std::left << std::setw(text_width) << 0.0;
          std::cout << std::left << std::setw(text_width) << it.second[0];
          std::cout << std::left << std::setw(text_width) << it.second[1];
          std::cout << std::left << std::setw(text_width) << it.second[2];
        }
        std::cout << "\n\n";
      }
      }
      if (path_decomposition <= 1){
        save_info.clear();
        for (int i=0; i<list_size; i++){
          list[i]->set_vol_costs(save_info);
        }
        std::cout << std::left << std::setw(decomp_text_width+1) << "Volume avg:";
        print_cost_model_tracker_header();
        for (auto& it : save_info){
          std::cout << "\n";
          std::cout << std::left << std::setw(decomp_text_width+1) << it.first;
          for (size_t j=0; j<num_decomp_vol_measures; j++){
            std::cout << std::left << std::setw(text_width) << it.second[j];
          }
        }
        std::cout << "\n";
      }
    }
  }
  if ((temp_mode) && (path_decomposition==2) && (path_count>0) && (symbol_timers.size()>0)){
    if (is_world_root){
      for (auto z=0; z<path_count; z++){
        std::cout << "***********************************************************************************************************************";
        std::vector<std::pair<std::string,std::array<float,6>>> sort_info(symbol_timers.size());
        for (int i=path_measure_index.size()-1; i>=0; i--){
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          // Reset symbol timers and sort
          size_t j=0;
          float cp_ref,pp_ref,vol_ref;
          for (auto& it : symbol_timers){
            if (it.second.start_timer.size() != 0) { std::cout << "Symbol " << it.first << " is not handled properly\n"; assert(it.second.start_timer.size() == 0); }
              sort_info[j++] = std::make_pair(it.second.name,
                                              std::array<float,6>{it.second.cp_numcalls[z][0],
                                                                  it.second.cp_excl_measure[z][i],	//TODO: why not replace with cp_exclusive_measures
                                                                  *it.second.pp_numcalls,
                                                                  it.second.pp_excl_measure[path_measure_index[i]],
                                                                  *it.second.vol_numcalls,
                                                                  it.second.vol_excl_measure[path_measure_index[i]]});
          }
          std::sort(sort_info.begin(),sort_info.end(),
                    [](std::pair<std::string,std::array<float,6>>& vec1,
                       std::pair<std::string,std::array<float,6>>& vec2){
                    return vec1.second[1] > vec2.second[1];});
          // These conditionals due to annoying fact that we don't consider CompKernelTime as among the path measures (so 6 total, not 7).
          cp_ref = cp_costs[path_measure_index[i]];
          pp_ref = max_pp_costs[path_measure_index[i]];
          vol_ref = vol_costs[path_measure_index[i]];

          // Exclusive
          std::cout << "\n\n\n\n" << std::left << std::setw(decomp_text_width+1) << get_measure_title(path_measure_index[i]);
          std::cout << std::left << std::setw(text_width) << "cp-#calls";
          std::cout << std::left << std::setw(text_width) << "pp-#calls";
          std::cout << std::left << std::setw(text_width) << "vol-#calls";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "cp-excl (s)";
          else std::cout << std::left << std::setw(text_width) << "cp-excl";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "pp-excl (s)";
          else std::cout << std::left << std::setw(text_width) << "pp-excl";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "vol-excl (s)";
          else std::cout << std::left << std::setw(text_width) << "vol-excl";
          std::cout << std::left << std::setw(text_width) << "cp-excl (%)";
          std::cout << std::left << std::setw(text_width) << "pp-excl (%)";
          std::cout << std::left << std::setw(text_width) << "vol-excl (%)";
          float cp_total_exclusive = 0.;
          float pp_total_exclusive = 0.;
          float vol_total_exclusive = 0.;
          float cp_total_invocations_exclusive = 0.;
          float pp_total_invocations_exclusive = 0.;
          float vol_total_invocations_exclusive = 0.;
          for (auto& it : sort_info){
            std::cout << "\n" << std::left << std::setw(decomp_text_width+1) << it.first;
            std::cout << std::left << std::setw(text_width) << it.second[0];
            std::cout << std::left << std::setw(text_width) << it.second[2];
            std::cout << std::left << std::setw(text_width) << it.second[4];
            std::cout << std::left << std::setw(text_width) << it.second[1];
            std::cout << std::left << std::setw(text_width) << it.second[3];
            std::cout << std::left << std::setw(text_width) << it.second[5];
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_exclusive += it.second[1];
            pp_total_exclusive += it.second[3];
            vol_total_exclusive += it.second[5];
            cp_total_invocations_exclusive += it.second[0];
            pp_total_invocations_exclusive += it.second[2];
            vol_total_invocations_exclusive += it.second[4];
          }
          std::cout << "\n" << std::left << std::setw(decomp_text_width+1) << "total";
          std::cout << std::left << std::setw(text_width) << cp_total_invocations_exclusive;
          std::cout << std::left << std::setw(text_width) << pp_total_invocations_exclusive;
          std::cout << std::left << std::setw(text_width) << vol_total_invocations_exclusive;
          std::cout << std::left << std::setw(text_width) << cp_total_exclusive;
          std::cout << std::left << std::setw(text_width) << pp_total_exclusive;
          std::cout << std::left << std::setw(text_width) << vol_total_exclusive;
          std::cout << std::left << std::setw(text_width) << 100.*cp_total_exclusive/cp_ref;
          std::cout << std::left << std::setw(text_width) << 100.*pp_total_exclusive/pp_ref;
          std::cout << std::left << std::setw(text_width) << 100.*vol_total_exclusive/vol_ref;
          std::cout << "\n";

          // Reset symbol timers and sort
          sort_info.clear(); sort_info.resize(symbol_timers.size());
          j=0;
          for (auto& it : symbol_timers){
            assert(it.second.start_timer.size() == 0);
            sort_info[j++] = std::make_pair(it.second.name,
                                            std::array<float,6>{it.second.cp_numcalls[z][0],
                                                                it.second.cp_incl_measure[z][i],
                                                                *it.second.pp_numcalls,
                                                                it.second.pp_incl_measure[path_measure_index[i]],
                                                                *it.second.vol_numcalls,
                                                                it.second.vol_incl_measure[path_measure_index[i]]});
          }
          std::sort(sort_info.begin(),sort_info.end(),
                    [](std::pair<std::string,std::array<float,6>>& vec1,
                       std::pair<std::string,std::array<float,6>>& vec2){
                    return vec1.second[1] > vec2.second[1];});
          cp_ref = cp_costs[path_measure_index[i]];
          pp_ref = max_pp_costs[path_measure_index[i]];
          vol_ref = vol_costs[path_measure_index[i]];

          // Inclusive
          std::cout << "\n" << std::left << std::setw(decomp_text_width+1) << get_measure_title(path_measure_index[i]);
          std::cout << std::left << std::setw(text_width) << "cp-#calls";
          std::cout << std::left << std::setw(text_width) << "pp-#calls";
          std::cout << std::left << std::setw(text_width) << "vol-#calls";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "cp-incl (s)";
          else std::cout << std::left << std::setw(text_width) << "cp-incl";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "pp-incl (s)";
          else std::cout << std::left << std::setw(text_width) << "pp-incl";
          if (i>=3) std::cout << std::left << std::setw(text_width) << "vol-incl (s)";
          else std::cout << std::left << std::setw(text_width) << "vol-incl";
          std::cout << std::left << std::setw(text_width) << "cp-incl (%)";
          std::cout << std::left << std::setw(text_width) << "pp-incl (%)";
          std::cout << std::left << std::setw(text_width) << "vol-incl (%)";
          float cp_total_inclusive = 0.;
          float pp_total_inclusive = 0.;
          float vol_total_inclusive = 0.;
          for (auto& it : sort_info){
            std::cout << "\n" << std::left << std::setw(decomp_text_width+1) << it.first;
            std::cout << std::left << std::setw(text_width) << it.second[0];
            std::cout << std::left << std::setw(text_width) << it.second[2];
            std::cout << std::left << std::setw(text_width) << it.second[4];
            std::cout << std::left << std::setw(text_width) << it.second[1];
            std::cout << std::left << std::setw(text_width) << it.second[3];
            std::cout << std::left << std::setw(text_width) << it.second[5];
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[1]/cp_ref;
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[3]/pp_ref;
            std::cout << std::left << std::setw(text_width) << std::setprecision(4) << 100.*it.second[5]/vol_ref;
            cp_total_inclusive = std::max(it.second[1],cp_total_inclusive);
            pp_total_inclusive = std::max(it.second[3],pp_total_inclusive);
            vol_total_inclusive = std::max(it.second[5],vol_total_inclusive);
          }
          std::cout << "\n" << std::left << std::setw(decomp_text_width+1) << "total";
          std::cout << std::left << std::setw(text_width) << "";
          std::cout << std::left << std::setw(text_width) << "";
          std::cout << std::left << std::setw(text_width) << "";
          std::cout << std::left << std::setw(text_width) << cp_total_inclusive;
          std::cout << std::left << std::setw(text_width) << pp_total_inclusive;
          std::cout << std::left << std::setw(text_width) << vol_total_inclusive;
          std::cout << std::left << std::setw(text_width) << 100.*cp_total_inclusive/cp_ref;
          std::cout << std::left << std::setw(text_width) << 100.*pp_total_inclusive/pp_ref;
          std::cout << std::left << std::setw(text_width) << 100.*vol_total_inclusive/vol_ref;
          std::cout << "\n";
        }
      }
    }
  }
}

}
}
}
