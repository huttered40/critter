#include "symbol_tracker.h"
#include "../util/util.h"

namespace critter{
namespace internal{
namespace profile{

std::unordered_map<std::string,symbol_tracker> symbol_timers;

symbol_tracker::symbol_tracker(std::string name_){
  if (path_count==0) return;
  this->name = std::move(name_);
  this->cp_exclusive_contributions.resize(path_count,nullptr);
  this->cp_exclusive_measure.resize(path_count,nullptr);
  this->cp_numcalls.resize(path_count,nullptr);
  this->cp_incl_measure.resize(path_count,nullptr);
  this->cp_excl_measure.resize(path_count,nullptr);

  //TODO: Can likely reduce from magic number of '4' to '3'
  size_t cp_path_select_offset = 4*num_decomp_cp_measures+1;
  size_t pp_path_select_offset = 4*num_decomp_pp_measures+1;
  size_t vol_path_select_offset = 4*num_decomp_vol_measures+1;
  size_t cp_offset = num_cp_measures;
  cp_offset += (symbol_timers.size()-1)*path_count*cp_path_select_offset;
  size_t pp_offset = num_pp_measures;
  pp_offset += (symbol_timers.size()-1)*path_count*pp_path_select_offset;
  size_t vol_offset = num_vol_measures;
  vol_offset += (symbol_timers.size()-1)*vol_path_select_offset;

  this->pp_numcalls = &max_pp_costs[pp_offset];
  this->pp_incl_measure = &max_pp_costs[pp_offset+1];
  this->pp_excl_measure = &max_pp_costs[pp_offset+num_decomp_pp_measures+1];
  this->pp_exclusive_contributions = &max_pp_costs[pp_offset+2*num_decomp_pp_measures+1];
  this->pp_exclusive_measure = &max_pp_costs[pp_offset+3*num_decomp_pp_measures+1];

  this->vol_numcalls = &vol_costs[vol_offset];
  this->vol_incl_measure = &vol_costs[vol_offset+1];
  this->vol_excl_measure = &vol_costs[vol_offset+num_decomp_vol_measures+1];

  for (auto i=0; i<path_count; i++){
    this->cp_numcalls[i] = &cp_costs[cp_offset+cp_path_select_offset*i];
    this->cp_incl_measure[i] = &cp_costs[cp_offset+1+cp_path_select_offset*i];
    this->cp_excl_measure[i] = &cp_costs[cp_offset+num_decomp_cp_measures+1+cp_path_select_offset*i];
    this->cp_exclusive_contributions[i] = &cp_costs[cp_offset+2*num_decomp_cp_measures+1+cp_path_select_offset*i];
    this->cp_exclusive_measure[i] = &cp_costs[cp_offset+3*num_decomp_cp_measures+1+cp_path_select_offset*i];
  }
  // memset not necessary because these members point into global arrays
}

void symbol_tracker::start(double save_time){
  if (path_count==0) return;
  if (symbol_stack.size()>0){
    // Add up exclusive time for most-recent symbol (the one about to be supplanted for this new one).
    auto last_symbol_time = save_time-symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<path_count; i++){
      for (auto j=0; j<path_measure_index.size(); j++){
        if (path_measure_index[j] == path_measure_select.size()-1){
          symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][j] += last_symbol_time;
          symbol_timers[symbol_stack.top()].cp_excl_measure[i][j] += last_symbol_time;
        }
      }
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  }
  auto last_comp_time = save_time - computation_timer;
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  symbol_stack.push(this->name);
  computation_timer = MPI_Wtime();
  this->start_timer.push((double)computation_timer);
}

void symbol_tracker::stop(double save_time){
  if (path_count==0) return;
  assert(this->start_timer.size()>0);
  auto last_symbol_time = save_time-this->start_timer.top();
  // Add up exclusive time for most-recent symbol (the one about to be removed from top-of-stack).
  for (auto i=0; i<path_count; i++){
    for (auto j=0; j<path_measure_index.size(); j++){
      if (path_measure_index[j] == path_measure_select.size()-1){
        this->cp_exclusive_measure[i][j] += last_symbol_time;
        this->cp_excl_measure[i][j] += last_symbol_time;
      }
    }
  }
  this->pp_exclusive_measure[num_decomp_pp_measures-1] += last_symbol_time;
  this->pp_excl_measure[num_decomp_pp_measures-1] += last_symbol_time;
  *this->pp_numcalls += 1.; *this->vol_numcalls += 1.;
  for (auto i=0; i<num_decomp_pp_measures; i++){
    this->pp_exclusive_contributions[i] += this->pp_exclusive_measure[i];
  }
  for (auto j=0; j<path_count; j++){
    *this->cp_numcalls[j] += 1.;
    for (auto i=0; i<num_decomp_cp_measures; i++){
      this->cp_exclusive_contributions[j][i] += this->cp_exclusive_measure[j][i];
    }
  }
  for (auto j=0; j<path_count; j++){
    memset(this->cp_exclusive_measure[j],0,sizeof(float)*num_decomp_cp_measures);
  }
  memset(this->pp_exclusive_measure,0,sizeof(float)*num_decomp_pp_measures);
  auto save_symbol = symbol_stack.top();
  this->start_timer.pop(); symbol_stack.pop();
  if (symbol_stack.size() > 0 && (save_symbol != symbol_stack.top())){
    for (auto j=0; j<path_count; j++){
      for (auto i=0; i<num_decomp_cp_measures; i++){
        symbol_timers[symbol_stack.top()].cp_exclusive_contributions[j][i] += this->cp_exclusive_contributions[j][i];
        this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
      }
    }
    for (auto i=0; i<num_decomp_pp_measures; i++){
      symbol_timers[symbol_stack.top()].pp_exclusive_contributions[i] += this->pp_exclusive_contributions[i];
      this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
    }
    for (auto j=0; j<path_count; j++){
      memset(this->cp_exclusive_contributions[j],0,sizeof(float)*num_decomp_cp_measures);
    }
    memset(this->pp_exclusive_contributions,0,sizeof(float)*num_decomp_pp_measures);
  } else if (symbol_stack.size() == 0){
    for (auto j=0; j<path_count; j++){
      for (auto i=0; i<num_decomp_cp_measures; i++){
        this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
      }
    }
    for (auto i=0; i<num_decomp_pp_measures; i++){
      this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
    }
    for (auto j=0; j<path_count; j++){
      memset(this->cp_exclusive_contributions[j],0,sizeof(float)*num_decomp_cp_measures);
    }
    memset(this->pp_exclusive_contributions,0,sizeof(float)*num_decomp_pp_measures);
  }
  auto last_comp_time = save_time - computation_timer;
  cp_costs[num_cp_measures-1] += last_comp_time;
  vol_costs[num_vol_measures-1] += last_comp_time;
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

}
}
}
