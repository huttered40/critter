#include "symbol_tracker.h"
#include "../util.h"

namespace critter{
namespace internal{

std::unordered_map<std::string,symbol_tracker> symbol_timers;

symbol_tracker::symbol_tracker(std::string name_){
  assert(name_.size() <= max_timer_name_length);
  assert(symbol_timers.size() < max_num_symbols);

  this->cp_exclusive_contributions.resize(num_critical_path_measures,0);
  this->pp_exclusive_contributions.resize(num_per_process_measures,0);
  this->cp_exclusive_measure.resize(num_critical_path_measures,0);
  this->pp_exclusive_measure.resize(num_per_process_measures,0);
 
  this->cp_incl_measure.resize(num_critical_path_measures);
  this->cp_excl_measure.resize(num_critical_path_measures);
  this->pp_incl_measure.resize(num_per_process_measures);
  this->pp_excl_measure.resize(num_per_process_measures);
  this->vol_incl_measure.resize(num_volume_measures);
  this->vol_excl_measure.resize(num_volume_measures);
 
  this->name = std::move(name_);
  this->cp_numcalls = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)]; *this->cp_numcalls = 0;
  this->pp_numcalls = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)]; *this->pp_numcalls = 0;
  this->vol_numcalls = &symbol_timer_pad_local_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)]; *this->vol_numcalls = 0;
  for (auto i=0; i<num_critical_path_measures; i++){
    this->cp_incl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*i+1]; *cp_incl_measure[i] = 0.;
    this->cp_excl_measure[i] = &symbol_timer_pad_local_cp[(symbol_timers.size()-1)*(num_ftimer_measures*num_critical_path_measures+1)+2*(i+1)]; *cp_excl_measure[i] = 0.;
  }
  for (auto i=0; i<num_per_process_measures; i++){
    this->pp_incl_measure[i] = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)+2*i+1]; *pp_incl_measure[i] = 0.;
    this->pp_excl_measure[i] = &symbol_timer_pad_local_pp[(symbol_timers.size()-1)*(num_ftimer_measures*num_per_process_measures+1)+2*(i+1)]; *pp_excl_measure[i] = 0.;
  }
  for (auto i=0; i<num_volume_measures; i++){
    this->vol_incl_measure[i] = &symbol_timer_pad_local_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*i+1]; *vol_incl_measure[i] = 0.;
    this->vol_excl_measure[i] = &symbol_timer_pad_local_vol[(symbol_timers.size()-1)*(num_ftimer_measures*num_volume_measures+1)+2*(i+1)]; *vol_excl_measure[i] = 0.;
  }
  this->has_been_processed = false;
}

void symbol_tracker::start(double save_time){
  if (symbol_stack.size()>0){
    auto last_symbol_time = save_time-symbol_timers[symbol_stack.top()].start_timer.top();
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
    *symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  }
  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer); }
  symbol_stack.push(this->name);
  computation_timer = MPI_Wtime();
  this->start_timer.push(computation_timer);
}

void symbol_tracker::stop(double save_time){
  assert(this->start_timer.size()>0);
  auto last_symbol_time = save_time-this->start_timer.top();
  this->cp_exclusive_measure[num_critical_path_measures-1] += last_symbol_time;
  this->cp_exclusive_measure[num_critical_path_measures-2] += last_symbol_time;
  this->pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
  this->pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
  *this->cp_excl_measure[num_critical_path_measures-1] += last_symbol_time;
  *this->cp_excl_measure[num_critical_path_measures-2] += last_symbol_time;
  *this->pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  *this->pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  *this->cp_numcalls += 1.; *this->pp_numcalls += 1.; *this->vol_numcalls += 1.;

  for (auto i=0; i<num_critical_path_measures; i++){ this->cp_exclusive_contributions[i] += this->cp_exclusive_measure[i]; }
  for (auto i=0; i<this->cp_exclusive_measure.size(); i++){ this->cp_exclusive_measure[i]=0.; }
  for (auto i=0; i<num_per_process_measures; i++){ this->pp_exclusive_contributions[i] += this->pp_exclusive_measure[i]; }
  for (auto i=0; i<this->pp_exclusive_measure.size(); i++){ this->pp_exclusive_measure[i]=0.; }

  auto save_symbol = symbol_stack.top();
  this->start_timer.pop(); symbol_stack.pop();
  if (symbol_stack.size() > 0 && (save_symbol != symbol_stack.top())){
    for (auto i=0; i<num_critical_path_measures; i++){ symbol_timers[symbol_stack.top()].cp_exclusive_contributions[i] += this->cp_exclusive_contributions[i];
                                                       *this->cp_incl_measure[i] += this->cp_exclusive_contributions[i]; }
    for (auto i=0; i<num_per_process_measures; i++){ symbol_timers[symbol_stack.top()].pp_exclusive_contributions[i] += this->pp_exclusive_contributions[i];
                                                     *this->pp_incl_measure[i] += this->pp_exclusive_contributions[i]; }
    for (auto i=0; i<this->cp_exclusive_contributions.size(); i++){ this->cp_exclusive_contributions[i]=0.; }
    for (auto i=0; i<this->pp_exclusive_contributions.size(); i++){ this->pp_exclusive_contributions[i]=0.; }
  } else if (symbol_stack.size() == 0){
    for (auto i=0; i<num_critical_path_measures; i++){ *this->cp_incl_measure[i] += this->cp_exclusive_contributions[i]; }
    for (auto i=0; i<num_per_process_measures; i++){ *this->pp_incl_measure[i] += this->pp_exclusive_contributions[i]; }
    for (auto i=0; i<this->cp_exclusive_contributions.size(); i++){ this->cp_exclusive_contributions[i]=0.; }
    for (auto i=0; i<this->pp_exclusive_contributions.size(); i++){ this->pp_exclusive_contributions[i]=0.; }
  }

  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<breakdown_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer); }
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

}
}
