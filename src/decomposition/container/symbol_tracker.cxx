#include "symbol_tracker.h"

namespace critter{
namespace internal{
namespace decomposition{

// Global namespace variable 'symbol_timers' must be defined here, rather than in src/util.cxx with the rest, to avoid circular dependence between this file and src/util.h
std::unordered_map<std::string,symbol_tracker> symbol_timers;

symbol_tracker::symbol_tracker(std::string name_){
  assert(name_.size() <= max_timer_name_length);
  assert(symbol_timers.size() < max_num_symbols);
  this->name = std::move(name_);
  this->cp_exclusive_contributions.resize(symbol_path_select_size,nullptr);
  this->cp_exclusive_measure.resize(symbol_path_select_size,nullptr);
  this->cp_numcalls.resize(symbol_path_select_size,nullptr);
  this->cp_incl_measure.resize(symbol_path_select_size,nullptr);
  this->cp_excl_measure.resize(symbol_path_select_size,nullptr);

  // Note that we use 'num_per_process_measures' instead of 'num_critical_path_measures' because we want to record the idle time along a path.
  //   A path being decomposed is not necessarily the critical path, thus idle time is possible. We don't set 'num_critical_path_measures'=='num_per_process_measures'
  //   because the latter specifies the number of global critical path metrics, none of which include idle time (wouldn't make sense).
  size_t cp_offset = (symbol_timers.size()-1)*symbol_path_select_size*(cp_symbol_class_count*num_per_process_measures+1);
  size_t pp_offset = (symbol_timers.size()-1)*(pp_symbol_class_count*num_per_process_measures+1);
  size_t vol_offset = (symbol_timers.size()-1)*(vol_symbol_class_count*num_volume_measures+1);

  this->pp_numcalls = &symbol_timer_pad_local_pp[pp_offset];
  this->pp_incl_measure = &symbol_timer_pad_local_pp[pp_offset+1];
  this->pp_excl_measure = &symbol_timer_pad_local_pp[pp_offset+num_per_process_measures+1];
  this->pp_exclusive_contributions = &symbol_timer_pad_local_pp[pp_offset+2*num_per_process_measures+1];
  this->pp_exclusive_measure = &symbol_timer_pad_local_pp[pp_offset+3*num_per_process_measures+1];

  this->vol_numcalls = &symbol_timer_pad_local_vol[vol_offset];
  this->vol_incl_measure = &symbol_timer_pad_local_vol[vol_offset+1];
  this->vol_excl_measure = &symbol_timer_pad_local_vol[vol_offset+num_volume_measures+1];

  size_t path_select_offset = cp_symbol_class_count*num_per_process_measures+1;// amount of memory allocated to each symbol, same for cp/pp/vol
  for (auto i=0; i<symbol_path_select_size; i++){
    this->cp_numcalls[i] = &symbol_timer_pad_local_cp[cp_offset+path_select_offset*i];
    this->cp_incl_measure[i] = &symbol_timer_pad_local_cp[cp_offset+1+path_select_offset*i];
    this->cp_excl_measure[i] = &symbol_timer_pad_local_cp[cp_offset+num_per_process_measures+1+path_select_offset*i];
    this->cp_exclusive_contributions[i] = &symbol_timer_pad_local_cp[cp_offset+2*num_per_process_measures+1+path_select_offset*i];
    this->cp_exclusive_measure[i] = &symbol_timer_pad_local_cp[cp_offset+3*num_per_process_measures+1+path_select_offset*i];
  }
  memset(&symbol_timer_pad_local_cp[cp_offset],0,sizeof(double)*symbol_path_select_size*path_select_offset);
  memset(&symbol_timer_pad_local_pp[cp_offset],0,sizeof(double)*(pp_symbol_class_count*num_per_process_measures+1));
  memset(&symbol_timer_pad_local_vol[vol_offset],0,sizeof(double)*(vol_symbol_class_count*num_volume_measures+1));
  this->has_been_processed = false;
}

void symbol_tracker::start(double save_time){
  if (symbol_stack.size()>0){
    auto last_symbol_time = save_time-symbol_timers[symbol_stack.top()].start_timer.top();
    for (auto i=0; i<symbol_path_select_size; i++){
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-1] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_exclusive_measure[i][num_per_process_measures-2] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-2] += last_symbol_time;
      symbol_timers[symbol_stack.top()].cp_excl_measure[i][num_per_process_measures-1] += last_symbol_time;
    }
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
    symbol_timers[symbol_stack.top()].pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  }
  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){
    critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer);
  }
  symbol_stack.push(this->name);
  computation_timer = MPI_Wtime();
  this->start_timer.push(computation_timer);
}

void symbol_tracker::stop(double save_time){
  assert(this->start_timer.size()>0);
  auto last_symbol_time = save_time-this->start_timer.top();
  for (auto j=0; j<symbol_path_select_size; j++){
    this->cp_exclusive_measure[j][num_per_process_measures-1] += last_symbol_time;
    this->cp_exclusive_measure[j][num_per_process_measures-2] += last_symbol_time;
    this->cp_excl_measure[j][num_per_process_measures-1] += last_symbol_time;
    this->cp_excl_measure[j][num_per_process_measures-2] += last_symbol_time;
  }
  this->pp_exclusive_measure[num_per_process_measures-1] += last_symbol_time;
  this->pp_exclusive_measure[num_per_process_measures-2] += last_symbol_time;
  this->pp_excl_measure[num_per_process_measures-1] += last_symbol_time;
  this->pp_excl_measure[num_per_process_measures-2] += last_symbol_time;
  *this->pp_numcalls += 1.; *this->vol_numcalls += 1.;
  for (auto i=0; i<num_per_process_measures; i++){
    this->pp_exclusive_contributions[i] += this->pp_exclusive_measure[i];
  }

  for (auto j=0; j<symbol_path_select_size; j++){
    this->cp_numcalls[j][0] += 1.;
    for (auto i=0; i<num_per_process_measures; i++){
      this->cp_exclusive_contributions[j][i] += this->cp_exclusive_measure[j][i];
    }
  }
  for (auto j=0; j<symbol_path_select_size; j++){
    memset(this->cp_exclusive_measure[j],0,sizeof(double)*num_per_process_measures);
  }
  memset(this->pp_exclusive_measure,0,sizeof(double)*num_per_process_measures);
  auto save_symbol = symbol_stack.top();
  this->start_timer.pop(); symbol_stack.pop();
  if (symbol_stack.size() > 0 && (save_symbol != symbol_stack.top())){
    for (auto j=0; j<symbol_path_select_size; j++){
      for (auto i=0; i<num_per_process_measures; i++){
        symbol_timers[symbol_stack.top()].cp_exclusive_contributions[j][i] += this->cp_exclusive_contributions[j][i];
        this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
      }
    }
    for (auto i=0; i<num_per_process_measures; i++){
      symbol_timers[symbol_stack.top()].pp_exclusive_contributions[i] += this->pp_exclusive_contributions[i];
      this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
    }
    for (auto j=0; j<symbol_path_select_size; j++){
      memset(this->cp_exclusive_contributions[j],0,sizeof(double)*num_per_process_measures);
    }
    memset(this->pp_exclusive_contributions,0,sizeof(double)*num_per_process_measures);
  } else if (symbol_stack.size() == 0){
    for (auto j=0; j<symbol_path_select_size; j++){
      for (auto i=0; i<num_per_process_measures; i++){
        this->cp_incl_measure[j][i] += this->cp_exclusive_contributions[j][i];
      }
    }
    for (auto i=0; i<num_per_process_measures; i++){
      this->pp_incl_measure[i] += this->pp_exclusive_contributions[i];
    }
    for (auto j=0; j<symbol_path_select_size; j++){
      memset(this->cp_exclusive_contributions[j],0,sizeof(double)*num_per_process_measures);
    }
    memset(this->pp_exclusive_contributions,0,sizeof(double)*num_per_process_measures);
  }
  critical_path_costs[num_critical_path_measures-2] += (save_time - computation_timer);		// update critical path computation time
  critical_path_costs[num_critical_path_measures-1] += (save_time - computation_timer);		// update critical path runtime
  volume_costs[num_volume_measures-2]        += (save_time - computation_timer);		// update local computation time
  volume_costs[num_volume_measures-1]        += (save_time - computation_timer);		// update local runtime
  for (size_t i=0; i<comm_path_select_size; i++){ critical_path_costs[critical_path_costs_size-1-i] += (save_time - computation_timer); }
  computation_timer = MPI_Wtime();
  if (symbol_stack.size()>0){ symbol_timers[symbol_stack.top()].start_timer.top() = computation_timer; }
}

}
}
}
